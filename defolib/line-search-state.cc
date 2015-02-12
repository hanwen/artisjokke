#include <string.h>
#include <math.h>
#include <stdio.h>

#include "node.hh"
#include "convergence-statistics.hh"
#include "mesh-feature.hh"
#include "big-vector.hh"
#include "mesh-connectivity.hh"
#include "line-search-state.hh"
#include "deformation-state.hh"
#include "setting.hh"

// #define PRINT_BETA_CHOICE
// #define PRINT_NEWTON
// #define PARANOIA



Line_search_state::Line_search_state ()
{
  no_more_improvement_b_ = false;
  restart_b_ = true;
  iter_count_ = 0;

  residual_tolerance_ = get_number_setting ("outer-loop-tolerance");

  
  char const * line_search = get_string_setting ("line-search");
  char const *selected = "";
  if (!strcmp (line_search, "secant"))
    {
      selected = "secant";
      step_length_function_ = &Line_search_state::secant_method;
    }
  else if (!strcmp (line_search, "linear"))
    {
      step_length_function_ = &Line_search_state::linear_alpha;
      selected = "linear";
    }
  else if (!strcmp (line_search, "semi-secant"))
    {
      selected = "semi-secant";
      step_length_function_ = &Line_search_state::semi_secant_newton_method;
    }
  else
    {
      selected = "newton";
      step_length_function_ = &Line_search_state::newton_method;
    }

  newton_tolerance_ = get_number_setting ("nonlinear-cg-newton-tolerance");
  printf ("Line search: %s, tol=%lf\n", selected, newton_tolerance_); 
}


void
Line_search_state::set_new_search_dir (Deformation_state*def)
{
  /*
    Steepest descent:
   */
  search_dir_ = def->constrained_residual_; 
  
}

  

/*
  COST:

  8*n flops
 */
void
Line_search_state::iteration_step (Deformation_state  *def) 
{
  int n = def->dimension ();

  global_iteration_count ++;
  
  set_new_search_dir (def);


  Real *dir = search_dir_.unsafe_access_array ();
  Real *disp =  def->displacements_.unsafe_access_array ();
  
#ifdef PARANOIA
  Real *residu = def->constrained_residual_.unsafe_access_array ();
  Real descent_dir = big_vector_ip (residu, dir, n);
  assert (descent_dir > 0);
#endif

  Real alpha;
  
#ifdef PRINT_BETA_CHOICE 
  if (polak_ribiere)
    {
      printf ("Polak-Ribiere\n");
    }
  else
    {
      printf ("Fletcher-Reeves\n");
    }
#endif

  alpha = (*step_length_function_)(this, def);

#if 0
  Real dirlen = big_vector_length (dir, n);
  printf ("Alpha: %lf, dirlen %lf\n", alpha, dirlen);
#endif
    
  if (alpha)
    {
      big_vector_axpy (disp, alpha, dir, disp, n);
      big_vector_copy (def->elastic_force_.unsafe_access_array(),
		       try_elastic_force_.access_array(), n);

      big_vector_copy (prev_residual_.unsafe_access_array (),
		       def->constrained_residual_.access_array (), n);
      big_vector_copy (def->constrained_residual_.unsafe_access_array(),
		       try_residual_.access_array(), n);

      
      def->set_residual_len (big_vector_length_squared (try_residual_.access_array(), n));
    }
  else
    {
      no_more_improvement_b_ = true;
    }
}

int
Line_search_state::dimension ()const
{
  return try_residual_.size ();
}

/*
  Calculate residual and 2nd derivative for SEARCH_DIR, leaving results in

  TRY_DISPLACEMENT
  MAT_TIMES_SEARCH_DIR
  TRY_RESIDUAL_

  and  lengths in

  DERIV, RESLEN
 */
void
Line_search_state::calculate_try_residual_and_derivative (Real alpha, Deformation_state const*def)
{
  Real *matdir = mat_times_search_dir_.unsafe_access_array ();
  Real *try_disp =  try_displacement_.unsafe_access_array ();
  Real * try_residual = try_residual_.unsafe_access_array ();
  Real * try_elastic_force = try_elastic_force_.unsafe_access_array ();  
  
  Real * dir = search_dir_.unsafe_access_array ();
  Real const * disp = def->displacements_.access_array ();
  int n =dimension ();

  if (alpha == 0.0)
    {
      big_vector_copy (try_disp, def->displacements_.access_array (), n);
      big_vector_copy (try_elastic_force, def->elastic_force_.access_array(), n);
      big_vector_copy (try_residual, def->constrained_residual_.access_array(), n );
      def->apply_force_derivative_and_residual (matdir, 0, 0, try_disp, dir);
    }
  else
    {
      big_vector_axpy  (try_disp, alpha, dir, disp, n );
      big_vector_nullify (matdir, n);
      big_vector_nullify (try_residual, n);
      def->apply_force_derivative_and_residual (matdir, try_elastic_force, try_residual, try_disp, dir);
    }
}

/*
  Interestingly, nonlinear CG doesn't make the residual length go down
  monotonously. We should look at the energy measure instead, -- that
  should go down.
  
 */
Real
Line_search_state::newton_method (Line_search_state *me,
				   Deformation_state const *def)
{
  int j_max  = 4;
  int j = 0;
  int n =me->dimension ();
  Real alpha = 0.0;
  Real const *d = me->search_dir_.access_array ();
  Real delta_0 = 0.0;

#ifdef PRINT_NEWTON
  Real dirlen = big_vector_length (d,n);
#endif
  
  while (j < j_max)
    {
      me->calculate_try_residual_and_derivative (alpha, def);
       
      Real delta = big_vector_ip (me->try_residual_.access_array (), d, n)
	/ big_vector_ip (d, me->mat_times_search_dir_.access_array (), n);
      if (j == 0)
	delta_0 =fabs(delta);
      
      
#ifdef PRINT_NEWTON
      printf ("iter %d alpha %lf delta %lf dirlen %lf, normalized alpha %lf\n",
	      j, alpha, delta, dirlen, alpha / dirlen);
#endif
      

      if (fabs (delta)  < me->newton_tolerance_ * delta_0)
	break;

      alpha += - delta;      
      
      j++;
    }

  return alpha;
}

/*
  Use Newton step for the first iteration, so we have a scale free measure for alpha_0 
*/
Real
Line_search_state::semi_secant_newton_method (Line_search_state *me,
					      Deformation_state const *def)
{
  int j_max  = 4;
  int j = 0;
  int n =me->dimension ();
  Real alpha = 0.0;
  Real const *dir = me->search_dir_.access_array ();
  Real delta_0 = 0.0;
  Real last_eta =0.0;
  Real last_alpha =0.0;
  Real const * try_res = me->try_residual_.access_array();

#ifdef PRINT_NEWTON
  Real dirlen = big_vector_length (dir, n);
#endif
  
  while (j < j_max)
    {
      Real delta =0.0;
      Real eta =0.0;      
      if (j == 0)
	{
	  me->calculate_try_residual_and_derivative (alpha, def);
	  eta = big_vector_ip (try_res, dir, n);
	  delta = - eta / big_vector_ip (dir, me->mat_times_search_dir_.access_array (), n);
	}
      else
	{
	  me->calculate_try_residual (alpha, def);
	  eta = big_vector_ip (try_res, dir, n);
	  delta = - (eta * (last_alpha - alpha)) / (last_eta - eta);
	}
      
      if (j == 0)
	delta_0 = fabs(delta);
	  
      
      
#ifdef PRINT_NEWTON
      printf ("iter %d, eta %lf, leta %lf\n", j, eta, last_eta);
      printf ("alpha %lf, lalpha %lf, delta %lf\n", alpha, last_alpha, delta);
#endif

      
      if (fabs (delta)  < me->newton_tolerance_ * delta_0)
	break;

      last_alpha = alpha;
      last_eta = eta;
      
      alpha +=delta;      
      
      j++;
    }

  return alpha;
}


/*
  We're doing Newton iteration, approximating the derivative with a
  differential. For finding the first differential, we should set
  alpha_0 , alpha_1 = small. Unfortunately, "small" is not yet
  scale-free, but set at 0.001, which is smaller than alpha values
  observed in practice.
 */
Real
Line_search_state::secant_method (Line_search_state *me,
				   Deformation_state const*def)
{
  int n = me->dimension();

  Real *dir = me->search_dir_.unsafe_access_array ();
  Real *res = me->try_residual_.unsafe_access_array ();

  Real last_alpha = -0.0;
  Real last_eta = big_vector_ip (def->constrained_residual_.access_array (),  dir, n);


  /*
    Small, but not scale free. (urg.)
   */
  Real alpha = 0.001;

  Real delta = 1.0;
  int j_max = 5;
  int j = 0;

  Real delta_0 = 0.0; 
  while (j < j_max )
    {
      me->calculate_try_residual (alpha, def);      
      Real eta =  big_vector_ip (res, dir,n);

      delta = (eta * (last_alpha - alpha)) / (last_eta - eta);

#ifdef PRINT_NEWTON
      printf ("iter %d, eta %lf, leta %lf\n", j, eta, last_eta);
      printf ("alpha %lf, lalpha %lf, delta %lf\n", alpha, last_alpha, delta);
#endif
      
      last_alpha = alpha;

      if (j == 0)
	{
	  delta_0 = fabs (alpha -delta); 
	}
      else if (fabs (delta) < delta_0 * me->newton_tolerance_ || j == j_max)
	{
	  break;
	}

      alpha = alpha - delta;
      last_eta = eta;

      j++;
    }

  return alpha;
}


void
Line_search_state::calculate_try_residual (Real alpha, Deformation_state const *def)
{
  Real *try_disp =  try_displacement_.unsafe_access_array ();
  Real * try_residual = try_residual_.unsafe_access_array ();
  Real * try_u_residual = try_elastic_force_.unsafe_access_array ();  
  Real * dir = search_dir_.unsafe_access_array ();
  Real const * disp = def->displacements_.access_array ();

  int n =dimension ();
  big_vector_axpy  (try_disp, alpha, dir, disp, n );
  big_vector_nullify (try_residual, n);  
  def->compute_elastic_force (try_u_residual, try_disp);
  def->constrained_residual_from_elastic_force (try_residual, try_u_residual);
}

Real
Line_search_state::linear_alpha (Line_search_state *me, Deformation_state const*def)
{
  Real * matdir = me->mat_times_search_dir_.unsafe_access_array();
  Real const *dir  = me->search_dir_.access_array();

  me->calculate_try_residual_and_derivative (0, def);

  /*
    todo: update me for reusing the above result and putting it back in DEF.
   */
  int n =me->dimension();
  Real alpha = - big_vector_ip (def->constrained_residual_.access_array (), dir, n) / big_vector_ip (matdir, dir, me->dimension( ));

  return alpha;
}



void
Line_search_state::signal_residual_change (Deformation_state *def)
{
  restart_b_ = true;
  no_more_improvement_b_ = false;
  iter_count_ = 0;      
  /*
    It is not necessary to clear SEARCH_DIR_ or
    MAT_TIMES_SEARCH_DIR_. This is implied by ITER_COUNT_ == 0
   */
}

void
Line_search_state::print_nodes () const
{
}

void
Line_search_state::print () const
{

}


extern bool contains (Array<int> const * is, int j);

Link_array< Array<Real> >
Line_search_state::vector_variables ()
{
  Array<Real> *vecs[] = {
    &mat_times_search_dir_, &search_dir_,
    &try_displacement_, &try_residual_,
    &try_elastic_force_, 
    &prev_residual_, 0,
  };

  Link_array< Array<Real> > vvs;
  for (int i= 0; vecs[i]; i++)
    vvs.push (vecs[i]);
  return vvs;
}

Line_search_state::~Line_search_state ()
{
}


/****************************************************************/

#if 0

void
Line_search_state::graph_alpha (Real start_step, Deformation_state const*def)
{
  printf ("start step %lf\n", start_step );
  int n = def->dimension();
  static int graphcount ;
  char filename[1024];
  
  printf("graph %d\n", graphcount);
  sprintf (filename, "tab%d", graphcount);

  FILE *gr =fopen (filename, "w");
  
  // plot graph.
  const int STEP_LIMIT= 40;
  Real start_graph= -0.5*start_step;
  Real stop_graph = 1.1*start_step;
  Real alpha = start_graph;
  Real delta = ((stop_graph - start_graph) / STEP_LIMIT); 
  int j =0;
  Real last_res = 0;
  do
    {
      Real res, deriv;
      calculate_try_residual_and_derivative (&res, &deriv, alpha, def);

      Real new_residual = def->compute_constrained_residual_force (try_residual_.unsafe_access_array (),
								      try_displacement_.unsafe_access_array ());

      if (j)
	{
	  Real deriv_est =  (sqrt (new_residual) - sqrt (last_res)) / delta;
	  fprintf (gr, "%f %f\n", alpha,  deriv_est);
	}
      last_res = new_residual;
      alpha += delta;
    }
  while (j ++ < STEP_LIMIT);
  fclose (gr);

  char s[1024];
  sprintf (s, "gnuplot %s", filename);
  system (s);

  graphcount++;
}
#endif

bool
Line_search_state::iteration_updates_residual() const
{
  return true;
}



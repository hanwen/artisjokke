#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "vector-io.hh"
#include "deformation-hook.hh"
#include "deformation-state.hh"
#include "mechanics.hh"
#include "mesh-connectivity.hh"
#include "mesh-feature.hh"
#include "big-vector.hh"
#include "setting.hh"
#include "convergence-statistics.hh"
#include "geometry.hh"
#include "element-state.hh"

const Real MEGA_REDUX = 1e6;
const int FRAME_RATE = 25;

Deformation_state::Deformation_state (int dim)
  : constraints_ (this)
{
  ignore_degenerate_ = false;
  iteration_frequency_ = 0;
  last_node_count_ = 0;
  external_force_len_sq_ =0.0;
  displacements_changed_ = true;
  iteration_updates_residual_ = false;
  tolerance_ = 0.0;
  diameter_ = 0.0;
  spatial_dimension_ = dim;
  forces_changed_ = false;

  elasticity_ = 0;

  constrained_residual_len_sq_ =0.0;
  prev_residual_len_sq_ = 0.0;
  tolerance_ = get_number_setting ("outer-loop-tolerance");
  char const * s = get_string_setting ("elasticity");  
  if (!strcmp (s, "nonlinear")
      || !strcmp (s, "venantkirchoff")
      )
    {
      elasticity_ = get_picinbono_elasticity (dim); 
    }
  else if (!strcmp (s, "picinbono2"))
    {
      elasticity_ = get_picinbono2_elasticity (dim);
    }
  else if (!strcmp (s, "linear"))
    {
      elasticity_ = get_linear_elasticity (dim);
      elasticity_->linear_ = true;
    }
  else if (!strcmp (s, "neohooke"))
    {
      elasticity_ = get_neo_hookean_elasticity (dim);
    }
  else if (!strcmp (s, "veronda"))
    {
      elasticity_ = get_veronda_elasticity (dim);
    }
  else if (!strcmp (s, "consistent-veronda"))
    {
      elasticity_ = get_consistent_veronda_elasticity (dim);
    }
  else
    {
      fprintf (stderr, "Unknown elasticity type");
      exit (2);
    }

  live_ = true;
}


/*
  Sync the vectors with the node structure.
 */
void
Deformation_state::update_topology (set<Element*> *changes)
{
  if (!changes)
    return ;
  
  for (iterof(i,*changes); i != changes->end(); i++)
    {
      Element *e = *i;
      
      if (!e->valid ())
	states_.erase (e);
      else
	{
	  Element_state *st = 0; 
	  if (!has_key (states_,e))
	    {
 	      st = get_element_state (spatial_dimension_);
	      states_[e] = st;
	    }
	  else
	    st = states_[e];

	  st->precompute (e, this);
	}
    }

  state_array_.set_size (states_.size());
  int j= 0;
  for (iterof(tp, states_); tp != states_.end() ; tp++)
    state_array_[j++] = tp->second;
}


Link_array< Array<Real> >
Deformation_state::interpolated_vector_variables ()
{
  Array<Real> *vecs[] = {
    &displacements_, 0
  };

  Link_array< Array<Real> > vvs;
  for (int i =0; vecs[i]; i++)
    vvs.push (vecs[i]);
  
  return vvs;
}

Link_array< Array<Real> >
Deformation_state::other_vector_variables ()
{
  Array<Real> *vecs[] = {
    &external_force_, 
    &reference_locations_,
    &statistic_vector1_, &statistic_vector2_,
    &elastic_force_, &constrained_residual_,    0
  };

  Link_array< Array<Real> > vvs;
  for (int i =0; vecs[i]; i++)
    vvs.push (vecs[i]);
  
  return vvs;
}


Link_array< Array<Real> >
Deformation_state::vector_variables ()
{
  Link_array< Array<Real> > vvs;

  vvs = other_vector_variables( );
  vvs.concat (interpolated_vector_variables ());
  
  return vvs;
}

void
Deformation_state::completize_nodal_arrays (int new_count)
{
  int n = spatial_dimension_ * new_count;
  Link_array<Array<Real> > vecs = vector_variables ();
  for (int i= 0; i < vecs.size (); i++)
    {
      completize_big_vector (vecs[i],n);
    }

  last_node_count_=  new_count;  
}

void
Deformation_state::set_residual_len (Real l)
{
  prev_residual_len_sq_ = constrained_residual_len_sq_;
  constrained_residual_len_sq_ = l;
}

/*
  update FORCE 
 */
void
Deformation_state::update_forces ()
{
  Real * ef = external_force_.unsafe_access_array();
  Real *ucr =elastic_force_.unsafe_access_array();
  Real *cr = constrained_residual_.unsafe_access_array(); 
    
  int n = dimension();


  bool give_signal = false;
  if (forces_changed_)
    {
      constraints_.apply_to_movement (ef, ef);
      constrained_residual_from_elastic_force (cr, ucr);
      constrained_residual_len_sq_ = big_vector_length_squared (cr, n);
      
      external_force_len_sq_ = big_vector_length_squared (ef,n);
      forces_changed_ = false;

      give_signal = true;
    }

  /*
    ugh.

    Should have a more waterproof way of monitoring changes requiring
    restart.
   */
  bool const_ch = constraints_.changed ();

  /*
    We assume that the routines below adapt to the latest settings for
    constraints.
   */
  constraints_.changed_ = false;
      

  if (!iteration_updates_residual_
      || displacements_changed_
      || const_ch)
    {
      compute_elastic_force (ucr, displacements_.access_array ());
      constrained_residual_from_elastic_force (cr, ucr);
      constrained_residual_len_sq_ = big_vector_length_squared (cr, n);
      displacements_changed_ = false;
      give_signal = true;      
    }

  if (give_signal)
    signal_residual_change ();
}

void
Deformation_state::constrained_residual_from_elastic_force (Real*cr, Real*ucr) const
{
  int n = dimension();
  big_vector_add (cr, ucr,  external_force_.access_array(), n);
  constraints_.apply_to_movement (cr, cr);
}

void
Deformation_state::set_reference_location (Node*nod, Real const *newpos)
{
  int n = spatial_dimension_ * (nod->number() + 1);
  if (reference_locations_.size ()  < n)
    {
      /*
	This is slightly hairy. We don't to spend too much time on
	resizing this vector over and over during mesh generation, so
	we resize it one big chunk. Note that it is vital that
	dimension() doesn't look at reference_locations_ for the
	problem size.
       */
      completize_big_vector (&reference_locations_,  2 * n + 1);
    }
  
  set_in_big_vector (nod, reference_locations_.unsafe_access_array(),
		     newpos, spatial_dimension_);
}

void
set_in_big_vector (Node *nod, Real * bv, Real const *v, int spatial_dimension)
{
  int idx = nod->number () * spatial_dimension;  
  for (int j=0; j < spatial_dimension; j++)
    {
      bv[idx + j] = v[j];
    }
}

void
Deformation_state::set_node_deformed_location (Node*nod, Real const * v)
{
  check_for_topology_update ();
  Real rv[3];
  int k = nod->number () * spatial_dimension_;
  for (int j = spatial_dimension_ ; j--;)
    rv[j] = v[j] - reference_locations_[k + j];
  
  set_in_big_vector (nod, displacements_.unsafe_access_array (),
		     rv, spatial_dimension_ ); 
  displacements_changed_ = true;
}


void
Deformation_state::add_force (Node * nod, Real const * df)
{
  int idx = nod->number()* spatial_dimension_;
  check_for_topology_update();


  Real *ef_arr = external_force_.unsafe_access_array();
  for (int j = spatial_dimension_; j--;)
    ef_arr[idx + j] += df[j];

  forces_changed_ = true; 
  constraints_.apply_to_node_movement (ef_arr + idx, nod, ef_arr + idx);
}

void
Deformation_state::set_force (Node * nod, Real const * df)
{
  int idx = nod->number()* spatial_dimension_;
  check_for_topology_update();


  Real *ef_arr = external_force_.unsafe_access_array();
  for (int j = spatial_dimension_; j--;)
    ef_arr[idx + j] = df[j];

  forces_changed_ = true; 
  constraints_.apply_to_node_movement (ef_arr + idx, nod, ef_arr);
}


int
Deformation_state::dimension () const
{
  return external_force_.size ();
}


void
Deformation_state::do_one_iteration()
{
}


void
Deformation_state::calibrate_deformation_speed ()
{
  Real calib_time = get_number_setting ("calibration-time");
  if (calib_time < 0)
    {
      iteration_frequency_ = -1;
      return ; 
    }
  else if (calib_time == 0)
    {
      iteration_frequency_ = 1;
      return ;
    }
  
  Real orig_tol = tolerance_;
  tolerance_ = 1e-18;		// want to get a good reading of the speed.

  
  iteration_frequency_ = 0;
  clock_t start = clock();
  clock_t now = start;


  clock_t delta_clock = clock_t (calib_time * CLOCKS_PER_SEC);
  long long fc = flop_count ();
  long long efc = element_flop_count;
  long long vfc = element_flop_count;
  
  do
    {
      bool finish = do_simulation_body ();

      if (finish)
	break ;
      iteration_frequency_++;
      now = clock();
    }
  while (now < start + delta_clock);

  fc = flop_count() - fc;
  vfc = vector_flop_count - vfc;
  efc = element_flop_count - efc;
  
  if (iteration_frequency_)
    {
      Real secs = Real (now - start) / Real (CLOCKS_PER_SEC);
      printf ("Calibrated %4.1lf MFLOPS (e %lf, v %4.1lf),  %d iters in %5.2f seconds\n",
	      fc *  1e-6,
	      efc * 1e-6,
	      vfc *1e-6,
	      iteration_frequency_,
	      secs);

      printf ("%4.1lf MFLOP/s %4.1lf iter/s %4.1lf MFLOP/iter\n",
	      fc / (1e6* secs),	      
	      iteration_frequency_/ secs,
	      (1e-6 *fc) / iteration_frequency_);
    }

  tolerance_ = orig_tol;
}

void
Deformation_state::flop_account_tick ()
{
  /*
    Flop accounting.
   */
  static int max_mflops = -1;
  if (max_mflops < 0)
    {
      max_mflops = (int) get_number_setting ("maximum-mflops");
    }
  if (max_mflops)
    {
      static Deformation_state * ref;
      if (!ref)
	ref = read_reference_state (mesh_);
      
      flop_stats(this, ref );

      if ( (element_flop_count * 1e-6) > max_mflops)
	{
#ifdef OPENGL
	  live_ = false;
#else  
	  exit (0);
#endif
	}
    }
}


/*
   Return value: solution reached? 
*/
bool
Deformation_state::do_simulation_body ()
{
  check_for_topology_update ();
  update_forces ();

  if (!live_)
    return false;

  bool good = good_solution ();


  if (!good)
    {
      do_one_iteration();
      for (int j = hooks_.size(); j--;)
	hooks_[j]->signal_iteration_done ();
    }
  
  flop_account_tick ();

  validate_numerical_sanity();

  return good;
}

void
Deformation_state::signal_residual_change ()
{
} 


Real
Deformation_state::scale_free_force_comparison () const
{
  Real t1= sqrt (external_force_len_sq_ +

		 constraints_.reaction_length_sq (elastic_force_.access_array()));


#if 0
  /*
    Ugh.

    We need this -- otherwise entry and exit take very much time:
    due to rounding errors, we get impossible-to-reach
    residual forces and tolerance 
   */
  t1  += lame_mu  * pow (diameter_, spatial_dimension_-1)
    * dimension () * 1e-10; 

  /*this is moved to the bottom....
   */
#endif

  return t1; 
}

bool
Deformation_state::good_solution ()const
{
  assert (!forces_changed_
	  && !displacements_changed_
	  && !constraints_.changed()
	  );

  static  Real eps_round;
  if  (!eps_round)
    eps_round = get_number_setting ("residual-rounding-tolerance");
  return sqrt (constrained_residual_len_sq_) <=
    (tolerance_ * scale_free_force_comparison ()
     + lame_mu  * pow (diameter_, spatial_dimension_-1) * sqrt (dimension ()) * eps_round);
}

void
Deformation_state::simulation_body ()
{
  if (!live_)
    return ; 
  
  if (!iteration_frequency_)
    calibrate_deformation_speed ();

  if (iteration_frequency_ < 0)
    {
      while (!do_simulation_body())
	;
      return;
    }
  
  for (int i = 1+ (iteration_frequency_ / FRAME_RATE); i--;)
    {
      if (do_simulation_body ())
	return ;
    }
}


void
Deformation_state::check_for_topology_update ()
{
  /*
    Check if vector sizes are up to date
  */
  int current_count = mesh_->node_count();

  bool top_change = 0;
  if (last_node_count_ != current_count)
    {
      completize_nodal_arrays (current_count);
      top_change = true;
    }

  if (set<Element*> *chs = get_changed_elements())
    {
      update_topology (chs);

      diameter_ = 1.0; 		// TODO
      
      delete chs;
      top_change = true;

      /*
	If the topology changes, residual changes too.  The following
	will force a restart before the next iteration.
       */
      forces_changed_ = true; 
    }

  if (constraints_.changed())
    {
      signal_boundary_condition_change();
    }

  if (top_change)
    for (int j = hooks_.size(); j--;)
      hooks_[j]->signal_topology_change ();
}

void
Deformation_state::signal_boundary_condition_change()
{
}

Deformation_state::~Deformation_state ()
{
  
}

long long element_flop_count = 0LL ;
int global_iteration_count;

void
Deformation_state::apply_force_function (Real *dest,  Real const * src) const
{
  big_vector_nullify (dest, dimension ());
  assert (elasticity_);
  
  for (int i = state_array_.size(); i--;)
    {
      Element_state * ts = state_array_[i];
      if (ignore_degenerate_ && ts->degenerate_b_)
	continue;

      element_flop_count += (ts)->elastic_force (elasticity_, dest, src);
    }
}

/*
  Calculate the derivative K(displ) * dir
*/
void
Deformation_state::apply_force_derivative (Real *deriv,
					   Real *elastic_force_vec,
					   Real *residual_vec,
					   Real const*displ, Real const *dir) const
{
  int n = dimension ();
  big_vector_nullify (deriv, n);
  if (residual_vec)
    {
      big_vector_nullify (residual_vec,n);
      big_vector_nullify (elastic_force_vec,n);
    }

  if (elasticity_->linear_)
    {
      /* linear case */
      apply_force_function (deriv, dir);
      if (residual_vec)
	apply_force_function (elastic_force_vec, displ);
    }
  else
    {
      for (int i = state_array_.size(); i--;)
	{
	  Element_state * ts = state_array_[i];
	  if (ignore_degenerate_ && ts->degenerate_b_)
	    continue;
	    
	  element_flop_count +=
	    ts->elastic_force_derivative (elasticity_, deriv,
					  elastic_force_vec,
					  displ, dir);
	}
    }    
}

void
Deformation_state::apply_force_derivative_and_residual (Real *deriv,
							Real *elastic_force_vec,
							Real * residual_vec,
							Real const*displ, Real const *dir) const
{
  /*
    elasticity
   */
  apply_force_derivative (deriv,
			  elastic_force_vec,
			  residual_vec, displ, dir);
  if (residual_vec)
    {
      constrained_residual_from_elastic_force (residual_vec, elastic_force_vec);
    }
  
  /*
    for consistency: we clear the deriv of fixed nodes as well. 
   */
  constraints_.apply_to_movement (deriv, deriv);
}

/*
  Should move into static-deformation-state?

  This gives the forces assuming an unconstrained situation.
*/
void
Deformation_state::compute_elastic_force (Real *residual_dest,
					  Real const *displacement) const
{
  /*
    Residu is  (stiffness  * displ )
  */
  apply_force_function (residual_dest, displacement); 
}



Real
Deformation_state::compute_constrained_residual_force (Real *residual_dest,
						       Real const *displacement) const
{
  compute_elastic_force (residual_dest, displacement);
  constrained_residual_from_elastic_force (residual_dest, residual_dest);

  return big_vector_length_squared (residual_dest, dimension());
}

void
Deformation_state::validate_numerical_sanity ()
{
  Real const *d = displacements_.access_array ();
  int n = displacements_.size();
  bool b = big_vector_sane (d, n);

  if (external_force_len_sq_)
    b = b &&  sqrt(constrained_residual_len_sq_ / external_force_len_sq_ ) < 1e10;

  if (!b)
    {
      fprintf(stderr, "Numerical La-La land reached. Kaboom.\n");
      exit (2);
    }
}

Real
Deformation_state::virtual_energy (Real const *disp) const
{
  Real e = 0.0;

  e = -big_vector_ip (external_force_.access_array (),
		      disp, dimension());
  
  for (int i = state_array_.size(); i--;)
    {
      Element_state * ts = state_array_[i];
      if (ignore_degenerate_ && ts->degenerate_b_)
	continue;

      e += ts->elastic_energy (elasticity_, disp);
      
      assert (!isnan (e));
    }
  return e;
}


Array<Real> 
Deformation_state::interpolate_deformation_variables (Real const * loc,
						      Element *e)
{
  check_for_topology_update ();
  
  Array<Real> vars;
  Real coef[MAX_DIMENSION + 1];
  int idx[MAX_DIMENSION + 1];

  for (int i = spatial_dimension_+ 1; i--;)
    {
      idx[i] = e->node(i)->number() *spatial_dimension_;
    }

  if (spatial_dimension_ ==3)
    {
      Vector3 v (loc[0], loc[1], loc[2]);
      Vector3 ps[MAX_DIMENSION + 1];
      for (int i = spatial_dimension_+ 1; i--;)
	ps[i] = reference_location3 (e->node(i), this);

      barycentric_coordinates3 (ps, coef, v);
    }
  else if (spatial_dimension_ == 2)
    {
      Vector2 v (loc[0], loc[1]);
      Vector2 ps[spatial_dimension_ + 1];
      for (int i = spatial_dimension_+ 1; i--;)
	ps[i] = reference_location2 (e->node(i), this);
      
      barycentric_coordinates2 (ps, coef, v);
    }
  
  Link_array< Array<Real> > bigvecs = vector_variables ();

  vars.set_size (spatial_dimension_ * bigvecs.size());
  for (int i = 0; i < bigvecs.size(); i++)
    {
      Vector3 v ;
      for (int j = 0; j < spatial_dimension_ + 1; j++)
	for (int k = 0; k < spatial_dimension_; k++)
	    v(k) += coef[j] * bigvecs[i]->elem(idx[j] + k);

      for (int j  = 0; j < spatial_dimension_; j++)
	vars[i*spatial_dimension_ + j] = v(j);
    }

  return vars;
}


void
Deformation_state::set_deformation_variables (Node* nod,
					      Array<Real> const &vars)
{
  check_for_topology_update ();
  
  int idx = nod->number() * spatial_dimension_;
  Link_array< Array<Real> > bigvecs = vector_variables ();
  for (int i = 0; i < bigvecs.size()  ; i++)
    {
      for (int k = spatial_dimension_; k--;)
	bigvecs[i]->elem_ref(idx + k) = vars[spatial_dimension_ * i + k];
    }
}

Deformation_state::Deformation_state( Deformation_state const &src)
  : Element_watcher (src),
    constraints_ (src.constraints_)
{
  assert(false);
}


void
Deformation_state::add_hook (Deformation_hook *h)
{
  hooks_.push (h);
  /*
    Should make special init function ? 
   */
  h->signal_topology_change ();
}

bool
Deformation_state::force_changed ()const
{
  return forces_changed_;
}

/*****************************************************************/

#include "static-deformation-state.hh"
#include "dynamic-deformation-state.hh"

Deformation_state*
make_new_deformation_state (Mesh_connectivity*top)
{
  Deformation_state *fem = 0;
  char const * pt = get_string_setting ("problem-type");
  if (!strcmp(pt , "static"))
    fem = new Static_deformation_state (top->dimension());
  else if (!strcmp( pt, "dynamic"))    
    fem = new Dynamic_deformation_state (top->dimension());
  else
    {
      if (pt)
	printf ("unknown deformation state type `%s'. Options: static and dynamic\n", pt);
      abort ();
    }

  top->add_watcher (fem);
  
  return fem;
}

void
Deformation_state::suicide()
{
  live_ = false;
}

/*
  GDB.
 */
extern "C" {
void
print_element (Element*e, Deformation_state*def)
{
  printf ("Element: ");
  for (int i = 0 ; i < e->simplex().count(); i++)
    {
      reference_location2(  e->node(i), def).print();
    }
}
}



#if 0
void
Deformation_state::update_incremental_step ()
{
  /*
    TODO: compute stiffness matrix for each edge and node.
   */
}
#endif

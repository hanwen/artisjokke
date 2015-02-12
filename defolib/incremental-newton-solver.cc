#include <stdio.h>
#include <math.h>

#include "incremental-newton-solver.hh"
#include "linear-cg-state.hh"
#include "deformation-state.hh"
#include "big-vector.hh"
#include "setting.hh"

int MAX_ITER;

Incremental_newton_solver::Incremental_newton_solver(Deformation_state *ps)
{
  if (!MAX_ITER)
    {
      MAX_ITER = (int) get_number_setting ("max-cg-iterations");
    }
  
  def_ = ps;
  cg_state_ = new Linear_cg_state;
  cg_residual_tolerance_ = get_number_setting("incremental-tolerance");
}

static void
compute_incremental_elastic_force (void *info, Real * force, Real const *dis)
{
  Incremental_newton_solver * me = (Incremental_newton_solver*)info;
  
  me->def_->apply_force_derivative_and_residual (force,
						    0, 0,
						    me->def_->displacements_.access_array(),
						    dis);

  /*
    We must negate since  (K(x) d) = - (ELASTICFORCE)
   */
  big_vector_negate (force,force, me->dimension());
  me->def_->constraints_.apply_to_movement (force,force);
}


Link_array< Array<Real> >
Incremental_newton_solver::vector_variables ()
{
  Array<Real> *vecs[]
    = { &direction_, &cg_elastic_force_ , & cg_constrained_residual_,
	&cg_right_hand_side_,	0};

  Link_array< Array<Real> > vvs;
  for (int i= 0; vecs[i]; i++)
    vvs.push (vecs[i]);

  vvs.concat (cg_state_->vector_variables());
  return vvs;
}


void
Incremental_newton_solver::iteration_step (Deformation_state *p)
{
  int n = def_->dimension();
  Real * cg_res = cg_constrained_residual_.unsafe_access_array();
  Real * cg_k_times_x = cg_elastic_force_.unsafe_access_array();
  
  Real * dir = direction_.unsafe_access_array ();

  Real reslen = big_vector_length_squared (p->constrained_residual_.access_array(), n);

  // def_->update_incremental_step ();

  printf ("Residual len : %lg\n", sqrt (reslen));
  
  cg_state_->cg_reinit();
  Real *rhs = cg_right_hand_side_ .unsafe_access_array();
  big_vector_nullify (dir, n);
  big_vector_add (rhs,
		  def_->external_force_.access_array(),
		  def_->elastic_force_ .access_array(), n);

  big_vector_nullify (cg_k_times_x, n);
  big_vector_copy (cg_res, rhs, n);
  def_->constraints_.apply_to_movement (cg_res,cg_res);
  
  
  Real cg_residual_len = reslen;
  Real prev_residual_len = 0.0;
  while (!MAX_ITER || cg_state_->iter_count_ < MAX_ITER)
    {
      cg_state_->cg_iteration_step (dir, cg_k_times_x, cg_res,
				    rhs,
				    cg_residual_len,
				    prev_residual_len,
				    &def_->constraints_,
				    (void*) this,
				    &compute_incremental_elastic_force);

      prev_residual_len = cg_residual_len;
      cg_residual_len = big_vector_length_squared (cg_res,  n);
      //      printf ("new reslen %lf\n", cg_residual_len);
      
      cg_state_->iter_count_ ++;

      bool good = cg_residual_len <= reslen * cg_residual_tolerance_;

      if (good)
	break ;
    }

  Real step = 1.0;
  
  if (big_vector_ip (dir,  p->constrained_residual_.access_array(), n) <0)
    {
      printf ("Huh? Not a descent direction.\n");
    }
  big_vector_axpy (p->displacements_.unsafe_access_array (),
		   step, dir,
		   p->displacements_.access_array (),
		   n);
}

/*
  Ugh, this is also more hairy than you'd think. This starts looping
  at the final step, I guess since roundoff errors mess up the
  comparison of energy values.

*/
Real
Incremental_newton_solver::find_descent_step (Real step)
{
  Real *curr_disp  = def_->displacements_.unsafe_access_array ();
  Real current_energy = def_->virtual_energy (curr_disp);
  Real *try_disp = cg_right_hand_side_.unsafe_access_array();
  Real * dir = direction_.unsafe_access_array();
  int n = dimension ();

  do 
    {
      big_vector_axpy (try_disp, step, dir, curr_disp, n);
      /*
      Real new_residual def_->compute_constrained_residual_force (try_res, try_disp);

      */
      Real new_energy = def_->virtual_energy (try_disp);
      if (new_energy  <= current_energy)
	return step;

      step *= 0.75;
      printf  ("Trying step %lf\n", step);
    }
  while (1);
  
}




Incremental_newton_solver::~Incremental_newton_solver()
{
  
}

bool
Incremental_newton_solver::iteration_updates_residual ()const
{
  return false;
}

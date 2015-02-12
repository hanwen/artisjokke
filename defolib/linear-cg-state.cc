#include <string.h>
#include <stdio.h>

#include <math.h>

#include "setting.hh"
#include "node.hh"
#include "mesh-feature.hh"
#include "big-vector.hh"
#include "debug.hh"
#include "mesh-connectivity.hh"
#include "linear-cg-state.hh"
#include "deformation-state.hh"
#include "convergence-statistics.hh"

// #define PARANOIA

Linear_cg_state::Linear_cg_state ()
{
  iter_count_ = 0;
}

void
deformation_gradient_function (void *ps, Real * dest,  Real const*src)
{
  Deformation_state *def = (Deformation_state*) ps;
  def->apply_force_function (dest, src);

  big_vector_negate (dest, dest, def->dimension());
}

void
Linear_cg_state::iteration_step (Deformation_state *def)
{
  int n = def->dimension();
  Real *rp = def->constrained_residual_.unsafe_access_array ();
  Real *elastic =def->elastic_force_.unsafe_access_array();

  cg_iteration_step (def->displacements_.unsafe_access_array (),
		     elastic,
		     rp, def->external_force_.access_array(),
		     def->get_residual_len_sq (),
		     def->prev_residual_len_sq_,
		     &def->constraints_,
		     (void*)def,
		     &deformation_gradient_function);

  Real rl =big_vector_length_squared (rp, n);
  def->set_residual_len (rl);
}


/*
  This code solves

  POSITIVE_DEFINITE_OPERATOR x = RHS

  where A is a POSITIVE_DEFINITE_OPERATOR.
  
 */
void
Linear_cg_state::cg_iteration_step (Real *displacement,
				    Real * elastic,
				    Real * c_residu,
				    Real const * rhs,
				    Real residual_len_sq, Real prev_len_sq,
				    Deformation_constraints* constraints,
				    void * optor_info,
				    void (*positive_definite_operator)(void*, Real* ,Real const*))
{
  global_iteration_count ++;
  
  int n = search_dir_.size ();
  
  if (iter_count_==0)
    log_message ("Restart!\n");
  else if (iter_count_ > n)
    {
      log_message ("Restart after %d iters!\n", iter_count_);      
      iter_count_ = 0;
    }
    
  Real beta =  (iter_count_ == 0)
    ?  0.0 : residual_len_sq / prev_len_sq;

  iter_count_ ++;
  Real *dir = search_dir_.unsafe_access_array ();
  Real *matdir = mat_times_search_dir_.unsafe_access_array ();

  big_vector_axpy (dir, beta, dir, c_residu, n);

#if 0
  assert (def->constraints_.satisfies_constraints (dir));
#endif

  (*positive_definite_operator) (optor_info, matdir, dir);
  Real alpha = residual_len_sq / big_vector_ip (matdir, dir, n);

#ifdef PARANOIA
  assert (!isnan (alpha) && !isinf (alpha));
#endif
  
  big_vector_axpy (displacement, alpha, dir, displacement, n);
  big_vector_axpy (elastic, - alpha, matdir, elastic, n);

#if 0
  log_message ("alpha, beta %lf, %lf\n", alpha, beta);
#endif
  
  big_vector_add (c_residu, rhs, elastic, n);
  if (constraints)
    constraints->apply_to_movement (c_residu, c_residu);
  
  /*
    updating residulen is done from a different routine; it also has
    to be done for "correct" solutions */
}

void
Linear_cg_state::signal_residual_change (Deformation_state *)
{
  iter_count_ = 0;
}


void
Linear_cg_state::cg_reinit ()
{
  if (iter_count_)
    log_message ("Reinit after %d iters\n", iter_count_);
  iter_count_ =0;
}


Link_array< Array<Real> >
Linear_cg_state::vector_variables ()
{
  Array<Real> *vecs[]
    = {
      &mat_times_search_dir_, &search_dir_,
      0};
  
  Link_array< Array<Real> > vvs;
  for (int i= 0; vecs[i]; i++)
    vvs.push (vecs[i]);
  return vvs;
}


bool
Linear_cg_state::iteration_updates_residual()const
{
  return true;
}

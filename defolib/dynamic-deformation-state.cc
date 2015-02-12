/*   
  dynamic-deformation-state.cc --  implement Dynamic_deformation_state

  (c) 2001 Han-Wen Nienhuys <hanwen@cs.uu.nl>
*/
#include <algorithm>

#include <sys/time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "setting.hh"
#include "big-vector.hh"
#include "dynamic-deformation-state.hh"
#include "element-state.hh"
#include "mechanics.hh"
#include "node.hh"
#include "simplex.hh"
#include "mesh-geometry.hh"
#include "convergence-statistics.hh"

using std::min;

/*
  Compute acceleration assuming viscosity_ * lumped_mass * velocity
  for damping.
 */
void
Dynamic_deformation_state::compute_acceleration (Real * acc,
						 Real const* velocity,
						 Real const* displacement)
{
  int n = dimension ();
  big_vector_nullify (acc, n);
  compute_constrained_residual_force (acc, displacement);

  Real * m = scratch4_.unsafe_access_array ();

  big_vector_pointwise_multiply (m, lumped_masses_.access_array(), velocity, n);
  big_vector_axpy (acc, -viscosity_, m, acc, n);

  
  Real * ms = inverse_masses_.unsafe_access_array ();
  big_vector_pointwise_multiply (acc, ms, acc, n);

  constraints_.apply_to_movement (acc, acc);

}

void
Dynamic_deformation_state::compute_acceleration_with_given_residual (Real * residual,
								     Real const* velocity,
								     Real const* displacement)
{
  Real *acc = residual;
  int n = dimension ();

  {
    Real * m = scratch4_.unsafe_access_array ();

    big_vector_pointwise_multiply (m, lumped_masses_.access_array(), velocity, n);
    big_vector_axpy (acc, -viscosity_, m, acc, n);
  }
  
  Real * m = inverse_masses_.unsafe_access_array ();
  big_vector_pointwise_multiply (acc, m, acc, n);

  constraints_.apply_to_movement (acc, acc);

  /*
    displacement is needed for "springs" with spring-damping.
   */
  (void)displacement;
}


/****************************************************************
   TIME INTEGRATION
 ****************************************************************/

/*
  Simple euler forward.
 */
void
Dynamic_deformation_state::euler_forward (Dynamic_deformation_state*me,
					  Real delta_t)
{
  int n = me->dimension ();

  Real *v = me->velocity_.unsafe_access_array ();
  Real *f = me->scratch1_.unsafe_access_array ();
  Real *d = me->displacements_.unsafe_access_array ();
  
  me->compute_acceleration_with_given_residual (me->constrained_residual_.unsafe_access_array(), v, d);
  big_vector_axpy (d, delta_t, v, d, n);
  big_vector_axpy (v, delta_t, f, v, n);

  me->constraints_.apply_to_movement (v, v);
}

void
Dynamic_deformation_state::classical_runge_kutta (Dynamic_deformation_state*me,Real delta_t)
{
  int n = me->dimension ();

  Real * a1 = me->scratch0_.unsafe_access_array ();
  Real * v = me->velocity_.unsafe_access_array ();
  Real * d = me->displacements_.unsafe_access_array ();
  Real * v1 =me->scratch1_.unsafe_access_array ();
  Real * d1 = me->scratch2_.unsafe_access_array ();
  Real * kv = me->scratch3_.unsafe_access_array ();
  Real * ka = me->scratch4_.unsafe_access_array ();

  me->constraints_.apply_to_movement (v,v);
  big_vector_copy (v1, v, n);
  big_vector_copy (d1, d, n);
  Real * cr = me->constrained_residual_.unsafe_access_array();
  me->compute_acceleration_with_given_residual (cr, v1, d1);
  big_vector_copy (ka, cr, n);
  big_vector_copy (kv, v1, n);

  big_vector_axpy (d1, delta_t* 0.5, v1, d1, n);
  big_vector_axpy (v1, delta_t* 0.5, a1, v1, n);  
  me->compute_acceleration (a1, v1, d1);
  big_vector_axpy (ka, 2.0, a1, ka, n);
  big_vector_axpy (kv, 2.0, v1, kv, n);
  
  big_vector_axpy (d1, delta_t* 0.5, v1, d1, n);
  big_vector_axpy (v1, delta_t* 0.5, a1, v1, n);  
  me->compute_acceleration (a1, v1, d1);
  big_vector_axpy (ka, 2.0, a1, ka, n);
  big_vector_axpy (kv, 2.0, v1, kv, n);

  big_vector_axpy (d1, delta_t, v1, d1, n);
  big_vector_axpy (v1, delta_t, a1, v1, n);  
  me->compute_acceleration (a1, v1, d1);
  big_vector_add (ka, ka, a1, n);
  big_vector_add (kv, kv, v1, n);

  big_vector_axpy (d, 1.0/6.0 * delta_t, kv, d, n);
  big_vector_axpy (v, 1.0/6.0 * delta_t, ka, v, n);  

#ifndef NDEBUG
  Real vlen = big_vector_ip (v, v, n);
  assert (!isinf (vlen) && !isnan (vlen));
#endif

  global_iteration_count ++;
}  


/*
  Explicit scheme from the Picinbono ICRA 2001 paper.

  This scheme is a SS22 multistep order 2 with theta_1=1/2 and theta_2 = 0.

  (see Zienkewicz vol. 1, 5th ed, page 530)

  Conditionally stable for time-step delta_t < h * sqrt(density/stiffness).

  We have h = Oh(n^{1/3})
  
  
  Note that the viscosity in their paper (gamma_i) is negative, while
  ours is positive.

  COST: 14 * n
*/
void
Dynamic_deformation_state::explicit_ss22_multistep (Dynamic_deformation_state*me,Real delta_t)
{
  int n = me->dimension ();
  Real * inv_factor = me->scratch0_.unsafe_access_array ();
  Real * this_t_contribution = me->scratch1_.unsafe_access_array ();
  Real * prev_t_contribution = me->scratch2_.unsafe_access_array ();
  Real * force_contribution = me->constrained_residual_ .unsafe_access_array ();
  Real * prev_t_displacements = me->velocity_.unsafe_access_array ();
  Real * d = me->displacements_.unsafe_access_array ();
  Real * m = me->lumped_masses_.unsafe_access_array ();

  Real * residual = me->constrained_residual_.unsafe_access_array();

  Real * elastic =me->elastic_force_.unsafe_access_array();
  me->compute_elastic_force (elastic, d);
  me->constrained_residual_from_elastic_force (residual, elastic);
  me->set_residual_len (big_vector_length_squared (residual, n));

  
  big_vector_pointwise_multiply (prev_t_contribution, m, prev_t_displacements, n);
  big_vector_scale (prev_t_contribution, (2 - me->viscosity_ * delta_t)/(2*sqr (delta_t)),
		    prev_t_contribution, n);

  big_vector_pointwise_multiply (this_t_contribution, m, d, n);
  big_vector_scale (this_t_contribution, 2.0 / sqr(delta_t), this_t_contribution, n);
  
  big_vector_add (this_t_contribution, this_t_contribution, force_contribution,n );
  big_vector_subtract (this_t_contribution, this_t_contribution, prev_t_contribution, n);

  big_vector_scale (inv_factor, (2.0 + me->viscosity_ * delta_t) /(2*sqr(delta_t)), m, n );
  big_vector_pointwise_inverse (inv_factor, inv_factor, n);

  big_vector_pointwise_multiply (this_t_contribution, this_t_contribution, inv_factor,n);
  big_vector_copy (prev_t_displacements, d,n);
  big_vector_copy (d, this_t_contribution, n);
}


/*
  Single step SS22, equivalent to GN22.

  Theta_1=1/2 is required to make the stability limit

  dt_crit proportional to h.  For Theta_1 > 1/2, we have dt_crit =
  Oh(h^2), which is a disadvantage.
  
 */
void
Dynamic_deformation_state::explicit_ss22_singlestep (Dynamic_deformation_state*me,
						     Real delta_t)
{
  Real theta_1 = 0.5;

  int n = me->dimension();
  Real * displ = me->displacements_.unsafe_access_array();
  Real * residual = me->constrained_residual_.unsafe_access_array();
  Real * uc_residual =me->elastic_force_.unsafe_access_array();
  Real * mean_displ = me->scratch0_.unsafe_access_array();
  Real * force_contribution = me->scratch1_.unsafe_access_array ();
  Real * velocity = me->velocity_. unsafe_access_array();
  Real * inverse_factor = me->scratch4_.unsafe_access_array();
  Real * masses = me->lumped_masses_.unsafe_access_array();

  big_vector_axpy (mean_displ, theta_1* delta_t, velocity, displ, n);
  
  me->compute_elastic_force (uc_residual, mean_displ);
  me->constrained_residual_from_elastic_force (residual, uc_residual);
  me->set_residual_len (big_vector_length_squared (residual, n));

  big_vector_scale (force_contribution, me->viscosity_, masses, n);
  big_vector_pointwise_multiply (force_contribution, force_contribution, velocity, n);
  big_vector_subtract (force_contribution, force_contribution, residual, n);
  
  big_vector_scale (inverse_factor, -(theta_1 * delta_t * me->viscosity_ + 1.0), masses, n);
  big_vector_pointwise_inverse (inverse_factor, inverse_factor, n);
  big_vector_pointwise_multiply (force_contribution, inverse_factor, force_contribution, n);

  me->constraints_.apply_to_movement (velocity,velocity);
  me->constraints_.apply_to_movement (force_contribution,force_contribution);
  

  
  big_vector_axpy (displ, 0.5 * sqr (delta_t),  force_contribution, displ, n);
  big_vector_axpy (displ, delta_t,  velocity, displ, n);

  big_vector_axpy (velocity, delta_t, force_contribution, velocity, n);
}




/****************************************************************/


Link_array< Array<Real> >
Dynamic_deformation_state::other_vector_variables ()
{
  Array<Real> *vecs[] = {
      &lumped_masses_, &inverse_masses_,
      &scratch0_, &scratch1_, &scratch2_, &scratch3_, &scratch4_,
      0 };

  Link_array< Array<Real> > vvs (Deformation_state::other_vector_variables ());
  for (int i =0; vecs[i]; i++)
    vvs.push (vecs[i]);
  
  return vvs;
}

Link_array< Array<Real> >
Dynamic_deformation_state::interpolated_vector_variables ()
{
  Array<Real> *vecs[] = {
    &velocity_, 0 };

  Link_array< Array<Real> > vvs (Deformation_state::interpolated_vector_variables ());
  for (int i =0; vecs[i]; i++)
    vvs.push (vecs[i]);
  
  return vvs;
}


/*
  We don't use CPU time for the time-integration, because it is not
  constant, certainly not on our platform.
*/

void
Dynamic_deformation_state::do_one_iteration ()
{
  (*time_integration_function_)(this, time_step_);

  time_ += time_step_;
}

Real
Dynamic_deformation_state::now()const
{
  return time_;
}

Dynamic_deformation_state::Dynamic_deformation_state (int d)
  : Deformation_state (d)
{
  viscosity_ = 1.0;
  density_ = 1.0;
  time_integration_function_ = 0;
  time_step_ = 1.0;
  time_ = 0.0;

  viscosity_ = get_number_setting ("viscosity");
  density_ = get_number_setting ("density");

  char const * vis = get_string_setting ("viscosity");
  char const * elas = get_string_setting ("elasticity");
  char const * integr = get_string_setting ("time-integration");
  if( !strcmp(integr, "RK4"))
    time_integration_function_ = &classical_runge_kutta;
  else if (!strcmp(integr, "eulerfw"))
    time_integration_function_ = &euler_forward;
  else if (!strcmp (integr, "ms22"))
    {
      iteration_updates_residual_ = true;            
      time_integration_function_ = &explicit_ss22_multistep;
    }
  else
    {
      if (!strcmp (integr, "ss22"))
	{
	  printf ("No integration specified, using SS22");
	}
      
      iteration_updates_residual_ = true;      
      time_integration_function_ = &explicit_ss22_singlestep;
    }
  
  
  printf (" using: viscosity = %s, elasticity=  %s,\ntime-integration = %s\n",
	   vis, elas, integr);


  init_lame_parameters();


  time_step_ = 0.0;
  
  
#if 0
  /*
    Fix this up:
   */
  
  Real crit_damping  = sqrt ( get_number_setting ("young") * get_number_setting ("density"));
  printf ("Damping relative to critical = %f\n", get_number_setting ("viscosity") / crit_damping);
#endif
}



Real
Dynamic_deformation_state::wave_speed ()const
{
  Real young = get_number_setting ("young");

  return sqrt (young / get_number_setting ("density"));
}

Real
Dynamic_deformation_state::critical_time_step () const
{
  if (time_integration_function_ == &explicit_ss22_multistep ||
      time_integration_function_ == &explicit_ss22_singlestep
      )
    return explicit_ss22_critical_time_step ();
  else if (time_integration_function_ == &euler_forward)
    return euler_forward_critical_time_step ();
  else if (time_integration_function_ == &classical_runge_kutta)
    {
    /*
      I think the stability is equal to euler FW stability for decay
      equations.  */
      assert (false);
      return 1.0;
    }
  
  return 0.0;
}

Real
Dynamic_deformation_state::explicit_ss22_critical_time_step () const
{
  /*
    Courant Friedrichs Lewy condition
    
    assuming lumped masses (see Zienkewicz vol I p 520).
    
    return sqrt (2) * minimum_element_length (fem_) / wave_speed ();

    Probably there is some issue with lumping, or something else that
    I don't quite understand, which makes me want to use the formula
    including sqrt(2).

    http://scienceworld.wolfram.com/physics/P-Wave.html
    
   */
  return sqrt(2) * minimum_edge_length_ *  sqrt (get_number_setting ("density")/ (lame_lambda + 2 * lame_mu));  
}

Real
Dynamic_deformation_state::euler_forward_critical_time_step ()const
{
  /*
   Zienkewicz p 498, 501, 502

   2D Euler FW is a single step Taylor serie collocation  of Euler FW 1D.

   Critical time step is Oh(h^2), which isn't all that good for FE computations.


   (? weird: lame parameters are no part of this?? )

   --note: p 498--502 apply to 1st order systems. The C and K relate to  to
   heat conductivity stuff.
   

  */

  return sqr (minimum_edge_length_) * viscosity_ / density_ ;
}

bool
Dynamic_deformation_state::good_solution ()
{
  return false;
}


void
Dynamic_deformation_state::update_topology (set<Element*> *changes)
{
  Deformation_state::update_topology (changes);

  minimum_edge_length_ = infty;

  /*
    This is slightly ugh, but hey, a relaxation step is Oh(N) anyway.
   */
  bool ign_degen = get_bool_setting ("ignore-degenerate-elements");
  
  for (int i = state_array_.size() ; i--; )
    {
      if (!state_array_[i]->degenerate_b_
	  || !ign_degen)
	minimum_edge_length_ = min (minimum_edge_length_,
				    state_array_[i]->minimum_edge_length_);
      
      assert (minimum_edge_length_ > 0);
    }
  time_step_ = get_number_setting ("time-step") * critical_time_step ();
  printf ("New timestep %lg\n", time_step_); 

  
  /*
    the new nodes should have meaningful positions for MS methods.
  */
  if (time_integration_function_ == &explicit_ss22_multistep)
    {
      /*
	UGH.  This is

	 Wrong, only newly added nodes should have copied displacements_,
       */
      velocity_ = displacements_;
    }
 
  int n = lumped_masses_.size();
  big_vector_nullify (lumped_masses_.unsafe_access_array (),n);

  Array<Real> lv (lumped_volumes (mesh_, this));

  //  total_mass / fem_->node_count ();
  for (int i = lv.size () ; i--;)
    for (int j = 3;j--;)
      {
	lumped_masses_[3*i + j] =  lv[i] * density_; 
      }

  big_vector_pointwise_inverse (inverse_masses_.unsafe_access_array (),
				lumped_masses_.access_array (),
				n);
}


Vector2
Dynamic_deformation_state::node_velocity2 (Node* n)const
{
  assert (time_integration_function_ == &explicit_ss22_singlestep);
  return extract_vector2 (n, velocity_.access_array());
}

void
Dynamic_deformation_state::set_node_velocity (Node* n , Real const *v)
{
  assert (time_integration_function_ == &explicit_ss22_singlestep);
  for (int i = 0; i < 3; i++)
    velocity_[i+spatial_dimension_*n->number()] = v[i];
}

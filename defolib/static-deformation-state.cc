#include <math.h>
#include <stdio.h>

#include "deformation-hook.hh"
#include "static-deformation-state.hh"
#include "mesh-connectivity.hh"
#include "mesh-feature.hh"
#include "linear-cg-state.hh"
#include "big-vector.hh"
#include "nonlinear-cg-state.hh"
#include "node.hh"
#include "setting.hh"
#include "incremental-newton-solver.hh"
#include "debug.hh"
#include "mechanics.hh"
#include "convergence-statistics.hh"
#include "element-state.hh"

Static_deformation_state::Static_deformation_state (int d)
  : Deformation_state (d)
{
  iter_state_ = 0;
  iteration_updates_residual_ = false;
  switch_to_newton_ = get_bool_setting ("compute-exact-solution");
  
  char const * relt = get_string_setting ("relaxation-type");
  if (!strcmp (relt, "linear"))
    {
      if (strcmp (get_string_setting ("elasticity"), "linear"))
	fprintf (stderr, " **** Requested non-linear elasticity and linear CG. Using nonlinear CG instead.;\n");
      else
	{
	  set_iteration_state ( new Linear_cg_state ());

	  iteration_updates_residual_ = !get_bool_setting ("exact-residual-update");
	}
    }
  else if (!strcmp (relt, "incremental-newton"))
    {
      set_iteration_state (new Incremental_newton_solver (this));
    }
  else if (strcmp (relt, "nonlinear"))
    {
      printf ("No such solver: using nonlinear CG");
    }
  
  if (!iter_state_)
    {
      set_iteration_state ( new Nonlinear_cg_state ());
    }
}

void
Static_deformation_state::update_topology (set<Element*> *changes)
{
  Deformation_state::update_topology (changes);
  
  /*
    In a static simulation, every component must be fixed to ascertain a
    unique solution.
  */
  
  //   fix_loose_components(this, top);
}

void
Static_deformation_state::signal_residual_change ()
{
  iter_state_->signal_residual_change (this);
}


void
Static_deformation_state::set_iteration_state (Iteration_state * iter)
{
  assert (iter_state_ != iter);
    
  delete iter_state_;
  iter_state_ = iter;
  
  iteration_updates_residual_ = iter_state_->iteration_updates_residual();

  completize_nodal_arrays (dimension()/spatial_dimension_);
}

void
Static_deformation_state::signal_boundary_condition_change ()
{
  iter_state_->signal_residual_change (this);
}

int global_max_iteration_count;

bool
Static_deformation_state::good_solution ()const
{
  bool b = Deformation_state::good_solution ();
  
  if (b && iter_state_->iter_count_
       && elasticity_->linear_)
    {
      /*
	UGH UGH. should be provided for by Deformation_state member
	function.
       */
      Deformation_state * def = (Deformation_state*) this; 
      Real *ucr = def->elastic_force_.unsafe_access_array();
      Real *cr  = def->constrained_residual_.unsafe_access_array ();
      compute_elastic_force (ucr, def->displacements_.access_array ());
      constrained_residual_from_elastic_force (cr, ucr);
      def->set_residual_len (big_vector_length_squared (cr, dimension ()));

      log_message ("Restarting against roundoff error after %d iters.\n", iter_state_->iter_count_);
      log_message ("reslen: %g  rhs:  %g\n", get_residual_len_sq(),
		   scale_free_force_comparison());

      b = Deformation_state::good_solution ();
      if (!b)
	{
	  log_message ("Not good enough, %d iters\n", iter_state_->iter_count_);
	  global_max_iteration_count = global_max_iteration_count >? iter_state_->iter_count_;
	  
	  iter_state_->iter_count_ = 0;
	}
    }

  if (b && iter_state_->iter_count_)
    {
      log_message ("Solution found in %d iterations\n ",  iter_state_->iter_count_);

      global_max_iteration_count = global_max_iteration_count >? iter_state_->iter_count_;
      
      log_message ("Scale free force comparison: %lg", scale_free_force_comparison ());
	
      iter_state_->iter_count_ = 0;

      for (int j = hooks_.size(); j--;)
	hooks_[j]->signal_converged();       
    }
  return b;
}

void
Static_deformation_state::do_one_iteration ()
{
  iter_state_->iteration_step (this);


  /*
    Find exact solution.
   */
  if (switch_to_newton_
      && get_residual_len_sq()
      && get_residual_len_sq() < 0.1  * external_force_len_sq_)
    {
      printf ("switching to Newton.\n");
      switch_to_newton_ = false;
      set_iteration_state (new Incremental_newton_solver (this));
    }
}

Link_array< Array<Real> >
Static_deformation_state::other_vector_variables ()
{
  Link_array<Array <Real> > vvs =  Deformation_state::other_vector_variables ();
  vvs.concat (iter_state_->vector_variables ());
  return  vvs;
}

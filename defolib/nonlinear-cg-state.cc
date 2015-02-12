/*
  Non linear CG .
  
*/

#include <string.h>
#include <stdio.h>
#include <math.h>

#include "debug.hh"
#include "setting.hh"
#include "big-vector.hh"
#include "nonlinear-cg-state.hh"
#include "deformation-state.hh"


Nonlinear_cg_state::Nonlinear_cg_state ()
  :  Line_search_state()
{
  use_polak_ribiere_ = use_hestenes_stiefel_ = use_fletcher_reeves_ = false;
  
  char const * beta_func = get_string_setting ("nonlinear-cg-betafunc");
  if (!strcmp(beta_func, "polak-ribiere"))
    {
      printf ("Using PR beta selection.");    
      use_polak_ribiere_= true;
      use_fletcher_reeves_ = false;
    }
  if (!strcmp(beta_func, "hestenes-stiefel"))
    {
      printf ("Using HS beta selection.");    
      use_hestenes_stiefel_ = true;
      use_polak_ribiere_= false;
      use_fletcher_reeves_ = false;
    }
  else if (!strcmp (beta_func, "fletcher-reeves"))
    {
      printf ("Using FR beta selection.");
      use_fletcher_reeves_ = true;
      use_polak_ribiere_ = false;
    }
  else
    {
      printf ("Using both PR and FR.");
      use_polak_ribiere_ = use_fletcher_reeves_ = true;
    }
}

Link_array< Array<Real> >
Nonlinear_cg_state::vector_variables ()
{
  Array<Real> *vecs[] = {
    &prev_residual_, 0,
  };

  Link_array< Array<Real> > vvs = Line_search_state::vector_variables();
  for (int i= 0; vecs[i]; i++)
    vvs.push (vecs[i]);
  return vvs;
}

void
Nonlinear_cg_state::set_new_search_dir (Deformation_state *def)
{
  static int restart_count = (int) get_number_setting ("nonlinear-cg-restart-count");
  if (restart_count < 0)
    restart_count = dimension();

  if (restart_count && (iter_count_ -last_restart_ >  restart_count))
     {
      restart_b_ = true;
       no_more_improvement_b_ = false;
     }

  
  Real beta =  0.0;

  /*
    In the past we had some troubles with Polak-Ribiere where
    Fletcher-Reeves performed better. They turned out to be related to
    programming errors -- optimizing an `incorrect' function gives
    interesting results.
  */
  if (no_more_improvement_b_)
    restart_b_ = true;

  int n = dimension();
  bool polak_ribiere = false;
  if (!restart_b_)
    {
      if (use_polak_ribiere_)
	{
	  /*
	    Polak-Ribiere method.
	  */
	  beta = (def->get_residual_len_sq ()  -
		  big_vector_ip (def->constrained_residual_.access_array (),
				 prev_residual_.access_array (), n)) / def->prev_residual_len_sq_;
	  polak_ribiere = true;
	}

      if (use_hestenes_stiefel_)
	{
	  Real const* cr = def->constrained_residual_.access_array ();
	  Real const* pr = prev_residual_.access_array ();
	  beta=
	    (def->get_residual_len_sq() - big_vector_ip (cr, pr, n) )/
	    (big_vector_ip(search_dir_.access_array (), pr, n)
	     - big_vector_ip (search_dir_.access_array(), cr , n ));
	     
	}
      
      if (use_fletcher_reeves_)
	{
	  /*
	    Fletcher-Reeves method -- can never result in negative beta.
	  */
	  if (beta <= 0.0)
	    {
	      beta = def->get_residual_len_sq() / def->prev_residual_len_sq_;
	      polak_ribiere = false;
	    }
	}

      if (beta < 0)
	{
	  beta = 0.0; // restarting
	  printf( "backing up!\n");
	}
    }
  bool restarted = ! beta;
  if (restarted)
    {
      log_message (" *** restarting. Residual force = %lf\n", sqrt (def->get_residual_len_sq()));
      last_restart_ = iter_count_;
    }
  
  restart_b_ = false;
  iter_count_ ++;

  Real *dir = search_dir_.unsafe_access_array ();
  Real *residu = def->constrained_residual_.unsafe_access_array ();
  
  big_vector_axpy (dir, beta, dir, residu,n);
  def->constraints_.apply_to_movement (dir, dir);  
}

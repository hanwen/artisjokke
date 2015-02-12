#include <time.h>
#include <math.h>

#include "deformation-state.hh"
#include "auto-inserter3.hh"
#include "needle-inserter3.hh"
#include "setting.hh"

Auto_needle_inserter3::Auto_needle_inserter3()
{
  Real h = pow (2.0, - get_number_setting ("refinement-level") / 3.0) * 0.1;
  speed_ = get_number_setting ("auto-insert-speed") * h;
  
  handle_ = Vector3(-0.01, get_number_setting ("auto-insert-y"),
		    get_number_setting ("auto-insert-z"));

  angle_ = get_number_setting ("auto-insert-angle");
  max_depth_ = get_number_setting ("auto-insert-depth");
  depth_ = 0.0;
  dir_ = Vector3 (1,0,0); // dir_vector (M_PI * get_number_setting ("auto-insert-angle")/ 180.0);
}

/*
   should move out of Needle_inserter.
 */

clock_t start_insert_time;   

void
Needle_inserter3::auto_needle_insert ()
{
  last_handle_ = auto_insert_->handle_;
  while (state_ == OUTSIDE)
    {
      auto_insert_->depth_ += auto_insert_->speed_;
      
      Vector3 t = auto_insert_->handle_
	+ auto_insert_->depth_ * auto_insert_->dir_;

      move_needle (last_handle_, t);


      if (state_ == INSIDE )
	start_insert_time = clock ();
    }
  
  do
    {
      if (state_ == INSIDE
	  && auto_insert_->depth_ < auto_insert_->max_depth_)
	{
	  auto_insert_->depth_ += auto_insert_->speed_;
	  Vector3 t = auto_insert_->handle_
	    + auto_insert_->depth_ * auto_insert_->dir_;

	  move_needle (last_handle_, t);
	}

      /*
	Something must have changed. Check if we should start
	calculating again.
      */
      deformation ()->update_forces();

      /*
	Ok, we've really changed something. Move control to either
	visualization or relaxation.
       */
      if (!deformation ()->good_solution ())
	break ; 

      /*
	It's still OK, then we're finished with the auto-insert.
       */
      if (auto_insert_->depth_ >= auto_insert_->max_depth_)
	{
	  // finish_auto_insert(mesh_, deformation (), this);

	  clock_t dt = clock ()  - start_insert_time;
	  Real secs = dt / (1.0* CLOCKS_PER_SEC);
	  printf ("Insertion time: %5.2lf seconds\n", secs );
	  printf ("Updates %d, avg update frequency: %lf\n", rearrange_count_,
		  rearrange_count_/secs);

#ifndef OPENGL
	  exit (0);
#else
	  auto_insert_ = 0;
#endif
	  
	  break ; 
	}
      
    } while (1);
}


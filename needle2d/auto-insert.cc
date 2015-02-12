#include <math.h>

#include "deformation-state.hh"
#include "auto-insert.hh"
#include "needle-inserter.hh"
#include "setting.hh"

Auto_needle_insert::Auto_needle_insert()
{
  speed_ = get_number_setting ("auto-insert-speed")
    * get_number_setting ("refinement-h");
  
  handle_ = Vector2(get_number_setting ("auto-insert-x"),
		    get_number_setting ("auto-insert-y")
		    );
  angle_ = get_number_setting ("auto-insert-angle");
  max_depth_ =get_number_setting ("auto-insert-depth");
  depth_ = 0.0;
  dir_ = dir_vector (M_PI * get_number_setting ("auto-insert-angle")/ 180.0);
}

/*
   should move out of Needle_inserter.
 */
void
Needle_inserter::auto_needle_insert ()
{
  last_handle_ = auto_insert_->handle_;
  while (state_ == OUTSIDE)
    {
      auto_insert_->depth_ += auto_insert_->speed_;
      
      Vector2 t = auto_insert_->handle_
	+ auto_insert_->depth_ * auto_insert_->dir_;

      move_needle (last_handle_, t);
    }
  
  do
    {
      if (state_ == INSIDE
	  && auto_insert_->depth_ < auto_insert_->max_depth_)
	{
	  auto_insert_->depth_ += auto_insert_->speed_;
	  Vector2 t = auto_insert_->handle_
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
	  finish_auto_insert(mesh_, deformation (), this);
#ifndef OPENGL
	  exit (0);
#else
	  auto_insert_ = 0;
#endif
	  break ; 
	}
      
    } while (1);

}


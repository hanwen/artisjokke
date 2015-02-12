#include <math.h>
#include <string.h>

#include "geometry.hh"
#include "mesh-feature.hh"
#include "setting.hh"
#include "bottom-fix-hook.hh"
#include "mesh-connectivity.hh"
#include "deformation-state.hh"

Bottom_fix_watcher::Bottom_fix_watcher(Deformation_state*def)
  : Deformation_hook (def)
{
  object_size_ = 0.1;		// ugh, should read from maubach-tree
  last_node_count_ = 0;
  def->get_mesh ()->add_watcher (this);
  //fix_state_ = RIGHTFIX;

  
  fix_state_ = 0;
  char const * ch[]= {"left", "right", "up", "down", "back", "front", 0};
  char const * fix = get_string_setting ("fix-planes");
  for (int i = 0; ch[i]; i++)
    {
      if (strstr (fix, ch[i]))
	fix_state_ = fix_state_ | (1<<i);
    }

  if (!fix_state_)
    fix_state_ = DOWNFIX;
}

Real EPS = 1e-14;

void
Bottom_fix_watcher::signal_topology_change ()
{
  if (set<Face*> *chf = get_changed_boundary())
    {
      for (iterof (i,*chf); i!=chf->end (); i++)
	{
	  Face * f = *i;
	  if (!f->valid ()
	      && f->mate ())
	    continue;

	  Vector3 n;
	  if (deformation ()->spatial_dimension () == 2)
	    n = simplex_normal2 (f->simplex(),
				 &reference_location2, deformation ());
	  else
	    n = simplex_normal3 (f->simplex(),
				 &reference_location3, deformation ());
	  bool fix = false;
	  fix = fix ||  (n(X_AXIS) > 0.5 && (fix_state_ & RIGHTFIX));
	  fix = fix ||  (n(X_AXIS) < -0.5 && (fix_state_ & LEFTFIX));
	  fix = fix ||  (n(Y_AXIS) > 0.5 && (fix_state_ & UPFIX));
	  fix = fix ||  (n(Y_AXIS) < -0.5 && (fix_state_ & DOWNFIX));
	  fix = fix ||  (n(Z_AXIS) < 0.5 && (fix_state_ & FRONTFIX));
	  fix = fix ||  (n(Z_AXIS) > -0.5 && (fix_state_ & BACKFIX));

	  if (fix)
	    for (int  i = f->simplex().count (); i--;)
	      deformation ()->constraints_.fix_node (f->node(i), true);
	}
    }
}


#ifndef BOTTOM_FIX_WATCHER_HH
#define BOTTOM_FIX_WATCHER_HH

#include "mesh-connectivity-watcher.hh"
#include "deformation-hook.hh"

class Bottom_fix_watcher :
  public Boundary_watcher,  public Deformation_hook
{
public:
  Bottom_fix_watcher (Deformation_state*);
  int last_node_count_;
  Real object_size_;
  enum { LEFTFIX = 0x1, RIGHTFIX=0x2, UPFIX=0x4, DOWNFIX=0x8, BACKFIX=0x10, FRONTFIX=0x20 };
  int fix_state_;   
  
  virtual void signal_topology_change ();
};

#endif

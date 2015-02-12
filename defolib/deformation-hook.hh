#ifndef DEFORMATION_HOOK_HH
#define DEFORMATION_HOOK_HH

#include "defo-proto.hh"

/*
Class to allow actions  when deformation routines are working.
Three types:

 -- when the deformation routines have processed a topology change
     (node/element count change), signal_topology_change() is called

 -- when a static iteration reaches the solution, signal_converged() is called

 -- when any iteration does a single step, signal_iteration_done() is called
 (note that this is called often, take care for efficiency.)

*/
class Deformation_hook
{

  Deformation_state * deformation_;

public:
  Deformation_state * deformation () const { return deformation_; }
  Deformation_hook (Deformation_state* def);
  
  virtual void signal_topology_change ();
  virtual void signal_converged ();  
  virtual void signal_iteration_done ();
  virtual ~Deformation_hook();
};

#endif

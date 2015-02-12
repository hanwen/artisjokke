#include "deformation-hook.hh"
#include "deformation-state.hh"

void
Deformation_hook::signal_iteration_done()
{
}

void
Deformation_hook::signal_topology_change ()
{

}

void
Deformation_hook::signal_converged()
{
}

Deformation_hook::~Deformation_hook()
{
}

Deformation_hook::Deformation_hook (Deformation_state*def)
{
  deformation_ = def;
  def->add_hook (this);
}

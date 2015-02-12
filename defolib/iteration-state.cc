#include <string.h>

#include <math.h>

/*
  #include "globals.hh"
#include "fem-model.hh"

*/
#include "node.hh"
#include "mesh-feature.hh"
#include "big-vector.hh"
#include "mesh-connectivity.hh"
#include "iteration-state.hh"
#include "deformation-state.hh"


Link_array<Array<Real> >
Iteration_state::vector_variables ()
{
  Link_array<Array<Real> > a;
  return a;
}

void
Iteration_state::iteration_step (Deformation_state  * ) 
{
}

void
Iteration_state::signal_residual_change (Deformation_state*)
{
  
}

Iteration_state::Iteration_state ()
{
  iter_count_ = 0;
}

bool
Iteration_state::good_solution (Deformation_state const * def)const
{
  return true;
}

void
Iteration_state::print_nodes () const
{
}

void
Iteration_state::print () const
{


}

 
extern bool contains (Array<int> const * is, int j);



Iteration_state::~Iteration_state ()
{
  
}

bool
Iteration_state::iteration_updates_residual()const
{
  return false;
}




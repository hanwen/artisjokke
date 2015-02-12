
#ifndef FRICTION_HH
#define FRICTION_HH

#include "array.hh"
#include "proto.hh" 
/*
  Should make virtual for more diverse friction mechanisms.
 */
struct Friction_definition
{
  Real peak_maximum_param;
  Real peak_length;
  Real peak_extra_force;
  Real force_density;
  Real entry_decrease_distance; // force at needle entry.
  
  Friction_definition();
};



void graph_friction_forces(Friction_definition const *fr);
Array<Real>
compute_friction_forces (Array<Real> needle_params,
			 Friction_definition const* friction);
			 
#endif

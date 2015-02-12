#ifndef MESH_GEOMETRY_HH
#define MESH_GEOMETRY_HH

#include "array.hh"
#include "mesh-connectivity.hh"

Array<Real> lumped_volumes (Mesh_connectivity*fem, Deformation_state*);

#endif

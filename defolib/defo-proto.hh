#ifndef DEFO_PROTO_HH
#define DEFO_PROTO_HH

#include "mesh-proto.hh"

#ifdef SINGLE_PRECISION
typedef float Real;
#else
typedef double Real;
#endif

const Real infty = 1e9;
const Real eps = 1e-9;

class Deformation_hook;
class Deformation_state;
class Deformation_constraints;
class Dynamic_deformation_state;
class Iteration_state;
class Line_search_state;
class Linear_cg_state;
class Matrix2;
class Matrix3;
class Node;
class Nonlinear_cg_state;
class Simplex;
class Static_deformation_state;
class Vector2;
class Vector3;
class Elasticity_functions;
class Element_state;

typedef Vector2 (* Node_vector_func2)(Node*, Deformation_state const*);
typedef Vector3 (* Node_vector_func3)(Node*, Deformation_state const*);
void init_deformation_module ();
void init_after_settings ();
#endif

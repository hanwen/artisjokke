#include "vector.hh"

#define ELEMENT_STATE Element_state2
#define elastic_energy consistent_veronda_elastic_energy2
#define elastic_force consistent_veronda_elastic_force2
#define elastic_force_derivative consistent_veronda_elastic_force_derivative2
#define MATRIX Matrix2
#define DIM 2

#include "consistent-veronda.icc"

#undef ELEMENT_STATE
#undef elastic_force
#undef elastic_force_derivative
#undef elastic_energy
#undef DIM
#undef MATRIX


#define ELEMENT_STATE Element_state3
#define elastic_energy consistent_veronda_elastic_energy3
#define elastic_force consistent_veronda_elastic_force3
#define elastic_force_derivative consistent_veronda_elastic_force_derivative3
#define MATRIX Matrix3
#define DIM 3

#include "consistent-veronda.icc"

IMPLEMENT_GET_ELASTICITY(consistent_veronda);

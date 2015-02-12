#ifndef VECTOR_HH
#define VECTOR_HH


#include "vector2.hh"
#include "vector3.hh"

#define DIM 3
#define VECTOR Vector3
#include "vector.icc"
#undef VECTOR
#undef DIM

#define DIM 2
#define VECTOR Vector2
#include "vector.icc"
#undef DIM
#undef VECTOR

#endif


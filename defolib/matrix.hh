#ifndef MATRIX_HH
#define MATRIX_HH

#include "vector.hh"

#include "matrix2.hh"
#include "matrix3.hh"

#define DIM 3
#define VECTOR Vector3
#define MATRIX Matrix3
#include "matrix.icc"
#undef VECTOR
#undef DIM
#undef MATRIX

#define DIM 2
#define VECTOR Vector2
#define MATRIX Matrix2
#include "matrix.icc"
#undef DIM
#undef VECTOR
#undef MATRIX

#endif



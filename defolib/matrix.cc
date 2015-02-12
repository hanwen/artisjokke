#define MATRIX_INLINE

#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "matrix.hh"
#include "vector.hh"

#define DIM 3
#define VECTOR Vector3
#define MATRIX Matrix3
#include "slow-matrix.icc"

#undef VECTOR
#undef DIM
#undef MATRIX

#define DIM 2
#define VECTOR Vector2
#define MATRIX Matrix2
#include "slow-matrix.icc"
#undef DIM
#undef VECTOR
#undef MATRIX





#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "matrix.hh"

#define DIM 2
// #define PARANOIA
/*
  COST:
 */
void
Matrix2::set_column_vectors (Vector2 const &v1, Vector2 const& v2)
{
  Vector2 const *vs[2] = {&v1, &v2};
  for (int j=0;  j < 2; j++)
    for (int k=0; k < 2; k++)
      {
	elts_[j][k] = vs[k]->elts_[j]; 
      }
}

void
Matrix2::set_row_vectors (Vector2 const &v1, Vector2 const& v2)
{
  Vector2 vs[2] = {v1, v2};
  set_row_vectors (vs);
}


void
Matrix2::set_row_vectors (Vector2 const *vs)
{
  for (int j=0;  j < 2; j++)
    for (int k=0; k < 2; k++)
      {
	elts_[j][k] = vs[j].elts_[k]; 
      }
}

Real
Matrix2::determinant() const
{
  return elts_[0][0]  * elts_[1][1] -  elts_[0][1]  * elts_[1][0];
}


void
Matrix2::cofactor_matrix(Matrix2 &dest, Matrix2 const &src)
{
  assert (&dest != &src);
  dest.elts_[0][0] = src.elts_[1][1];
  dest.elts_[1][1] = src.elts_[0][0];
  dest.elts_[1][0] = -src.elts_[1][0];
  dest.elts_[0][1] = -src.elts_[0][1];
}



/*
  can be more efficient: src.determinant () uses src.sub_determinant()

  PRECONDITION

  DEST != SRC


  (setq Matrix2::invert_to (+ Matrix2::cofactor_matrix Matrix2::scale Matrix2::determinant))

  COST: 48 flops

*/
void
Matrix2::invert_to (Matrix2 &dest , Matrix2 const &src)
{
  assert (&dest != &src);
  cofactor_matrix (dest,src);
  
  Matrix2::scale (dest, 1/ src.determinant (), dest);

#ifdef PARANOIA
  Matrix2 r, u;
  multiply_mm (r, src, dest);
  u.diag (1);
  u.subtract (r, u, r);
  if (r.max_norm () >  (1e-8) * src.max_norm ())
    {
      fprintf (stderr, "%s: norm = %lf\n", __PRETTY_FUNCTION__, r.max_norm());
    }
#endif
}

/*
  (setq Matrix2::invert_to_with_det (+ Matrix2::cofactor_matrix Matrix2::scale))
 */
void
Matrix2::invert_to_with_det (Matrix2 &dest , Matrix2 const &src, Real det)
{
  cofactor_matrix(dest, src);
  Matrix2::scale (dest, 1/ det, dest);
}


/*
  COST: 2 gonio ops = ? flops.  

  (setq Matrix2::set_rotate 9)   ? 
 */
void
Matrix2::set_rotate (int axis, Real angle)
{
  int a1 = (axis + 1)%DIM;
  int a2 = (axis + 2)%DIM;

  Real c = cos (angle);
  Real s = sin (angle);

  fill (0.0);

  elts_[a1][a1] = c;
  elts_[a1][a2] = s;
  elts_[a2][a2] = c;
  elts_[a2][a1] = -s;
  elts_[axis][axis] = 1.0;
}


void
Matrix2::to_opengl_matrix (Real *om, Matrix2 const &m)
{
  memset (om, 0, sizeof (Real)*15);
  om[15] = 1.0;			// homogenous coords.
  for (int i=0; i < DIM; i++)
    for (int j=0; j < DIM; j++)
      om[i + j * 4] = m.elts_[i][j];
}

void
Matrix2::from_opengl_matrix (Matrix2 *m, Real const *om)
{
  m->fill(0.0);

  for (int i=0; i < DIM; i++)
    for (int j=0; j < DIM; j++)
      m->elts_[i][j] = om[i + j * 4];
}

void
Matrix2::from_opengl_matrixf (Matrix2 *m, float const *om)
{
  m->fill(0.0);

  for (int i=0; i < DIM; i++)
    for (int j=0; j < DIM; j++)
      m->elts_[i][j] = om[i + j * 4];
}




Matrix2
gram_matrix2 (Vector2 v1, Vector2 v2)
{
  Real d1d1 = v1 * v1;
  Real d1d2 = v1 * v2;
  Real d2d2 = v2 * v2;            

  Matrix2 m;
  m.elts_[0][0] = d1d1;
  m.elts_[1][1] = d2d2;
  m.elts_[1][0] = d1d2;
  m.elts_[0][1] = d1d2;

  return m;
}



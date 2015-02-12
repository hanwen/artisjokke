
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "matrix.hh"

#define DIM 3

// #define PARANOIA
/*
  COST:
 */
void
Matrix3::set_column_vectors (Vector3 const &v1, Vector3 const& v2, Vector3 const& v3)
{
  Vector3 const *vs[3] = {&v1, &v2, &v3};
  for (int j=0;  j < 3; j++)
    for (int k=0; k < 3; k++)
      {
	elts_[j][k] = vs[k]->elts_[j]; 
      }
}

void
Matrix3::set_row_vectors (Vector3 const &v1, Vector3 const& v2, Vector3 const& v3)
{
  Vector3 vs[3] = {v1, v2, v3};
  set_row_vectors (vs);
}



void
Matrix3::set_row_vectors (Vector3 const *vs)
{
  for (int j=0;  j < 3; j++)
    for (int k=0; k < 3; k++)
      {
	elts_[j][k] = vs[j].elts_[k]; 
      }
}



/*
  cost: 3 flops

  (setq  Matrix3::sub_determinant 3)
  
*/
Real
Matrix3::sub_determinant (int i, int j) const
{
  Real a = elts_[(i+1)%DIM][(j+1)%DIM];
  Real b = elts_[(i+1)%DIM][(j+2)%DIM];
  Real c = elts_[(i+2)%DIM][(j+1)%DIM];
  Real d = elts_[(i+2)%DIM][(j+2)%DIM];

  int parity =  1 - 2 * ((i+j) % 2) ;
  return (a * d - b * c)*parity;
}

/*
  cost: 13 flops

  (setq  Matrix3::determinant (* (+ Matrix3::sub_determinant 1) 3))
 */
Real
Matrix3::determinant() const
{
  Real d = 0.0;
  for (int i=0; i < DIM; i++)
    {
      int parity = 1 - 2 * (i % 2);
      
      d += elts_[i][0]* sub_determinant (i, 0) * parity;
    }
  return d;
}


/*
  (setq Matrix3::cofactor_matrix (* Matrix3::sub_determinant 9))

  cost: 37
*/
void
Matrix3::cofactor_matrix(Matrix3 &dest, Matrix3 const &src)
{
  assert (&dest != &src);
  for (int i=0; i < DIM; i++)
    for (int j=0; j < DIM; j++)
      {
	int parity = 1- 2 * ((i+j) % 2);
	dest (i,j) = parity * src.sub_determinant (j,i);
      }
}



/*
  can be more efficient: src.determinant () uses src.sub_determinant()

  PRECONDITION

  DEST != SRC


  (setq Matrix3::invert_to (+ Matrix3::cofactor_matrix Matrix3::scale Matrix3::determinant))

  COST: 48 flops


*/
void
Matrix3::invert_to (Matrix3 &dest , Matrix3 const &src)
{
  cofactor_matrix(dest, src);
  Matrix3::scale (dest, 1/ src.determinant (), dest);

  // #define PARANOIA
  
#ifdef PARANOIA
  Matrix3 r, u;
  multiply_mm (r, src, dest);
  u.diag (1);
  u.subtract (r, u, r);
  if (r.max_norm () >  (1e-8) * src.max_norm ())
    {
      fprintf (stderr, "%s : norm = %lf\n", __PRETTY_FUNCTION__, r.max_norm());
    }
#endif
}

/*
  (setq Matrix3::invert_to_with_det (+ Matrix3::cofactor_matrix Matrix3::scale))
 */
void
Matrix3::invert_to_with_det (Matrix3 &dest , Matrix3 const &src, Real det)
{
  cofactor_matrix(dest, src);
  Matrix3::scale (dest, 1/ det, dest);
}



/*
  COST: 3 gonio ops = ? flops.  

  (setq Matrix3::set_rotate 9)   ? 
 */
void
Matrix3::set_rotate (int axis, Real angle)
{
  int a1 = (axis + 1)%3;
  int a3 = (axis + 3)%3;

  Real c = cos (angle);
  Real s = sin (angle);

  fill (0.0);

  elts_[a1][a1] = c;
  elts_[a1][a3] = s;
  elts_[a3][a3] = c;
  elts_[a3][a1] = -s;
  elts_[axis][axis] = 1.0;
}

  
Matrix3
Matrix3::rotate_about (Vector3 v, Real theta)
{
  v/= v.length ();
  Matrix3 r;
  r.set_rotate (0, theta);

  Vector3 n1, n3;
  n1(0) = 1;
  n1 = Vector3::cross (n1, v);
  Real l =n1.length ();
  if (l > 1e-8)
    {
      n1 /= l;
    }
  else
    {
      n1(1) = 1;
      n1 = Vector3::cross (n1, v);
      n1 /= n1.length ();
    }

  n3 = Vector3::cross (v, n1);

  Matrix3 B, Binv;
  B.set_column_vectors (v, n1,n3);

  invert_to (Binv, B);

  Matrix3 rotmat (r * Binv);
  rotmat = B * rotmat;

  return rotmat;
}

void
Matrix3::to_opengl_matrix (Real *om, Matrix3 const &m)
{
  memset (om, 0, sizeof (Real)*15);
  om[15] = 1.0;			// homogenous coords.
  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      om[i + j * 4] = m.elts_[i][j];
}

void
Matrix3::from_opengl_matrix (Matrix3 *m, Real const *om)
{
  m->fill(0.0);

  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      m->elts_[i][j] = om[i + j * 4];
}

void
Matrix3::from_opengl_matrixf (Matrix3 *m, float const *om)
{
  m->fill(0.0);

  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      m->elts_[i][j] = om[i + j * 4];
}


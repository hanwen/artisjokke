/*   
  matrix.hh -- declare Matrix
  
  (c) 1999 Han-Wen Nienhuys <hanwen@cs.uu.nl>
  
 */

#ifndef MATRIX2_HH
#define MATRIX2_HH

#include "vector.hh"

/*
  WARNING:

  In my infinite wisdom, I decided that this class does not have a
  default ctor , meaning that a Matrix2 will be filled with junk at
  instantiation.  */

struct Matrix2
{
  Real elts_[2][2];

  void fill (Real);
  void diag (Real);
  /**
     Make the matrix [VA|VB|VC]
   */

  void set_column_vectors (Vector2 const &va, Vector2 const &vb);
  void set_row_vectors (Vector2 const &va, Vector2 const &vb);
  void set_row_vectors (Vector2 const *);
  
  Real operator ()(int i, int j) const { return elts_[i][j]; } 
  Real &operator ()(int i, int j) { return elts_[i][j]; }
  Real trace () const;
  Real determinant () const;
  Real sub_determinant (int, int ) const;
  void print () const;
  Real max_norm () const;
  Real hilbert_schmidt_norm () const;
  Real max_row_sum () const;
  Real max_col_sum () const;

  void add_diag (Real);
  static Real inner_product (Matrix2 const &,Matrix2 const &);
  static void negate (Matrix2 &dest, Matrix2 const&);
  static void square (Matrix2 & dest, Matrix2 const &);
  static void multiply_mtm (Matrix2 &dest, Matrix2 const&m1, Matrix2 const&m2);
  static void multiply_mm (Matrix2 &dest, Matrix2 const&m1, Matrix2 const&m2);
  static void scale (Matrix2 &dest, Real r, Matrix2 const&m2);
  static void subtract (Matrix2 & dest, Matrix2 const &, Matrix2 const&);
  static void add (Matrix2 & dest, Matrix2 const &, Matrix2 const&);
  static void transpose_to (Matrix2 &dest, Matrix2 const&);
  static void invert_to (Matrix2 &dest, Matrix2 const&);
  static void invert_to_with_det (Matrix2 &dest, Matrix2 const&,Real);
  static void cofactor_matrix (Matrix2 &dest, Matrix2 const&);
  static void multiply_mv (Vector2 &dest, Matrix2 const &m, Vector2 const &rhs);
  static void multiply_vm (Vector2 &dest, Vector2 const &lhs, Matrix2 const &m);
  static void axpy (Matrix2 &dest, Real a,  Matrix2 const&, Matrix2 const&);
  /**
     Set rotation matrices; angle is in radians.
   */

  static void symmetrify (Matrix2 &dest, Matrix2 const &src);
  static void double_symmetrify (Matrix2 &dest, Matrix2 const &src);  
  
  void set_rotate (int axis,  Real angle);
  static void to_opengl_matrix (Real *rp, Matrix2 const&m);
  static void from_opengl_matrix (Matrix2 *m, Real const *rp );
  static void from_opengl_matrixf (Matrix2 *m, float const *rp );  
  void operator *= (Real r);
  static Matrix2 rotate_about (Vector2, Real);
  bool kosher () const;
};

inline Matrix2
operator * (Matrix2 const &m1, Matrix2 const &m2)
{
  Matrix2 ret;
  Matrix2::multiply_mm (ret, m1, m2);
  return ret;
}


inline Vector2
operator * (Matrix2 const &m1, Vector2 const &v)
{
  Vector2 w;
  Matrix2::multiply_mv (w, m1, v);
  return w;
}

inline Matrix2
operator +(Matrix2 const &m1, Matrix2 const &m2) 
{
  Matrix2 r;  
  Matrix2::add (r, m1,m2);
  return r;
}


inline Matrix2
operator -(Matrix2 const &m1, Matrix2 const &m2) 
{
  Matrix2 r;  
  Matrix2::subtract (r, m1,m2);
  return r;
}

inline Matrix2
operator * (Real r, Matrix2 const & m)
{
  Matrix2 ret;
  Matrix2::scale (ret, r, m);
  return ret;
}

Matrix2 gram_matrix2 (Vector2, Vector2);
#endif /* MATRIX_HH */


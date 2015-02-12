/*   
  matrix.hh -- declare Matrix
  
  (c) 1999 Han-Wen Nienhuys <hanwen@cs.uu.nl>
  
 */

#ifndef MATRIX3_HH
#define MATRIX3_HH

#include "vector.hh"

/*
  WARNING:

  In my infinite wisdom, I decided that this class does not have a
  default ctor , meaning that a Matrix3 will be filled with junk at
  instantiation.  */

struct Matrix3
{
  Real elts_[3][3];

  void fill (Real);
  void diag (Real);
  /**
     Make the matrix [VA|VB|VC]
   */

  void set_column_vectors (Vector3 const &va, Vector3 const &vb ,Vector3 const &vc);
  void set_row_vectors (Vector3 const &va, Vector3 const &vb ,Vector3 const &vc);
  void set_row_vectors (Vector3 const *);
  
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
  static Real inner_product (Matrix3 const &,Matrix3 const &);
  static void negate (Matrix3 &dest, Matrix3 const&);
  static void square (Matrix3 & dest, Matrix3 const &);
  static void multiply_mtm (Matrix3 &dest, Matrix3 const&m1, Matrix3 const&m3);
  static void multiply_mm (Matrix3 &dest, Matrix3 const&m1, Matrix3 const&m3);
  static void scale (Matrix3 &dest, Real r, Matrix3 const&m3);
  static void subtract (Matrix3 & dest, Matrix3 const &, Matrix3 const&);
  static void add (Matrix3 & dest, Matrix3 const &, Matrix3 const&);
  static void transpose_to (Matrix3 &dest, Matrix3 const&);
  static void invert_to (Matrix3 &dest, Matrix3 const&);
  static void invert_to_with_det (Matrix3 &dest, Matrix3 const&,Real);
  static void cofactor_matrix (Matrix3 &dest, Matrix3 const&);
  static void multiply_mv (Vector3 &dest, Matrix3 const &m, Vector3 const &rhs);
  static void multiply_vm (Vector3 &dest, Vector3 const &lhs, Matrix3 const &m);
  static void axpy (Matrix3 &dest, Real a,  Matrix3 const&, Matrix3 const&);
  /**
     Set rotation matrices; angle is in radians.
   */

  static void symmetrify (Matrix3 &dest, Matrix3 const &src);
  static void double_symmetrify (Matrix3 &dest, Matrix3 const &src);  
  
  void set_rotate (int axis,  Real angle);
  static void to_opengl_matrix (Real *rp, Matrix3 const&m);
  static void from_opengl_matrix (Matrix3 *m, Real const *rp );
  static void from_opengl_matrixf (Matrix3 *m, float const *rp );  
  void operator *= (Real r);
  static Matrix3 rotate_about (Vector3, Real);
  bool kosher () const;
};

inline Matrix3
operator * (Matrix3 const &m1, Matrix3 const &m3)
{
  Matrix3 ret;
  Matrix3::multiply_mm (ret, m1, m3);
  return ret;
}


inline Vector3
operator * (Matrix3 const &m1, Vector3 const &v)
{
  Vector3 w;
  Matrix3::multiply_mv (w, m1, v);
  return w;
}

inline Matrix3
operator +(Matrix3 const &m1, Matrix3 const &m3) 
{
  Matrix3 r;  
  Matrix3::add (r, m1,m3);
  return r;
}


inline Matrix3
operator -(Matrix3 const &m1, Matrix3 const &m3) 
{
  Matrix3 r;  
  Matrix3::subtract (r, m1,m3);
  return r;
}

inline Matrix3
operator * (Real r, Matrix3 const & m)
{
  Matrix3 ret;
  Matrix3::scale (ret, r, m);
  return ret;
}

#endif /* MATRIX_HH */


/*   
  big-vector.hh -- declare Big_vector.

  (c) 2001 Han-Wen Nienhuys <hanwen@cs.uu.nl>

 */

#ifndef BIG_VECTOR_HH
#define BIG_VECTOR_HH

#include "array.hh"
#include "defo-proto.hh"

/*
  ugh.  GDB 5.0 suckiness

  TODO:

  This is all very nice, but this should really be done through (for
  instance) GSL (GNU scientific library). Lowlevel routines in BLAS /
  LAPACK can be 2 to 3 times as fast as the normal ones.
  
 */
extern "C" {
  Real big_vector_1_norm (Real const * v, int n);
  Real big_vector_distance (Real const * v1, Real const *v2, int);
  Real big_vector_max_distance (Real const * v1, Real const *v2, int);  
  Real big_vector_ip (Real const * p1, Real const * p2, int n);
  Real big_vector_length (Real const *s, int n);
  Real big_vector_length_squared (Real const *s, int n);  
  void big_vector_add (Real  *dest, Real const * p1, Real const * p2, int n);
  void big_vector_copy (Real *dest, Real const * src, int n);
  void big_vector_fill (Real *dest, Real val, int n);
  void big_vector_negate (Real *dest, Real const*src, int);
  void big_vector_nullify (Real *ar, int);
  bool big_vector_is_null (Real const *ar, int);  
  void big_vector_partial_copy (Real *d, bool const*, Real const *src,int n); 
  void big_vector_partial_nullify (Real * dest, Real const *src, bool const * onoff,  int n);
  void big_vector_partial_inv_nullify (Real * dest, Real const *src, bool const * onoff,  int n);
  void big_vector_pointwise_inverse (Real *dest, Real const *src, int);
  void big_vector_pointwise_multiply (Real *dest, Real const * a, Real const *b, int );
  void big_vector_print (Real const *ar, int);
  void big_vector_axpy (Real  * dest, Real a, Real const *x, Real const *y, int n);
  void big_vector_scale (Real * d, Real fact, Real *s, int n); 
  void big_vector_subtract (Real *dest, Real const *a, Real const *b, int n);
  bool big_vector_sane (Real const * , int n);
  
}

void completize_big_vector (Array<Real> *a, int n);

#endif /* BIG_VECTOR_HH */


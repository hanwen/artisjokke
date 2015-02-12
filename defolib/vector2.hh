/*   
  vector.hh -- declare Vector2, Matrix2, Real
  
  (c) 1999 Han-Wen Nienhuys <hanwen@cs.uu.nl>
  
 */

#ifndef VECTOR2_HH
#define VECTOR2_HH

#include "defo-proto.hh"

enum Axes {
  X_AXIS = 0,
  Y_AXIS,
  Z_AXIS,
} ;


struct Vector2
{
  Real elts_[2];

  
  Real operator ()(int i) const { return elts_[i]; } 
  Real &operator ()(int i) { return elts_[i]; }  
  Vector2 ();
  Vector2 (Vector3);
  Vector2 (Real x,Real y);
  Vector2 (Real *);
  
  void setnull();
  void fill (Real);
  void set (Real x, Real y);
  operator Real const * () { return elts_; } 
  static void scalar_multiply (Vector2 &dest, Real s, Vector2 const& v);
  static void add (Vector2 &dest, Vector2 const &w, Vector2 const &v);
  static void sub (Vector2 &dest, Vector2 const &w, Vector2 const &v);
  char * str()const;
  static Real ip (Vector2 const &, Vector2 const&);
  static Real static_length_squared (Vector2 const&);
  static void cross (Vector2 &, Vector2 const&,Vector2 const&);
  void normalize ();
  Vector2 normalized () const; 
  Vector2 operator - ()
  {
    Vector2 v(*this);
    v*=-1;
    return v;
  }
  
  Vector2 &operator *= (Real r)
    {
      scalar_multiply (*this, r, *this);
      return *this;
    }
  Vector2 &operator /= (Real r)
    {
      return (*this) *= 1/r;
    }
  Vector2 &operator += (Vector2 const &r)
    {
      add (*this, r, *this);
      return *this;
    }
  Vector2 &operator -= (Vector2 const &r)
    {
      sub (*this, *this, r);
      return *this;
    }

  Real length_squared () const { return static_length_squared (*this); }
  Real length () const ;
  static Real cross (Vector2 const& v1, Vector2 const &v2);
  static Vector2 combination (Real lp, Vector2 p, Real lq, Vector2 q);  
  void print ()const;
};

inline Vector2
rotate90 (Vector2 z) {
  return Vector2 (-z(Y_AXIS), z(X_AXIS));
}

Real euclidian_distance (Vector2 const &,Vector2 const &);
Real vector_angle (Vector2 n1 ,Vector2 n2);
Vector2 dir_vector (Real angle);

Vector2 operator + (Vector2 v1, Vector2 const &v2);
Vector2 operator - (Vector2 v1, Vector2 const &v2);

Vector2 operator * (Real r, Vector2 v1);
Vector2 operator * (Vector2 v1,Real r);

Vector2 operator / (Vector2 v1,Real r);


Real operator *(Vector2 const &v, Vector2 const& w);

inline Real sqr (Real r)
{
  return r*r;
}

/*fsking solaris.
 */
inline double
my_fabs (double d)
{
  if (d<0)
    return -d;
  return d;
}

inline int
sign (double d)
{
  if (d> 0)
    return 1;
  else if (d < 0)
    return -1 ;
  return 0;  
}


inline
void
Vector2::set (Real x, Real y)
{
  elts_[0]  = x;
  elts_[1]  = y;
}



inline
Vector2::Vector2 (Real x, Real y)
{
  elts_[0]  = x;
  elts_[1]  = y;
}


#endif /* VECTOR2_HH */


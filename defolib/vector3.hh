/*   
  vector.hh -- declare Vector3, Matrix3, Real
  
  (c) 1999 Han-Wen Nienhuys <hanwen@cs.uu.nl>
  
 */

#ifndef VECTOR3_HH
#define VECTOR3_HH

#include "defo-proto.hh"

struct Vector3
{
  Real elts_[3];

  Real operator ()(int i) const { return elts_[i]; } 
  Real &operator ()(int i) { return elts_[i]; }  
  Vector3 ();
  Vector3 (Real x,Real y, Real z){ elts_[0] =x;elts_[1] =y;elts_[2] =z;}
  Vector3 (Real *xs);
  Vector3 (Vector2 x){ elts_[0] =x(0);elts_[1] =x(1);elts_[2] =0;}
  void fill (Real);
  void setnull ();
  void set (Real x ,Real y,Real z) { elts_[0] =x;elts_[1] =y;elts_[2] =z;}
  operator Real const * () { return elts_; } 

  
  static void scalar_multiply (Vector3 &dest, Real s, Vector3 const& v);
  static void add (Vector3 &dest, Vector3 const &w, Vector3 const &v);
  static void sub (Vector3 &dest, Vector3 const &w, Vector3 const &v);
  char * str()const;
  static Real ip (Vector3 const &, Vector3 const&);
  static Real static_length_squared (Vector3 const&);
  static void cross (Vector3 &, Vector3 const&,Vector3 const&);
  void normalize ();
  Vector3 normalized () const; 
  Vector3 operator - ()
  {
    Vector3 v(*this);
    v*=-1;
    return v;
  }
  
  Vector3 &operator *= (Real r)
    {
      scalar_multiply (*this, r, *this);
      return *this;
    }
  Vector3 &operator /= (Real r)
    {
      return (*this) *= 1/r;
    }
  Vector3 &operator += (Vector3 const &r)
    {
      add (*this, r, *this);
      return *this;
    }
  Vector3 &operator -= (Vector3 const &r)
    {
      sub (*this, *this, r);
      return *this;
    }

  Real length_squared () const { return static_length_squared (*this); }
  Real length () const ;
  static Vector3 cross (Vector3 v1, Vector3 v3);
  static Vector3 combination (Real lp, Vector3 p, Real lq, Vector3 q);  
  void print ()const;
};

Real euclidian_distance (Vector3 const &,Vector3 const &);
Real vector_angle (Vector3 n1 ,Vector3 n3);

Vector3 random_unit_vector ();

Vector3 operator + (Vector3 v1, Vector3 const &v3);
Vector3 operator - (Vector3 v1, Vector3 const &v3);

Vector2 operator / (Vector2 v1,Real r);


Vector3 operator * (Real r, Vector3 v1);
Vector3 operator * ( Vector3 v1,Real r);

Real operator *(Vector3 const &v, Vector3 const& w);



#endif /* VECTOR_HH */


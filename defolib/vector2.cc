

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include "matrix.hh"
#include "vector.hh"



Vector2::Vector2 (Vector3 v)
{
  elts_[0]  = v(0);
  elts_[1]  = v(1);
}

Real
Vector2::cross (Vector2 const &v1, Vector2 const &v2)
{
  return (v1(0) * v2(1) - v2(0) * v1(1));
}



// damn GDB
Real
inner_product (Vector2 a, Vector2 b)
{
  return a * b;
}

Real
Vector2::length () const
{
  return sqrt (length_squared ());
}


void
Vector2::normalize ()
{
  *this /= length();
}

Real
euclidian_distance (Vector2 const & a, Vector2 const &b)
{
  Vector2 c (a);
  c -= b;
  return c.length ();
}

Vector2
Vector2::combination (Real lp, Vector2 p, Real lq, Vector2 q)
{
  return lp* p  + lq * q;
}

void
Vector2::print ()const
{
  Vector2 v = *this;
  fprintf (stdout, "(%lf,%lf)\n", v(0), v(1));
}


char * 
Vector2::str ()const
{
  char s [1024];
  Vector2 v = *this;
  sprintf (s, "(%lf, %lf)", v(0), v(1));

  return strdup (s);
}



/*
  Compute the angle between T->triangle (I) and T->triangle (J)
 */
Real
vector_angle (Vector2 n1 ,Vector2 n2)
{
  Real cosine = (n1*n2) / (n1.length  ()  * n2.length ());

  /*
    Safety check: rounding errors  may put COSINE slightly above 1.0, causing
    NaN results.
   */
  if (fabs (cosine) > 1)
    {
      return cosine > 0 ? M_PI : 0.0;
    }
  
  Real a = M_PI - acos( cosine);


  assert (!isnan (a));

  return a;
}


Vector2
Vector2::normalized () const
{
  Vector2 rv;
  rv = *this;
  rv.normalize();
  return rv;
}

Vector2
dir_vector (Real angle)
{
  return Vector2 (cos (angle),
		  sin (angle));
}

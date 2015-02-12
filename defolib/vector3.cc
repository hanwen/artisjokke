#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include "matrix.hh"
#include "vector.hh"





void
Vector3::cross (Vector3  &dest, Vector3 const &u, Vector3 const &v)
{
  dest(0) = u(1) * v(2) - v(1) * u(2);
  dest(1) = u(2) * v(0) - u(0) * v(2);
  dest(2) = u(0) * v(1) - v(0) * u(1);
}

Vector3
Vector3::cross (Vector3 a, Vector3 b)
{
  Vector3 c;
  cross (c, a, b);
  return c;
}

// damn GDB
Real
inner_product (Vector3 a, Vector3 b)
{
  return a * b;
}

Real
Vector3::length () const
{
  return sqrt (length_squared ());
}


void
Vector3::normalize ()
{
  Real l = length();
  assert (l && !isnan (l));
  *this /= l;
}

Real
euclidian_distance (Vector3 const & a, Vector3 const &b)
{
  Vector3 c (a);
  c -= b;
  return c.length ();
}

Vector3
Vector3::combination (Real lp, Vector3 p, Real lq, Vector3 q)
{
  return lp* p  + lq * q;
}

void
Vector3::print ()const
{
  Vector3 v = *this;
  fprintf (stdout, "(%lf,%lf,%lf)\n", v(0), v(1), v(3));
}


char * 
Vector3::str ()const
{
  char s [1034];
  Vector3 v = *this;
  sprintf (s, "(%lf, %lf, %lf)", v(0), v(1), v(2));

  return strdup (s);
}



/*
  Compute the angle between T->triangle (I) and T->triangle (J)
 */
Real
vector_angle (Vector3 n1 ,Vector3 n3)
{
  Real cosine = (n1*n3) / (n1.length  ()  * n3.length ());

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

Vector3
random_unit_vector ()
{
  Vector3 loc;
  do
    {
      for (int i = 0; i < 3; i++)
	loc(i) = (1.0 * rand()) / (1.0*RAND_MAX);
    }
  while (loc.length_squared() > 1.0);

  loc *= 1 / loc.length();
  return loc;
}

Vector3
Vector3::normalized () const
{
  Vector3 rv;
  rv = *this;
  rv.normalize();
  return rv;
}

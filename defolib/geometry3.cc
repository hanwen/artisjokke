#include <algorithm>
#include <math.h>
#include <assert.h>

#include "vector.hh"
#include "matrix.hh"
#include "geometry.hh"

using std::min;

/*
  Compute the normal of a 2-simplex, using its parity.
 */
Vector3
simplex_normal3 (Simplex const &p, Node_vector_func3 func, Deformation_state *def)
{
  assert (p.dimension()== TRIANGLE);
  
  Vector3 vs[3];

  for (int i = 0; i < 3 ; i++)
    {
      vs[i] = (*func) (p.node(i), def);
    }

  Vector3 n = Vector3::cross (vs[1] -vs[0], vs[2] - vs[0]);
  n.normalize();

  if (!p.parity())
    n = -n;

#if 0
  static enum {
    UNSET, ONEDIR, OTHERDIR
  } ccw;

  if (ccw == UNSET)
    {
      if (strcmp (get_string_setting ("normal-orientation"), "ccw"))
        ccw = OTHERDIR;
      else
        ccw = ONEDIR;
    }

  if  (ccw == ONEDIR)
    n *= -1;
#endif

  return n;
}

Real
tetrahedral_volume (Vector3 *vs)
{
  return robust_orient3d (vs[0].elts_,
                          vs[1].elts_,
                          vs[2].elts_,
                          vs[3].elts_) / 6.0;

}

Real
simplex_volume3 (Simplex const &ts, Node_vector_func3 func, Deformation_state *def)
{
  assert (ts.dimension() == TETRAHEDRON);
  
  Vector3 vs[4];
  for (int i = 4; i-- ;)
    vs[i] = (*func)(ts.node(i),def);

  Real v = tetrahedral_volume (vs);
  if (ts.parity())
    v *= -1;
  return v;
}





Vector3
simplex_centroid3 (Simplex const&s,
		  Node_vector_func3 func,
		  Deformation_state*def)
{
  Vector3 v ;

  for (int i = s.count() ; i--;)
    v += (*func) (s.node(i),def);

  v *= 1/double (s.count());
  return v;
}


/*
  return barycentric coords for the tet with POINTS and origin as
  verrtices.
 */
Vector3
barycentric_origin_tetrahedron_coordinates (Vector3 *points,
				     Vector3 point)
{
  Matrix3 m;
  m.set_column_vectors (points[0], points[1], points[2]);

  Matrix3 inv;
  Matrix3::invert_to(inv,m);

  return inv * point;
};

void
barycentric_coordinates3 (Vector3 const *points,
			  Real *coef,
			  Vector3 point)
{
  Vector3 vs[3];
  for (int i = 0; i < 3; i++)
    vs[i] = points[i] - points[3];

  Vector3 bar = barycentric_origin_tetrahedron_coordinates (vs,
							    point - points[3]);

  coef[3] = 1.0;
  for (int i = 0; i < 3;i++)
    {
      coef[i] = bar(i);
      coef[3] -= bar(i); 
    }
};

/*
  This is the volume of F extended with V. This is called virtual
  since V is not a vertex of the complex, and the simplex does not
  exist, strictly spoken.
 */
Real
virtual_deformed_simplex_volume3 (Simplex const & fplex, 
				  Vector3 const &v,Node_vector_func3 func,
				  Deformation_state * def)
{
  Vector3 vs[MAX_DIMENSION];
  for (int i = 0;i <fplex.count() ; i++)
    vs[i] = (*func) (fplex.node(i), def);
 
  int sign = fplex.parity() ? 1 : -1; 
  return sign * robust_orient3d (vs[0].elts_,
				 vs[1].elts_,
				 vs[2].elts_,
				 v.elts_);
}

/****************************************************************/

bool
point_in_triangle (Vector3 point,
		   Vector3 * vs)
{
  Real coef[2];
  
  Vector3 v1  = vs[0] -vs[2];
  Vector3 v2  = vs[1] -vs[2];
  Vector3 p = point -vs[2];    

  Matrix2 m;
  m(0,0) = v1*  v1;
  m(1,1) = v2*  v2;
  m(1,0) = m(0,1) = v1*  v2;    
  Matrix2 inv;
  Matrix2::invert_to (inv, m);
  Vector2 rhs (v1 * p, v2 * p);

  Vector2 lm = inv *rhs ;
  coef[0] = lm(0);
  coef[1] = lm(1);
  

  return (coef[0] >= 0 && coef[1] >= 0 && coef[1] + coef[0] <= 1.0);  
}

bool
intersect_triangle_plane (Real * param, Vector3 *point,
			  Vector3 *vs,
			  Vector3 p, Vector3 q)
{
  Vector3 v0 = vs[0] - vs[2];
  Vector3 v1 = vs[1] - vs[2];
  
  Vector3 normal = Vector3::cross (v0, v1);

  Real off = vs[2] * normal;

  Vector3 dir = q -p ;
  Vector3 org = p;
  Real dot = normal * dir;
  if (fabs (dot) < 1e-8)	// ugh
    return false;
      
  Real dist = ( off -  org * normal)  / dot;

  *param = dist;
  *point = dir * dist + org;

  return true;
}

bool
intersect_line_simplex3 (Real *param,
			 Simplex const &pl,
			 Node_vector_func3 func,
			 Deformation_state *def,
			 Vector3 p, Vector3 q) 
{
  Vector3 t;
  Vector3 vs[3];
  for (int i = 0; i < 3; i++)
    vs[i] = (*func) (pl.node (i), def);
  
  bool s = intersect_triangle_plane (param, &t, vs, p, q);
  if (!s)
    return false;
    
  return point_in_triangle(t, vs);
}


/****************************************************************/

Vector3
project_point_on_line3 (Vector3 p, Vector3 l1, Vector3 l2)
{
  Vector3 dir = l2 - l1;
  dir /= dir.length();
  p = p -  l1;
  return (dir * p)* dir + l1;
}

Real
point_to_line_distance3 (Vector3 p, Vector3 l1, Vector3 l2)
{
  Vector3 dir = l2 - l1;
  dir /= dir.length();
  p = p - l1;
  p = p - (dir * p)* dir;
  return p.length();
}

/*
  distance is a 1-param function d(l), where l in [0,1]. Either the
  minimum is in the interior (achieved by projecting P on l1 l2), or
  it is achieved at the end points of the segment.
 */
Real
point_to_line_segment_distance3 (Vector3 p, Vector3 l1, Vector3 l2)
{
  Vector3 dir = l2 - l1;

  p = p -  l1;

  Real dot = (dir * p) / (dir * dir);
  if (dot >= 0 && dot < 1.0)
    {
      p = p -  dot * dir;
      return p.length();
    }
  else
    {
      return min(euclidian_distance (l1, p),
		 euclidian_distance (l2, p));
    }
}

/****************************************************************/


/*
  return distance to plane, if the projected point falls in the
  triangle. Otherwise return infty.
 */
Real
point_to_triangle_distance3 (Vector3 p, Vector3 l1, Vector3 l2, Vector3 l3)
{
  Vector3 vs[3] = {l1, l2 , l3} ;

  Vector3 n = Vector3::cross (l1 - l3, l2-l3);
  
  
  Vector3 in_tri = p - (((p - l3) * n ) / (n*n)) *n ;

  if (point_in_triangle (in_tri, vs))
    {
      return euclidian_distance (in_tri, p);
    }
  else
    return 1e18;
}


Real point_to_simplex_distance3 (Simplex const &plex, Node_vector_func3 func,
				 Deformation_state *def, Vector3 point)
{
  Vector3 vs[MAX_DIMENSION +1];

  for (int i = 0; i < plex.count () ; i++)
    {
      vs[i] = (*func) (plex.node (i), def);
    }

  Real dist = 1e18;
  
  for (int i = 0 ; i < plex.count();  i++)
    for (int j = i+1 ; j < plex.count();  j++)
      {
	dist = min (dist, 
		    point_to_line_segment_distance3 (point, vs[i], vs[j]));
      }

  for (int i = 0 ; i < plex.count();  i ++ )
    for (int j = i+1 ; j < plex.count();  j ++ )
      for (int k = j+1 ; k < plex.count();  k ++ )          
	{
	  dist = min (dist,
		      point_to_triangle_distance3 (point, vs[i], vs[j], vs[k]));
	}
  return dist;
}


/****************************************************************/


/*
  Find L, M, b such that

   L P + (1-L) Q + b n = M A + (1-M) B

   n = (P- Q)x (A-B)

   RETURN:

   L P + (1-L) Q
   
 */
bool 
cross_lines3 (Vector3 *point, Real *pparam, Real *aparam,
	     Vector3 p, Vector3 q, Vector3 a, Vector3 b)
{
  Matrix3 m;

  Vector3 n =Vector3 ::cross (p-q, a-b);
  m.set_column_vectors (p-q, b-a, n);

  /*
    ugh.
   */
  if (m.determinant ()== 0.0)
    {

      return false;
    }
  
  Vector3 rhs = b-q;

  Matrix3 i;
  Matrix3::invert_to (i,m);
  
  
  
  Vector3 sol =i * rhs;

  *pparam = sol(0);
  *aparam = sol(1);

  *point = Vector3::combination (*pparam, p, 1-*pparam, q);
  return true;
}


/*
  distance for two line segments is  a 2 -param function,

  d(lambda, mu).  The minimum may be attained for

  I = [0,1], dI = {0,1}, Iint = (0,1)
  
  lambda in  dI, or lambda in Iint and
  mu in dI, or mu in Iint.
  
  The case that both are in Iint is done by cross_lines().  The other
  cases are covered by the 4 point_to_line_segment_distance() calls.

*/
Real
line_segment_to_line_segment_distance (Vector3 p, Vector3 q, Vector3 l1, Vector3 l2)
{
  Real a, b;
  
  Vector3 cp1;
  bool succ = cross_lines3 (&cp1, &a, &b, p, q, l1, l2);
  
  
  

  if (succ && 0.0 <= a && a <= 1.0 && a >= 0.0
      && 0.0 <= b && b <= 1.0)
    {
      Vector3 cp2 = Vector3::combination (b, l1, 1-b, l2);
      return euclidian_distance(cp1,cp2);
    }
  else
    {
      /*
	We do some work doubly, since the endpoints are compared in
	every call of point_to_line_segment_distance
       */
      return
	min(point_to_line_segment_distance3 (p, l1, l2),
	    min (point_to_line_segment_distance3 (q, l1, l2),
		 min (point_to_line_segment_distance3 (l1, p, q),
		      point_to_line_segment_distance3 (l2, p, q)))); 
    }
}



Real
triangle_to_line_segment_distance3 (Simplex const & plex,Node_vector_func3 func,
				  Deformation_state *def, Vector3 p1, Vector3 p2)
{
  Vector3 vs[3];
  Real dist = 1e18;

  assert (plex.count () == 3);
  for (int i = 0; i < 3; i++)
    {
      vs[i] = (*func) (plex.node (i), def);
    }

  
  for (int i = 0; i <  3; i++)
    {
      dist = min (dist, 
		  line_segment_to_line_segment_distance (vs[(i+1)%3], vs[(i+2)%3],
							 p1, p2));
    }

  dist = min (dist, point_to_triangle_distance3 (p1, vs[0], vs[1], vs[2]));
  dist = min (dist, point_to_triangle_distance3 (p2, vs[0], vs[1], vs[2]));


  Real param ;
  Vector3 inter ;
  bool b = intersect_triangle_plane (&param, &inter,
				     vs, p1, p2);

  if (b && param >= 0.0 && param <= 1.0
      &&  point_in_triangle(inter, vs))
    return 0.0;

  return dist;
}


Real
tetrahedron_to_line_segment_distance3 (Simplex const & plex, Node_vector_func3 func,
				      Deformation_state *def, Vector3 p1, Vector3 p2)
{
  Real dist = 1e18;
  if (tetrahedron_contains_point3 (plex, p1, func, def , 0.0)
      || tetrahedron_contains_point3 (plex, p2, func, def , 0.0))
    return 0.0;
  
  for (int i = 0; i < 4; i++)
    {
      dist = min (dist, triangle_to_line_segment_distance3 (plex.get_subset (i),
							    func, def,
							    p1, p2));
    }

  return dist;
}

/****************************************************************/

bool
tetrahedron_contains_point3 (Simplex const & plex, Vector3 loc,
			    Node_vector_func3 func, Deformation_state *def, Real eps)
{
  Vector3 vs[4];
  for (int i = 0; i < plex.count (); i++)
    vs[i]= (*func) (plex.node(i), def);

  Real b[4];
  barycentric_coordinates3 (vs, b, loc);

  for (int j = 4; j--;)
    if (b[j] < - eps)
      return false;
  return true;
}


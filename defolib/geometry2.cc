
#include <assert.h>

#include "simplex.hh"
#include "geometry2.hh"
#include "matrix.hh"
#include "deformation-state.hh"

Real
simplex_area (Simplex const &s,
	      Node_vector_func2 func, 
	      Deformation_state *def)
{
  Vector2 v1, v2;

  v1 = (*func) (s.node (0), def) - (*func) (s.node (2), def);
  v2 = (*func) (s.node (1), def) - (*func) (s.node (2), def);
  
  Real a = Vector2::cross (v1, v2);

  if (s.parity())
    a *= -1;

  a /= 2;
  
  return a;
}

Real
simplex_volume2 (Simplex const &s, Node_vector_func2 f, Deformation_state*d )
{
  if (s.count () == 3)
    {
      return simplex_area (s, f,d);
    }
  else
    return 0.0;
}


/*
  Return in COEF X expressed in barycentric coords. X should be on the
  triangle.
*/
void
barycentric_coordinates2 (Vector2 const * vertices,
			 Real *coef, Vector2 ax) 
{
  Vector2 v1 = vertices[1] - vertices[0];
  Vector2 v2 = vertices[2] - vertices[0];
  Vector2 x = ax - vertices[0];

  Matrix2 mat = gram_matrix2 (v1,v2);
  Vector2 rhs (x  * v1, x * v2);
  Matrix2 inv;
  Matrix2::invert_to (inv, mat);

  Vector2 coefs;
  coefs =  inv *rhs  ;
  coef[1] = coefs(0); 
  coef[2] = coefs(1); 
  coef[0] = 1- coefs(1)- coefs(0);

  // #define PARANOIA
  
#ifdef PARANOIA
  Vector2 y;
  for (int i = 3; i--;)
    y+= coef[i] * vertices[i];

  Real defect = (ax-y).length  () / v1.length ();
  if (defect > 1e-6)
    {
      printf ("defect %lf\n", defect);
      assert (false);
    }
#endif
}



Vector2
simplex_centroid2 (Simplex const&s,
		  Node_vector_func2 func,
		  Deformation_state*def)
{
  Vector2 v ;

  for (int i = s.count() ; i--;)
    v += (*func) (s.node(i),def);

  v *= 1/double (s.count());
  return v;
}

Real
edge_length2 (Simplex const&e, Node_vector_func2 func,
	     Deformation_state const*def)
{
  return euclidian_distance ((*func) (e.node(0), def), (*func)(e.node(1), def));
}


Real
edge_length3 (Simplex const & pl,
	      Node_vector_func3 func, Deformation_state const *def)
{
  return  euclidian_distance ((*func)(pl.node(0), def),
			      (*func)(pl.node(1), def));
	
}




bool
intersect_lines (Vector2 *point, Real *l , Real *m, 
		 Vector2 p1, Vector2 p2, Vector2 q1, Vector2 q2)
{/*
	"""Return (L,M,X), st.
	
	X = L p1 + (1-L) p2 = M q1 + (1-M) q2

	"""
 */
  Matrix2  mat;

  mat.set_column_vectors (p1-p2, q2-q1);
  if (mat.determinant() == 0.0)
    return  false;
	

  Vector2 rhs = q2-p2;

  Matrix2 inv ;
  Matrix2::invert_to (inv, mat);
  Vector2 sol = inv * rhs;

  if (point)
    *point = (sol(0) * p1 + (1-sol(0))* p2);
  if (l)
    *l = sol(0);
  if (m)
    *m = sol(1);
  
  return true;
}

bool intersect_segments2 (Vector2 * point,
			 Real *l, Real *m,
			 Vector2 p1, Vector2 p2,
			 Vector2 q1, Vector2 q2)
{
  bool b = intersect_lines (point, l,m,
			 p1,p2,q1,q2);

  if (b && -eps <= *l
      && *l <= 1 +eps
      &&  -eps <= *m &&  *m < 1+ eps)
    return true;
  else
    return false;
}

Vector2
project_point_on_line2 (Vector2 p, Vector2 l1, Vector2 l2)
{
  Vector2 dir = l2 - l1;
  dir /= dir.length();
  p = p -  l1;
  return (dir * p)* dir + l1;
}

Vector2
simplex_normal2 (Simplex const &s, Node_vector_func2 func , Deformation_state* def)
{
  Vector2 v1 = (*func) (s.node(0) , def);
  Vector2 v2 = (*func) (s.node(1) , def);
  Vector2 n = rotate90 (v1- v2);
  if (s.parity())
    n *= -1;

  
  return n.normalized(); 
}

/*
  ugh ugh.

  this code can be condensed by a factor 3, I guess.
 */
Vector2
transform_back2 (Simplex const &s , Vector2 x, Deformation_state * def)
{
  Real const  * displacement_vec = def->displacements_.access_array ();
  int spatial_dim = 2;
  Matrix2 ux_inv ;

  Matrix2 displ_mat;

  Real const *d2 = displacement_vec + s.node (spatial_dim)->number () *  spatial_dim;
  for (int i=spatial_dim; i--;)
    {
      Real const *di = displacement_vec + s.node (i)->number () * spatial_dim;
      for (int j=spatial_dim; j--;)
        displ_mat(j,i) =  di[j] -  d2[j];
    }
  Matrix2 locmat, invlocmat;
  
  for (int i=0; i < spatial_dim; i++)
    for (int j=0; j < spatial_dim; j++)
      {
        locmat(j,i) = reference_location2 (s.node (i), def) (j)
	  -  reference_location2 (s.node (spatial_dim), def)(j);
      }
  
  Matrix2::invert_to (invlocmat, locmat);
  Matrix2::multiply_mm (ux_inv, displ_mat, invlocmat);

  // ugh. Use disp. vector here. 
  Vector2 x2 = deformed_location2 (s.node (spatial_dim), def);
  Vector2 X2 = reference_location2 (s.node (spatial_dim), def);

  x += ux_inv * X2 - (x2 - X2); 

  for (int i=0; i < spatial_dim; i++)
    ux_inv(i,i)+= 1.0;
  
  Matrix2 xu_inv;
  Matrix2::invert_to (xu_inv, ux_inv);
  
  return xu_inv * x;
}




/*
  This is the volume of F extended with V. This is called virtual
  since V is not a vertex of the complex, and the simplex does not
  exist, strictly spoken.
 */
Real
virtual_deformed_simplex_volume2 (Simplex const & fplex, 
				  Vector2 const &v,Node_vector_func2 func,
				  Deformation_state * def)
{
  Vector3 vs[fplex.count ()];
  for (int i = 0;i <fplex.count() ; i++)
    vs[i] = (*func) (fplex.node(i), def);
 
  int sign = fplex.parity() ? -1 : 1; 
  return sign * robust_orient2d (vs[0].elts_,
				 vs[1].elts_,
				 v.elts_);
}



#ifndef GEOMETRY2_HH
#define GEOMETRY2_HH

#include "defo-proto.hh"
#include "simplex.hh"


Real simplex_area (Simplex const&, Node_vector_func2, Deformation_state*);
Real simplex_volume2 (Simplex const&,Node_vector_func2, Deformation_state*);
void barycentric_coordinates2 (Vector2 const * vertices,
			       Real *coef, Vector2 x);


Vector2 simplex_centroid2 (Simplex const& s,Node_vector_func2, Deformation_state*);
Real edge_length2 (Simplex const&e, Node_vector_func2 func, Deformation_state const *def);

Vector2 simplex_normal2 (Simplex const & s, Node_vector_func2 , Deformation_state* );




bool intersect_lines2 (Vector2 *point,
		       Real *l, Real *m, 
		       Vector2 p1, Vector2 p2,
		       Vector2 q1, Vector2 q2);

bool intersect_segments2 (Vector2 * point,
			 Real *l, Real *m,
			 Vector2 p1, Vector2 p2,
			 Vector2 q1, Vector2 q2);

Vector2 project_point_on_line2 (Vector2 p, Vector2 l1, Vector2 l2);

Real robust_orient2d (Real const*,Real const*,Real const*);

Real
virtual_deformed_simplex_volume2 (Simplex const & fplex, 
				  Vector2 const &v,Node_vector_func2 func,
				  Deformation_state * def);

/*
  mesh geo.
 */
Element* locate_element_containing_point (Mesh_connectivity* mesh, Element *entry, Vector2 v, 
					  Node_vector_func2 func, Deformation_state * def);



#endif

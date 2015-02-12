#ifndef GEOMETRY3_HH
#define GEOMETRY3_HH

Vector3 simplex_normal3 (Simplex const & s, Node_vector_func3, Deformation_state* );

Real simplex_volume3 (Simplex const &ts, Node_vector_func3 func, Deformation_state *def);
Real tetrahedral_volume (Vector3 *vs);
Real robust_orient3d(Real const *pa, Real const *pb, Real const *pc, Real const *pd);
Real edge_length3 (Simplex const & pl,
		   Node_vector_func3 func, Deformation_state const *def);
Real virtual_deformed_simplex_volume3 (Simplex const & fplex, 
				       Vector3 const &v,Node_vector_func3 func,
				       Deformation_state * def);

void barycentric_coordinates3 (Vector3 const *points, Real *coef, Vector3 point);


Vector3
simplex_centroid3 (Simplex const&s,
		  Node_vector_func3 func,
		   Deformation_state*def);
bool
intersect_line_simplex3 (Real *param,
		    Simplex const &pl,
		    Node_vector_func3 func,
		    Deformation_state *def,
		    Vector3 p, Vector3 q) ;


Real point_to_line_distance3 (Vector3 p, Vector3 l1, Vector3 l2);
Real point_to_line_segment_distance3 (Vector3 p, Vector3 l1, Vector3 l2);
Vector3 project_point_on_line3 (Vector3 p, Vector3 l1, Vector3 l2);

Real point_to_simplex_distance3 (Simplex const & , Node_vector_func3 func,
				 Deformation_state* ,
				 Vector3 point);

Real line_segment_to_line_segment_distance3 (Vector3 p, Vector3 q, Vector3 l1, Vector3 l2);
Real tetrahedron_to_line_segment_distance3 (Simplex const & plex, Node_vector_func3 func,
					    Deformation_state *def, Vector3 p1, Vector3 p2);
bool tetrahedron_contains_point3 (Simplex const & plex, Vector3 loc,
				 Node_vector_func3 func, Deformation_state *def, Real eps);


#endif

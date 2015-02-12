#ifndef GEOMETRY_HH
#define GEOMETRY_HH


#include "geometry2.hh"
#include "geometry3.hh"



Element*
locate_element_containing_point2 (Mesh_connectivity *mesh,
				  Element *entry,
				  Vector2 v, 
				  Node_vector_func2 func,
				  Deformation_state * def);


Element*
locate_element_containing_point3 (Mesh_connectivity *mesh,
				  Element *entry,
				  Vector3 v, 
				  Node_vector_func3 func,
				  Deformation_state * def);
void minmax_edge_length2 (Real *min, Real *max, Simplex const & pl,
			 Node_vector_func2 func, Deformation_state const *def
			 );
void minmax_edge_length3 (Real *min, Real *max, Simplex const & pl,
			 Node_vector_func3 func, Deformation_state const *def
			 );


#endif

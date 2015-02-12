#ifndef NEEDLE_2D_MISC
#define NEEDLE_2D_MISC

#include "proto.hh"

void plot_reference_state_distance (Maubach_tree * tree,
				    Deformation_state *def,
				    Deformation_state *ref,
				    char const *fn);



void
mesh2d_print_boundary (Mesh_connectivity* top,
	    char const * fn,
	    Node_vector_func2 func,
	    Deformation_state * def,
	    Needle_inserter * ni,
		       Real scale);
void mesh2d_print (Mesh_connectivity* top, char const * fn, Node_vector_func2 func, Deformation_state * def, Needle_inserter*, Real scale );

void mesh2d_dot_print (Mesh_connectivity* top, char const * fn, Node_vector_func2 func, Deformation_state * def, Real scale );

void mesh2d_view (Mesh_connectivity*);

#endif

#ifdef OPENGL

#ifndef MESH_VISUALIZE3_HH
#define MESH_VISUALIZE3_HH

#include "artisjokke-drawer.hh" 
#include "defo-proto.hh"
#include "opengl.hh"


void
opengl_draw_face (Artisjokke_drawer*draw, Face*t, bool solid,
		  Node_vector_func3 func, 
		  Deformation_state*def);


#endif
#endif

#ifdef OPENGL

#include <math.h>

#include "opengl.hh"
#include "artisjokke-drawer2.hh"
#include "mesh-connectivity.hh"
#include "deformation-state.hh"
#include "misc2d.hh"
#include "maubach-tree.hh"

float needle_color [] = {1,0,0};

void
Artisjokke_drawer2::process_key (char c)
{
  switch (c) {
  case 'p':
    mesh2d_print (mesh_,
		"deformed-mesh.eps",
		&deformed_location2,
		deformation_,
		needle_inserter_,
		10);
    mesh2d_print (mesh_,
		"reference-mesh.eps",
		&reference_location2,
		deformation_,
		needle_inserter_,
		10);
    mesh2d_dot_print (mesh_,
		    "dot-mesh.eps",
		    &deformed_location2,
		    deformation_,
		    10);

    mesh2d_print_boundary (mesh_,
			   "reference-boundary.eps",
			   &reference_location2,
			   deformation_,
			   needle_inserter_,
			   10);
    mesh2d_print_boundary (mesh_,
			   "deformed-boundary.eps",
			   &deformed_location2,
			   deformation_,
			   needle_inserter_,
			   10);

    break ;

  }
}

Artisjokke_drawer2::Artisjokke_drawer2(Mesh_connectivity *top, Deformation_state*def)
: Artisjokke_drawer (top,def)
{
  this->needle_inserter_ =0; 

}

Maubach_tree* 
Artisjokke_drawer2::tree( ) const
{
  return dynamic_cast <Maubach_tree*> (mesh_);
}

void
Artisjokke_drawer2::set_needle()
{
  if (!needle_inserter_->live())
    return;
  
  Vector3 needle_handle_;
  Vector3 needle_tip_;

  needle_handle_.set (-1.0,0.073,0.0);
  Vector3 needle_dir = Vector3 (0.99,0,0);

  float om[4][4];
  build_rotmatrix (om, needle_drag_.quat_);
  Matrix3 m;
  Matrix3::from_opengl_matrixf (&m, &om[0][0]);
  
  const Real FACTOR_LIMIT = 100.0;
  Real fact = pow (2, needle_drag_.scale_(0));

  if (fact < 1/FACTOR_LIMIT)
    fact = 1/FACTOR_LIMIT;
  else if (fact > FACTOR_LIMIT)
    fact = FACTOR_LIMIT;

  needle_dir *= fact;
  //  needle_dir = m * needle_dir;

  needle_handle_ = needle_handle_ + (needle_drag_. translate_)*0.10;
  needle_tip_ = needle_handle_ + needle_dir;

  needle_inserter_->move_needle (needle_handle_, needle_tip_);
}


void
Artisjokke_drawer2::draw_needle()
{
  set_color (needle_color);
  glLineWidth(4.0);
  glBegin (GL_LINES);

  Vector3 t= needle_inserter_->tip ();
  Vector3 h= needle_inserter_->handle ();
  
  enter_vertex (h.elts_);
  enter_vertex (t.elts_);  
  glEnd ();
}

void
Artisjokke_drawer2::opengl_visualize_mesh_connectivity()
{
  opengl_visualize_mesh_connectivity2 (this);
}
#endif

#ifdef OPENGL

#include <algorithm>

#include <math.h>

#include "geometry.hh"
#include "opengl.hh"
#include "proto.hh"
#include "scenario.hh"
#include "artisjokke-drawer3.hh"
#include "mesh-connectivity.hh"
#include "needle-inserter3.hh"
#include "deformation-state.hh"
#include "needle-inserter-visualize.hh"
#include "setting.hh"
#include "auto-inserter3.hh"

using std::max;
using std::min;

float needle_color [] = {1,0,0, 0.8};
float ground_color [] = {1.0,1.0,1.0, 0.0};


Artisjokke_drawer3::Artisjokke_drawer3 (Mesh_connectivity *top,
					Deformation_state*def,
					Needle_inserter3 *n)
  : Artisjokke_drawer(top,def)
{
  axes_ = true;
  needle_inserter_ = n;

  lighting_ = true;
  needle_surface_ = true;
  fog_ = true;
  if (!get_bool_setting ("side-view"))
    trackball(view_drag_.quat_, 0.0, 0.0, -0.3, -0.1 );

  point_weight_ = 2.5;
}

void
Artisjokke_drawer3::opengl_visualize_mesh_connectivity ()
{
  opengl_visualize_mesh_connectivity3 (this);
}

void
Artisjokke_drawer3::process_key (char c )
{
  switch (c) {
  case '!':
    {
      
      Real  mine = 1e6;
      for (iterof (i, *mesh_->elements());
	   i !=  mesh_->elements ()->end (); i++)
	{
	  Real e, E;
      
	  minmax_edge_length3 (&e, &E, (*i)->simplex(),
			       &reference_location3, deformation_);
	  mine = min (mine, e);
	}
      printf ("min edge length %lf\n", mine); 
    }
    break ;
    
  case 'N':
    needle_surface_ = !needle_surface_;
    lighting_ = needle_surface_ || solid_;
    set_lighting();
    break;

  case 'A':
    {
      Needle_inserter3 * ins = needle_inserter_ ;
      ins->auto_insert_ = new Auto_needle_inserter3 ();
      
      ins->auto_needle_insert();
    }
    break ;
    
  default: 
    printf ("Unknown key %c\n", c);
    break ;
    
  }
}



void
Artisjokke_drawer3::set_needle()
{
  if (!needle_inserter_->live() || needle_inserter_->auto_insert_)
    return;
  
  Vector3 needle_handle;
  Vector3 needle_tip;

  needle_handle.set (-1.0,0.073,0.05);
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

  needle_dir = m * needle_dir;

  needle_handle = needle_handle + (needle_drag_. translate_)*0.10;
  needle_tip = needle_handle + needle_dir;

  needle_inserter_->move_needle (needle_handle, needle_tip);
}

void
Artisjokke_drawer3::draw_ground ()
{
  Real eps = 0.02;
  Real dim = 0.1;
  Vector3 vs [] = {
    Vector3(-eps, 0, -eps ),
    Vector3(dim + eps, 0, -eps ),
    Vector3(dim + eps, 0, eps + dim ),
    Vector3(-eps, 0,eps + dim )
  };

  set_color (ground_color);  
  glBegin (GL_POLYGON);
  enter_normal (Vector3 (0,1,0));

  // for (int i = 0; i  < 4;  i++i)
  //  ;
  for (int i = 4;  i--; )
    enter_vertex (vs[i]);
  glEnd ();
}


void
Artisjokke_drawer3::draw_needle()
{
  draw_ground( );
  set_color (needle_color);
  glLineWidth(4.0);
  glBegin (GL_LINES);

  Vector3 t= needle_inserter_->tip ();
  Vector3 h= needle_inserter_->handle ();
  
  enter_vertex (h.elts_);
  enter_vertex (t.elts_);  
  glEnd ();

  if (!solid_ && needle_surface_)
    visualize_needle_inserter (this, needle_inserter_);
  
}


#endif


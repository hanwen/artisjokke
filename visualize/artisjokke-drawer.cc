#ifdef OPENGL

#include <algorithm>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>

#include "opengl.hh"
#include "mesh-connectivity.hh"
#include "artisjokke-drawer.hh"
#include "deformation-state.hh"
#include "edge-watcher.hh"

using std::max;

float background [] = {1,1.2,1};
float background_strength =0.9 ;

void
enter_vertex (Vector3 const &v)
{
#ifdef SINGLE_PRECISION
  glVertex3fv (v.elts_);  
#else
  glVertex3dv (v.elts_);  
#endif
}

void
enter_normal (Vector3 const &v)
{
#ifdef SINGLE_PRECISION
  glNormal3fv (v.elts_);  
#else
  glNormal3dv (v.elts_);  
#endif
}

void
Artisjokke_drawer::rotate_drag (Quaternion * current,
				int sx, int sy,
				int beginx, int beginy,
				int ex, int ey)
{
  /* drag in progress, simulate trackball */
  float spin_quat[4];
  trackball(spin_quat,
	    (2.0*beginx - sx) / sx,
	    (sy - 2.0*beginy) / sy,
	    (2.0*ex - sx) / sx,
	    (sy - 2.0*ey) / sy);

  
  add_quats (spin_quat, *current, *current);

  if(current == &view_drag_.quat_)
    {
      set_camera (0);
    }
  else if (current == & needle_drag_.quat_)
    {
      set_needle ();
    }
  else if(current == &light_drag_.quat_)
    {
      set_lighting ();
    }
}

void
Artisjokke_drawer::focus_on_model ()
{
  Real rad =0.15;

  distance_ = max ((rad * 10), 2.0); 
  Real  fov = 2.5 * atan (rad  / distance_) * 180.0 / M_PI ; 
  view_drag_.scale_(0) = log (fov) / log (2);
  view_drag_.translate_ = Vector3 (-0.01,-0.01,0);
}

void
Artisjokke_drawer::setup()
{
  glDisable (GL_LINE_STIPPLE);
  glClear (GL_COLOR_BUFFER_BIT);
  glEnable (GL_BLEND);
  glBlendFunc (GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
  
  focus_on_model();
  Real m = 0;
  for (int i=0; i<3; i++)
    m = max (m, Real(background[i]));

  if (m)
    {
      for (int i=0; i<3; i++)
	 background[i] *= background_strength / m;
    }

  glClearColor (background[0],background[1],background[2],background[3]);
  set_frustrum();
  set_lighting ();
  set_depth_cue ();
  //  set_needle();
}


void
Artisjokke_drawer::set_frustrum ()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity ();

  Real fov =  pow (2, view_drag_.scale_(0));
  gluPerspective(fov, 1.0, 0.01, 50);
  set_camera (0);
}

void
Artisjokke_drawer::set_camera (int stereo_direction)
{
  /* view */
  glMatrixMode(GL_MODELVIEW);

  /* transformations */
  GLfloat m[4][4];

  glLoadIdentity();

  // stereo_direction_ * 0.5 * get_number_setting ("eye-distance");  
  Real eyedist = 0.0;

  Real d =  distance_ *  pow (2, view_drag_.scale_ (1)); 
  glTranslated (eyedist, 0, - d);

  Vector3 space_trans = view_drag_.translate_ ;
  glTranslated (space_trans(0),
		space_trans (1),
		space_trans (2));
  Real phi = atan(eyedist / d);

  glRotated (phi, 0,1,0);

  build_rotmatrix( m, view_drag_.quat_ );
  glMultMatrixf( &m[0][0] );
}  


void
Artisjokke_drawer::visualize_axes ()
{
  set_color (black_col);
  
  glBegin (GL_LINES);
  Vector3 x[3];
  for (int  i=3;i--;)
      x[i](i) = 0.1;
  Vector3 n (-1,0 ,1) ;
  for (int  i=3;i--;)
    {
      enter_vertex (n.elts_);
      enter_vertex ((x[i]+n).elts_);
    }
  glEnd();

  set_camera(0);
}

void
Artisjokke_drawer::draw_frame()
{
  glClear (GL_COLOR_BUFFER_BIT);
  glClear (GL_DEPTH_BUFFER_BIT);

  if (edges_ && mesh_->dimension() == 3)
    {
      Array<Simplex> *plex = edge_watcher_->edge_array ();

      if (plex)
	{
	  delete edge_array_;
	  edge_array_ = plex;
	}
    }
  
  opengl_visualize_mesh_connectivity();

  this->draw_needle ();
  if (this->axes_)
    this->visualize_axes( );
  glFlush();  
  deformation_->simulation_body ();
}
void
Artisjokke_drawer::opengl_visualize_mesh_connectivity()
{
  
}

Artisjokke_drawer::Artisjokke_drawer(Mesh_connectivity *top, Deformation_state*def)
{
  this->mesh_ = top;
  this->original_ = false;
  this->solid_ = false;
  this->deformation_ = def;
  this->boundary_ = false;
  this->edges_ = false;
  this->nodes_ = true;
  this->distance_ = 2.0;
  this->axes_ = false;
  this->lighting_ = false;
  this->fog_ = false;
  this->point_weight_= 2.0;
  this->line_weight_ = 1.0;
  
  edge_watcher_ = 0;
  edge_array_ = 0;

  if (top->dimension() == 3)
    {
      edge_watcher_ = new Edge_watcher;
      top->add_watcher (edge_watcher_);
    }
}

void
Artisjokke_drawer::set_color(float *c)
{
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, c);
  glColor4fv (c);
}

void
Artisjokke_drawer::scale_drag(Vector3 *current,
			    int sx, int sy,
			    int beginx, int beginy,
			    int ex, int ey)
{
  Vector3 add = Vector3(  (Real (ex) - beginx ) / sx,
			  (beginy - Real (ey) ) / sy ,
			  0.0);

  *current += add;

  if (current == &view_drag_.scale_)
    {
      // this is unintuitive.
      //      view_drag_.translate_(2) = 2 * view_drag_.scale_(1);
      
      set_frustrum ();
    }
  else if (current == &needle_drag_.scale_)
    {
      set_needle ();
    }
  else if(current == &light_drag_.scale_)
    set_lighting ();
}

void
Artisjokke_drawer::translate_drag (Vector3 *current,
				 int sx, int sy,
				 int beginx, int beginy,
				 int ex, int ey)
{
  Vector3 add = Vector3((Real (ex) - beginx ) / sx,
			(beginy - Real (ey) ) / sy ,
			0.0);

  *current += add;

  if (current == &view_drag_.translate_)
    {
      set_camera (0);
    }
   else if (current == &needle_drag_.translate_)
    {
      set_needle ();
    }
  else if(current == &light_drag_.translate_)
    {
      set_lighting ();
    }
}

void
Artisjokke_drawer::set_needle ()
{
}

void
Artisjokke_drawer::draw_needle ()
{
}

void
Artisjokke_drawer::set_lighting ()
{
  if (lighting_)
    {
      float m[4][4];
      build_rotmatrix(m, light_drag_.quat_);

      Matrix3 lm;
      Matrix3::from_opengl_matrixf (&lm, &m[0][0]);
      
      Vector3 light_pos = Vector3 (5, 5, 10);
      light_pos = lm * light_pos;

      /*
	4th arg : 0 == directional.
       */
      float position[] = {light_pos(0),
			  light_pos(1),
			  light_pos(2),0};

      // use this to show "back" of model.
      //      glFrontFace(GL_CCW);
      glEnable(GL_NORMALIZE);
      glEnable(GL_DEPTH_TEST);
      glLightfv(GL_LIGHT0, GL_POSITION, position);
      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT0);
      glShadeModel(GL_FLAT);
      glDisable ( GL_CULL_FACE);      
    }
  else
    {
      glShadeModel(GL_FLAT);            
      glDisable ( GL_CULL_FACE);
      glDisable (GL_LIGHTING);      
    }
}

void
Artisjokke_drawer::print_stats()
{
  set<Element*> *tris =mesh_->elements ();
  iterof (i,*tris);

  set<Node*> used;

  int tri_count = 0;
  int bound = 0;
  for (;i != tris->end(); i++)
    {
      for (int j = 03; j--;)
	{
	  used.insert ((*i)->node(j));
	  if (!(*i)->face (j)->mate ())
	    bound++;
	}
      tri_count++;
    }

  int nc = 0;
  for (iterof(i,used); i!= used.end(); i++)
    {
      nc++;
    }
  fprintf (stdout, "Elements %d, B-Faces %d, nods %d\n", tri_count, bound, nc); 
}
void
Artisjokke_drawer::process_key (char)
{
  
}

Artisjokke_drawer::~Artisjokke_drawer()
{
  
}

void
Artisjokke_drawer::set_depth_cue ()
{
  if (fog_)
    {
      glHint (GL_FOG_HINT, GL_NICEST);
      glFogf(GL_FOG_START, distance_ - 2);
      glFogf(GL_FOG_END, distance_ + 2);
      glEnable (GL_FOG);
      glEnable (GL_DEPTH_TEST);
      glFogi (GL_FOG_MODE, GL_LINEAR); // no typo.

      float fog_color[] =  {0.0, 0., 0., 0.0};
      glFogfv(GL_FOG_COLOR, fog_color);
    }
  else
    glDisable (GL_FOG);


 glFlush ();
}


#endif


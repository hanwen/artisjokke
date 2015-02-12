#ifndef BAZZOEN_DRAWER_HH
#define BAZZOEN_DRAWER_HH

#ifdef OPENGL

#include <time.h>

#include "vector.hh"
#include "matrix.hh"
#include "mesh-feature.hh"
#include "matrix.hh"
#include "trackball.hh"
#include "defo-proto.hh"
#include "array.hh"

struct Drag_state {
  Vector3 scale_;
  Vector3 translate_;
  Quaternion quat_;

  Drag_state (){
    reset ();
  }
  void reset ()
  {
    trackball(quat_, 0.0, 0.0, 0.0, 0.0 );
    scale_ = Vector3 (0,0,0);
    translate_ = Vector3 (0,0,0);
  }
};


class Artisjokke_drawer
{
public:
  Drag_state view_drag_;
  Drag_state clip_drag_;
  Drag_state needle_drag_;
  Drag_state light_drag_;
  Mesh_connectivity * mesh_;
  Real distance_ ;
  Real line_weight_;
  Real point_weight_;
  bool solid_;
  bool nodes_;
  bool edges_;
  bool boundary_;
  bool original_;
  bool axes_;
  bool lighting_; 
  bool fog_;
  
  Deformation_state * deformation_;

  Edge_watcher * edge_watcher_;
  Array<Simplex> * edge_array_;

  
  Artisjokke_drawer (Mesh_connectivity *, Deformation_state*);

  void print_stats();
  void visualize_axes ();  
  void set_lighting();
  void draw_frame();
  void rotate_drag (Quaternion*,int,int,int,int, int, int);
  void translate_drag (Vector3*,int, int, int, int, int, int);
  void scale_drag (Vector3*,int, int, int, int, int, int);
  void set_camera(int);
  void set_frustrum();
  void setup();
  void set_color (float*);
  void focus_on_model();
  void set_depth_cue ();
  
  virtual void opengl_visualize_mesh_connectivity();
  virtual void set_needle ();
  virtual ~Artisjokke_drawer  ();
  virtual void process_key (char );
  virtual void draw_needle ();
};

void enter_vertex (Vector3 const & );
void enter_normal (Vector3 const & );
void opengl_visualize_mesh_connectivity2 (Artisjokke_drawer*);
void opengl_visualize_mesh_connectivity3 (Artisjokke_drawer* draw);


extern float black_col[];
extern float body_col[] ;
extern float original_col [] ;
extern float black_col [] ;
extern float force_col [] ;

#endif
#endif

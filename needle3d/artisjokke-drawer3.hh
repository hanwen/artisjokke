#ifdef OPENGL

#ifndef ARTISJOKKE_DRAWER3_HH
#define ARTISJOKKE_DRAWER3_HH

#include "artisjokke-drawer.hh"
#include "stl.hh"
#include "proto.hh"

class Artisjokke_drawer3 : public Artisjokke_drawer
{
  Needle_inserter3 * needle_inserter_;

  bool needle_surface_ ;
  
public:
  Artisjokke_drawer3(Mesh_connectivity *top, Deformation_state*def,
		     Needle_inserter3*);
  
protected:
  void draw_ground ();
  virtual void opengl_visualize_mesh_connectivity();
  virtual void set_needle ();
  virtual void process_key (char );
  virtual void draw_needle ();
};


#endif
#endif

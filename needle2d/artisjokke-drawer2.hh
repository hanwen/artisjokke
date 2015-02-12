#ifdef OPENGL
#ifndef ARTISJOKKE_DRAWER2_HH
#define ARTISJOKKE_DRAWER2_HH

#include "artisjokke-drawer.hh"
#include "needle-inserter.hh"

class Artisjokke_drawer2 : public Artisjokke_drawer
{
public:
  Maubach_tree * tree  () const ;
  Needle_inserter * needle_inserter_;

  Artisjokke_drawer2(Mesh_connectivity *top, Deformation_state*def);

protected:
  virtual void opengl_visualize_mesh_connectivity();
  virtual void set_needle ();
  virtual void process_key (char );
  virtual void draw_needle ();

};


#endif
#endif

#ifndef NEEDLE_INSERTER3_HH
#define NEEDLE_INSERTER3_HH


#include "proto.hh"
#include "vector.hh"
#include "deformation-hook.hh"
#include "mesh-connectivity-watcher.hh"

class Needle_inserter3 : public Deformation_hook,  public Element_watcher
{
  enum Needle_state {
    INSIDE,
    OUTSIDE
  };
  Needle_state state_;
  bool live_;

  set<Element*> needle_elements_;
  map<Node*,Real> needle_node_radii_;

  
  bool filter_needle_rotations_ ;
  Real friction_density_ ; 
  Real dynamic_friction_factor_ ; 
  Real needle_radius_;
  Vector3 last_handle_;
  Vector3 last_tip_;

  
  Maubach_tree3 *mesh_;

  int refinement_level_;
  int rearrange_count_ ; 

  void drag_nodes (Vector3,Vector3);
  void rearrange_boundary_conditions ();
  Real friction_force (Real *,  Simplex const&, Vector3, Vector3, Real);
  Real force_distribution (Real,Real);
  void select_elements_around_needle (Element*, Vector3,Vector3);
  void update_boundary_nodes ();
  void update_needle_elements ();
  Real get_entry_parameter (Vector3,Vector3)const;
  void handle_tip_changes (Vector3, Vector3);
  void accrue_needle_elements (Vector3, Vector3);    
  void try_insert_element (Element*);
  bool is_needle_element (Element*, Vector3 , Vector3, Real ) const;
  virtual void signal_converged ();
public:
  Auto_needle_inserter3 * auto_insert_;
  void auto_needle_insert ();

  void refine_around_needle_tip (Vector3, Vector3);
  
  Needle_inserter3 (Maubach_tree3*, Deformation_state*);
  bool live ()const { return live_ ; }
  Vector3 tip () const;
  Vector3 handle () const;
  void move_needle (Vector3,Vector3);
  set<Face*> get_needle_surface ();
};


extern float needle_color [];

#endif

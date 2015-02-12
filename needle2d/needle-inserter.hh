/*
  declare Needle_inserter
 */

#ifndef NEEDLE_INSERTER
#define NEEDLE_INSERTER

#include "vector.hh"
#include "deformation-hook.hh"
#include "proto.hh"
#include "parray.hh"

class Needle_inserter : public Deformation_hook
{
  Maubach_tree * mesh_;
  Real edge_length_threshold_;
  Real dynamic_factor_;
  
  bool live_;
  bool relocate_;
  
  /*
    From TIP to  HANDLE.
   */
  Vector2 last_dir_;
  Vector2 last_tip_;
  Vector2 last_handle_;
  
  Link_array<Node> needle_;
  Friction_definition * friction_;
  int rearrange_count_ ;
  
  /*
    X-coordinate: along axis, Y coord: perpendicular.
  */
  Array<Vector2> needle_params_;
  enum
    {
      OUTSIDE,
      INSIDE,
    }
  state_;
public:
  Auto_needle_insert * auto_insert_;
  
  bool live () const;
  void suicide();
  int rearrange_count() const { return rearrange_count_ ; }
  void dump_needle_forces (const char *);
  
  Vector2 tip () const { return last_tip_; }
  Vector2 handle () const { return last_handle_; }
  
  Needle_inserter (Maubach_tree*, Deformation_state*);
  void move_needle (Vector2 h, Vector2 t);
  Vector2 tip_reference_location() const;

  void auto_needle_insert ();
protected:
  void rearrange_boundary_conditions ();
  void refine_around_needle_tip (Vector2,Vector2);
  void rearrange_boundary_conditions (Vector2,Vector2);
  void drag_fixed_nodes (Vector2,Vector2);
  void print_needle ();  
  void signal_converged ();
};


struct Element_intersection
{
  Element * elt_;
  Face * entry_face_;
  Face * exit_face_;

  Real entry_param_;
  Real exit_param_;

  Element_intersection();
};


Array<Element_intersection>
track_line_through_deformed_mesh (Maubach_tree * mesh,
				  Deformation_state * def,
				  Vector2 tip,
				  Vector2 handle
				  );
#endif

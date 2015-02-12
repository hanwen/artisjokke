#ifndef DEFORMATION_CONSTRAINT_HH
#define DEFORMATION_CONSTRAINT_HH

#include "defo-proto.hh"
#include "array.hh"
#include "vector.hh"
#include "stl.hh"

/*
  This specifies constraints on movements (displacements,
  velocity, acceleration, forces, etc.).

  TODO: we should have constraints of the form

  (Node => Vector)[3]

  (0 : fixed, 1 : line constraint, 2 : plane normal constraint)

  The fixed array allows for all these constraints, but only on
  coordinate aligned axes.
*/
class Deformation_constraints
{
  Deformation_state * parent_;
  /*
    These nodes can move along a line only.
   */
  map<int, Real * > line_constraints_;

  /*
    Move perpendicular to a line.
   */
  map<int, Real * > ortho_line_constraints_; 
  set<int> position_constraints_;
  
  bool changed_;

protected:
  friend class Deformation_state;

public:
  Deformation_constraints (Deformation_state *);
  void write_to_file (char const *fn) const;

  /*
    Change  constraints
   */
  void change_ortho_movement_constraint (Node *add, Real const *addv,
					 Node *remove);
  void change_linear_movement_constraint (Node *add, Real const *addv,
					  Node *remove);
  void fix_node (Node* nod, bool fix);

  void remove_all_node_constraints (Node* nod);
  /*
    Use constraints
   */
  void find_reactions (Real *dest, Real const * src)const;
  void apply_to_movement (Real * dest, Real const *mov)const;
  void apply_to_node_movement (Real *dest, Node * dest, Real const* mov)const;
  void node_reaction (Real*, Node*, Real const* ) const;
  Real reaction_length_sq (Real const *) const;
  bool satisfies_constraints (Real const*)const;

  /*
    inspect constraints
  */
  bool is_fixed (Node*)const;
  bool has_linear_movement_constraint (Node*)const;
  bool has_ortho_movement_constraint (Node*)const;  
  bool read_from_file (char const *fn);
  bool changed () const; 
  int dimension() const;
  
};

void fix_loose_components (Deformation_constraints *cons, Mesh_connectivity *top);


#endif

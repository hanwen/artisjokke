#ifndef MAUBACH_TREE_HH
#define MAUBACH_TREE_HH

#include "mesh-connectivity.hh"
#include "proto.hh"
#include "mesh-feature.hh"

Node * right_angle_node (Element*);
Face* long_edge (Element* e);

class Maubach_element : public Element 
{
public:
  Maubach_element * children_[2];
  Maubach_element (Simplex p);
};

class Maubach_tree :  public Mesh_connectivity 
{
  Maubach_element * top_nodes_[2];

protected:
  virtual Element* construct_element (Simplex);

public:
  Element * locate (Maubach_element* , Vector2,
		    Node_vector_func2 func, Deformation_state *);
  Element * locate (Vector2, Node_vector_func2 func,Deformation_state *);  

  Node *simple_bisect(Face *, Deformation_state*, Deformation_state*);
  Node *bisect_edge (Face *, Deformation_state*);
  void set_geometry (Deformation_state*, Real);
  Maubach_tree ();
};

void refine_around2 (Maubach_tree * mb,  Real h, Vector2 x,
		    Node_vector_func2, Deformation_state*);
void refine_uniformly2 (Maubach_tree *tree, Deformation_state *, Real factor);

#endif

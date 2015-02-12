#ifndef MAUBACH_TREE3_HH
#define MAUBACH_TREE3_HH

#include "mesh-connectivity.hh"
#include "mesh-feature.hh"
#include "defo-proto.hh"

class Maubach_element3 : public Element 
{
  int level_;
  friend class Maubach_tree3;
  
public:
  Maubach_element3 * children_[2];
  Node * unordered_simplex_[MAX_DIMENSION + 1];
  int level() const { return level_; }
  
  Simplex bisection_edge ();
  Maubach_element3 (Simplex p);
};


class Maubach_tree3 :  public Mesh_connectivity 
{
  Maubach_element3 ** top_nodes_;
  int dim_factorial_; 
protected:
  virtual Element* construct_element (Simplex);

public:

  Element * locate (Maubach_element3*, Vector3,
		    Node_vector_func3 func, Deformation_state *);
  Element * locate (Vector3, Node_vector_func3 func,Deformation_state *);  

  Node *simple_bisect(Element *, Deformation_state*, Deformation_state*);
  Node *bisect_edge (Element *, Deformation_state*);
  void set_geometry (Deformation_state*, Real);
  Maubach_tree3 ();
};

void refine_around3 (Maubach_tree3 * mb,  int level, Vector3 x,
	       Node_vector_func3 func, Deformation_state*def);

void subdivide_element_to (Maubach_tree3* mb, int level ,
			    Maubach_element3 * e);


void refine_uniformly3 (Maubach_tree3 * mb ,  int level, Deformation_state *ref);
void get_sub_elements (Element* e, set<Element*> * es);
#endif

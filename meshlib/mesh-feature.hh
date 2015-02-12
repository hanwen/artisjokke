#ifndef TOPOLOGICAL_OBJECT
#define TOPOLOGICAL_OBJECT
#include "simplex.hh"

class Face;
class Element;

class Mesh_feature
{
protected:
  Simplex simplex_;
  bool valid_;
  virtual ~Mesh_feature();
  
public:
  void print ()const;
  Mesh_feature(Simplex);
  Simplex simplex () const { return simplex_; }
  Simplex & simplex () { return simplex_; }  
  bool valid () { return valid_; }
  Node* node (int i) const { return simplex_.node(i); }
  bool parity ()  const{ return simplex_.parity(); }
  bool contains (Node*n) const { return simplex_.contains (n); }

  friend class Mesh_connectivity;
};

class Element : public Mesh_feature
{
  Face*faces_[MAX_DIMENSION+1];
public:
  Element (Simplex);
  bool has_face (Face*e) const;
  Face*face (int i) const { return faces_[i]; }
  Face * opposite_face (Node*) const;
  Node * opposite_node (Face*) const;
  Face * incident_face (Node*, int) const;
  Simplex opposite_simplex (Simplex super, Simplex sub) ;
  void set_incident_faces (Face**, Node*);
  
  friend class Mesh_connectivity;
  friend class Maubach_tree;

  static int compare (Element  * const&,Element  * const&);
};

class Face : public Mesh_feature
{
  Element* element_;
  Face * mate_;
public:
  
  Face(Simplex);
  Element* element() const { return element_;}
  Face* mate () const { return mate_;}

  friend class Mesh_connectivity;
};

Element * to_element (Mesh_feature *obj);
Face * to_face (Mesh_feature *obj);


#endif

#ifndef MESH_CONNECTIVITY_HH
#define MESH_CONNECTIVITY_HH

#include "mesh-proto.hh"
#include "parray.hh"
#include "stl.hh"
#include "simplex.hh"
#include "stl.hh"

typedef map<Simplex,Face*,Simplex_less> Face_map;

class Mesh_connectivity
{
  set<Element*> elements_;
  Face_map boundary_;
  int dimension_;

  Link_array<Node> nodes_; 
  Link_array<Mesh_connectivity_watcher> watchers_;
  
  Face *make_face (Simplex);
  void remove_face_pair (Face*);
  Face* make_face_pair(Simplex);
  void remove_face (Face*);
  void add_changed_element(Element*);
  void add_changed_boundary (Face*);
  void validate () const;
  void add_boundary (Face*);
  Face * remove_boundary (Simplex const &);

public:
  void add_watcher (Mesh_connectivity_watcher*); 
  Link_array<Node> const* node_array() const;
  Mesh_connectivity(int element_dim);
  int dimension () const { return dimension_ ; }
  int node_count () const;
  set<Face*> *faces () const;

  // duh
  set<Element*> * elements ()  { return  & elements_; }
  set<Face*> boundary() const;
  Face* find_boundary (Simplex)const;
  virtual Element* construct_element (Simplex);
  virtual void change_elements (map <Element*,  map<Node * , Node *> *> const &tri_map);
  virtual Link_array<Element> replace_elements (Link_array<Element>, Array<Simplex>);
  virtual ~Mesh_connectivity();
  virtual Node *make_node();
  void print  () const;
};


Link_array<Element> find_edgestar_elements (Face * entry, int k);
Link_array<Element> find_edgestar_elements_both (Face * entry, int k);


#if 0

/*
  Stuff from older prototypes
 */
void surface_flip_face_substitution (Link_array<Element> *dest,
				     Array<Simplex> *sdest,
				     Face* e1, Node * opp_surf_nod);
void flip_edge_substitution (Link_array<Element> *tets,
			     Array<Simplex> *sdest,
			     Face * entry, int i);
void flip_face_substitution (Link_array<Element> *tets,
			     Array<Simplex> *sdest,
			     Face * entry);
Link_array<Element> flip_face (Mesh_connectivity*, Face*);
Link_array<Element> flip_edge (Mesh_connectivity * mesh, Face * entry, int i);
Link_array<Element> surface_flip_face (Mesh_connectivity *m, Face* e1, Node * opp_surf_nod);

void remove_elements_containing (Mesh_connectivity*, Link_array<Node>); 
bool remove_element (Mesh_connectivity *, Element *);
void contract_face_to (Mesh_connectivity *top, Face*e, Node *dest);
Link_array<Element> insert_node_in_face (Mesh_connectivity*, Face*, Node*);
Link_array<Element> insert_node_in_single_face (Mesh_connectivity*, Face*, Node*);
Link_array<Element> insert_node_in_edge (Mesh_connectivity*m , Face* entry, int idx,  Node* n);
Element* remove_node_in_single_element (Mesh_connectivity *t,
					Element *contains, Node *victim);
Link_array<Element> insert_node_in_single_element  (Mesh_connectivity*m, Element*t, Node*n);

set<Element*> find_star_elements (Node*,Face*entry);
set<Element*> find_star_elements (Node *node, Element *entry);

Link_array<Face> inefficient_find_boundary_face_star (Element * e, Node *nod);
Link_array<Face> find_boundary_face_star (Face* entry, Node * nod);
Face* find_neighboring_boundary_face (Face * outer, Node *opp);

bool is_boundary_node (Node* node, Element * entry);
bool is_boundary_node (Node* nod, Face * e);
bool topologically_flippable_face (Face*f);
bool topologically_flippable_edge (Face * entry, int);

map<Node*, Element*> * find_all_nodes (Mesh_connectivity*m);
Element*find_element_containing_node (Mesh_connectivity*m,Node* v);
int element_count (Mesh_connectivity*);
int edge_count (Mesh_connectivity* top);


// dimension hopping.
map<Node*, Node*> *import_nodes (Mesh_connectivity*to, Mesh_connectivity*from);
Simplex import_simplex (map<Node*,Node*>  *, Simplex) ;
Mesh_connectivity* make_boundary_complex (Mesh_connectivity*) ;
Mesh_connectivity * reorder_mesh_connectivity (Mesh_connectivity* from);  

// io
Mesh_connectivity* read_element_set (char const *fn, Simplex_dimension);
void write_element_set (set<Element*> *tets, char const *fn);

map<int, Node*> *read_node_file (Mesh_connectivity *top, char const *fn);
void read_element_file (Mesh_connectivity*top, char const *fn,
			map<int,Node*> *nodmap);
void write_mesh_boundary (Mesh_connectivity *top, char const *fn);
void write_mesh_connectivity (Mesh_connectivity*mesh, char const *fn);

// inspection.
void print_element_set (set<Element*>);

#endif

#endif


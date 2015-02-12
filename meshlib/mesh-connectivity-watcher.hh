
#ifndef MESH_CONNECTIVITY_WATCHER
#define MESH_CONNECTIVITY_WATCHER

#include "mesh-proto.hh"
#include "stl.hh"
#include "mesh-connectivity.hh"	// urg.

class Mesh_connectivity_watcher
{
protected:
  Mesh_connectivity * mesh_;
  friend class Mesh_connectivity;

  Mesh_connectivity_watcher ();
  virtual ~Mesh_connectivity_watcher (){}
  virtual void add_changed_element (Element*){}
  virtual void add_changed_boundary (Face*){ }
  virtual void init_elements (set<Element*> const*);
  virtual void init_boundary (Face_map const*);
public:
  Mesh_connectivity * get_mesh () const
  {
    return mesh_;
  }
};

class Element_watcher : public virtual Mesh_connectivity_watcher
{
  set<Element*> * element_changes_;

protected:
  set<Element*> * get_changed_elements ();
  virtual void init_elements (set<Element*> const*);
  Element_watcher ();
  ~Element_watcher();

public:
  virtual void add_changed_element (Element*);
};

class Boundary_watcher : public virtual Mesh_connectivity_watcher
{
  set<Face*> * face_changes_;
protected:
  ~Boundary_watcher();
  Boundary_watcher();
  set<Face*> * get_changed_boundary ();
  virtual void init_boundary (Face_map const*);
public:
  virtual void add_changed_boundary (Face*);  
};

class Element_boundary_watcher : public Boundary_watcher, public Element_watcher{
public:
  Element_boundary_watcher();
protected:
  virtual void add_changed_element (Element*e);
  virtual void add_changed_boundary (Face*);
  virtual void init_elements (set<Element*> const*);
  virtual void init_boundary (Face_map const*);
};
				 
#endif

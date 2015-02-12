#include "mesh-connectivity-watcher.hh"


Mesh_connectivity_watcher::Mesh_connectivity_watcher()
{
  mesh_ = 0;
}

void
Element_watcher::add_changed_element (Element*e)
{
  if (!element_changes_)
    element_changes_ = new set<Element*>;
  element_changes_->insert (e);
}

void
Boundary_watcher::add_changed_boundary (Face*e)
{
  if (!face_changes_)
    face_changes_ = new set<Face*>;
  face_changes_->insert (e);
}

set<Element*> *
Element_watcher::get_changed_elements ()
{
  set<Element*> *r  = element_changes_;
  if (element_changes_)
    {
      element_changes_ = 0;
    }
  
  return r ;
}

Element_watcher::Element_watcher ()
{
  element_changes_ = 0;
}

Boundary_watcher::Boundary_watcher()
{
  face_changes_ = 0;
}

set<Face*> *
Boundary_watcher::get_changed_boundary ()
{
  set<Face*> *r  = face_changes_;
  if (face_changes_)
    {
      face_changes_ = 0;
    }
  
  return r ;
}

Element_watcher::~Element_watcher ()
{
  delete element_changes_; 
}

Boundary_watcher::~Boundary_watcher()
{
  delete face_changes_; 
}

void
Element_watcher::init_elements (set<Element*> const *elts)
{
  assert (!element_changes_);
  element_changes_ = new set<Element*> (*elts);
}

void
Boundary_watcher::init_boundary (Face_map const *bdry)
{
  for (iterof(i, *bdry); i != bdry->end (); i++)
    {
      add_changed_boundary (i->second);
    }
}

void
Mesh_connectivity_watcher::init_boundary (Face_map const * ) {}

void
Mesh_connectivity_watcher::init_elements (set<Element*> const * ) {}

void Element_boundary_watcher::add_changed_element (Element*e)
{
  Element_watcher::add_changed_element (e);
}

void
Element_boundary_watcher::add_changed_boundary (Face*f)
{
  Boundary_watcher::add_changed_boundary (f);
}

void
Element_boundary_watcher::init_elements (set<Element*> const*s)
{
  Element_watcher::init_elements (s);
}

void
Element_boundary_watcher::init_boundary (Face_map const*s)
{
  Boundary_watcher::init_boundary (s);  
}
						
Element_boundary_watcher::Element_boundary_watcher(){

}

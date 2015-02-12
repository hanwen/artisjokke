#include <assert.h>
#include <string.h>

#include "mesh-feature.hh"


Mesh_feature::Mesh_feature (Simplex s)
  : simplex_(s)
{
  valid_ = 0;
}

Element::Element(Simplex s)
  : Mesh_feature(s)
{
  for (int i = s.count (); i--;)
    faces_[i] = 0;
}

Face*
Element::opposite_face (Node*n) const
{
  for (int i = 0; i < simplex_.count(); i ++)
    if (node(i) == n)
      return face(i);

  return 0;
}

Node*
Element::opposite_node (Face*n) const
{
  for (int i = 0; i < simplex_.count(); i ++)
    if (face(i) == n)
      return node(i);

  assert(false);
  return 0;  
}

Face::Face( Simplex s)
  : Mesh_feature(s)
{
  mate_ = 0;
  element_ = 0;
}

void
Element::set_incident_faces (Face**es, Node*n)
{
  int i = simplex_.index (n);
  int c = simplex_.count ();
  for (int j = 1; j  < c; j++)
    {
      es[j] = faces_[(i+j)%c];
    }
}

void
Mesh_feature::print ()const
{
  simplex_.print();
}

bool
Element::has_face (Face*e) const
{
  for (int i = simplex_.count(); i--; )
    if (this->faces_[i] == e)
      return 1;
  return 0;
}

Face* 
Element::incident_face  (Node*n , int k)const
{
  int i = simplex_.index  (n);

  assert (k < simplex_.count () - 1);
  return faces_[(1 + i + k)%simplex_.count()];
}

int
Element::compare (Element  * const&a, Element  * const&b)
{
  return Simplex::compare (a->simplex(), b->simplex());
}


Element*
to_element (Mesh_feature *o )
{
  return dynamic_cast<Element*> (o);
}

Face*
to_face (Mesh_feature *o )
{
  return dynamic_cast<Face*> (o);
}

Mesh_feature::~Mesh_feature()
{

}

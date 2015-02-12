#include "deformation-state.hh"
#include "node.hh"


Vector2
deformed_location2 (Node* nod, Deformation_state const *def)
{
  int sd = def->spatial_dimension ();
  assert (sd == 2);
  Vector2 v = reference_location2 (nod, def);
  if (def)
    {
      ((Deformation_state*)def)->check_for_topology_update();
      int idx = sd * nod->number ();
      for (int j = sd; j--;)
	v(j) += def->displacements_[idx + j];
    }

  return v;
}

Vector2
reference_location2 (Node *node, Deformation_state const* ref)
{
  assert (ref->spatial_dimension() == 2);
  return extract_vector2 (node, ref->reference_locations_.access_array());
}

Vector2
get_external_force2 (Node* nod, Deformation_state const *def) 
{
  assert (def->spatial_dimension() == 2);
  return extract_vector2 (nod, def->external_force_.access_array());
}


Vector2
extract_vector2 (Node* n , Real const * b )
{
  int dim = 2;
  int j = n->number()* dim;
  Vector2 v;
  for (int i = 0; i < dim; i++)
    v(i) = b[j+i];
  return v;
}



Vector3
deformed_location3 (Node* nod, Deformation_state const *def)
{
  int sd = def->spatial_dimension ();
  assert (sd==3);  
  Vector3 v = reference_location3 (nod, def);
  if (def)
    {
      ((Deformation_state*)def)->check_for_topology_update();
      int idx = sd * nod->number ();
      for (int j = sd; j--;)
	v(j) += def->displacements_[idx + j];
    }

  return v;
}

Vector3
reference_location3 (Node *node, Deformation_state const* ref)
{
  assert (ref->spatial_dimension()==3);  
  return extract_vector3 (node, ref->reference_locations_.access_array());
}

Vector3
get_external_force3 (Node* nod, Deformation_state const *def) 
{
  assert (def->spatial_dimension()==3);  
  return extract_vector3 (nod, def->external_force_.access_array());
}


Vector3
extract_vector3 (Node* n , Real const * b )
{
  int dim = 3;
  int j = n->number()* dim;
  Vector3 v;
  for (int i = 0; i < dim; i++)
    v(i) = b[j+i];
  return v;
}

#include "deformation-state.hh"
#include "mesh-connectivity.hh"
#include "geometry2.hh"
#include "mesh-geometry.hh"
#include "setting.hh"

void
apply_gravity (Deformation_state*def, int axis, int sign)
{
  Mesh_connectivity*top =def->get_mesh ();
  
  Array<Real> lv = lumped_volumes (top, def);
  Real f = 0.0;
  Link_array<Node> const*nods = top->node_array (); 

  Real d = get_number_setting ("density")* get_number_setting ("gravity");
  for (int i = lv.size(); i--;)
    {
      Real df = lv[i]*d;

      Vector2 dfv;
      dfv(axis) = -df* sign;

      def->add_force (nods->elem(i),dfv);
      f += df;
    }
}

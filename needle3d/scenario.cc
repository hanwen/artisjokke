#include <time.h>
#include <math.h>

#include "maubach-tree3.hh"
#include "deformation-state.hh"
#include "proto.hh"
#include "mesh-geometry.hh"
#include "setting.hh"
#include "auto-inserter3.hh"
#include "needle-inserter3.hh"

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

      Vector3 dfv;
      dfv(axis) = -df* sign;

      def->add_force (nods->elem(i),dfv);
      f += df;
    }
}


clock_t needle_timer;

void
auto_insert_scenario (Maubach_tree3 *,
		      Deformation_state *def,
		      Needle_inserter3 *ins)
{
  ins->auto_insert_ = new Auto_needle_inserter3 ();
  needle_timer = clock();

  ins->auto_needle_insert();
  
#ifndef OPENGL
  def->update_forces();
  while (!def->good_solution() || ins->auto_insert_)
    {
      def->simulation_body ();
    }
#endif


  
}


void
auto_select_scenario (Maubach_tree3 * tree,
		      Deformation_state *def,
		      Needle_inserter3 * ins )
{
  Real y = 0.0701;
  Vector3 handle = Vector3 (-0.1, y, 0.0501);
  ins->move_needle (handle, 
		    Vector3 (-0.05, y, 0.0501));

  Real h  = 0.1 * pow (2.0 , - get_number_setting ("refinement-level") / 3.0) * 0.5;
  
  for (Real x = 0.0; x < 0.05; x += h)
    ins->refine_around_needle_tip (handle, Vector3 (x, handle(Y_AXIS), handle (Z_AXIS)));
  ins->move_needle (handle,
		    Vector3 (0.05, y, 0.0501));
}

void
test_scenario (char const *ch,
	       Maubach_tree3 * tree,
	       Deformation_state *def,
	        Needle_inserter3 * ins )
{
  if (!ch[0])
    return;

  if (0)
    {
      ;
    }
  else if (!strcmp(ch, "gravity"))
    {
      apply_gravity (def, 0, -1);      
    }
  else if (!strcmp(ch, "auto-insert"))
    {
      auto_insert_scenario (tree, def, ins);
    }
  else if (!strcmp(ch, "auto-select"))
    {
      auto_select_scenario (tree, def, ins);
    }
  else if (!strcmp (ch, "refinement-demo"))
    {
      int rl = (int) get_number_setting ("refinement-level");
      refine_around3 (tree, rl, Vector3 (0.0852, 0.0701, 0.01111),
		      reference_location3, def);
    }
 else
    {
      printf ("Unknown scenario %s, try\n{gravity,auto-insert,auto-select}\n", ch);
      exit (2);
    }
      
}



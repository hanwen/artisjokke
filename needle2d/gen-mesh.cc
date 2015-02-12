#include <stdio.h>
#include <math.h>

#include "vector.hh"
#include "node.hh"
#include "mesh-connectivity.hh"
#include "geometry2.hh"
#include "matrix.hh"
#include "deformation-state.hh"

void
generate_square_mesh (Mesh_connectivity* m, Deformation_state * def, int w, int h)
{
  assert (h>= 2 && w >= 2);

  Link_array<Node> last_col;
  for (int x = 0 ; x < w; x++)
    {
      Link_array<Node> col;
      for (int y = 0 ; y < h; y++)
	{
	  Vector2 loc;

	  loc (0) = (1.0*x)/(w-1);
	  loc (1) = (1.0*y)/(h-1);

	  Node * n = m->make_node ();
	  def->set_reference_location (n, loc);
	  
	  col.push (n);
	}

      if ( last_col.size())
	{
	  for (int y = 1 ; y < h; y++)
	    {
	      Node *ns1[] = {last_col[y], col[y], col[y-1]};
	      Node *ns2[] = {last_col[y], col[y-1], last_col[y-1]};

	      Simplex s1 (TRIANGLE, ns1, 0);
	      Simplex s2 (TRIANGLE, ns2, 0);
	      if (simplex_area (s1, &reference_location2, def  ) < 0)
		s1 = s1.get_mate ();
	      if (simplex_area (s2, &reference_location2, def) < 0)
		s2 = s2.get_mate ();

	      Array<Simplex> ns;

	      ns.push (s1);
	      ns.push (s2);
	      Link_array<Element> ts;

	      m->replace_elements (ts, ns);
	    }
	}
      last_col = col;
    }
}


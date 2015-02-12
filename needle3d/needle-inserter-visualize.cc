#ifdef OPENGL
#include "artisjokke-drawer3.hh"
#include "opengl.hh"
#include "needle-inserter3.hh"
#include "deformation-state.hh"
#include "mesh-visualize3.hh"

void
visualize_needle_inserter (Artisjokke_drawer3 * draw,
			   Needle_inserter3 * needle)
{
  set<Face*> fs = needle->get_needle_surface ();

  draw->set_color (needle_color);
  for (iterof (i, fs); i != fs.end(); i++)
    {
      Face *  f =(*i);

      opengl_draw_face (draw, f, true,
			&deformed_location3, needle->deformation ());
    }
}




#endif

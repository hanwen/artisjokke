#include "mesh-connectivity.hh"
#include "mesh-feature.hh"

Link_array<Element>
find_edgestar_elements_both (Face * entry, int j )
{
  Link_array<Element> es = find_edgestar_elements (entry, j);

  if (entry->mate ()
      && entry->mate ()->element () != es.top ())
    {
      es.concat (find_edgestar_elements (entry->mate () , j));
    }

  return es;
}


Link_array<Element>
find_edgestar_elements (Face * entry, int j )
{
  Link_array<Element> tets;

  Face * f = entry;
  Node  *v = entry->node (j);
  while (f)
    {
      Element* tet = f->element();
      tets.push (tet);

      Face *exit = tet->opposite_face (v);
      v = tet->opposite_node (f);
      f = exit->mate ();

      if (f == entry)
        f = 0;
    }
  
  return tets;
}



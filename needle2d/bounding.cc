#include "node.hh"
#include "deformation-state.hh"
#include "geometry2.hh"
#include "bounding.hh"

void
nodes_bounding_sphere (Vector2 *center, Real * radius,
		       Link_array<Node> const &nodes, Deformation_state *def
		       )
{
  *center = Vector2 (0,0);
  
  for (int i = nodes.size();
       i--;)
    {
      *center += deformed_location2 (nodes[i], def);
    }

  *center *=  1.0 / nodes.size();
  *radius =0.0;  

  for (int i = nodes.size();
       i--;)
    {
      *radius = *radius >? (deformed_location2 (nodes[i], def) - *center).length();
    }
}


void nodes_bounding_box (Vector2 *min, Vector2 * max,
			 Link_array<Node> const &nodes, Deformation_state*def
			 )

{
  int DIM = def->spatial_dimension ();
  const Real inf = 1e8;
  for (int i= DIM; i--; )
    {
      (*min)(i) = inf;
      (*max)(i) = -inf;
    }
  
  for (int j = nodes.size();
       j--;)
    {
      Vector2 l (deformed_location2 (nodes[j],def));
      
      for (int i=DIM; i--; )
	{
	  (*min)(i) = (*min)(i) <? l(i);
	  (*max)(i) = (*max)(i) >? l(i);
	}
    }
}

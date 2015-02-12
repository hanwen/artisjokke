#include "geometry.hh"
#include "mesh-connectivity.hh"
#include "big-vector.hh"
#include "mesh-feature.hh"
#include "deformation-state.hh"

/*
  TODO: this is Oh (n log (n)). (ugh)
 */
Array<Real>
lumped_volumes (Mesh_connectivity*fem, Deformation_state*def)
{
  Array<Real> vols;
  vols.set_size (fem->node_count ());

  big_vector_nullify (vols.unsafe_access_array (), vols.size());

  int dim = fem->dimension();

  set<Element*> *es = fem->elements();
  for (iterof(i,*es); i != es->end(); i++)
    {
      /*
	TODO: dimension dependence!
       */
      Real vol = 0.0;

      if (dim == 2)
	vol = simplex_volume2 ((*i)->simplex(), &reference_location2, def);
      else
	vol = simplex_volume3 ((*i)->simplex (), &reference_location3, def);

      /*
	distribute over nodes.
       */
      vol /=  dim + 1;
      
	
      
      for (int j= fem->dimension() + 1; j--;)
	{
	  vols[(*i)->node (j)->number ()] +=  vol;
	}
    }
  return vols;
}


#define LOCATE_FUNC locate_element_containing_point2
#define VEC Vector2
#define FUNC Node_vector_func2
#define VOLUME simplex_volume2
#define DEFORMED_SIMPLEX_VOLUME virtual_deformed_simplex_volume2
#define MINMAX_EDGE_LENGTH minmax_edge_length2

#include "mesh-geometry.icc"

#undef DEFORMED_SIMPLEX_VOLUME
#undef LOCATE_FUNC
#undef VEC
#undef FUNC
#undef VOLUME
#undef MINMAX_EDGE_LENGTH 

#define MINMAX_EDGE_LENGTH minmax_edge_length3

#define LOCATE_FUNC locate_element_containing_point3
#define VEC Vector3
#define FUNC Node_vector_func3
#define VOLUME simplex_volume3
#define DEFORMED_SIMPLEX_VOLUME virtual_deformed_simplex_volume3

#include "mesh-geometry.icc"

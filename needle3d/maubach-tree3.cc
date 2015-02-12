#include <string.h>

#include "geometry.hh"
#include "maubach-tree3.hh"
#include "vector.hh"
#include "deformation-state.hh"

Maubach_element3::Maubach_element3 (Simplex p)
  : Element (p)
{
  level_ = 0;
  children_[0] = children_[1] = 0;
}

Element* 
Maubach_tree3::construct_element (Simplex p)
{
  return new Maubach_element3 (p); 
}


void
Maubach_tree3::set_geometry (Deformation_state * def, Real h)
{
  Vector3 vs[8] =
    {
      Vector3 (0,0,0),
      Vector3 (h,0,0),
      Vector3 (0,h,0),
      Vector3 (0,0,h),
      Vector3 (h,h,0),
      Vector3 (0,h,h),
      Vector3 (h,0,h),
      Vector3 (h,h,h)
    };

  assert (node_count() == 8);

  Link_array<Node> const *nods = node_array();
  for (int i = 8; i--;)
    def->set_reference_location (nods->elem (i), vs[i]);

  iterof (i,*elements());

  assert (simplex_volume3 ((*i)->simplex(), reference_location3, def) > 0);   
}

Maubach_tree3::Maubach_tree3 ()
  : Mesh_connectivity (3)
{
  dim_factorial_ = 6;
  top_nodes_ = new Maubach_element3 *[dim_factorial_];
  
  int plexes[][4]  ={
    {0,1,4,7},
    {0,1,6,7},
    {0,2,4,7},
    {0,2,5,7},
    {0,3,5,7},
    {0,3,6,7},
  };
  int parities[] = {
    0,
    1,
    1,
    0,
    1,
    0};
  
  Node * nods [8];
  for (int i = 0 ; i <8 ; i++)
    nods[i]=     make_node ();

  
  Link_array<Element> old;
  Array<Simplex> newelts;

  for (int i = 0 ; i <6 ; i++)
    {
      Node * plex_nods [4];
      for (int j = 4; j--;)
	{
	  plex_nods[j] = nods[plexes[i][j]];
	}
      Simplex p (TETRAHEDRON, plex_nods, ! parities[i]);

      newelts.push (p);
    }
  old= replace_elements(old, newelts);
  for (int i =  old.size () ; i--;)
    {
      Maubach_element3* me = dynamic_cast<Maubach_element3*> (old[i]);
      top_nodes_[i] =  me;
      
      for (int j = 4; j--;)
	me->unordered_simplex_[j] = nods[plexes[i][j]];
    }
}


Simplex
Maubach_element3::bisection_edge ()
{
  int n = simplex().dimension();
  int k = n - (level_ %  n);

  Node *x0 = unordered_simplex_[0];
  Node *xk = unordered_simplex_[k];
  Node *xs[] = {x0, xk}; 
  Simplex e (EDGE, xs, 0);
  if (e.parity())
    e = e.get_mate ();

  return e;
}

void
face_index_for_edge (Element* e,
		     Node *x0, Node *x1,
		     Face **fp, int *jp)
{
  Simplex other = e->simplex().get_subset (x0).get_subset (x1);
  
  *fp = e->opposite_face (other.node (0));
  *jp = (*fp)->simplex().index (other.node (1));
}
		     

Node*
Maubach_tree3::bisect_edge (Element* e,
			    Deformation_state * def)
{
  Maubach_element3 * me = dynamic_cast<Maubach_element3*> (e);

  Simplex edge = me->bisection_edge ();
  

  bool change = false;
  do
    {
      /*
	We keep looping until all elements around EDGE are compatibly divisible.
       */
      change = false;
      assert (e->valid ());
      
      Face * f;
      int j ;
      face_index_for_edge (me, edge.node(0), edge.node(1), &f, &j);

      Link_array<Element> elements = find_edgestar_elements_both (f,j);
      for (int i = 0; i < elements.size(); i++)
	{
	  Maubach_element3 * me = dynamic_cast<Maubach_element3*> (elements [i]);
	  if (!elements[i]->valid ())
	    continue;
	  
	  Simplex ee = me->bisection_edge ();
	  if (ee.parity())
	    ee = ee.get_mate();
	  
	  if (ee != edge)
	    {
	      change = true;
	      bisect_edge (elements[i], def);
	    }
	}
      
    }
  while (change);

  return  simple_bisect (e, def, def);
}

Node*
Maubach_tree3::simple_bisect (Element* e,
			      Deformation_state * ref,
			      Deformation_state * def)
{
  Maubach_element3 * me = dynamic_cast<Maubach_element3*> (e);
  int k = dimension () - (me->level_ %  dimension());

  Node *x0 = me->unordered_simplex_[0];
  Node *xk =  me->unordered_simplex_[k];
  Simplex other = e->simplex().get_subset (x0).get_subset (xk);
  
  Vector3 ctr = 0.5 * (reference_location3 (x0, ref) +
		       reference_location3 (xk, ref));

  Array<Real> vars; 
  if (def)
    vars = def->interpolate_deformation_variables (ctr, e);
  
  Node * znod = make_node();
  ref->set_reference_location (znod, ctr);
  if (def)
    {
      def->set_deformation_variables (znod, vars);
    }

  
  Face * f = e->opposite_face (other.node (0));
  int j = f->simplex().index (other.node( 1));
  
  Link_array<Element> elements = find_edgestar_elements_both (f,j);
  Array<Simplex> newplexes ;
 
  for (int i = 0; i < elements.size(); i++)
    {
      
      assert (dynamic_cast<Maubach_element3*>(elements[i])->level_ == me->level_);
      
      Simplex p1 = elements[i]->simplex().get_substituted (xk, znod);
      Simplex p2 = elements[i]->simplex().get_substituted (x0, znod);

      newplexes.push (p1);
      newplexes.push (p2);
    }

  Link_array<Element> newelts
    = replace_elements (elements, newplexes);

  for (int i = 0; i < elements.size(); i++)
    {
      Maubach_element3 *orig =dynamic_cast<Maubach_element3*> (elements[i]);
      Maubach_element3* me31 = dynamic_cast<Maubach_element3*> (newelts[2*i+0]);
      Maubach_element3* me32 = dynamic_cast<Maubach_element3*> (newelts[2*i+1]);

      memcpy (me31->unordered_simplex_ ,
	      orig->unordered_simplex_,
	      sizeof (orig->unordered_simplex_));
      memcpy (me32->unordered_simplex_ ,
	      orig->unordered_simplex_,
	      sizeof (orig->unordered_simplex_));
      
      me31->level_ = orig->level_+1;
      me32->level_ = orig->level_+1;

      orig->children_[0] = me31;
      orig->children_[1] = me32;
      
      me31->unordered_simplex_[k] = znod;
      for (int  i = 0; i < k ; i++)
	me32->unordered_simplex_[i] = me32->unordered_simplex_ [i+1];
      me32->unordered_simplex_[k] = znod;
    }


  return znod;
}


/*
  Due to deformations, the bounding volume hierarchy does not work
  well for point location in the deforming object. Solution: lookup
  point in reference location, and use that as initial result for
  jump-and-walk point location.
 */
Element* 
Maubach_tree3::locate  (Vector3 loc, Node_vector_func3 func,
		       Deformation_state*def)
{
  Element *e = 0;
  for (int j = 0; !e && j < dim_factorial_; j++)
    e=  locate (top_nodes_[j], loc, &reference_location3, def);

  if (func == &reference_location3)
    return e;

  if (!e)
    {
      e = * elements ()->begin ();
    }
  e = locate_element_containing_point3 (this, e, loc, func, def);
  
  return e;
}


/*
  This is very much ugh, but unfortunately, we have query nodes that
  are structurally on edges in the mesh.
 */
#define CONTAINMENT_EPS 1e-4

Element* 
Maubach_tree3::locate  (Maubach_element3 * tr, Vector3 loc,
			Node_vector_func3 func, Deformation_state* def)
{
  if (!tetrahedron_contains_point3 (tr->simplex(), loc, func, def, CONTAINMENT_EPS))
    return 0;

  if (!tr->children_[0])
    return tr;
  
  for (int j = 0; j < 2; j++)
    {
      Element* e  = locate (tr->children_[j], loc, func, def);
      if (e)
	return  e;
    }

  return 0;
}



void
refine_around3 (Maubach_tree3 * mb,  int level, Vector3 x,
	       Node_vector_func3 func, Deformation_state*def)
{

  do
    {
      Element *e = mb->locate (x, func,def);
      if (!e)
	return ;

      Maubach_element3* me = dynamic_cast<Maubach_element3 *> (e);
      if (me->level () >= level)
	break;
      
      mb->bisect_edge (e, def);
    }
  while (1);
}
     

void
refine_uniformly3 (Maubach_tree3 * mb ,  int level, Deformation_state *ref)
{
  for (int l = 0;  l <= level ; l++)
    {
      set<Element*> es = *mb->elements();
	
      for (iterof (i, es); i != es.end (); i++)
	{
	  Element * e = *i;

	  if (!e->valid ())
	    continue;

	  Maubach_element3 * mbe = dynamic_cast<Maubach_element3*> (e);
	  if (mbe->level () < l)
	    mb->simple_bisect (mbe, ref, 0);
	}
    }
}



void
get_sub_elements (Element* e, set<Element*> * es)
{  
  Maubach_element3* me = dynamic_cast <Maubach_element3*> (e); 
  
  if (me->children_[0])
    {
      get_sub_elements (me->children_[0], es);
      get_sub_elements (me->children_[1], es);
    }
  else
    {
      assert (e->valid());
      es->insert (e);
    }
}

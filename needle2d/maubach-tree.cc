#include <math.h>

#include "maubach-tree.hh"
#include "geometry.hh"
#include "mesh-feature.hh"
#include "deformation-state.hh"
#include "setting.hh"

/*
  By construction, the newest node in the element is always opposite
  the long edge. Hence the right angle is the last node in the simplex
  (assuming that simplexes are sorted in ascending order)
 */
Node *
right_angle_node (Element*e)
{
  return e->node (e->simplex().count ()-1);
}

Face*
long_edge (Element* e)
{
  return e->opposite_face (right_angle_node (e));
}



Maubach_tree::Maubach_tree ()
  : Mesh_connectivity (2)
{
  Node*ns [4]; 
  for (int i = 0; i < 4; i++)
    {
      ns[i] = make_node ();
    }
  Node * tris[2][3]= {
    {ns[0], ns [1], ns[2]},
    {ns[0], ns [1], ns[3]},
  };
  
  Simplex s1 (TRIANGLE, tris[0], 1);
  Simplex s2 (TRIANGLE, tris[1], 0);

  Array<Simplex> plexes;
  plexes.push (s1);
  plexes.push (s2);

  Link_array<Element> elts;

  elts = replace_elements (elts, plexes);

  top_nodes_ [0] = dynamic_cast< Maubach_element*> (elts[0]);
  top_nodes_ [1] = dynamic_cast< Maubach_element*> (elts[1]);  
}

void
Maubach_tree::set_geometry (Deformation_state * def, Real h)
{
  Vector2 vs[4] =
    {
      Vector2 (0,0),
      Vector2 (h,h),
      Vector2 (h,0),
      Vector2 (0,h)
    };

  assert (node_count() == 4);

  Link_array<Node> const *nods = node_array();
  for (int i = 4; i--;)
    def->set_reference_location (nods->elem (i), vs[i]);
}

void
refine_uniformly_once (Maubach_tree *tree, Deformation_state *def)
{
  set<Face * > to_refine;
  set<Element*> *es = tree->elements();
  for (iterof (i, *es); i != es->end (); i++)
    {
      Face * f = (*i)->opposite_face (right_angle_node (*i));
      to_refine.insert (f);
    }

  for (iterof (i, to_refine); i!=to_refine.end (); i++)
    {
      if ( (*i)->valid ())
	tree->simple_bisect (*i,  def, 0 );
    }
}

void
refine_uniformly2 (Maubach_tree *tree, Deformation_state *def, Real factor )
{
  int times = int (ceil (2*log (factor) /log (2.0)));

  while (times --)
    {
      refine_uniformly_once (tree, def);
    }
}


/*
  This is very much ugh, but unfortunately, we have query nodes that
  are structurally on edges in the mesh.
 */
#define CONTAINMENT_EPS 1e-4

bool
contains(Simplex const & plex, Vector2 loc, Node_vector_func2 func, Deformation_state *def)
{
  Vector2 vs[3];
  for (int i = 0; i < plex.count (); i++)
    vs[i]= (*func) (plex.node(i), def);

  Real b[3];
  barycentric_coordinates2 (vs, b, loc);

  for (int j = 3; j--;)
    if (b[j] < - CONTAINMENT_EPS)
      return false;
  return true;
}

Element* 
Maubach_tree::locate  (Maubach_element * tr, Vector2 loc,
		       Node_vector_func2 func, Deformation_state* def)
{
  if (!contains(tr->simplex(), loc, func, def))
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


/*
  Due to deformations, the bounding volume hierarchy does not work
  well for point location in the deforming object. Solution: lookup
  point in reference location, and use that as initial result for
  jump-and-walk point location.
 */
Element* 
Maubach_tree::locate  (Vector2 loc, Node_vector_func2 func,
		       Deformation_state*def)
{
  Element *e = 0;
  for (int j = 0; !e && j < 2; j++)
    e=  locate (top_nodes_[j], loc, &reference_location2, def);

  if (func == &reference_location2)
    return e;


  static  bool reloc_initted;
  static  bool reloc;
  if (!reloc_initted)
    {
      reloc_initted = true;
      reloc = get_bool_setting ("relocate-nodes");
    }
 
  if (reloc)
    {    
      iterof (i, *elements());

      /*
	Jump and walk: walk from element until we find it... 
      */
      while (!e && i != elements()->end())
	{
	  /*
	    Ugh. this ugly loop is needed for searching with relocation
	    code (flat elements cause jump & walk searches to fail.).

	    It causes strange behavior when inverting elements. UGh.
	  */
	  Element * start = *i;
	  e =locate_element_containing_point2 (this, start, loc, func, def);
	  i++;
	}
    }
  else
    e =locate_element_containing_point2 (this, *elements()->begin(),
					loc, func, def);
  
  return e;
}


Node*
Maubach_tree::simple_bisect (Face* f,
			     Deformation_state *ref,
			     Deformation_state *def)
{
  Link_array<Element> rm;
  Array<Simplex> ns;
  
  Element *e = f->element();

  assert (e->opposite_node (f) 
	  == right_angle_node ( e));
  
  rm.push (e);

  Vector2 ctr = simplex_centroid2 (f->simplex(), &reference_location2, ref);

  Array<Real> vars; 
  if (def)
    vars = def->interpolate_deformation_variables (ctr, e);

  Node* nod = make_node ();

  ref->set_reference_location (nod, ctr);
  if (def)
    {
      def->set_deformation_variables (nod, vars);
    }
  
  ns.push (e->simplex().get_substituted(f->node (0), nod));
  ns.push (e->simplex().get_substituted(f->node (1), nod));
  if (f->mate())
    {
      e = f->mate()->element();
      
      assert (e->opposite_node (f->mate ())
	      == right_angle_node (e));
      rm.push (e);
      ns.push (e->simplex().get_substituted(f->node (0), nod));
      ns.push (e->simplex().get_substituted(f->node (1), nod));
    }
  
  Link_array<Element> newes = replace_elements (rm, ns);
  int k = 0;
  for (int j = 0; j < rm.size(); j++)
    {
      for (int i = 0; i < 2 ; i++)
	{
	  Maubach_element * e = dynamic_cast<Maubach_element*> (rm[j]);
	  
	  e->children_[i] = dynamic_cast<Maubach_element*> (newes[k++]);
	}
    }
  return nod;
}

/*
  Return one of the bisected elements.
 */
Node*
Maubach_tree::bisect_edge (Face * f,
			   Deformation_state *def)
{
  assert (f->valid ());
  
  Element * e = f->element ();
  if (e->opposite_node (f) != right_angle_node (e))
    {
      bisect_edge (e->opposite_face (right_angle_node (e)), def);
      e = f->element();
    }
  if (f->mate ())
    {
      Element* e2 = f->mate()->element();
      if (e2->opposite_node (f->mate()) != right_angle_node (e2))
	{
	  bisect_edge (e2->opposite_face (right_angle_node (e2)), def);
	}
    }

  Node * n = simple_bisect (f, def, def);
  return n; 
}

void
refine_around2 (Maubach_tree * mb,  Real h, Vector2 x, Node_vector_func2 func, Deformation_state*def)
{
  int k = 0;
  int LIMIT = 10;
  do
    {
      Element *e = mb->locate (x, func,def);
      if (!e)
	return ;
      
      Face * edge = e ->opposite_face (right_angle_node (e));
      Real el = edge_length2 (edge->simplex(), &reference_location2, def);

      if (el < h)
	break;
      mb->bisect_edge (edge, def);
    }
  while (k++ < LIMIT);
}

Maubach_element::Maubach_element(Simplex p)
  : Element (p)
{
  children_[0] = children_[1] = 0;
}

Element*
Maubach_tree::construct_element (Simplex p)
{
  return new Maubach_element (p);
}

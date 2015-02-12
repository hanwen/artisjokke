#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "deformation-state.hh"
#include "maubach-tree.hh"
#include "mesh-feature.hh"
#include "geometry.hh"
#include "mesh-connectivity.hh"
#include "mesh-feature.hh"
#include "vector-io.hh"
#include "setting.hh"
#include "geometry2.hh"
#include "needle-inserter.hh"

void
node_print (FILE*f, Node *n, Node_vector_func2 func, Deformation_state* def)
{
  Vector2 v= (*func)(n,def);
  fprintf (f, "%lf %lf ", v(0), v(1));
}

char const * eps_header =
"%%!PS-Adobe \n"
"%%%%BoundingBox: -20 -20 120 120 \n"
"%%%%EndComments\n"
;
char const * eps_prologue =
" \n"
"/setbodycolor { 0.8 0.48 0.24 setrgbcolor } bind def \n"
"/setedgecolor { 0.1 0.1 0.1 setrgbcolor } bind def \n"
"/setbackcolor { 0.72 0.864 0.72 setrgbcolor } bind def \n"
" \n"
"100 100 scale \n"
"0 setlinecap \n"
"1 setlinejoin \n"
"/enternode { 0.0005 0 360 arc closepath stroke fill } bind def\n"
"/draw_refnode { 0.0005 0 360 arc closepath  fill } bind def\n"
"/draw_defnode { 0.0005 0 360 arc closepath stroke } bind def\n"
"/draw_fixnode { 0.00075 0 360 arc closepath fill } bind def\n"
"setbackcolor \n"
"newpath \n"
"-0.2 -0.2 moveto \n"
" 1.2 -0.2 lineto \n"
" 1.2 1.2 lineto  \n"
" -0.2 1.2 lineto \n"
" -0.2 -0.2 lineto \n"
" closepath  \n"
"fill \n"
" \n"
"-0.2 -0.2  1.4 1.4 rectclip\n"
"%lf dup scale\n"
"0.0003 setlinewidth \n"
;

void
tri_print (FILE *f, Element* e, Node_vector_func2 func, Deformation_state* def,
	   bool edges
	   )
{
  fprintf (f, " newpath ");
  node_print (f, e->node(0), func, def );
  fprintf (f, " moveto ");
  
  node_print (f, e->node(1), func, def );
  fprintf (f, " lineto ");
  node_print (f, e->node(2), func, def );
  fprintf (f, " lineto ");

  if (edges)
    {
      fprintf (f, " closepath "
	       "gsave "
	       " setbodycolor "
	       "fill "
	       "grestore\n");
  
      fprintf (f," setedgecolor  stroke ");
    }
  else
    {
      fprintf (f,
	       " setbodycolor "
	       " closepath "
	       " fill "
	       );
    }
}

void
edge_print (FILE *f, Face* fac, Node_vector_func2 func, Deformation_state* def)
{
  fprintf (f, " newpath ");
  node_print (f, fac->node(0), func, def );
  fprintf (f, " moveto ");
  
  node_print (f, fac->node(1), func, def );
  fprintf (f, " lineto ");
  fprintf (f," setedgecolor  stroke ");
}


void
mesh2d_dot_print (Mesh_connectivity* top,
		char const * fn,
		Node_vector_func2 func,
		Deformation_state * def,
		Real scale)
{
  FILE * f = xfopen (fn, "w");
  fprintf(f,eps_header);
  write_invocation_header (f, "%");
  fprintf(f,eps_prologue, scale);		 

  fprintf (f, "setedgecolor\n");
  Link_array<Node> const *nodarr = top->node_array();
  for (int i = 0; i< nodarr->size(); i++)
    {
      node_print (f, nodarr->elem (i), func, def);
      fprintf (f, " enternode\n");
    }
  fclose (f);
}

void
print_needle (FILE * f,
	      Mesh_connectivity* top,
	      Deformation_state* def,
	      Vector2 t, Vector2 h)
{
  Array<Element_intersection> eis
    = track_line_through_deformed_mesh (dynamic_cast<Maubach_tree*> (top),
					def, t, h);
  
  for (int i = 0; i < eis.size (); i++)
    {
      Element_intersection ei = eis[i];
      if (!ei.entry_face_ || ! ei.exit_face_)
	continue;
      
      Vector2 p1 = reference_location2 (ei.entry_face_->node(0), def);
      Vector2 p2 = reference_location2 (ei.entry_face_->node(1), def);
      Vector2 p = ei.entry_param_ * p1 + (1-ei.entry_param_) * p2;
      
      Vector2 q1 = reference_location2 (ei.exit_face_->node(0), def);
      Vector2 q2 = reference_location2 (ei.exit_face_->node(1), def);
      Vector2 q = ei.exit_param_ * q1 + (1-ei.exit_param_) * q2;

      fprintf (f, " %lf %lf moveto %lf %lf lineto stroke %% needle \n ",
	       p(0), p(1),
	       q(0), q(1));
    }
}

void
mesh2d_print (Mesh_connectivity* top,
	    char const * fn,
	    Node_vector_func2 func,
	    Deformation_state * def,
	    Needle_inserter * ni,
	    Real scale)
{
  FILE * f = xfopen (fn, "w");
  assert (f);
  
  fprintf(f,eps_header);
  write_invocation_header (f, "%");
  fprintf(f,eps_prologue, scale);		 


  set<Element*> *elts = top->elements();
  for (iterof (i, *elts); i != elts->end(); i++)
    {
      tri_print (f, *i, func, def, 1);
    }

  Link_array<Node> const*nodarr = top->node_array();
  for (int i = 0; i < nodarr->size(); i ++)
    {
      Node * n = nodarr->elem (i);

      if (def->constraints_.is_fixed (n))
	{
	  Vector2 loc = (*func) (n,def);
	  fprintf (f, "%lf %lf draw_fixnode\n", loc(X_AXIS), loc (Y_AXIS));
	}
    }

  if (ni)
    {
      Vector2 tip = ni->tip();
      Vector2 handle = ni->handle(); 
      
      if (func == deformed_location2)
	{
	  fprintf (f, "%lf %lf moveto %lf %lf lineto stroke ",
		   tip(0), tip(1),
		   handle (0), handle(1));
	}
      else if (func == &reference_location2)
	{
	  print_needle (f, top, def, tip, handle);
	}
  
    }      
  fclose (f);
}


void
mesh2d_print_boundary (Mesh_connectivity* top,
	    char const * fn,
	    Node_vector_func2 func,
	    Deformation_state * def,
	    Needle_inserter * ni,
	    Real scale)
{
  FILE * f = xfopen (fn, "w");
  assert (f);
  
  fprintf(f,eps_header);
  write_invocation_header (f, "%");
  fprintf(f,eps_prologue, scale);		 


  set<Element*> *elts = top->elements();
  for (iterof (i, *elts); i != elts->end(); i++)
    {
      tri_print (f, *i, func, def, false);
    }

  set<Face*> bdry = top->boundary ();
  for (iterof (i, bdry); i !=bdry.end(); i++)
    {
      edge_print (f, *i, func, def); 
    }
  
  if (ni)
    {
      Vector2 tip = ni->tip();
      Vector2 handle = ni->handle(); 
      
      if (func == deformed_location2)
	{
	  fprintf (f, "%lf %lf moveto %lf %lf lineto stroke ",
		   tip(0), tip(1),
		   handle (0), handle(1));
	}
      else if (func == &reference_location2)
	{
	  print_needle (f, top, def, tip, handle);
	}
  
    }      
  fclose (f);
}

void
mesh2d_view (Mesh_connectivity* t, Node_vector_func2 func, Deformation_state * def)
{
  mesh2d_print (t, "mesh.eps", func, def, 0, 5);
  system ("gs -g200x200 -q  mesh.eps");
}


/*
  This is clumsy: both def and ref are point based, but we only have a
  triangulation for def. It would be better to dump both deformation
  fields, and use the R triangle mesh facility to construct a Delaunay
  mesh of both fields, and use that to construct an interpolation.
 */
void
plot_reference_state_distance (Maubach_tree * tree,
			       Deformation_state *def,
			       Deformation_state *ref,
			       char const *fn)
{
#ifdef PS_ERROR
  /*
    Large files, R-plots look nicer too. 
  */ 
  char psfn [100];
  strcpy (psfn, fn);
  strcat (psfn, ".eps");
  FILE *ps = xfopen (psfn, "w");
  
  Real scale = 5;
  fprintf(ps,eps_header);
  write_invocation_header (ps, "%");
  fprintf(ps,eps_prologue, scale);		 
  fprintf (ps, " setedgecolor ");  
#endif


  char fnx[1000],fny[1000];
  strcpy (fnx,fn);
  strcpy (fny,fn);
  strcat (fnx, ".xerror");
  strcat (fny, ".yerror");  

  bool relocate = get_bool_setting ("relocate-nodes");
  FILE * f =xfopen (fn, "w");
  FILE * f_x =xfopen (fnx, "w");
  FILE * f_y =xfopen (fny, "w");    


  write_invocation_header (f, "#"); 

  Element * random_e = * tree->elements()->begin();
  
  int dim = def->spatial_dimension ();
  int nod_count = ref->dimension ()/ dim ;
  for (int i = 0; i < nod_count; i ++)
    {
      Vector2 z(ref->reference_locations_[dim *i],
		ref->reference_locations_[dim *i+ 1]);
      
      Vector2 u_ref(ref->displacements_[dim *i],
		    ref->displacements_[dim *i + 1]);

      Vector2 p_ref = u_ref + z;
      
      Element * e = relocate
	? locate_element_containing_point2 (tree, random_e, z, &reference_location2, def)
	: tree->locate (z, &reference_location2, def);


      if (!e)
	{
	  printf ("Can't find element; continuing.");
	  continue;
	}
      Vector2 points[dim + 1];
      for (int j = dim+ 1; j--;)
	points[j] = reference_location2 (e->node(j), def);

      Real coef[dim + 1];      
      barycentric_coordinates2 (points, coef, z);

      Vector2 p;
      for (int j = dim+1 ; j--;)
	p += coef[j]* deformed_location2(e->node(j), def);

      /*
	Should plot error vector field.
       */
      Vector2 errvec = p - p_ref;
      fprintf (f, "%lf %lf %le\n",  z(0), z(1), errvec.length ());
      fprintf (f_x, "%lf %lf %le\n",  z(0), z(1), errvec (X_AXIS));
      fprintf (f_y, "%lf %lf %le\n",  z(0), z(1), errvec (Y_AXIS)); 
      
#ifdef PS_ERROR
      fprintf (ps, "%lf %lf draw_refnode \n", p_ref(0), p_ref(1));
      fprintf (ps, "%lf %lf draw_defnode \n", p(0), p(1));
#endif
    }

#ifdef PS_ERROR
  fclose (ps);
#endif
  
  fclose (f);
  fclose (f_x);
  fclose (f_y);
}


#ifdef OPENGL

#include <stdio.h>

#include "opengl.hh"
#include "node.hh"
#include "artisjokke-drawer.hh"
#include "mesh-connectivity.hh"
#include "deformation-state.hh"
#include "geometry.hh"
#include "mesh-visualize3.hh"

void opengl_visualize_unclipped_mesh_topology(Artisjokke_drawer* draw, Node_vector_func3,
					      Deformation_state*);
void opengl_visualize_fixed_nodes (Artisjokke_drawer* draw, Deformation_state *def);

void opengl_visualize_nodes3  (Artisjokke_drawer *draw, Node_vector_func3 func);
#define opaque_black_col black_col

float edge_col[] = {
  0.6, 0.6, 0.6, 1.0
};


void
get_clip (Vector3 * n, Real * off, Artisjokke_drawer *draw)
{
  *n= Vector3  (1,0,0);

  float om[4][4];

  build_rotmatrix (om, draw->clip_drag_.quat_);
  Matrix3 m;
  Matrix3::from_opengl_matrixf (&m, &om[0][0]);

  *n = m* (*n) ;
  *off  = draw->clip_drag_.scale_(0);
}

void
enter_vertex_location (Node_vector_func3 func, Deformation_state * def, Node *nod)
{
  Vector3 v = (*func) (nod, def);
  enter_vertex (v);
}
		       
void
enter_edge_vertices (Node_vector_func3 func, Deformation_state*def, Simplex ss)
{
  enter_vertex_location (func, def, ss.node(1));
  enter_vertex_location (func, def, ss.node(0));
}


void
opengl_draw_face (Artisjokke_drawer*draw, Face*t, bool solid,
		  Node_vector_func3 func, 
		  Deformation_state*def)
{
  bool p = true;
  
  int e   = p ? 2 :  0;
  int b   = p ? 0 :  2;
  int dir = p ? 1 : -1;


  if (solid)
    {
      draw->set_color (body_col);
      glBegin(GL_TRIANGLES);
      
      Vector3 n = simplex_normal3 (t->simplex(), func, def);

      enter_normal (n.elts_);

      for (int j =b; j <= e; j = j + dir)
	{
	  enter_vertex_location (func, def, t->node (j));
	}
      glEnd();
    }

  if (draw->edges_)
    {
      glLineWidth(2.0 * draw->line_weight_);
      draw->set_color (edge_col); // opaque_
      glBegin (GL_LINES);
      for (int  k = 3; k--;)
	{
	  Simplex ss = t->simplex().get_subset(k);
	  if (ss.parity())
	    enter_edge_vertices (func, def, ss);
	}
      glEnd();
    }
}

#if 0 
void
opengl_visualize_clipped_mesh_topology (Artisjokke_drawer *draw, Deformation_state * def)
{
  Mesh_connectivity * mesh = draw->mesh_;

  assert  (mesh->dimension () == 3);
  
  set<Element*> * tets (mesh->elements());
  iterof (i, *tets);

  Vector3 clip_normal;
  Real clip_off = 0.0;


  get_clip (&clip_normal, &clip_off, draw);
  for (i = tets->begin(); i != tets->end(); i++)
    {
      Element *tet = *i;
      
      Vector3 c = simplex_centroid3 (tet->simplex());
      
      if (clip_normal *c < clip_off)
	continue;
      
      for (int j = 0; j < 4; j++)
	{
	  Face *t = tet->face (j);
	  if (!t->mate ())
	    opengl_draw_face (draw, t, draw->solid_,def);
	  else
	    {
	      Vector3 c = simplex_centroid (t->mate ()->element()->simplex());
	      if (clip_normal * c < clip_off)
		opengl_draw_face (draw, t, draw->solid_,def );
	    }
	}

      if (draw->inside_)
	{
	  draw->set_color (edge_col);
	  glLineWidth(1.0 * draw->line_weight_);
	  glBegin (GL_LINES);
	  for (int j = 0; j < 4 ; j++)
	    for (int k = 3; k--;)
	      {
		Simplex ss = tet->face (j)->simplex().get_subset(k);

		if (ss.parity ())
		  enter_edge_vertices (def, ss);
	      }
	  glEnd();
	  glBegin (GL_POINTS);
	  glPointSize (10 * draw->point_weight_);
	  for (int j = 0; j < 4 ; j++)
	    if (def->constraints_.is_fixed (tet->node(j)))
	      enter_vertex_location (def, tet->node(j));
	  glEnd();
	}
    }
}
#endif

void
opengl_visualize_mesh_connectivity3 (Artisjokke_drawer* draw)
{
  Deformation_state * def = draw->deformation_;

  if (draw->nodes_)
    {
      opengl_visualize_nodes3 (draw,  deformed_location3);
      if (draw->original_)
	opengl_visualize_nodes3 (draw,  reference_location3);
    }

  opengl_visualize_unclipped_mesh_topology (draw, &deformed_location3, def);

  
#if 0
  if (draw->show_degenerate_)
    {
      body_col = degenerate_body_col;
      opengl_visualize_degenerate_elements (draw, def);
    }
  if(draw->clip_tets_)
    {
      body_col = transparent_body_col;
      opengl_visualize_clipped_mesh_topology (draw, def);
    }
  else
    {
      if (!draw->show_degenerate_)
	{
	  body_col = opaque_body_col;
	}
      else
	body_col = transparent_body_col;
      
      opengl_visualize_unclipped_mesh_topology (draw, def);
    }
#endif
}

void
opengl_visualize_nodes3 (Artisjokke_drawer *draw, Node_vector_func3 func)
{
  draw->set_color (black_col);

  Link_array<Node> const *nodarr = draw->mesh_->node_array();
  Deformation_state * def = draw->deformation_;

  glPointSize (1.0 * draw->point_weight_);
  glBegin (GL_POINTS);
  for (int i = 0; i< nodarr->size(); i++)
    {
      Node * n = nodarr->elem (i);
      if (!def->constraints_.is_fixed (n))
	{
	  Vector3 v ((*func) (n, def));
	  enter_vertex (v.elts_);
	}
    }
  glEnd ();
  
  glPointSize (3.0 * draw->point_weight_);
  glBegin (GL_POINTS);
  for (int i = 0; i< nodarr->size(); i++)
    {
      Node * n = nodarr->elem (i);
      if (def->constraints_.is_fixed (n))
	{
	  Vector3 v ((*func) (n, def));
	  enter_vertex (v.elts_);
	}
    }
  glEnd ();
}

void
opengl_visualize_cached_edges (Artisjokke_drawer* draw,
			       Node_vector_func3 func,
			       Deformation_state *def)
{
  glLineWidth(1.0 * draw->line_weight_);
  draw->set_color (edge_col);
  
  glBegin (GL_LINES);
  for (int i = draw->edge_array_->size( );  i--; )
    {
      Simplex p (draw->edge_array_->elem (i));

      Vector3 v1 ((*func) (p.node (0), def));
      Vector3 v2 ((*func) (p.node (1), def));

      enter_vertex (v1.elts_);
      enter_vertex (v2.elts_);
    }
  glEnd ();
}


void
opengl_visualize_unclipped_mesh_topology(Artisjokke_drawer* draw,
					 Node_vector_func3 func,
					 Deformation_state * def)
{
  set<Element*> * tets (draw->mesh_->elements());
  iterof (i, *tets);

  if (draw->edges_ && !draw->solid_ && draw->mesh_->dimension() == 3)
    {
      opengl_visualize_cached_edges (draw, func, def);
    }
  else
    {
      if (draw->edges_ || draw->solid_)
	{
	  for (i = tets->begin(); i != tets->end(); i++)
	    {
	      Element *tet = *i;

	      for (int j = 0; j < 4 ; j++)
		{
		  Face *t = tet->face (j);
		  if (!t->mate ())
		    {
		      opengl_draw_face (draw, t, draw->solid_,
					func,
					def);
		    }
		}
	    }
	}
  
      if (draw->edges_ && !draw->boundary_)
	{
	  if (draw->original_ && !def)
	    draw->set_color (original_col);
	  else
	    draw->set_color (edge_col);
      
	  glLineWidth(1.0 * draw->line_weight_);
	  glBegin(GL_LINES);
	  {
	    i = tets->begin ();
	    for (; i != tets->end (); i++)
	      {
		Element *tet = *i;
		for (int j = 0; j < 4 ; j++)
		  for (int k = 3; k--;)
		    {
		      Simplex ss = tet->face (j)->simplex().get_subset(k);
		      enter_edge_vertices (func, def, ss);
		    }
	      }
	  }
	  glEnd ();

	}
    }
}
#endif

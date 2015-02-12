#ifdef OPENGL

#include "opengl.hh"
#include "node.hh"
#include "artisjokke-drawer.hh"
#include "mesh-connectivity.hh"
#include "deformation-state.hh"

float body_col[] = {.8, 0.48, 0.24, 1};
float original_col [] = {0.4,0.4,0.4, 0};
float black_col [] = {0,0,0,1};
float force_col [] = { 1, 0, 0, 1};

void
opengl_visualize_solid_tris (set<Element*> *tris,
			     float *col,
			     Node_vector_func2 func,
			     Artisjokke_drawer *draw
			     )
{
  Deformation_state *def =draw->deformation_;
  draw->set_color (col);

  glBegin(GL_TRIANGLES);
  for (iterof (i,*tris); i != tris->end(); i++)
    {
      int e = 2;
      int b = 0;
      int dir = 1;

      Vector3 n (0,0,- 1);
      enter_normal (n.elts_);

      for (int j =b; j <= e; j = j + dir)
	{
	  Node *nd = (*i)->node (j);
	  Vector3 v ((*func) (nd, def));

	  enter_vertex (v.elts_);
	}
    }
  glEnd();
}

void
opengl_visualize_wire_tris (set<Element*>*tris,
			    float* col,
			    Node_vector_func2 func,
			    Artisjokke_drawer *draw)
{
  Deformation_state *def =draw->deformation_;

  draw->set_color (col);
  glLineWidth (1.0 * draw->line_weight_);
  glPointSize (6.0 * draw->point_weight_);
  for (iterof (i,*tris); i != tris->end (); i++)
    {
      Element *t = *i;
	  
      for (int j = 0; j<  3 ; j++)
	{
	  Face * e = t->face (j);
	  if (e->mate ()
	      && e->parity())
	    continue;

	  if (draw->boundary_ && e->mate ())
	    continue;
	      
	  Node *nn =t->node ((j+1)%3);	  
	  Node *nd =t->node ((j+2)%3);

	  Vector3 n1 ((*func) (nn, def));
	  Vector3 n2 ((*func) (nd, def));
	      
	  if (draw->deformation_->constraints_.is_fixed (nn))
	    {
	      glBegin (GL_POINTS);
	      enter_vertex (n1.elts_);
	      glEnd();
	    }
	        
	  if (draw->deformation_->constraints_.is_fixed (nd))
	    {
	      glBegin (GL_POINTS);
	      enter_vertex (n2.elts_);
	      glEnd();
	    }
	      
	  glBegin(GL_LINES);
	  enter_vertex (n1.elts_);
	  enter_vertex (n2.elts_);
	  glEnd ();
	}
    }
}

void
opengl_visualize_nodes2 (Artisjokke_drawer *draw)
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
	  Vector3 v (deformed_location2 (n, def));
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
	  Vector3 v (deformed_location2 (n, def));
	  enter_vertex (v.elts_);
	}
    }
  glEnd ();
}

void
opengl_visualize_mesh_connectivity2 (Artisjokke_drawer* draw)
{
  set<Element*> *tris (draw->mesh_->elements());
  if (draw->nodes_)
    {
      opengl_visualize_nodes2 (draw);
    }
  
  if (draw->solid_)
    {
      opengl_visualize_solid_tris (tris, body_col, &deformed_location2, draw);
    }
  
  if (draw->edges_)
    {
      opengl_visualize_wire_tris (tris, black_col, &deformed_location2, draw);
    }
  if (draw->original_)
    {
      opengl_visualize_wire_tris (tris, original_col, &reference_location2, draw);
    }
}

#endif

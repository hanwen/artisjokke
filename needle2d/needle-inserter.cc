#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "main.hh"
#include "mesh-feature.hh"
#include "deformation-state.hh"
#include "needle-inserter.hh"
#include "maubach-tree.hh"
#include "geometry2.hh"
#include "big-vector.hh"
#include "setting.hh"
#include "debug.hh"
#include "vector-io.hh"
#include "friction.hh"

#define needle_assert(a) do { if (!a) { fprintf(stderr, "NEEDLE ASSERT FAILED: %s\n", #a); suicide(); return ; }}  while (0)


Link_array<Element>
find_node_star (Face *entry, Node * n)
{
  Face * f = entry;
  Link_array<Element> es;
  while (f)
    {
      Element* e = f->element();
      es.push (e);
      Node *nod = f->simplex().get_subset (n).node(0);
      Face *exit = e->opposite_face (nod);

      f = exit->mate ();
      if (f == entry)
	f = 0;
    }
  
  return es;
}


bool
Needle_inserter::live()const
{
  return live_;
}

Needle_inserter::Needle_inserter(Maubach_tree * tr, Deformation_state * st)
  : Deformation_hook (st)
{
  rearrange_count_ = 0;
  mesh_ = tr;
  live_ = true;
  last_handle_ = Vector2 (-0.1, 0.05);
  last_tip_ = Vector2 (-0.05, 0.05);  
  last_dir_ = Vector2 (1,0);

  
  /*
    1.01 is to avoid extra refinement due rounding,

    sqrt(2), is to control the inter needle-node spacing.
  */
  edge_length_threshold_ = get_number_setting ("refinement-h")
    * sqrt (2) * 1.01;
  state_ = OUTSIDE;
  friction_ = new Friction_definition ();
  auto_insert_ = 0;
  relocate_ = get_bool_setting ("relocate-nodes");
  dynamic_factor_ = get_number_setting ("dynamic-friction-factor");
}

void
Needle_inserter::suicide ()
{
  live_ = false;
  deformation ()->suicide();
}

void
Needle_inserter::refine_around_needle_tip (Vector2 handle , Vector2 tip)
{
  Vector2 dir = (handle-tip).normalized ();
  
  /*
    Make sure that there are enough elements.
   */
  Real tip_param = infty;
  if (needle_.size())
    {
      tip_param  = ((deformed_location2 (needle_.top(), deformation ()) - tip)*dir);
    }

  /*
    init  to length_threshold_  to refine in advance of the needle.
   */
  Real param = - edge_length_threshold_;
  
  while (1)
    {
      if (param > tip_param)
	break ;
      
      Vector2 l = tip  + dir * param;
      
      Element* e = mesh_->locate (l, &deformed_location2, deformation ());
      
      if (e && edge_length2 (long_edge (e)->simplex(), &reference_location2, deformation ()) > edge_length_threshold_)
	{
	  refine_around2 (mesh_, edge_length_threshold_, l, &deformed_location2, deformation ());

	  /*
	    This refines in a region around the needle.  Have not
	     verified if this is truly necessary.
	   */
#if 1
	  refine_around2 (mesh_, edge_length_threshold_,
			 l - rotate90 (dir)   * edge_length_threshold_ , &deformed_location2, deformation ());
	  refine_around2 (mesh_, edge_length_threshold_,
			 l + rotate90 (dir)   * edge_length_threshold_ , &deformed_location2, deformation ());
#endif

	  e = mesh_->locate (l, &deformed_location2, deformation ());
	}

      if (!e && ! needle_.size())
	{
	  /*
	    We're outside.
	   */
	  break ; 
	}

      param += edge_length_threshold_;
    }

  Vector2 l = tip;
  Element* e = mesh_->locate (l,  &deformed_location2, deformation ());


  /*
    Duh. this can hang if elements are inverted
   */
  Link_array<Node> add_nodes;
  Face * entry =0;

  /*
    Walk from tip to the first node that's already on the needle.
   */
  while (e)
    {
      Face * exit = 0;
      if (entry
	  && needle_.size()
	  && e->opposite_node (entry) == needle_.top ())
	break;
      
      Real lambda,mu ;
      Node * newnod  = 0; 
      for (int i = 0 ; !newnod && i < 3 ; i++)
	if (e->face (i) != entry
	    && intersect_segments2 (0, &lambda, &mu, tip, handle,
				   deformed_location2 (e->face(i)->node(0), deformation ()),
				   deformed_location2 (e->face(i)->node(1), deformation ())))
	  {
	    exit = e->face (i);
	    newnod  = exit->node((mu > 0.5 ) ? 0 : 1);
	  }

      if (needle_.find_l (newnod))
	{
	  break ; 
	}
      else if (newnod)
	add_nodes.push (newnod);

      if (!exit)
	{
	  log_message ("Whugh, we've been inverted.\n");

	  
	  Vector2 c = simplex_centroid2 (e->simplex(), &deformed_location2, deformation ());
	  Element *new_e  =  e;
	  while (new_e == e)
	    {
	      c += edge_length_threshold_ * dir;
	      new_e = mesh_->locate (c,  &deformed_location2, deformation ());
	    }
	  e = new_e;
	  entry = 0;
	}
      else
	{
	  entry = (exit->mate ());
	  e = (entry) ? entry->element() : 0;
	}
    }
  
  add_nodes.default_sort();
  add_nodes.uniq();
  if (add_nodes.size( ))
    {
      log_message ("Adding %d nodes\n", add_nodes.size ());
    }

  while (add_nodes.size())
    {
      Node *n = add_nodes.pop();

      assert (!needle_.find_l (n));
      Vector2 x = deformed_location2 (n, deformation ());
      Vector2 proj_loc = project_point_on_line2 (x, tip, handle);

      
      Vector2 param((x - tip) * dir, (x- proj_loc) * rotate90 (dir));

      /*
	We can't break out of the loop: if there are inverted elts,
	we're hosed.
       */
      
      if (param (X_AXIS) < 0)
	continue;
      

      /*
	Relocation:
       */
      if (relocate_)
	{
	  Element *e = mesh_->locate (proj_loc, &deformed_location2,
				      deformation ());

	  /*
	    Ugh, relocation sucks.
	   */
	  if (!e)
	    e = mesh_->locate (reference_location2 (n, deformation ()), &reference_location2, deformation ());
	  
	  Vector2 y = transform_back2 (e->simplex (), proj_loc, deformation ());
	  assert (e);

	  Face * f = e->incident_face (n,0);
	  Link_array<Element> star = find_node_star (f,n);
	  if (f->mate())
	    star.concat (find_node_star (f->mate(),n));

	  for (int i= star.size(); i--;)
	    deformation ()->add_changed_element (e);
	  deformation ()->set_reference_location (n, y);
	  param(Y_AXIS) = 0.0;
	}
      
      deformation ()->constraints_.fix_node  (n, true);
      needle_.push (n);
      needle_params_.push (param);
    }

  /*
    TODO: add nodes to needle_
  */
}

void
Needle_inserter::drag_fixed_nodes (Vector2 h, Vector2 t)
{
  Vector2 dir = (h - t).normalized(); 
  for (int i = 0; i < needle_.size (); i++)
    {
      Node * n = needle_[i];
      Vector2 param = needle_params_[i];
      if (!deformation ()->constraints_.is_fixed (n))
	{
	  Vector2 x = deformed_location2 (n, deformation ());
	  param(X_AXIS) = (x - t)*dir ;
	  deformation ()->constraints_.change_linear_movement_constraint (n, dir, 0);
	}
      
      deformation ()->set_node_deformed_location (n, t + dir * param(X_AXIS) +
						param (Y_AXIS) * rotate90 (dir)
						);
      
    }
}


void
Needle_inserter::rearrange_boundary_conditions ()
{
  if (!needle_.size())
    return;


  rearrange_count_ ++;
  Vector2 dir = last_dir_;

  /*
   */
  
  Real *elforce = deformation ()->elastic_force_.unsafe_access_array();

  Array<Real> friction_thresholds;
  for (int i =  needle_.size(); i--;)
    {
      Vector2 x =  deformed_location2 (needle_[i], deformation ());

      friction_thresholds.push ((x - last_tip_)*dir);
    }

  friction_thresholds = compute_friction_forces (friction_thresholds,
						 friction_);
  friction_thresholds.reverse();
  big_vector_scale (friction_thresholds.unsafe_access_array(),
		    1.0/dynamic_factor_,
		    friction_thresholds.unsafe_access_array(),
		    friction_thresholds.size());
    
  Real total_elastic_force =0.0; 
  for (int i = 0; i < needle_.size(); i++)
    {
      Node * n = needle_[i];
      bool fix = deformation ()->constraints_.is_fixed (n);

      Vector2 elastic_force = extract_vector2 (n, elforce);
      Real directed_elastic_force  = elastic_force * dir;

      if (fabs (directed_elastic_force) > friction_thresholds[i])
	{
	  if (fix)
	    deformation ()->constraints_.change_linear_movement_constraint (n, dir, 0);

	  deformation ()->set_force (n, 
				   -dynamic_factor_* sign (directed_elastic_force) * dir * friction_thresholds[i]);

	  total_elastic_force += sign (directed_elastic_force)  * friction_thresholds[i];
	}
      else
	{
	  if (!fix)
	    deformation ()->constraints_.fix_node (n, true);

	  Vector2 p = needle_params_[i];

	  needle_params_[i]
	    = Vector2( (deformed_location2 (n, deformation ()) - last_tip_)*  last_dir_, p(Y_AXIS));

	  total_elastic_force += directed_elastic_force;
	}
    }

  // printf ("Elastic force: %lf\n", total_elastic_force);
  int popped = 0;

  for (int i = needle_.size(); i--;)
    {
      Node * n = needle_[i];
      Vector2 x = deformed_location2 (n, deformation ());

      /*
	Should use smaller threshold? 
       */
      if ((x - last_tip_)  * last_dir_ < 0)
	{
	  if (debug_out)
	    {
	      log_message ("nod pos/tip pos:");
	      x.print();
	      last_tip_.print();
	      last_dir_.print();
	    }
	  deformation ()->constraints_.change_linear_movement_constraint (0, Vector2(0,0), n);
	  deformation ()->constraints_.fix_node (n, false);
	  needle_.del (i);
	  needle_params_.del(i);
	  popped ++;  
	}
    }
  
  if (popped)
    log_message ("Popped %d\n", popped);

  if (!needle_.size())
    {
      state_ = OUTSIDE;
      log_message ("Leaving mesh\n");
    }
}

void
Needle_inserter::print_needle ()
{
  log_message ("Needle: "); 
  for (int i = needle_.size(); i--;)
    {
      bool f = deformation ()->constraints_.is_fixed (needle_[i]);
      bool l = deformation ()->constraints_.has_linear_movement_constraint (needle_[i]);
      printf ("%d:  Fix %d, linear %d\n", i, f,l);
    }
}

void
Needle_inserter::signal_converged ()
{
  if (!live_)
    return ;

  rearrange_boundary_conditions ();
  refine_around_needle_tip (last_handle_, last_tip_ );

  if (auto_insert_)
    {
      auto_needle_insert ();
    }
}

void
Needle_inserter::move_needle (Vector2 h, Vector2 t)
{
  if (!live_)
    return ;

  last_dir_ =(h-t).normalized();
  last_tip_ = t ;
  last_handle_ =h; 

  if (state_ == OUTSIDE)
    {
      Element * e  = mesh_->locate (t,  &deformed_location2, deformation ());

      if (!e)
	return ;


      if (e)
	{
	  state_ = INSIDE;
	}
    }
  
  if (state_ == INSIDE)
    {
      drag_fixed_nodes (h, t);
      refine_around_needle_tip (h, t); 
    }
}

void
Needle_inserter::dump_needle_forces (const char *fn)
{
  FILE * f = xfopen (fn, "w");

  write_invocation_header (f, "# ");
#if 0
  fprintf (f,
	   "# FORMAT:\n"
	   "# NODE-COUNT\n"
	   "# NODE-COUNT*2 lines, with REF-LOC, EXTERNAL-FORCE\n");
  
  fprintf (f, "%d\n", needle_.size ());

#endif
  
  fprintf (f,
	   "# FORMAT:\n"
	   "# node-count lines of REF-LOC{X,Y} DEF-LOC{X,Y} EXT-FORCE{X,Y} FIXED \n");

  Real total_force  = 0.0;
  for (int i =0; i < needle_.size (); i++)
    {
      Node* n = needle_[i];

      Vector2 loc = reference_location2 (n, deformation ());
      Vector2 ext_force =extract_vector2 (n, deformation ()->external_force_.access_array ());
      Vector2 reaction_force =extract_vector2 (n, deformation ()->elastic_force_.access_array ());
      deformation ()->constraints_.node_reaction (reaction_force.elts_, n, reaction_force.elts_);
      deformation ()->constraints_.apply_to_node_movement (ext_force.elts_, n, ext_force.elts_);
      

      bool fix = deformation ()->constraints_.is_fixed (n);
      Vector2 force = ext_force - reaction_force;

      Vector2 defloc = deformed_location2 (n, deformation ());
      fprintf (f, "%.15g %.15g %.15g %.15g %.15g %.15g %d\n",
	       loc(0), loc(1),
	       defloc (0), defloc (1),
	       force (0), force (1),
	       ! (!fix)
	       
	       );
      
#if 0
      write_vector2 (f, loc);
      write_vector2 (f, force);
#endif
    }

  
  fclose (f);
  fprintf (stderr, "\n");
}
Element_intersection::Element_intersection()
{
  exit_param_ = entry_param_ =-1.0;
  elt_ = 0;
  entry_face_ = exit_face_  = 0;
}


Array<Element_intersection>
track_line_through_deformed_mesh (Maubach_tree * mesh,
				  Deformation_state * def,
				  Vector2 tip,
				  Vector2 handle
				  )
{
  Array<Element_intersection> retval;  
  Vector2 l = tip;

  Vector2 dir  = handle - tip;
  
  Element* e = mesh->locate (l,  &deformed_location2, def);
  Real edge_length_threshold = get_number_setting ("refinement-h");
  
  Face * entry =0;

  /*
    Walk from tip to the first node that's already on the needle.
   */
  while (e)
    {
      Face * exit = 0;
      
      Real lambda,mu ;

      Element_intersection einter;
      einter.elt_ = e;

      if (retval.size ())
	{
	  einter.entry_face_ = entry;

	  /*
	    don't reverse. The sequence is not inverted, only the
	    parity bit.
	  */
	  einter.entry_param_ = retval.top().exit_param_;
	}
      
      for (int i = 0 ; !exit && i < 3 ; i++)
	if (e->face (i) != entry
	    && intersect_segments2 (0, &lambda, &mu, tip, handle,
				   deformed_location2 (e->face(i)->node(0), def),
				   deformed_location2 (e->face(i)->node(1), def)))
	  {
	    exit = e->face (i);
	    einter.exit_face_ = exit;
	    einter.exit_param_ = mu;
	  }


      if (!exit)
	{
	  log_message ("Whugh, we've been inverted.\n");

	  
	  Vector2 c = simplex_centroid2 (e->simplex(), &deformed_location2, def);
	  Element *new_e  =  e;
	  while (new_e == e)
	    {
	      c += edge_length_threshold * dir;
	      new_e = mesh->locate (c,  &deformed_location2, def);
	    }
	  e = new_e;
	  entry = 0;
	}
      else
	{
	  entry = (exit->mate ());
	  e = (entry) ? entry->element() : 0;
	}

      retval.push (einter); 
    }

  return retval;
}

Vector2
Needle_inserter::tip_reference_location () const
{
  Vector2 t= tip();
  if (state_ == INSIDE)
    {
      Element * e = mesh_->locate (t, deformed_location2, deformation ());
      if (!e)
	fprintf (stderr, "Can't locate element of tip....\n");
      else
	{
	  t = transform_back2 (e->simplex(), t, deformation ());
	}
    }

  return t;
}

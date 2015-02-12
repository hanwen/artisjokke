#include <algorithm>

#include <math.h>
#include <string.h>

#include "deformation-state.hh"
#include "setting.hh"
#include "maubach-tree3.hh"
#include "needle-inserter3.hh"
#include "geometry.hh"

using std::max;
using std::min;

/*
  The needle is characterized by the set NEEDLE_ELEMENTS_.

  All nodes that are constrained (boundary of the NEEDLE_ELEMENTS_
  subcomplex) also have a needle radius (necessary for moving the needle.)

  The needle is described by vectors (H, T), where H is the handle
  (outside o f the complex) and T is the tip (inside). All
  `measurements' are done from the tip , using DIR , the normalized
  vector pointing from T to H.

  Handling NEEDLE_ELEMENTS_ is not done particularly
  efficient. However, the total program is typically bound by the FEM
  computations, so this is probably moot.
 */
Needle_inserter3::Needle_inserter3 (Maubach_tree3* tree,
				    Deformation_state*def)
  : Deformation_hook (def)
{
  mesh_ = tree;
  state_ = OUTSIDE;
  live_  = true;
  refinement_level_ = (int) get_number_setting ("refinement-level");
  needle_radius_ = get_number_setting ("needle-radius");
  auto_insert_ = 0;
  dynamic_friction_factor_ = get_number_setting ("dynamic-friction-factor");
  rearrange_count_ = 0;
  /*
    This only makes sense for exact strain formulations.
   */
  filter_needle_rotations_ = strcmp (get_string_setting ("elasticity"), "linear");
  mesh_->add_watcher (this);
  friction_density_ = get_number_setting ("friction-density");
}
				    
Vector3
Needle_inserter3::handle () const
{
  return last_handle_ ;
}

Vector3
Needle_inserter3::tip () const
{
  return last_tip_ ;
}

void
Needle_inserter3::move_needle (Vector3 h, Vector3 t)
{
  if (state_ == OUTSIDE)
    {
      Element * e = mesh_->locate (t, &reference_location3, deformation ());

      if (e)
	{
	  state_ = INSIDE;
	}
    }

  if (state_ == INSIDE)
    {
      drag_nodes (h , t);
    }

  last_tip_ = t;
  last_handle_ = h;

  if (state_ == INSIDE)
    handle_tip_changes (h, t);
}


/*
  All faces of NEEDLE_ELEMENTS_ whose mates are not in
  NEEDLE_ELEMENTS_.

  We only include the faces that do have a mate (i.e. are not in
  boundary of the large complex).
 */
set<Face*>
Needle_inserter3::get_needle_surface ()
{
  set<Face*> needle_surf;
  
  for (iterof (i, needle_elements_); i != needle_elements_.end(); i++)
    {
      for (int j= 4; j--;)
	{
	  Face * f = (*i)->face (j);

	  if (f->mate () && 
	      !has_elt(needle_elements_, f->mate ()->element()))
	    {
	      needle_surf .insert (f);
	    }
	}
    }
  return needle_surf;
}

/*
  Return friction force / area.
  
  */
Real
Needle_inserter3::force_distribution (Real tip_distance, Real entry_param)
{
  Real density = friction_density_ /* N/m  */
    / (2 * M_PI *  0.0005);	// DiMaio's measurement.
  
  Real entry_ramp_distance= 75e-3;
  if ((entry_param -tip_distance)  < entry_ramp_distance)
    {
      density *=  (entry_param -tip_distance) / entry_ramp_distance;
    }

  /*
    TODO: include the extra cutting force at the tip. Left to the
    reader for an excercise.
  */
  
  return density;
}



/*
  Compute the friction force for triangle PLEX, given needle
  coordinates (H,T). The signed area is also stored in AREA: if the
  triangle is facing away from the centerline, it is considered to
  have positive area, otherwise it has negative area.

  This should approximate the area of the needle relatively accurate.

  Printouts in the code below indicate that there are errors of 10 to
  30 % when ref-level is small (h large) in these computed net areas
  meaning that forces (and hence, displacements) will be off by some
  10 to 30 % too. I don't quite understand why, but there's not enough
  time to figure it out.

 */
Real
Needle_inserter3::friction_force (Real * area, Simplex const &plex, Vector3 h, Vector3 t,
				  Real entry_param
				  )
{
  Vector3 x[3];  
  Vector3 proj_x[3];
  Vector3 phi_z_coordinates [3];
  for (int j = 3; j--;)
    {
      x[j] = deformed_location3 (plex.node(j), deformation ());
      proj_x[j] = project_point_on_line3 (x[j],  t, h);
    }

  Vector3 normal = simplex_normal3 (plex, &deformed_location3, deformation ());
  Vector3 centroid = (x[0] + x[1] + x[2]) / 3.0;

  Vector3 dir = h -t;
  dir.normalize ();

  /*
    beyond the tip of the needle:
   */
  if ( (centroid - t) * dir  < 0)
    {
      *area = 0;
      return 0.0;
    }
  
  Real sign = ((x[0] - proj_x[0] ) * normal > 0) ?  1 : -1;


  /*
    We could discard tris if the triangle intersects the needle line, but
    that doesn't help if h_ref << radius.
  */
  Vector3 dir1 = x[1] - proj_x[1];
  if (dir1.length () < needle_radius_ / 100)
    {
      printf ("Short distance to needle, %lf\n", dir1.length());
    }
  
  dir1.normalize ();
  Vector3 dir2 = Vector3::cross (dir , dir1);
  
  for (int j = 3; j--;)
    {
      Vector3 dir0 = x[j] - proj_x[j];
      Vector2 xy  (dir0 * dir1 , dir0 * dir2);

      if  (xy.length () < needle_radius_/ 100)
	{
	  printf ("Short distance to needle : %lf\n", xy.length ());
	  continue; 
	}
      
      Real phi = atan2(xy(0), xy(1));
      
      phi_z_coordinates[j](Y_AXIS) = (x[j] -t)* dir ;
      phi_z_coordinates[j](X_AXIS) = phi;
    }

  *area = Vector3::cross (phi_z_coordinates[0] -phi_z_coordinates[2],
			  phi_z_coordinates[1] - phi_z_coordinates[2]).length () * 0.5 * needle_radius_ *sign;
  
  return force_distribution ((centroid -t) * dir, entry_param) * (*area);
}


/*
  Change the needle position: express the locations of all needle
  nodes (including internal nodes of NEEDLE_ELEMENTS_) in coordinates
  relative to the needle, and set their positions to by changing the
  coordinate system to the new needle position.
 */
void
Needle_inserter3::drag_nodes (Vector3 h, Vector3 t)
{
  Vector3 dir = last_handle_ -last_tip_;

  dir.normalize();
  Vector3 ez (0,0,1);
  Vector3 o1 = Vector3::cross (ez, dir);
  if (o1.length () < 1e-6)
    {
      o1 = Vector3 ::cross (dir, Vector3(0,1,0));
    }

  o1.normalize();
  Vector3 o2 = Vector3::cross (dir, o1);

  map<Node*, Vector3> parameterized_locs;
  for (iterof (i, needle_elements_); i != needle_elements_.end (); i++)
    {
      for (int k = 4; k--;)
	{
	  Node * n = (*i)->node (k);
	  if (!has_key (parameterized_locs, n))
	    {
	      Vector3 x = deformed_location3 (n, deformation ()) - last_tip_;  
	      parameterized_locs[n] = Vector3(x * o1,
				       x * o2,
				       x * dir);
	    }
	}
    }

  dir = h - t;
  dir.normalize ();

  ez = Vector3 (0,0,1);
  o1 = Vector3::cross (ez, dir);
  if (o1.length () < 1e-6)
    {
      o1 = Vector3 ::cross (dir, Vector3(0,1,0));
    }

  o1.normalize();
  o2 = Vector3::cross (dir, o1);

  for (iterof (i, parameterized_locs); i != parameterized_locs.end (); i++)
    {
      Vector3 p (i->second);
      Vector3 x = p(X_AXIS) * o1+  p(Y_AXIS) * o2 + p(Z_AXIS) * dir + t;

      deformation ()->set_node_deformed_location (i->first,  x);

      if (deformation ()->constraints_.has_ortho_movement_constraint (i->first))
	{
	  Vector3 d = p(X_AXIS) * o1 + p(Y_AXIS) * o2;
	  deformation ()->constraints_.change_ortho_movement_constraint (i->first, d.elts_, 0);
	}
    }
}

/*
  ENTRY_PARAM is where the needle enters the tissue.
*/
Real
Needle_inserter3::get_entry_parameter (Vector3 h, Vector3 t) const
{
  Vector3 dir = h - t ;
  dir.normalize();


  Real maxp = 0.0; 
  for (iterof (i, needle_elements_); i != needle_elements_.end () ; i++)
    {
      Vector3 c = simplex_centroid3 ((*i)->simplex(),
				     &deformed_location3, deformation());
      maxp = max (maxp, (c-t) * dir);
      
    }
  return maxp;
}

/*
  Compute friction thresholds, and change boundary conditions on the
  basis of these thresholds.
 */
void
Needle_inserter3::rearrange_boundary_conditions ()
{
  rearrange_count_ ++;
    
  map<Node*,Real> frictions;
  set<Face*> needle_surf = get_needle_surface ();
  for (iterof (i, needle_surf); i != needle_surf.end(); i++)
    {
      for (int j = 3; j--; )
	frictions[(*i)->node(j)] = 0.0;
    }
  
  Vector3 dir = (last_handle_ - last_tip_ );
  dir.normalize ();

  Real maxp = -1e6;
  Real minp = 1e6;
  Real tipp = (last_tip_  * dir);
  for (iterof (i, frictions); i !=  frictions.end (); i++)
    {
      Vector3 dx = deformed_location3  (i->first,  deformation());
      Real p = (dx * dir) - tipp ;
      maxp = max (maxp, p);
      minp = min (minp, p);
    }

  Real area = 0.0;

  int inverted = 0;
  for (iterof (i, needle_surf); i != needle_surf.end(); i++)
    {
      Real a = 0.0;
      Real frict = friction_force (&a ,(*i)->simplex () , last_handle_, last_tip_,
				   maxp
				   );

      if (a < 0)
	inverted ++;
      area += a;
      
      for (int j = 3; j--; )
	frictions[(*i)->node(j)] += frict / 3.0;
    }

#if 1
  Real theoretical_area = (maxp - minp) * needle_radius_ * 2 * M_PI;
  printf ("Approximated area %lf, theoretical: %lf, relative diff %4.1lf %% \n",
	  area,
	  theoretical_area,
	  (theoretical_area - area) / theoretical_area * 100.0);

  printf ("Surface triangles: %d, inverted: %d (%4.1f %%)\n",
	  needle_surf.size(),
	  inverted,
	  (100.0 * inverted )/  needle_surf.size());
#endif
  
  /*
    Use friction thresholds to
    fix, or unfix nodes.

    Fixing, Unfixing happens in 2 directions.
   */
  Real *elforce = deformation ()->elastic_force_.unsafe_access_array();

  for (iterof (i,  frictions); i != frictions.end (); i++)
    {
      Node *n = i->first;
      Real dynamic_nodal_friction = i->second;
      Real static_nodal_friction = dynamic_nodal_friction / dynamic_friction_factor_;
      
      assert (has_key (needle_node_radii_, n));
      Real node_radius = needle_node_radii_[n]; 

      Vector3 elastic_force = extract_vector3 (n,  elforce);

      Vector3 x = deformed_location3 (n, deformation ());
      Vector3 proj_x = project_point_on_line3 (x, last_handle_, last_tip_);
      Vector3 radial_dir =  x - proj_x;

      /*
	rotations introduce strains 
      */
      radial_dir.normalize(); 

      Vector3 radial_force = (elastic_force  *  radial_dir) * radial_dir;
      Vector3 planar_elastic_force = elastic_force - radial_force;


      if (filter_needle_rotations_)
	{
	  Vector3 new_location  = proj_x + radial_dir * node_radius;
	  deformation ()->set_node_deformed_location (n, new_location);
	}

      
      if (fabs (static_nodal_friction) < planar_elastic_force.length ())
	{
	  deformation ()->constraints_.change_ortho_movement_constraint (n, radial_dir, 0);
	  deformation ()->set_force (n,
				   - (planar_elastic_force / planar_elastic_force.length ())
				   * dynamic_nodal_friction);
	}
      else 
	{
	  deformation ()->constraints_.fix_node (n, true);
	}

#if 0      
      printf("Friction %lf, force %lf, fixed: %d \nel force: ",
	     nodal_friction, planar_elastic_force.length(),
	     deformation ()->constraints_.is_fixed (n)
	     ) ;
      elastic_force.print();
#endif
    }

  /*
    Pop off elements slid off the needle.
   */
  Link_array<Element> erase;
  for (iterof (i, needle_elements_); i != needle_elements_.end (); i++)
    {
      bool slid= true;
      for (int j = 4; slid && j--;)
	{
	  Node * n = (*i)->node(j);
	  Vector3 x = deformed_location3 (n, deformation ());
	  slid = slid && ((x - last_tip_) * dir < 0);
	}
      if (slid)
	erase.push (*i);
    }


  
  for (int i =  erase.size(); i--;)
    needle_elements_.erase (erase[i]);

  if (!needle_elements_.size ())
    {
      state_ = OUTSIDE;
    }
}

void
Needle_inserter3::handle_tip_changes (Vector3 h, Vector3 t)
{
  refine_around_needle_tip (h, t);

  update_needle_elements ();
  accrue_needle_elements (h, t); 
}

/*
  Refine the mesh around the needle tip.

  This happens in a circle perpendicular to the needle.
  
  This assumes that the needle doesn't move a big deal between
  BC rearrangements.
*/
void
Needle_inserter3::refine_around_needle_tip (Vector3 h, Vector3 t)
{
  Vector3 dir  = h -t  ;
  dir.normalize ();



  Vector3 ez (0,0,1);
  Vector3 o1 = Vector3::cross (ez, dir);
  if (o1.length () < 1e-6)
    {
      o1 = Vector3 ::cross (dir, Vector3(0,1,0));
    }

  o1.normalize();
  Vector3 o2 = Vector3::cross (dir, o1);
  Real refinement_h = 0.1 * pow (2, - (1.0/3.0) * refinement_level_);

  int steps =  int (needle_radius_ * 2 * M_PI / refinement_h);
  steps *= 2;
  for (int i = 0; i  < steps; i++)
    {
      Real angle = (2 * M_PI  / steps) * i;

      Vector3 rad = (o1 * sin (angle)  + o2 * cos (angle));
      Vector3 x = t + rad * needle_radius_;

      /*
	We want accurate sols in the vicinity too.
       */
      Vector3 farx = t + rad * (refinement_h +  needle_radius_);
      
      refine_around3 (mesh_, refinement_level_, x, & deformed_location3, deformation ());
      refine_around3 (mesh_, refinement_level_, farx, & deformed_location3, deformation ());      
    }

  Vector3 fart = t  - refinement_h * dir;
  refine_around3 (mesh_, refinement_level_, fart, & deformed_location3, deformation ());      
  

}


bool print_dist = false;

bool
Needle_inserter3::is_needle_element (Element *e, Vector3 h, Vector3 t,  Real entry_param) const
{
  Vector3 c = simplex_centroid3 (e->simplex(), &deformed_location3, deformation ());

  Vector3 dir = h - t;

  Real r = point_to_line_segment_distance3 (c, h, t);
  Real p = dir * (c - t);

  Real simp_dist = tetrahedron_to_line_segment_distance3 (e->simplex (), &deformed_location3,
							  deformation (), h, t);

  if (print_dist)
    printf ("Rad %lf, simplex dist %lf, param %lf, ep %lf\n",
	    r, simp_dist, p, entry_param);
  return (p >= 0.0 && p <= entry_param)  && (r < needle_radius_ || simp_dist < r /2 ) ;
}


/*
  Add elements to NEEDLE_ELEMENTS_. We start from the tip, and
  recurively add neighbors.
*/
void
Needle_inserter3::accrue_needle_elements (Vector3 h, Vector3 t)
{

  Element*  tip_e = mesh_->locate (t, &deformed_location3, deformation ());
  if (!tip_e)
    {
      fprintf (stderr, "No elts around tip??\n");  
      return;
    }

  Vector3 dir  = h -t  ;
  dir.normalize ();
  
  Link_array<Element> todo;
  todo.push (tip_e);

  /*
    Don't use entry_param for selection.
   */
  Real entry_param =   1e6;
  set<Element*> done ;
  while (todo.size ())
    {
      Element * e = todo.pop ();

      bool already_done = has_elt (done, e);
      if (already_done)
	continue;

      done.insert (e);
      assert (e->valid());

      if (has_elt (needle_elements_, e))
	continue;

      if (is_needle_element (e, last_handle_, last_tip_, entry_param))
	{
	  needle_elements_.insert (e);
      
	  for (int j = 4 ; j--;)
	    {
	      Face* f  = e->face (j);
	  
	      if (!f->mate ())
		continue;
	      todo.push (f->mate()->element());
	    }
	}
    }  
  update_boundary_nodes ();
}


/*
  Find new needle surface nodes by comparing NEEDLE_ELEMENTS_ and
  NEEDLE_NODE_RADII_, and set appropriate radii.

  Remove node that are no longer part of the needle, junking the
  constraints as well.
 */
void
Needle_inserter3::update_boundary_nodes ()
{
  set<Face*> surf = get_needle_surface ();
  set<Node*> nods;
  
  for (iterof (i, surf); i != surf.end(); i++)
    {
      for (int j = (*i)->simplex().count();j--;)
	nods.insert ((*i)->node (j)); 
    }

  for (iterof (i, nods); i != nods.end (); i++)
    {
      if (!has_key (needle_node_radii_, *i))
	{
	  Vector3 y  = deformed_location3 (*i, deformation ());
	  needle_node_radii_[*i] = point_to_line_distance3 (y, last_handle_, last_tip_);
	  deformation ()->constraints_.fix_node (*i, true);
	}
    }

  
  Link_array<Node> rm;
  for (iterof (i, needle_node_radii_); i != needle_node_radii_.end ();  i++)
    {
      if (!has_elt (nods, i->first))
	{
	  rm.push (i->first);
	}
    }
  
  for (int i = rm.size () ; i--; )
    {
      needle_node_radii_.erase (rm[i]);
      deformation()->constraints_.remove_all_node_constraints (rm[i]);
    }
}


/*

*/

void
Needle_inserter3::signal_converged ()
{
  if (!live_)
    return ;

  rearrange_boundary_conditions ();
  update_boundary_nodes ();
  
  handle_tip_changes (last_handle_, last_tip_ );
  if (auto_insert_)
    {
      auto_needle_insert ();
    }
}


/*
  Put topological changes in the NEEDLE_ELEMENTS_ set: subdivision may
  render existing elements of NEEDLE_ELEMENTS_ invalid. a
 */
void
Needle_inserter3::update_needle_elements ()
{
  Vector3 dir = last_handle_ - last_tip_;
  dir.normalize ();

  Real entry_param = get_entry_parameter (last_handle_, last_tip_);
    
  if (set<Element*> *chs = get_changed_elements ())
    {
      for (iterof (i ,*chs); i != chs->end (); i++)
	{
	  Element * e = *i;
	  bool v = (e->valid ());
	  bool h = has_elt (needle_elements_, e);
	    
	  if (!v && h)
	    {
	      needle_elements_.erase (*i);

	      set<Element*> es;
	      
	      get_sub_elements (e, &es);
	      for (iterof (j, es); j != es.end (); j++)
		{
		  if (is_needle_element (*j, last_handle_, last_tip_, entry_param))
		    needle_elements_.insert (*j);
		}
	    }
	  else if (v && !h)
	    {
	      if (is_needle_element (e, last_handle_, last_tip_, entry_param))
		needle_elements_.insert (e);
      	    }
	}

      update_boundary_nodes ();
    }

}

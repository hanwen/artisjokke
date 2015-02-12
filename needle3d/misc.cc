#if 0


/*
  This doesn't work: the split routine recurses a lot , and there are
  still a lot of smallish peaks around the needle.  
 */
void
Needle_inserter3::split_ortho_to_needle ()
{
  Vector3 dir  = last_handle_ - last_tip_ ;

  set<Face*> wrong ;
  dir.normalize ();

  set<Face*> surf = get_needle_surface ();
  for (iterof (i, surf); i != surf.end (); i++)
    {
      Face * f = *i;
      Vector3 c = simplex_centroid3 (f->simplex (), &deformed_location3, deformation());

      
      Vector3 n = simplex_normal3 (f->simplex(), &deformed_location3, deformation());

      Real angle = acos (dir * n);

      if (angle < M_PI * 10.0 / 180.0
	  && (c - last_tip_) * dir > 0)
	{
	  wrong.insert (f); 
	}
    }

  int k = 0;
  for (iterof (i, wrong); i != wrong.end (); i ++)
    {
      Face * f = *i;
      
      if (f->valid ())
	{
	  mesh_->bisect_edge (f->element(), deformation());
	  k ++;
	}
    }

  if (k)
    printf ("Subdivided %d elements to combat flat needle triangles.\n", k);
  
}



Link_array<Element>
track_line_through_mesh (Mesh_connectivity * mesh,
			 Element * start,
			 Vector3 handle,
			 Vector3 tip,
			 Deformation_state * def)
{
  Vector3 dir = handle-  tip;
  dir.normalize();

  Element * e  = start;
  Link_array<Element> elts ;
  Face *entry = 0;
  while (e)
    {
      elts.push (e);
      Face * exit = 0;
      for (int i = 4;  !exit && i --; )
	{
	  Face * f =  e->face (i);

	  if (f != entry)
	    {
	      Real p ;
	      bool succ = intersect_line_simplex3 (&p,  f->simplex(),
					     &deformed_location3,
					     def, handle, tip);

	      if (succ)
		exit = f;
	    }
	}

      assert  (!exit);
      
      entry = exit->mate();
      e = entry ? entry->element() : 0;
    }

  return elts;
}

/*
  unused.
 */
void
Needle_inserter3::select_elements_around_needle (Element * e ,
						 Vector3 h, Vector3 t)
{
  Link_array<Element> elts = track_line_through_mesh (mesh_, e, h, t,
						      deformation ());
  
  for (int i = elts.size(); i--;)
    {
      Element * e = elts[i];

      if (has_elt (needle_elements_, e))
	break;

      Link_array<Element> todo;
      todo.push (e);
      while (todo.size())
	{
	  Element * k = todo.pop ();
	  if (has_elt(needle_elements_, k))
	    continue;

	  Vector3 cen = simplex_centroid3 (k->simplex(),
					   &deformed_location3, deformation ());
	  
	  if (point_to_line_segment_distance3 (cen, h, t) > needle_radius_)
	    continue;
	  
	  needle_elements_.insert (k);
	  for (int j = 4 ; j--;)
	    if (Face *f = k->face (j)->mate ())
	      {
		todo.push (f->element ());
	      }
	}
    }
}

#endif

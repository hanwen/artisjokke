#include "edge-watcher.hh"
#include "mesh-feature.hh"

Edge_watcher::Edge_watcher()
{
  changed_ = false;
}

set<Simplex> const *
Edge_watcher::edges ()
{
  process_changes();
  return &edges_;
}

Array<Simplex> *
Edge_watcher::edge_array ()
{
  process_changes();
  if (!changed_)
    return 0;

  Array<Simplex> *plexes = new Array<Simplex>;
  plexes->set_size (edges_.size());

  int j = 0;
  for (iterof (i, edges_);
       i!= edges_.end (); 
       i++)
    {
      plexes->elem_ref (j++) = *i;
    }

  changed_ =false;
  return plexes;
}

void
Edge_watcher::process_feature(Mesh_feature * f, bool old)
{
  if (f->valid () == old)
    return;

  Simplex const &p = f->simplex ();
  if (p.count () == 4)
    {
      for( int i = p.count (); i--;)
	for (int j = i ; j-- >0 ; )
	  {
	    Simplex e (p.get_subset (i). get_subset (j));
	    if (e.parity())
	      continue ;

	    if (old)
	      edges_.erase (e);
	    else
	      edges_.insert (e);
	  }
    }
  else if (p.count () == 3)
    {
      for( int i = p.count (); i--;)
	{
	  Simplex e (p.get_subset (i));
	  if (e.parity())
	    continue ;

	  if (old)
	    edges_.erase (e);
	  else
	    edges_.insert (e);
	}
    }
}

void
Edge_watcher::process_changes()
{
  set<Face*> * bdry = get_changed_boundary ();
  set<Element*> * elts = get_changed_elements ();

  if (bdry)  
    for (iterof (i,*bdry);
	 i != bdry->end (); i++)
      process_feature (*i, true);
  
  if (elts)  
    for (iterof (i,*elts);
	 i != elts->end (); i++)
      process_feature (*i, true);

  if (bdry)  
    for (iterof (i,*bdry);
	 i != bdry->end (); i++)
      process_feature (*i, false);
  
  if (elts)  
    for (iterof (i,*elts);
	 i != elts->end (); i++)
      process_feature (*i, false);

  if (elts || bdry)
    {
      changed_ = true;
      delete elts;
      delete bdry;
    }
}


Edge_watcher::~Edge_watcher()
{
}

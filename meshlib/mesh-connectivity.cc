#include <stdio.h>
#include <math.h>

#include "mesh-feature.hh"
#include "mesh-connectivity.hh"
#include "node.hh"
#include "mesh-connectivity-watcher.hh"

// #define PARANOIA



Mesh_connectivity::Mesh_connectivity(int spatial_dimension)
{
  dimension_ = spatial_dimension;
}

set<Face*> * 
Mesh_connectivity::faces () const
{
  set<Face*> *es  = new set<Face*>;

  for (iterof (i,elements_); i != elements_.end (); i++)
    {
      for (int j = dimension_ + 1; j--;)
	es->insert ((*i)->face(j));
    }
  return es;  
}

set<Face*> 
Mesh_connectivity::boundary () const
{
  set<Face*> es;
  for (iterof (i, boundary_); i != boundary_.end(); i++)
    es.insert (i->second);

  return es;
}

Face*
Mesh_connectivity::find_boundary (Simplex s) const
{
  map<Simplex,Face*>::const_iterator i = boundary_.find (s);
  if (i != this->boundary_.end ())
    {
      return (*i).second;
    }

  i = this->boundary_.find (s.get_mate ());
  if (i != this->boundary_.end ())
    {
      return (*i).second;
    }

  return 0; 
}

Face* 
Mesh_connectivity::make_face (Simplex s)
{
  Face * e = new Face (s);

  Face* mate = 0;
  map<Simplex, Face*> ::const_iterator i =this->boundary_.find (s.get_mate());
  if ( i != this->boundary_.end ())
    {
      mate = (*i).second;
    }

  if (mate)
    {
      remove_boundary (s.get_mate ());
      
      assert (!mate->mate_);
      
      mate->mate_ = e;
      e->mate_ = mate;

    }
  else
    {
      add_boundary (e); 
    }
  
  return e;
}

void
Mesh_connectivity::remove_face (Face*e)
{
  if (e->mate_)
    {
      e->mate_->mate_ = 0;
      add_boundary (e->mate_);
    }
  else
    {
      remove_boundary (e->simplex());
    }
}

void
Mesh_connectivity::remove_face_pair (Face* e)
{
  e->valid_ = e->mate_->valid_ = 0;
}

Face* 
Mesh_connectivity::make_face_pair (Simplex s)
{
  Face * e = new Face(s);
  Face * m = new Face (s.get_mate ());

  m->mate_ = e;
  e->mate_ = m;
  e->valid_ = m->valid_ = 1;
  
  return e;
}
  

/*
  Notes: the order of the return value corresponds with the order of NEWPLEXES, i.e.

  RETVAL[j].simplex() =  NEWPLEXES[j]
*/
Link_array<Element>
Mesh_connectivity::replace_elements (Link_array<Element> old_elts,
				  Array<Simplex> newplexes)
{
#ifdef PARANOIA
  this->validate ();
#endif
    
  Face_map old_faces;
  Face_map::const_iterator ei;


  for (int i=0;  i < old_elts.size (); i++)
    {
      for (int j = 0; j < dimension_ + 1;j++)
	{
	  Face*e = old_elts[i]->face(j);
	  assert (!has_key (old_faces, e->simplex_));
	  
	  old_faces[e->simplex_] = e;
	  e->valid_ = 0;

	  old_elts[i]->faces_[j] = 0;
	}
      old_elts[i]->valid_ = 0;

      add_changed_element (old_elts[i]);
    }


  /*
    Save all the faces that can be reused.
  */
  set<Simplex, Simplex_less> todo_faces;
  Face_map new_face_dict;
  for (int  i = newplexes.size();i--;)
    for (int j = dimension_ + 1; j--;)
      {
	Simplex eplex = newplexes[i].get_subset(j);
	ei = old_faces.find (eplex);
	if (ei != old_faces.end ())
	  {
	    assert (!has_key (new_face_dict, eplex));
	    new_face_dict[eplex] = (*ei).second;
	    old_faces.erase (eplex);
	  }
	else
	  {
	    assert (!has_elt (todo_faces, eplex));

	    todo_faces.insert(eplex);
	  }
      }

  /*
    Make newly required faces. Req'd pairs are created 2 at a time,
    saving boundary ops.
  */
  for (iterof(i, todo_faces); i != todo_faces.end (); i++)
    {
      Simplex splex = *i;
      if (has_elt (todo_faces, splex.get_mate()))
	{
	  if ((*i).parity())
	    {
	      Face *e = this->make_face_pair (splex);
	      new_face_dict[splex] = e;
	      new_face_dict[splex.get_mate()] = e->mate_;
	    }
	}
      else
	{
	  new_face_dict[splex] = this->make_face (splex);
	}
    }

  Link_array<Element> new_elts;
  for (int i= 0; i< newplexes.size(); i++)
    {
      Simplex plex =  newplexes[i];
      Element *t = 0;

      /*
	if set to 1: reuse old elements.
       */
      if (0 && old_elts.size())
	{
	  t = old_elts.pop();
	  t->simplex_ = plex;
	}
      else
	{
	  t = construct_element (plex);
	  this->elements_.insert (t);
	}

      for (int j = dimension_ + 1; j-- ;)
	{
	  Simplex eplex = plex.get_subset (j);

	  Face * e = new_face_dict [eplex];
	  t->faces_[j] = e;
	  e->valid_ = 1;
	  e->element_ = t;
	}
      t->valid_ = 1;
      new_elts.push (t);
    }

  for (map<Simplex,Face*,Simplex_less>::const_iterator i (old_faces.begin ());
       i != old_faces.end () ; i++)
    {
      Simplex p = (*i).first;
      Face *e = (*i).second;
      if (has_key (old_faces, p.get_mate()))
	{
	  if (p.parity())
	    this->remove_face_pair (e);
	}
      else
	{
	  this->remove_face (e);
	}
    }
  
  for (int i = old_elts.size(); i--;)
    elements_.erase (old_elts[i]);

#ifdef PARANOIA
  this->validate ();
#endif

  for (int j = new_elts.size(); j--;)
    add_changed_element (new_elts[j]); 
  
  return new_elts;
}

Element*
Mesh_connectivity::construct_element (Simplex p)
{
  return new Element (p);
}

void
Mesh_connectivity::validate () const
{
#ifndef NDEBUG
  printf ("Validating topology... ");
  fflush (stdout);

  int face_count = 0;   
  set<Face*> *faces =  this->faces();
  for (set<Face*> ::const_iterator i (faces->begin ());
       i != faces->end (); i++)
    {
      Face * e= *i;

      face_count ++;
      assert (dimension_  == e->simplex().count());

      assert (e->valid_);
      if (e->mate_)
	{
	  assert(!has_key (boundary_, e->simplex_));
	  assert(e->mate_->simplex_ == e->simplex_.get_mate ());
	}
      else
	{
	  assert (has_key (this->boundary_, e->simplex_));
	}
    }


  int elt_count = 0;  
  for (set<Element*>::const_iterator i (this->elements_.begin());
       i != this->elements_.end (); i++)
    {
      Element* t = *i;
      assert (dimension_ + 1 == t->simplex().count());
	      
      assert (t->valid_);
      for (int j = 0 ; j < dimension_ + 1; j++)
	{
	  Face * f = t->faces_[j];
	  assert (f->simplex_ == t->simplex_.get_subset(j));
	  assert (f->element_ == t);

	  if (f->mate_)
	    {
	      Element * mate_elt = f->mate_->element_ ;

	      assert (mate_elt != t);

	      /*
		This is a element flipped onto itself.
	       */
	      assert (mate_elt->simplex_ != t->simplex_.get_mate ());
	    }
	}
      elt_count ++;
    }

  assert ((dimension_ + 1)*elt_count == face_count);
  
  for (map<Simplex,Face*, Simplex_less>::const_iterator i (this->boundary_.begin());
       i != this->boundary_.end (); i++)
    {
      Face* e = (*i).second;
      Simplex s = (*i).first;
      assert (e->simplex_ == s);
      assert (e->valid_);
    }
  printf ("done\n");
#endif
}

/*
  the number is for comparisons. The address of the object would be
  sufficient, but it is sensitive to runtime details, making debugging
  a true pain.

  See also node.hh
*/
Node*
Mesh_connectivity::make_node ()
{
  int k = this->nodes_.size();

  Node* n = new Node ();
  n->number_ = k;  

  this->nodes_.push (n);

  return n;
}

Link_array<Node> const*
Mesh_connectivity::node_array()const
{
  return &this->nodes_;
}

void
Mesh_connectivity::print ()const
{
  for (iterof (i,this->elements_);
       i != this->elements_.end (); i++)
    {
      Element* t = *i;
      t->print ();
    }
}

/*
  Assume that a node subsitution is performed on a set of elements,
  such that the old and new nodes are related somehow, then the
  following function ensures that equivalent faces are persistent.

  elt->face (old_simplex.index (old_node)) ==
  elt->face (new_simplex.index (new_node)),

  where new_simplex is substution (old_simplex).

*/


/*

WARNING: this routine hasn't been tested after adding Boundary_watcher
and add_boundary()/remove_boundary ().

Anyway, this code should give you the rough idea.

*/
void
Mesh_connectivity::change_elements (map<Element *,
				 map< Node * , Node * > * > const &elt_map)
{
#ifdef PARANOIA
  validate();
#endif
  
  map<Simplex, Simplex, Simplex_less> simplex_map;
  map<Simplex, Face *, Simplex_less> old_faces;
  map<Simplex, Face *, Simplex_less> new_faces;
  map<Simplex, Element *, Simplex_less> new_elts;  
  map<Simplex, Face*, Simplex_less> surrounding_faces ;
  map<Simplex, Face*, Simplex_less> old_boundary_faces;
  
  for (iterof (t, elt_map); t != elt_map.end(); t++)
    {
      Element *tr =  t->first;

      add_changed_element (tr);
      
      assert (tr->simplex().count () == dimension_ + 1);

      Simplex old_plex =tr->simplex ();
      Simplex new_plex = old_plex;
      
      map<Node*, Node*> *nodmap = t->second;
      int k = 0;
      for (iterof (n, *nodmap); n!= nodmap->end (); n++)
	{
	  if (n->first != n->second)
	    new_plex = new_plex.get_substituted (n->first, n->second);
	  k++ ;
	}

      assert (k <= dimension_ + 1);
      
      for (int i = dimension_+1; i --; )
	{
	  Face * e = tr->face (i);


	  if (has_key (old_faces, e->simplex()))
	    assert (false);
	  old_faces[e->simplex()] = e;
	  
	  Node *nn = tr->node(i);
	  if (has_key (*nodmap, nn))
	    {
	      nn = (*nodmap)[nn];
	    }
	  
	  Simplex new_eplex = new_plex.get_subset (new_plex.index (nn));

	  assert (!has_key (new_faces,new_eplex));
	  new_faces[new_eplex] = e;
	}

      assert (!has_key (simplex_map, old_plex));
      simplex_map[old_plex] = new_plex;
    }

  /*
    Disconnect faces.
  */
  for (iterof (i, old_faces); i != old_faces.end(); i++)
    {
      Face *e = i->second;

      if (e->mate ())
	{
	  if (!has_key (old_faces, e->simplex().get_mate ()))
	    {
	      assert (!has_key (surrounding_faces, e->simplex().get_mate()));
	     
	      surrounding_faces[e->mate ()->simplex()] = e->mate();
	    }


	  e->mate_->mate_ = 0;
	  e->mate_ = 0;
	}
      else
	{
	  assert (!has_key (old_boundary_faces, e->simplex()));
	  old_boundary_faces [e->simplex()] = e;
	}
    }

  /*
    Change simplexes.
  */
  for (iterof (t, elt_map); t != elt_map.end(); t++)
    {
      Element *elt = t->first;
      Simplex old_elt_simplex = elt->simplex_;
      elt->simplex_ = simplex_map[old_elt_simplex];
     
      for (int i = dimension_+1; i--; )
	{
	  Simplex new_face_simplex = elt->simplex_.get_subset (i);
	  Face*f = new_faces [new_face_simplex];
	  elt->faces_[i] = f;

	  Simplex old_face_simplex = f->simplex_;
	  f->simplex_ = new_face_simplex;

	  Simplex mateplex= new_face_simplex.get_mate();
	  Face *mate = 0;

	  if (has_key (new_faces, mateplex))
	    mate = new_faces[mateplex];
	  else if (has_key (surrounding_faces, mateplex))
	    mate = surrounding_faces[mateplex];
	  else if (has_key (boundary_, mateplex))
	    {
	      mate = remove_boundary (mateplex);
	    }

	  if (mate)
	    {
	      if (mate->mate_)
		{
		  assert (mate->mate_ == f); 
		}
	      else
		{
		  f->mate_  =mate;
		  mate->mate_ = f;
		}

	      if (has_key (old_boundary_faces, old_face_simplex))
		{
		  remove_boundary (old_face_simplex); 
		}
	    }
	  else
	    {
	      if (has_key (old_boundary_faces, old_face_simplex))
		{
		  if (old_face_simplex != f->simplex_)
		    {
		      remove_boundary (old_face_simplex);
		      add_boundary (f);
		    }
		}
	      else
		{
		  assert (!has_key (boundary_, f->simplex_));
		  add_boundary (f);
		}
	    }
	}
    }

  /*
    Fixup surrounding faces.  All the surrounding faces were originally
    connected to the faces of the victim elements. If their mate_
    fields are still not connected to new faces, then they are a new
    part of the boundary.
  */
  for (iterof(i, surrounding_faces); i != surrounding_faces.end(); i++)
    {
      if (!i->second->mate_)
	{
	  add_boundary (i->second);
	}
    }

  /*
    Old boundary faces may have left the boundary
  */
  for (iterof (i, old_boundary_faces); i != old_boundary_faces.end (); i++)
    {
      if (i->second->mate_  )
	{
	  remove_boundary (i->first);
	}
    }
#ifdef PARANOIA
  validate();
#endif
}

int
Mesh_connectivity::node_count() const
{
  return nodes_.size();
}


void
Mesh_connectivity::add_changed_boundary (Face * f)
{
  for (int i = 0; i < watchers_.size ();   i++)
    {
      watchers_[i]->add_changed_boundary (f); 
    }
}


void
Mesh_connectivity::add_changed_element (Element * e)
{
  for (int i = 0; i < watchers_.size ();   i++)
    {
      watchers_[i]->add_changed_element (e); 
    }
  
}


/*
  TODO: this can be done more generically by superclassing
  Deformation_state.
 */
void
Mesh_connectivity::add_watcher (Mesh_connectivity_watcher*w)
{
  watchers_.push (w);
  w->mesh_ = this;
  w->init_elements (&elements_);
  w->init_boundary (&boundary_);  
}

Mesh_connectivity::~Mesh_connectivity()
{
  
}

void
Mesh_connectivity::add_boundary (Face * f)
{
  assert (!has_key (boundary_, f->simplex()));
  boundary_[f->simplex()] = f;
  add_changed_boundary (f);	  
}

Face*
Mesh_connectivity::remove_boundary (Simplex const& p)
{
  map<Simplex,Face*>::iterator i = boundary_.find (p);
  assert (i!=boundary_.end());

  Face * f = i->second;
  
  boundary_.erase (i);

  add_changed_boundary (f);
  return f;
}

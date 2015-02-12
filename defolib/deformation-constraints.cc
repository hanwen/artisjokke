#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defo-proto.hh"
#include "deformation-constraints.hh"
#include "deformation-state.hh"
#include "node.hh"
#include "big-vector.hh"
#include "mesh-feature.hh"
#include "convergence-statistics.hh"


/*
  remove constraint on REMOVE, and add or overwrite the constraint on ADD. 
 */
void
Deformation_constraints::change_linear_movement_constraint (Node *add,
							    Real const *addv,
							    Node *remove)
{
  assert (add != remove);
  int sd = parent_->spatial_dimension();
  if (add)
    {
      Real * v =new Real[sd];
      Real l = 0.0; 
      for (int i = 0; i < sd ; i++)
	l += sqr (addv[i]);

      l = 1/sqrt (l);
      for (int i = 0; i < sd ; i++)
	v[i]  = addv[i] * l;

      int k = add->number();
      line_constraints_[k] = v ; 

      position_constraints_.erase(k);
    }
  
  if (remove)
    {
      int k = remove->number();
      line_constraints_.erase(k);      
    }

  changed_ = changed_ || remove || add;
}

/*
  remove constraint on REMOVE, and add or overwrite the constraint on ADD. 
 */
void
Deformation_constraints::change_ortho_movement_constraint (Node *add,
							   Real const *addv,
							   Node *remove)
{
  assert (add != remove);
  int sd = parent_->spatial_dimension();
  if (add)
    {
      Real * v =new Real[sd];
      Real l = 0.0; 
      for (int i = 0; i < sd ; i++)
	l += sqr (addv[i]);

      l = 1/sqrt (l);
      for (int i = 0; i < sd ; i++)
	v[i]  = addv[i] * l;

      int k = add->number();

      ortho_line_constraints_[k] = v ; 
      line_constraints_ .erase (k);
      position_constraints_.erase(k);
    }
  
  if (remove)
    {
      int k = remove->number();
      line_constraints_.erase(k);      
    }

  changed_ = changed_ || remove || add;
}

/*
  Find the reaction forces; these are perpendicular on the movement
  constraints.
*/
void
Deformation_constraints::find_reactions (Real *dest, Real const * src) const
{
  assert (src != dest);
  int n = dimension();
  big_vector_nullify (dest, n);
  
  int spatial_dim = parent_->spatial_dimension (); 

  for (iterof (i,line_constraints_); i != line_constraints_.end (); i++)
    {
      int k = (*i).first;
      
      Real ip = 0.0;
      for (int j = spatial_dim; j--;)
	ip += src[spatial_dim*k+j] * (*i).second[j];
      for (int j = spatial_dim; j--;)
	dest[spatial_dim*k+j] = src[spatial_dim*k+j] - ip * (*i).second[j];
    }

  for (iterof (i, ortho_line_constraints_); i != ortho_line_constraints_.end (); i++)
    {
      int k = (*i).first;
      Real *dir = (*i).second;
      Real ip = 0.0;
      for (int j = spatial_dim; j--;)
	ip += src[spatial_dim*k+j] * dir[j];
      for (int j = spatial_dim; j--;)
	dest[spatial_dim*k+j] = ip * dir[j];
    }

  for (iterof (i, position_constraints_); i != position_constraints_.end (); i++)
    {
      int k = spatial_dim * (*i);
      for (int j = spatial_dim; j--;)
	{
	  dest[k] = src[k];
	  k++;
	}
    }
}

/*
  Find the reaction forces; these are perpendicular on the movement
  constraints.
*/
bool
Deformation_constraints::satisfies_constraints (Real const *src) const
{
  assert (false);
  
#if 0
  /*
    TODO make spatial dimension independent.
   */
  int spatial_dim = parent_->spatial_dimension (); 

  for (iterof (i, position_constraints_); i != position_constraints_.end (); i++)
    {
      int k = spatial_dim * *i;
      for (int j = spatial_dim; j--;)
	{
	  if (src[k])
	    return false;
	  k++;
	}
    }
    
  for (iterof (i,line_constraints_); i != line_constraints_.end (); i++)
    {
      int k = (*i).first;
      
      Vector2 src_disp;
      for (int j = spatial_dim; j--;)
	src_disp(j) = src[spatial_dim*k+j];


      Vector2 proj = src_disp -  (src_disp * (*i).second) * (*i).second;
      if (proj.length () / src_disp.length () > 1e-6)
	return false;
    }
#endif
  
  return true;
}


/*
  Project MOV along the direction constraints for a single node NOD.

  MOV and DEST are have length SPATIAL_DIMENSION
 */
void
Deformation_constraints::apply_to_node_movement (Real * dest,
						 Node*nod, Real const * mov)const
{
  int k = nod->number();
  int spatial_dim = parent_->spatial_dimension();
  if (has_key (line_constraints_, k))
    {
      Real * dir = line_constraints_.find (k)->second;

      Real ip = 0.0;
      for (int j = spatial_dim; j--;)
	ip += mov[j] * dir[j];
      for (int j = spatial_dim; j--;)
	dest[j] = ip * dir[j];
    }
  else if (has_key (ortho_line_constraints_, k))
    {
      Real * dir = ortho_line_constraints_.find (k)->second;

      Real ip = 0.0;
      for (int j = spatial_dim; j--;)
	ip += mov[j] * dir[j];
      for (int j = spatial_dim; j--;)
	dest[j] = mov[j] - ip * dir[j];
    }
  else if (has_elt (position_constraints_, k))
    {
      for (int j = spatial_dim; j--;)
	dest[j] = 0.0;
    }

}

/*
  Project MOV along the direction constraints for a single node NOD.

  MOV and DEST are have length SPATIAL_DIMENSION
 */
void
Deformation_constraints::node_reaction (Real * dest,
					Node*nod, Real const * mov)const
{
  int k = nod->number();
  int spatial_dim = parent_->spatial_dimension();
  if (has_key (line_constraints_, k))
    {
      Real * dir = line_constraints_.find (k)->second;

      Real ip = 0.0;
      for (int j = spatial_dim; j--;)
	ip += mov[j] * dir[j];
      for (int j = spatial_dim; j--;)
	dest[j] = mov[j] - ip * dir[j];
    }
  else if (has_key (ortho_line_constraints_, k))
    {
      Real * dir = ortho_line_constraints_.find (k)->second;

      Real ip = 0.0;
      for (int j = spatial_dim; j--;)
	ip += mov[j] * dir[j];
      for (int j = spatial_dim; j--;)
	dest[j] = ip * dir[j];
    }
  else  if (has_elt (position_constraints_, k))
    {
      for (int j = spatial_dim; j--;)
	dest[j] = mov[j];
    }
  else
    {
      for (int j = spatial_dim; j--;)
	dest[j] = 0.0;
    }
}


/*
  Project MOV along the direction constraints.
 */
void
Deformation_constraints::apply_to_movement (Real * dest, Real const *mov)const
{
  int n = dimension();
  int spatial_dim = parent_->spatial_dimension ();
  
  big_vector_copy (dest, mov, n);

  for (iterof(i,line_constraints_); i != line_constraints_.end (); i++)
    {
      int k = (*i).first;
      Real nod_move[3];
      for (int j = spatial_dim; j--;)
	nod_move[j] = dest[spatial_dim*k+j];
      
      Real ip = 0.0;
      for (int j = spatial_dim; j--;)
	ip += nod_move[j] * (*i).second[j];
      
      for (int j = spatial_dim; j--;)
	dest[spatial_dim*k+j] = ip * (*i).second[j];
    }

  for (iterof(i,ortho_line_constraints_); i != ortho_line_constraints_.end (); i++)
    {
      int k = (*i).first;
      Real nod_move[3];
      for (int j = spatial_dim; j--;)
	nod_move[j] = dest[spatial_dim*k+j];
      
      Real ip = 0.0;
      Real *dir = (*i).second;
      for (int j = spatial_dim; j--;)
	ip += nod_move[j] * dir[j];
      
      for (int j = spatial_dim; j--;)
	dest[spatial_dim*k+j] = nod_move[j] - ip * dir[j];
    }

  vector_flop_count += line_constraints_.size() * 8  +
    ortho_line_constraints_.size () *8  
    ;
  for (iterof (i, position_constraints_); i != position_constraints_.end (); i++)
    {
      int k = spatial_dim * (*i);
      for (int j = spatial_dim; j--;)
	{
	  dest[k] = 0;
	  k++;
	}
    }
}

bool
Deformation_constraints::has_linear_movement_constraint (Node*n)const
{
  return has_key (line_constraints_, n->number());
}

bool
Deformation_constraints::has_ortho_movement_constraint (Node*n) const
{
  return has_key (ortho_line_constraints_, n->number());
}

bool
Deformation_constraints::is_fixed (Node *n)const
{
  parent_->check_for_topology_update ();

  return has_elt (position_constraints_, n->number());
}

/*
  Also removes line constraint if the node is completely fixed.
 */
void
Deformation_constraints::fix_node (Node* nod, bool fix)
{
  parent_->check_for_topology_update ();
  
  int idx = nod->number () ;
  if (fix)
    {
      position_constraints_.insert (idx); 
      line_constraints_.erase (idx);
      ortho_line_constraints_ .erase (idx);
    }
  else
    {
      position_constraints_.erase (idx);
    }
  
  changed_ = true;		// todo: check if really changed.
}


bool
Deformation_constraints::changed () const
{
  return changed_;
}


Deformation_constraints::Deformation_constraints (Deformation_state * def)
{
  parent_ = def;
  changed_ = false;
}

int
Deformation_constraints::dimension () const
{
  return parent_->dimension();
}

void
Deformation_constraints::remove_all_node_constraints (Node * nod)
{
  int idx = nod->number () ;
  position_constraints_.erase (idx); 
  line_constraints_.erase (idx);
  ortho_line_constraints_ .erase (idx);
}

Real
Deformation_constraints::reaction_length_sq (Real const *src) const
{
  int spatial_dim = parent_->spatial_dimension ();
  Real reaction_len = 0.0;

  for (iterof (i,line_constraints_); i != line_constraints_.end (); i++)
    {
      int k = (*i).first;
      
      Real ip = 0.0;
      for (int j = spatial_dim; j--;)
	ip += src[spatial_dim*k+j] * (*i).second[j];
      
      for (int j = spatial_dim; j--;)
	reaction_len += sqr (src[spatial_dim*k+j] - ip * (*i).second[j]);
    }

  for (iterof (i, position_constraints_);
       i != position_constraints_.end (); i++)
    {
      int k = spatial_dim * (*i);
      for (int j = spatial_dim; j--;)
	{
	  reaction_len += sqr (src[k]);
	  k++;
	}
    }

  return reaction_len;
}

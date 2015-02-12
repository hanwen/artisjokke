#ifndef AUTO_INSERT_HH
#define AUTO_INSERT_HH

#include "proto.hh"
#include "vector.hh"

struct Auto_needle_insert
{
  Real speed_;

  Real angle_;
  Real max_depth_;
  Real depth_;
  Vector2 dir_;
  Vector2 handle_;
  
  Auto_needle_insert();
};

void finish_auto_insert (Maubach_tree * tree,
	       Deformation_state *def,
	       Needle_inserter * ins);


#endif

#include <stdio.h>

#include "node.hh"

Node::Node()
{
  number_= -1;
}

void
Node::print () const
{
  fprintf (stderr,  "%d ", number_);
}



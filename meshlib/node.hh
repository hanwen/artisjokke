
#ifndef NODE_HH
#define NODE_HH


/*
  Node is a relic of a time when it also contained a location. Right
  now, the only thing over an int is the fact it can be distinguished
  from ints in type-checks.

  Other than that, the extra indirection adds a little overheid.
*/
class Node
{
  int number_;

  Node();
public:

  bool valid () const { return true ; }	// todo
  int number () { return number_ ; }
  void print ()const;
  friend class Mesh_connectivity;
};
#endif

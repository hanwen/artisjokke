/*   
topol-set.hh -- declare  Simplex_dimension, Simplex

(c) 1999 Han-Wen Nienhuys <hanwen@cs.uu.nl>

 */

#ifndef Simplex_HH
#define Simplex_HH

#include <stdio.h>

#include "node.hh"

enum Simplex_dimension {
  NODE = 0,
  EDGE = 1,
  TRIANGLE = 2,
  TETRAHEDRON = 3,
  MAX_DIMENSION=3,

};

struct Simplex
{
protected:
  Node * elements_[MAX_DIMENSION+1];
  char count_;
  char parity_;
public:
  char parity () const { return parity_; } 
  void print_on (FILE*)const;
  void short_print_on (FILE*)const;  
  void print () const;
  void short_print () const;
  Node * node (int i)const { return elements_[i]; }
  int count () const { return count_; } 
  
  Simplex_dimension dimension () const { return Simplex_dimension (count_ -1); }
  Simplex ();
  Simplex (Simplex_dimension, Node **, int par);
  Simplex get_superset (Node*) const;
  Simplex get_subset (Node*n) const { return get_subset (index (n)); }
  Simplex get_subset (int) const;
  Simplex get_mate () const;
  Simplex get_opposite (Simplex sub) const;

  Simplex intersection_nodes (Simplex other) const;
  Node** node_array() const {
    /*
      Warning: this function yields R/W access, but is only meant for
      read-only applications
    */
    return (Node **) elements_ ;

  }
  void sort ();
  
  bool contains (Node*) const;
  bool has_subset (Simplex) const;
  bool valid_nodes () const;
  
  Simplex get_substituted (Node*o, Node*d) const;
  void swap (int,int);
  static unsigned int ordered_hash (Simplex g);
  int index (Node *) const;
  static int compare (Simplex const & t1 , Simplex const & t2);
};

struct Simplex_less {
  bool operator () (Simplex const &t1, Simplex const &t2) const
  {
    return Simplex::compare (t1, t2) < 0;
  }
};


bool operator == ( Simplex const &g1,Simplex const &g2);
bool operator != ( Simplex const &g1,Simplex const &g2);
inline
bool operator < (Simplex const &g1,Simplex const &g2)
{
  return Simplex::compare (g1, g2) < 0;
}


inline int
Simplex::compare (Simplex const &t1, Simplex const & t2)
{
  if (t1.count_ - t2.count_)
    return t1.count_ - t2.count_;

  for (int i =0; i < t1.count_; i++)
    {
      int d = t1.node (i)->number ()  - t2.node (i)->number () ;
      if (d)
	return d;
    }

  return t1.parity_ - t2.parity_;
}


#endif /* Simplex_HH */

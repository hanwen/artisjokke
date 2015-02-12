/*   
topological-set.cc --  implement Simplex

(c) 1999 Han-Wen Nienhuys <hanwen@cs.uu.nl>

 */

#include <assert.h>

#include "simplex.hh"
#include "node.hh"


Simplex
Simplex::get_mate () const
{
  Simplex ret (*this);
  ret.parity_ = ! parity_;
  
  return ret;
}

void
Simplex::sort ()
{
  // bubble sort.
  for (int i = 0;  i < count_; i++)
    for (int j = i +1; j < count_ ; j++)
      {
	if (elements_[i]->number ()  > elements_[j]->number () )
	  {
	    Node * t = elements_[i];
	    elements_[i] = elements_[j];
	    elements_[j] = t;
	    parity_ = ! parity_;
	  }
	else
	  /*
	    All elts must be unique.
	   */
	  assert (elements_[i]->number ()  < elements_[j]->number () ); 
      }
}

void
Simplex::print () const
{
  print_on (stdout);
}

void
Simplex::short_print () const
{
  short_print_on (stdout);
}

void
Simplex::print_on (FILE *out)const
{
  if (!count_)
    fprintf (out, "()");
  for (int i =0; i < count_; i++)
    {
      Node * n = elements_[i];
      fprintf (out, "%d ", n->number());
      if (i < count_-1)
	fprintf (out, ", ");
    }

  fprintf (out, " %seven\n", (parity_) ? "un" : "") ;
}

void
Simplex::short_print_on (FILE *out)const
{
  if (!count_)
    fprintf (out, "()");
  for (int i =0; i < count_; i++)
    {
      Node * n = elements_[i];
      fprintf (out, "%d", n->number());
      if (i < count_-1)
	fprintf (out, ", ");
    }

  fprintf (out, " %seven\n", (parity_) ? "un" : "") ;
}

unsigned int 
Simplex::ordered_hash (Simplex g)
{
  unsigned int h = 0;
  for (int i=0; i < g.count_; i++)
    {
      h = h  + (i+1) * g.elements_[i]->number ()  ; 
    }
  h = 2  * h  + g.parity_ ? 1 :0;
  return h;
}


Simplex::Simplex ()
{
  count_ = 0;
  parity_ = 0;
}

Simplex::Simplex (Simplex_dimension n, Node **gs, int parity)
{
  count_ = n+ 1;		// face = 1-dim, 2 elts.
  for (int i=0; n && i < count_; i++)
    elements_[i] = gs[i];

  parity_ = parity;
  sort ();
}

bool
Simplex::valid_nodes ()const
{
  bool  v = true;
  for (int i= count (); i--;)
    v = v && node(i)->valid ();
  
  return v; 
}

bool
operator==(Simplex const &g1,Simplex const &g2)
{
  return 0 == Simplex::compare (g1,g2);
}

bool
operator!=(Simplex const &g1,Simplex const &g2)
{
  return ! (g1 == g2);
}




bool
Simplex::contains (Node *np) const
{
  for (int i=0; i < count_; i++)
    if (elements_[i] == np)
      return true;
  return false;
}


/*
  TODO:

  provide proof of this routine. It is troublesome that a bug still was found
  years after it was first written.
  
 */
bool
Simplex::has_subset (Simplex s) const
{
  int not_found_idx  = -1;

  Node ** a = s.elements_;
  Node **aend = s.elements_ + s.count_;
  Node *const * b = elements_;
  Node *const * bend = elements_ + count_;

  /*

    Assuming a and b are both monotone, we have
    
    a sub b <= (x:a sub x:b)

    a sub b <= (a sub x:b)

    [] sub b <= true

    a sub [] => false
    
   */

  while (b < bend && a < aend)
    {
      if (*a == *b)
	{
	  a++;
	  b++;
	}
      else
	{
	  not_found_idx = b - elements_ ;
	  b++;
	}
    }

  if (a < aend && b == bend)
    return false;

  
  
  /*
    we have a == aend, b < bend
   */
  if (s.count_ == count_-1)
    {
      /*
	In case not_found_idx is the last one.
       */
      if (not_found_idx < 0)
	not_found_idx = count_ - 1;
      return s.parity_ == (parity_ ^ (not_found_idx % 2));	
    }
  else
    {
      return true;
    }
}

Simplex
Simplex::get_subset (int k) const
{
  Simplex ret;
  int j =0;
  ret.count_ =  count_ - 1;
  
  for (int i=0; i  < count_;  i++)
    if (i != k)
      ret.elements_[j++] = elements_[i];

  ret.parity_ = (k%2) ? !parity_  : parity_;

  return ret;
}

Simplex
Simplex::get_superset (Node *n) const
{
  Simplex ret;
  ret.count_ = count_ + 1;
  ret.elements_[0] = n;
  for (int i = 1; i < ret.count_; i++)
    ret.elements_[i] = elements_[i-1];

  ret.parity_ = parity_;
  ret.sort();
  return ret;
}


Simplex
Simplex::get_substituted (Node *o, Node *d) const
{
  if (o == d)
    return *this;
  
  Simplex t(*this);
  
  for (int k = 0;  k < t.count_; k++)
    if (t.elements_[k] == o)
      t.elements_[k] = d;
  
  t.sort ();
  return t;
}


int
Simplex::index (Node*n)const
{
  for (int  i= count_; i--;)
    if (n == elements_[i])
      return i;

  return -1;
}


  
Simplex
Simplex::get_opposite (Simplex sub)const
{
  Simplex s (*this);
  
  for (int i = sub.count_; i--;)
    s = s.get_subset (s.index (sub.node(i)));
  return s;
}

  
Simplex
Simplex::intersection_nodes (Simplex other) const
{
  Simplex r;
  for (int i = 0 ;  i < count_;i++)
    {
      if (other.contains(elements_[i]))
	r.elements_[r.count_++] = elements_[i];
    }
  

  return r;
}

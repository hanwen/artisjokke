#ifndef STL_HH
#define STL_HH
#include <map>
#include <set>

using std::set;
using std::map;

#define iterof(i,s) typeof((s).begin()) i((s).begin())

template<class T, class L> 
inline bool
has_elt (set<T, L> const &s, T const &k)
{
  return s.find (k) != s.end ();
}

template<class T, class V, class L> 
inline bool
has_key (map<T, V, L> const &s , T const &k)
{
  return s.find (k) != s.end ();
}



template<class T, class L> 
inline bool
is_empty (set<T, L> const &s)
{
  return s.begin () == s.end ();
}

#endif

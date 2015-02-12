#ifndef EDGE_WATCHER
#define EDGE_WATCHER

#include "stl.hh"
#include "mesh-connectivity-watcher.hh"

class Edge_watcher : public Element_boundary_watcher
{
  set<Simplex> edges_;
  bool changed_;
  
  void process_changes ();
  void process_feature (Mesh_feature*, bool);
public:
  Array<Simplex> * edge_array ();
  set<Simplex> const * edges();
  Edge_watcher();

  virtual ~Edge_watcher();
};



#endif 

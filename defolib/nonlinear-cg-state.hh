/*   
  linear-cg-state.hh -- declare Mesh_connectivity_state
  
  (c) 1999 Han-Wen Nienhuys <hanwen@cs.uu.nl>
  
 */

#ifndef NONLINEAR_CG_STATE_HH
#define NONLINEAR_CG_STATE_HH

#include "line-search-state.hh"

struct Nonlinear_cg_state
  : public Line_search_state
{
  bool use_hestenes_stiefel_;
  bool use_polak_ribiere_;
  bool use_fletcher_reeves_;

  Array<Real> prev_residual_;  
  int last_restart_; 
  
  Nonlinear_cg_state ();
  virtual void set_new_search_dir(Deformation_state*);
protected:
  Link_array< Array<Real> > vector_variables ();
};

#endif /* LINEAR-CG_STATE_HH */


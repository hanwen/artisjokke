/*   
  linear-cg-state.hh -- declare Mesh_connectivity_state
  
  (c) 1999 Han-Wen Nienhuys <hanwen@cs.uu.nl>
  
 */

#ifndef LINEAR_CG_STATE_HH
#define LINEAR_CG_STATE_HH

#include "defo-proto.hh"
#include "parray.hh"
#include "iteration-state.hh"

struct Linear_cg_state : public Iteration_state
{
  Array<Real> search_dir_;
  Array<Real> mat_times_search_dir_;

  Linear_cg_state ();

  void cg_iteration_step (Real*,Real*,Real*,
			  Real const*, Real, Real,
			  Deformation_constraints *,
			  void*,
			  void (*)(void*, Real*,Real const*));
  void cg_reinit ();
  
  Link_array< Array<Real> > vector_variables ();
protected:
  virtual void iteration_step(Deformation_state *);
  virtual void signal_residual_change (Deformation_state *);

  virtual bool iteration_updates_residual () const;

};

#endif /* LINEAR-CG_STATE_HH */


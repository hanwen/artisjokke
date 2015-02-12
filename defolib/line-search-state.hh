/*   
  linear-cg-state.hh -- declare Mesh_connectivity_state
  
  (c) 1999 Han-Wen Nienhuys <hanwen@cs.uu.nl>
  
 */

#ifndef LINESEARCH_STATE_HH
#define LINESEARCH_STATE_HH

#include "defo-proto.hh"
#include "parray.hh"
#include "iteration-state.hh"


struct Line_search_state
  : public Iteration_state
{
  Real residual_tolerance_;
  Array<Real> search_dir_;
  Array<Real> mat_times_search_dir_;
  
  Array<Real> try_displacement_;
  Array<Real> try_residual_;
  Array<Real> try_elastic_force_;
  Array<Real> prev_residual_;  
  
  bool no_more_improvement_b_;	// set to true, if CG can't improve current sol.
  bool restart_b_;

  Real (*step_length_function_)(Line_search_state *state,
				Deformation_state const*);
  Real newton_tolerance_;

  Line_search_state ();
  
  virtual void signal_residual_change (Deformation_state*);
  virtual void set_new_search_dir (Deformation_state*);  
  virtual void iteration_step(Deformation_state *);

  virtual void print_nodes ()const;
  virtual void print () const;
  virtual void calculate_try_residual_and_derivative (Real, Deformation_state const*);
  virtual ~Line_search_state ();
protected:
  Link_array< Array<Real> > vector_variables ();

  static Real newton_method (Line_search_state*, Deformation_state const*);
  static Real secant_method (Line_search_state*, Deformation_state const*);
  static Real semi_secant_newton_method (Line_search_state*, Deformation_state const*);
  static Real linear_alpha (Line_search_state*, Deformation_state const*);
  virtual bool iteration_updates_residual () const;

  void calculate_try_residual (Real, Deformation_state const*);
  
  void graph_alpha (Real, Deformation_state const*);
  int dimension () const;
};



#endif /* LINEAR-CG_STATE_HH */


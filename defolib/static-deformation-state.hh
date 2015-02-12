/*   
static-deformation-state.hh -- declare Static_deformation_state

(c) 2001 Han-Wen Nienhuys <hanwen@cs.uu.nl>

 */

#ifndef STATIC_DEFORMATION_STATE_HH
#define STATIC_DEFORMATION_STATE_HH


#include "deformation-state.hh"

class Static_deformation_state : public Deformation_state 
{
  Iteration_state * iter_state_;
  bool switch_to_newton_;

public:
  Static_deformation_state (int);
  bool good_solution () const;
protected:
  virtual void do_one_iteration ();
  
  virtual void signal_boundary_condition_change (); 
  virtual void update_topology (set<Element*> *);
  virtual void signal_residual_change ();

  void set_iteration_state (Iteration_state*);
  Link_array< Array<Real> > other_vector_variables ();

};

#endif /* STATIC_DEFORMATION_STATE_HH */



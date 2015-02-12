/*   
iteration-state.hh -- declare Iteration_state

(c) 2001 Han-Wen Nienhuys <hanwen@cs.uu.nl>

 */

#ifndef ITERATION_STATE_HH
#define ITERATION_STATE_HH

#include "defo-proto.hh"
#include "parray.hh"

class Iteration_state
{
public:
  
  int iter_count_;

  Iteration_state ();


  virtual Link_array< Array<Real>  > vector_variables ();
  virtual void iteration_step (Deformation_state*);
  virtual void signal_residual_change (Deformation_state*);
  virtual void print_nodes () const;
  virtual void print () const;
  virtual bool good_solution (Deformation_state const*) const;


  virtual bool iteration_updates_residual () const;
  
  virtual ~Iteration_state ();
};

#endif /* ITERATION_STATE_HH */


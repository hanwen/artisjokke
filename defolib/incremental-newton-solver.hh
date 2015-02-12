/*   
  nonlinear-newton.hh -- declare  Incremental_newton_solver

  source file of the GNU LilyPond music typesetter

  (c) 2001 Han-Wen Nienhuys <hanwen@cs.uu.nl>

 */

#ifndef NONLINEAR_NEWTON_HH
#define NONLINEAR_NEWTON_HH

#include "iteration-state.hh"

class Incremental_newton_solver : public Iteration_state
{
public:
  Array<Real> cg_constrained_residual_;
  Array<Real> cg_elastic_force_;
  Array<Real> cg_right_hand_side_;
  
  Real cg_residual_tolerance_;

  Array<Real> direction_;
  
  Linear_cg_state * cg_state_;
  Deformation_state * def_; 
  
  virtual void iteration_step (Deformation_state*p);


  virtual Real find_descent_step (Real);
  virtual bool iteration_updates_residual ()const;
  
  virtual Link_array< Array<Real > > vector_variables();


  int dimension () const { return cg_constrained_residual_.size(); }

  Incremental_newton_solver (Deformation_state*);
  ~Incremental_newton_solver();
};

#endif /* NONLINEAR_NEWTON_HH */


/*   
  Dynamic-deformation-state.hh -- declare Dynamic_deformation_state

  (c) 2001 Han-Wen Nienhuys <hanwen@cs.uu.nl>

 */

#ifndef DYNAMIC_DEFORMATION_STATE_HH
#define DYNAMIC_DEFORMATION_STATE_HH


#include "deformation-state.hh"

class Dynamic_deformation_state :public Deformation_state
{
  Array<Real> velocity_;
  Array<Real> lumped_masses_;
  Array<Real> inverse_masses_;
  Array<Real> scratch0_, scratch1_, scratch2_, scratch3_, scratch4_;

  Real time_step_;
public:
  Real viscosity_, density_;
  Real minimum_edge_length_;
  Real time_;

  Real get_time_step () { return time_step_; }
  
  void (*time_integration_function_)(Dynamic_deformation_state*, Real);

  Dynamic_deformation_state (int);
  void relaxation_step (Real);
  
  void compute_acceleration (Real*,Real const *, Real const*);
  void compute_acceleration_with_given_residual (Real * , Real const* , Real const* );
  static void explicit_ss22_multistep (Dynamic_deformation_state*,Real);
  static void explicit_ss22_singlestep (Dynamic_deformation_state*,Real);  
  static void classical_runge_kutta (Dynamic_deformation_state*,Real);
  static void euler_forward (Dynamic_deformation_state*,Real);
  virtual bool good_solution ();
  Real euler_forward_critical_time_step () const;
  Real explicit_ss22_critical_time_step () const;
  Real wave_speed () const;
  Real critical_time_step () const;
  
  Vector2 node_velocity2 (Node*) const;
  void set_node_velocity (Node*, Real const*);
  Real now () const; 
  
  virtual void do_one_iteration ();


  //  virtual void signal_boundary_condition_change ();

  virtual Link_array< Array<Real> > other_vector_variables ();
  virtual Link_array< Array<Real> > interpolated_vector_variables ();
public:
  
  virtual void update_topology (set<Element*> *);
};

#endif /* DYNAMIC_DEFORMATION_STATE_HH */


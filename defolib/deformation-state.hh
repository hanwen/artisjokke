/*   
  problem-state.hh -- declare deformation_state

  
  (c) 2001 Han-Wen Nienhuys <hanwen@cs.uu.nl>
 */

#ifndef DEFORMATION_STATE_HH
#define DEFORMATION_STATE_HH

#include "defo-proto.hh"
#include "parray.hh"
#include "deformation-constraints.hh"
#include "mesh-connectivity-watcher.hh"



class Deformation_state : public Element_watcher
{
  Real constrained_residual_len_sq_;
  
protected:
  int iteration_frequency_;
  Link_array<Deformation_hook> hooks_;
  /*
  properties of the iteration type:
  */
  bool iteration_updates_residual_;


  //
  bool need_reinit_b_;
  bool live_;
  bool ignore_degenerate_;
  bool displacements_changed_;
  bool forces_changed_;
  
  int last_node_count_;
  int spatial_dimension_;
public:
  int spatial_dimension () const { return spatial_dimension_; }

  Array<Real> displacements_;	// naming!

  /*
    TODO: should move to base class Mesh_geometry
   */
  Array<Real> reference_locations_;
  Array<Real> external_force_;
  Array<Real> constrained_residual_;
  Array<Real> elastic_force_;

  Real prev_residual_len_sq_;
  Deformation_constraints constraints_;

  void suicide ();
  void add_hook (Deformation_hook * hook_ptr);
  
  /*
    run-time state:
   */
  Real external_force_len_sq_;
  Real tolerance_;

  /*
    Needed for calibrating the stopping criterion.
   */
  Real diameter_;


  /*
    for generating flop stats.
   */
  Array<Real> statistic_vector1_, statistic_vector2_;
  
  map<Element*, Element_state*> states_;
  Link_array<Element_state> state_array_;

  Deformation_state * reference_state_;

protected:
  Elasticity_functions * elasticity_;

  void apply_force_derivative (Real *der, Real *force, Real*, Real const*displ, Real const *ir) const;

public:
  

  virtual bool good_solution()const;
  int dimension () const;
  void apply_force_function (Real *dest, Real const*src) const;
  void validate_numerical_sanity ();

  void apply_force_derivative_and_residual (Real *der, Real *force, Real*, Real const*displ, Real const *ir) const;
  Real compute_constrained_residual_force (Real *dest, Real const * displ) const;
  void compute_elastic_force (Real *dest, Real const * displ) const;
  void constrained_residual_from_elastic_force (Real*cr, Real*ucr) const;
  void add_force (Node*, Real  const *);
  void set_node_deformed_location (Node*, Real const* );
  void set_reference_location (Node*, Real const*);

  void calibrate_deformation_speed (); 
  bool force_changed () const;

  Deformation_state (int);
  Real virtual_energy (Real const*) const;
  void completize_nodal_arrays (int);
  void check_for_topology_update ();
  void simulation_body ();

  void set_force (Node*, Real const*);

  Array<Real> interpolate_deformation_variables (Real const *, Element *e);
  void set_deformation_variables (Node*, Array<Real> const &);
  void flop_account_tick ();
  void set_residual_len(Real);
  Real get_residual_len_sq () const { return constrained_residual_len_sq_; }
  Real scale_free_force_comparison ()const;
  virtual ~Deformation_state() ;
  virtual void update_forces ();

protected:
  virtual void do_one_iteration();
  virtual void update_topology (set<Element*> *);
  virtual void signal_boundary_condition_change ();
  virtual void signal_residual_change ();
  
  Link_array< Array<Real> > vector_variables ();
  Link_array< Array<bool> > bool_variables ();
  
  virtual Link_array< Array<Real> > other_vector_variables ();
  virtual Link_array< Array<Real> > interpolated_vector_variables ();  

private:
  bool do_simulation_body ();

  /*
    Don't copy.
   */
  Deformation_state (Deformation_state const & );
};




void apply_standard_force_configuration (Deformation_state *def,
					 Mesh_connectivity *top);
bool only_positive_volumes (Mesh_connectivity* mesh, Deformation_state *def);
void set_in_big_vector (Node* nod, Real * big_array, Real const *, int);

Deformation_state *make_new_deformation_state (Mesh_connectivity*top);

void write_deformation_state (Deformation_state const*def, char const *fn);
bool read_deformation_state (Deformation_state * def, char const * fn);


void apply_gravity (Deformation_state*def, int axis, int sign);


/*
  2D vectors.
 */
Vector2 get_external_force2 (Node*, Deformation_state const * def) ;
Vector2 reference_location2 (Node*,Deformation_state const * ref) ;
Vector2 deformed_location2  (Node* nod, Deformation_state const *def);
Vector2 extract_vector2 (Node* nod, Real const * big_array);
Vector2 transform_back2 (Simplex const &, Vector2 x, Deformation_state * def);

/*
  3D
 */
Vector3 get_external_force3 (Node*, Deformation_state const * def) ;
Vector3 reference_location3 (Node*,Deformation_state const * ref) ;
Vector3 deformed_location3  (Node* nod, Deformation_state const *def);
Vector3 extract_vector3 (Node* nod, Real const * big_array);
Vector3 transform_back3 (Simplex const &, Vector3 x, Deformation_state * def);




#endif /* DEFORMATION_STATE_HH */


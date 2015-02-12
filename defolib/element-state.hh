/*   
  element-state.hh -- declare Element_state
  
  (c) 1999 Han-Wen Nienhuys <hanwen@cs.uu.nl>
 */

#ifndef ELEMENT_STATE_HH
#define ELEMENT_STATE_HH

#include "defo-proto.hh"
#include "matrix3.hh"
#include "matrix2.hh"

struct Element_state
{
  bool degenerate_b_;
  Real volume_;

  Real minimum_edge_length_;
  
  virtual void precompute (Element*, Deformation_state const*); 
  virtual Real linear_elastic_energy (Real const *) const;  
  virtual Real elastic_energy (Elasticity_functions const*, Real const *) const;
  virtual int elastic_force (Elasticity_functions const*, Real*,Real const*) const;
  virtual int elastic_force_derivative (Elasticity_functions const*, Real*,Real*, Real const*,Real const*) const;
  virtual ~Element_state() ;
};

class Element_state;
typedef Real (*Element_energy) (Real,Real,Real);

typedef int (*Element_force_function2) (Element_state const*,Matrix2*,Matrix2 const &);
typedef int (*Element_dforce_function2) (Element_state const*,Matrix2*,Matrix2*, Matrix2 const &,Matrix2 const &);

struct Elasticity_functions
{
  int dim_;
  bool linear_;
  virtual ~Elasticity_functions(){}
  Elasticity_functions();
};

/****************/


struct Elasticity_functions2 : Elasticity_functions
{
  Element_force_function2 force_function_;
  Element_dforce_function2 dforce_function_;
  Element_energy energy_function_;
};

struct Element_state2 : public Element_state
{
  int incident_nodes_[3];
  
  /*
    Entry i,j -> force on i-th entry of node number j.

    Node no. 3 (4th node) is
    
   */
  Matrix2 inverse_location_;
  Matrix2 invloc_transpose_;

  
  Element_state2 ();

  void precompute (Element*, Deformation_state const*); 
  
  void big_vector_to_matrix (Matrix2 * dest_mat, Real const*src_vec) const; 
  void matrix_to_big_vector (Real * dest_vec, Matrix2 const &src_mat) const;
  
  virtual Real linear_elastic_energy (Real const *) const;  
  virtual Real elastic_energy (Elasticity_functions const*, Real const *) const;
  virtual int elastic_force (Elasticity_functions const*, Real*,Real const*) const;
  virtual int elastic_force_derivative (Elasticity_functions const*, Real*,Real*, Real const*,Real const*) const;
};

/****************/


typedef int (*Element_force_function3) (Element_state const*,Matrix3*,Matrix3 const &);
typedef int (*Element_dforce_function3) (Element_state const*,Matrix3*,Matrix3*, Matrix3 const &,Matrix3 const &);

struct Element_state3 : public Element_state
{
  int incident_nodes_[4];
  
  /*
    Entry i,j -> force on i-th entry of node number j.

    Node no. 3 (4th node) is
    
   */
  Matrix3 inverse_location_;
  Matrix3 invloc_transpose_;

  Element_state3 ();

  void precompute (Element*, Deformation_state const*); 
  
  void big_vector_to_matrix (Matrix3 * dest_mat, Real const*src_vec) const; 
  void matrix_to_big_vector (Real * dest_vec, Matrix3 const &src_mat) const;
virtual Real linear_elastic_energy (Real const *) const;  
  virtual Real elastic_energy (Elasticity_functions const*, Real const *) const;
  virtual int elastic_force (Elasticity_functions const*, Real*,Real const*) const;
  virtual int elastic_force_derivative (Elasticity_functions const*, Real*,Real*, Real const*,Real const*) const;
 
};

struct Elasticity_functions3 : Elasticity_functions
{
  Element_force_function3 force_function_;
  Element_dforce_function3 dforce_function_;
  Element_energy energy_function_;
};


Element_state * get_element_state (int dim);

#endif /* TETRAHEDRON_STATE_HH */



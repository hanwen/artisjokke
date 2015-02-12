/*
   Veronda & Westmann proposed a function with exponential
  nonlinearity, leaving the exact I3 dependency unspecified.

  This function is based on the paper "An implicit FE method for 
  elastic solids in contact by Gentaro Hirota et.al."

  However, it is not consistent with the linear model.
    
  W(I1,I2, I3) = 0.5 (mu/gamma (exp (gamma(I1 -3)) - 1) - mu/2 (I2 - 3) +
    lambda (I3 - ln (I3)))


  S =  [mu exp(gamma (I1-3)) I - mu/2 (I1 I - C) + lambda [I3 - 1] C^-1]
  
  force = - F S X^-T
  
  For the derivative we get

  Dforce = - (DF S  + F DS) X^-T

  where

  DS = mu exp(gamma (I1-3)) gamma tr(DC) Id - mu/2 (tr(DC) Id. - DH)
       + lambda (I3 C^-1 : DC) C^-1 - lambda (I3-1) C^-1 DC C^-1

     = (exp (..) gamma - 1/2) mu tr(DC) Id.+ mu/2 DH
       + lambda ( I3 C^-1 : DC - (I3-1) C^-1 DC) C^-1
  
 */

#include <math.h>

#include "element-state.hh"
#include "matrix.hh"
#include "mechanics.hh"



#include "vector.hh"

#define ELEMENT_STATE Element_state2
#define elastic_energy veronda_elastic_energy2
#define elastic_force veronda_elastic_force2
#define elastic_force_derivative veronda_elastic_force_derivative2
#define MATRIX Matrix2
#define DIM 2

#include "veronda-material.icc"

#undef ELEMENT_STATE
#undef elastic_force
#undef elastic_force_derivative
#undef elastic_energy
#undef DIM
#undef MATRIX


#define ELEMENT_STATE Element_state3
#define elastic_energy veronda_elastic_energy3
#define elastic_force veronda_elastic_force3
#define elastic_force_derivative veronda_elastic_force_derivative3
#define MATRIX Matrix3
#define DIM 3

#include "veronda-material.icc"

IMPLEMENT_GET_ELASTICITY(veronda);

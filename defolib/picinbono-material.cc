/*
  Compute elastic reactions using

    W(C) = lambda tr(C)^2 + mu/2 tr(C^2)
   D_C(W) = mu (C : . )  + lambda tr(C) (I  :  . )
      =   A :  (.)

   with

      A = (mu C + lambda I1 Id)
     
   
  this gives
  
    Sigma = lambda  tr(E) + 2 mu E
    E = 1/2 ((I + Grad U )^T (I + Grad U) - I)

  This is linear material combined with nonlinear geometry.


  For small deformations, you get more numerical precision since
  F = I + grad U, where grad U = Oh(Eps), and I = Oh(1).

  Working with E, E = Oh(Eps), and you loose less significant digits. 
  
 */

#include <stdio.h>
#include <assert.h>

#include "element-state.hh"
#include "mechanics.hh"




#define ELEMENT_STATE Element_state2
#define elastic_energy picinbono_elastic_energy2
#define elastic_force picinbono_elastic_force2
#define elastic_force_derivative picinbono_elastic_force_derivative2
#define MATRIX Matrix2
#define DIM 2

#include "st-venant-kirchoff.icc"

#undef ELEMENT_STATE
#undef elastic_force
#undef elastic_force_derivative
#undef elastic_energy
#undef DIM
#undef MATRIX


#define ELEMENT_STATE Element_state3
#define elastic_energy picinbono_elastic_energy3
#define elastic_force picinbono_elastic_force3
#define elastic_force_derivative picinbono_elastic_force_derivative3
#define MATRIX Matrix3
#define DIM 3

#include "st-venant-kirchoff.icc"

IMPLEMENT_GET_ELASTICITY(picinbono);

/*
  This file contains the continuum mechanics formulas.  In terms of
  flop count, it is suboptimal; also -- it uses lots of local 3x3
  matrix variables, leading to bad cache behavior ?

  Possible optimizations: lots of the tensorfields are symmetric, so we
  can cut back on both space (cache!) and flop requirements for those.

  Possible optimizations: we can cut back on transposed matrices, if we
  use multiplication routines that do transposition internally.
  
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mesh-feature.hh"
#include "element-state.hh"
#include "setting.hh"
#include "node.hh"
#include "mechanics.hh"
#include "geometry.hh"
#include "deformation-state.hh"


/*
  
  We read lambda and mu only once from get_number_setting () to save costly
  calls to get_number_setting ().

  TODO: should put these in Element_state, or better yet, make
  them tunable from outside?
  
*/

Real lame_lambda = 0.0;   
Real lame_mu = 0.0;
Real lame_gamma =0.0;

void
init_lame_parameters ()
{
  Real poisson =get_number_setting ("poisson");
  Real young = get_number_setting ("young");      

  lame_mu =  young / (2 * (1 + poisson));
  lame_lambda =  poisson * young /((1  +poisson) *(1 - 2*poisson));
  lame_gamma = get_number_setting ("material-nonlinearity");

  if (get_bool_setting ("plane-stress"))
    {
      lame_lambda = 2 * lame_lambda * lame_mu/ (2*lame_mu  +lame_lambda);
    }

  printf ("Elasticity: lambda =%lf, mu = %lf, gamma=%lf ", lame_lambda, lame_mu, lame_gamma );
}

Element_state*
get_element_state (int dim)
{
  if (dim == 3)
    return new Element_state3;
  else if (dim == 2)
    return new Element_state2;
  else
    assert (false);
}



#define DIM 2
#define MATRIX Matrix2
#define ELEMENT_STATE Element_state2
#define ELASTICITY_FUNCS Elasticity_functions2
#define REFERENCE_LOCATION reference_location2
#define MINMAX_EDGE_LENGTH minmax_edge_length2
#include "element-state.icc"

#undef DIM
#undef MATRIX 
#undef ELEMENT_STATE
#undef ELASTICITY_FUNCS
#undef REFERENCE_LOCATION
#undef MINMAX_EDGE_LENGTH

#define DIM 3
#define MATRIX Matrix3
#define ELEMENT_STATE Element_state3
#define ELASTICITY_FUNCS Elasticity_functions3
#define REFERENCE_LOCATION reference_location3
#define MINMAX_EDGE_LENGTH minmax_edge_length3

#include "element-state.icc"


Elasticity_functions::Elasticity_functions()
{
  dim_ = -1;			// todo
  linear_ = false;
}


void
Element_state::precompute (Element*, Deformation_state const*){} 
Real
Element_state::linear_elastic_energy (Real const *) const{ return -1;}  
Real
Element_state::elastic_energy (Elasticity_functions const*, Real const *) const{ return -1;}
int
Element_state::elastic_force (Elasticity_functions const*, Real*,Real const*) const{return 0;}
int
Element_state::elastic_force_derivative (Elasticity_functions const*, Real*,Real*, Real const*,Real const*) const
{return 0;}

Element_state::~Element_state()
{
}

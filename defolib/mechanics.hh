#ifndef MECHANICS_HH
#define MECHANICS_HH

Real linear_elastic_energy (Real i1, Real i2, Real i3);
int linear_elastic_force (Element_state const*, Matrix2* dest, Matrix2 const & src);


#define IMPLEMENT_GET_ELASTICITY(NAME) \
static Elasticity_functions **elasticity_functions;			\
Elasticity_functions*							\
get_ ##  NAME ## _elasticity (int d)					\
{									\
  if (!elasticity_functions)						\
    {									\
      elasticity_functions = new (Elasticity_functions*)[2];		\
									\
      Elasticity_functions3 *e3 = new Elasticity_functions3;		\
      e3->force_function_ = &NAME ## _elastic_force3;			\
      e3->dforce_function_ = &NAME ## _elastic_force_derivative3;	\
      e3->energy_function_ = &NAME ## _elastic_energy3;			\
									\
      Elasticity_functions2 *e2 = new Elasticity_functions2;		\
      e2->force_function_ = &NAME ## _elastic_force2;			\
      e2->dforce_function_ = &NAME ## _elastic_force_derivative2;	\
      e2->energy_function_ = &NAME ## _elastic_energy2;			\
      									\
      elasticity_functions[1] = e3;					\
      elasticity_functions[0] = e2;					\
    }									\
									\
  return elasticity_functions[d -2];					\
}

Elasticity_functions * get_linear_elasticity (int dim);
Elasticity_functions * get_picinbono2_elasticity (int dim);
Elasticity_functions * get_picinbono_elasticity (int dim);
Elasticity_functions * get_neo_hookean_elasticity (int dim);
Elasticity_functions * get_neohookean_elasticity (int dim);
Elasticity_functions * get_veronda_elasticity (int dim);
Elasticity_functions * get_consistent_veronda_elasticity (int dim);


extern Real lame_lambda, lame_mu, lame_alpha, lame_gamma;

void init_lame_parameters();

#endif

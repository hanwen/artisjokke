#include <gsl/gsl_integration.h>
#include <stdio.h>
#include "friction.hh"
#include "big-vector.hh"
#include "vector-io.hh"


  /*
    All these hard coded numbers are inspired by DiMaio.
  */
Friction_definition::Friction_definition()
{
  peak_maximum_param = 0.0075; 
  peak_length = 0.02;
  peak_extra_force = 14.0/48.0;
  force_density = 48;
  entry_decrease_distance = 0.0075; // force at needle entry.
}


struct Friction_params {
  Real min_p;
  Real max_p;
  Real a,b;
  
  Friction_definition const * friction;

  Real friction_density (Real);
};



/*
  Density function , indepedent of discretization
 */
Real
Friction_params::friction_density (Real x)
{
  Real f = friction->force_density;
  if (x < min_p + friction->peak_maximum_param)
    {
      f *= (1.0 + friction-> peak_extra_force
	    * (x - min_p)  / (friction->peak_maximum_param  - min_p));
    }
  else if (x < min_p + friction-> peak_length)
    {
      f *= (1.0 +  friction->peak_extra_force
	    * ((min_p + friction-> peak_length) - x) / ( friction->peak_length -  friction->peak_maximum_param));
    }
  
  if (x > max_p - friction->entry_decrease_distance)
    {
      f *= (max_p - x ) / friction->entry_decrease_distance;
    }

  return f;
}


/*
      Dump forces, so we can verify continuity.
*/
void
graph_friction_forces(Friction_definition const *fr)
{
  Friction_params p ;
  p.friction = fr;
  p.min_p =0.0;
  p.max_p = 0.07;
  
  FILE * graph = xfopen ("force-graph.txt", "w");
  for (int i= 0; i < 100; i++)
    {
      Real x = (i / 99.0) * 0.07;
      fprintf (graph, "%lg %lg\n", x, p.friction_density (x));
    }
  fclose (graph);
}

/*
  left func: interval left of peak
 */
double
left_friction_func (double x, void * params)
{
  Friction_params * p = (Friction_params*)params;
  assert (x <= p->b &&   x >= p->a);
  Real unit_func =((x- p->a) / (p->b - p->a));
  return p->friction_density (x) * unit_func; 
}

/*
  right func: interval right of peak
 */
double
right_friction_func (double x, void * params)
{
  Friction_params * p = (Friction_params*)params;
  assert (x <= p->b && x >= p->a);
  Real unit_func = ((p->b - x) / (p->b - p->a));
  return p-> friction_density (x) * unit_func;
}

/*
  Integrate friction function (given by min_p, max_p, friction)
  multiplied with a linear interpolation from either (a,0) to (b,1)
  (dir = LEFT = -1), or (a,1) to (b,0) (dir = RIGHT = 1), over [a,b].
 */

Real
integrate_friction_interval (Friction_params const * p, int dir)
{
  
  if (p->b < p->max_p - p->friction->entry_decrease_distance
      && p->a > p->min_p + p->friction->peak_length)
    {
      return (p->b - p->a)*  p->friction->force_density * 0.5;
    }
  else
    {
      /*
	This is overkill, but we're lazy, and getting the details of
	the integration over the boundaries right is the kind of
	trickyness I don't want to get involved with on a sunday
	afternoon.
      */
      gsl_function F;
      F.function = (dir == 1) ? right_friction_func :left_friction_func;
      F.params = (void*)p;

      size_t neval;
      double err;
      double res;
      double epsrel = 0.01;	// what does this mean? ?
      double epsabs = 0.01 *  (p->b-p->a) * p->friction->force_density;
#if 0 
      int succ = gsl_integration_qng (&F, p->a, p->b,
				      epsabs, epsrel,
				      &res, &err, &neval);
#else
      static gsl_integration_workspace *wsp;

      int limit = 14;
      if (!wsp)
	wsp = gsl_integration_workspace_alloc (limit);
      int succ = gsl_integration_qag (&F, p->a, p->b,
				      epsabs, epsrel, limit,
				      GSL_INTEG_GAUSS21,
				      wsp,
				      &res, &err);
      
#endif
      assert (!succ);
      return res;
    }
}

Array<Real>
compute_friction_forces (Array<Real> needle_params,
			 Friction_definition const* friction)
{
  Array<Real> forces;
  completize_big_vector (&forces, needle_params.size());


  /*
  NEEDLE_PARAMS is ascending.
  */
  Friction_params p;
  p.friction = friction;

  p.min_p = needle_params[0];
  p.max_p = needle_params.top();
  
  for (int i = 0; i <  needle_params.size()-1; i++)
    {
      p.a = needle_params[i];
      p.b = needle_params[i+1];

      Real sign = 1.;
      if (p.a > p.b)
	{
	  p.a = needle_params[i+1];
	  p.b = needle_params[i];
	  sign = -1;
	}
      forces[i] += sign * integrate_friction_interval (&p,  1);
      forces[i+1] += sign * integrate_friction_interval (&p,  -1);      
    }
  
  return forces;
}




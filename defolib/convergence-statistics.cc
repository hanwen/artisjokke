#include <math.h>

#include "mesh-connectivity.hh"
#include "vector-io.hh"
#include "deformation-state.hh"
#include "setting.hh"
#include "convergence-statistics.hh"
#include "big-vector.hh"

/*
really, too ugh for words.
*/

/*
  TODO: should flush file.
 */
static FILE * stat_file = 0;

Real reference_state_distance (Deformation_state * def, Deformation_state * ref);
Real reference_state_energy_distance (Deformation_state * def , Deformation_state *ref);


Deformation_state* 
read_reference_state (Mesh_connectivity*top) 
{
  char const *fn = get_string_setting ("reference-state-file");

  Deformation_state *ref = new Deformation_state(top->dimension ());

  bool succ = read_deformation_state (ref, fn);
  if (!succ)
    {
      ref->completize_nodal_arrays (top->dimension()*top->node_count () +  1);
    }
  else
    {
      top->add_watcher (ref);	  
    }
  return ref;
}

/*
note: all flop counts are off for 2D.
*/

void
flop_stats (Deformation_state * def, Deformation_state *ref) 
{
  long long remember_e = element_flop_count ;
  long long remember_v = vector_flop_count ;

  if (!stat_file)
    {
      
      stat_file = xfopen ("convergence-statistics.txt", "w");
      if (!stat_file)
	fprintf (stderr ,"Can't open statistics file. Exiting.\n");

      write_invocation_header (stat_file, "# ");
      fprintf (stat_file, "# FORMAT: MFLOPs\t RES(rel)\tERR(max)\tERR(energy)\tENERGY\n");
    }

  long long fc = flop_count();
  printf ("MFLOPS: %f",  fc * 1e-6);


  Real compare = def->scale_free_force_comparison();

  Real res =sqrt (def->get_residual_len_sq ());
  Real relative_res = res / compare;
  Real dsd = reference_state_distance (def, ref);
  Real esd = reference_state_energy_distance (def, ref);  
  printf(" REL_FORCE_LEN: %g (abs: %lf) SOL_DIFFERENCE: %g (energy norm %g)", relative_res, res, dsd, esd);
  
  Real e  = def->virtual_energy (def->displacements_.access_array ());

  printf( " ENERGY: %f\n", e);

  if (stat_file)
    {
      fprintf (stat_file, "%f\t%g\t%g\t%g\t%g\n", fc * 1.0e-6, relative_res, dsd, esd, e);
      fflush (stat_file);
    }
  
  element_flop_count = remember_e;
  vector_flop_count = remember_v;
}

Real
reference_state_energy_distance (Deformation_state * def , Deformation_state *ref)
{
  if (ref->dimension () != def->dimension ())
    {
      return -1.0;
    }
  else
    {
      Real *stat1 = def->statistic_vector1_.unsafe_access_array();
      Real *stat2 = def->statistic_vector2_.unsafe_access_array();

      Real const *exact = ref->displacements_.access_array ();

      int n = def->dimension();
      
      big_vector_subtract (stat1, def->displacements_.access_array (), exact, n);
      def->apply_force_derivative_and_residual (stat2, 0, 0, exact, stat1);

      Real esd = - big_vector_ip (stat1, stat2, n);
      return esd;
    }
}


/*
  Compare PSTATE to reference state, returning the euclidian distance. 
 */
Real
reference_state_distance (Deformation_state * def, Deformation_state * ref)
{
  if (ref->dimension () != def->dimension ())
    {
      return -1.0;
    }
  else

    /*
      Used to have distance in 2-norm, but infty norm is dimension
      independent.
     */
    return big_vector_max_distance (ref->displacements_.access_array (),
				    def->displacements_.access_array (),
				    def->dimension ());
}

long long
flop_count ()
{
  return element_flop_count + vector_flop_count;
}

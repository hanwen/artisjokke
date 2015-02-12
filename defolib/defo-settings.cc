#include "setting.hh"
#include "mechanics.hh"

static
Number_option_list_entry init_alist[] = {
  {"calibration-time", "Time to use for calibration in the start. Negative means run till solution. \n0.0 one iteration per frame.", -1},
  {"density", "material density for mass-lumping", 1e3},
  {"incremental-tolerance", "Tolerance for the inner loop of an incremental method.", 0.25},
  {"material-nonlinearity", "Exponent factor for non-linear material", 8.0},
  {"max-cg-iterations", "Maximum number of CG iterations (0 = unlimited)", 0.0},
  {"maximum-iterations", "maximum number of iterations",0.0},
  {"maximum-mflops", "Print speed statistics per MFLOP, and exit after X mflops", 0.0},
  {"outer-loop-tolerance", "Tolerance for the outer loop of an CG or Newton iteration. We stop of |residual| < TOL * |force|", 0.01},
  {"residual-rounding-tolerance", "Make bigger to prevent excessive computations for small external forces.",  1e-8},    
  {"poisson","Poisson's ratio", 0.34},
  {"time-step", "Time step relative to critical TS.", 0.8},
  {"viscosity", "Spring constant: viscosity, relative to mass", 0.02},
  {"young", "Young modulus (times slab thickness)", 10e3},
  {"nonlinear-cg-newton-tolerance", "Error tolerance for newton solution in nonlinear CG", 1e-2},
  {"nonlinear-cg-restart-count", "Restart after how many iterations? (negative = N, 0 = no restart)", 0},
  {0,0,0},
};
    

static
String_option_list_entry init_str_alist[] = {
  {"problem-type", "{dynamic,static}", "static"},
  {"time-integration", "Type of time-integration {eulerfw,RK4,ss22,ms22}", "ss22"},
  {"elasticity", "Which kind of elasticity:\n"
   "\n"
   "\tlinear -- standard Hooke's law\n"
   "\tnonlinear -- linear stress/strain, exact geometry\n"
   "\tsprings -- simple springs\n"
   "\tneohooke -- compressible neo-hookean\n"
   "\tveronda -- compressible nonlinear Veronda material\n"
   "\tconsistent-veronda -- compressible nonlinear Veronda material, with consistent linearization\n", "linear"},
  {"viscosity", "Which velocity to use for viscosity calculation. {spring,node}", "node"},
  {"relaxation-type", "Which optimization routines? {linear,nonlinear,incremental-newton}", "linear"},
  {"line-search", "Which method for 1D optimization? {newton,secant}", "secant"},
  {"nonlinear-cg-betafunc", "Which method for selecting BETA? {fletcher-reeves,polak-ribiere,both}", "fletcher-reeves"},
  {"reference-state-file", "file with solution state", "reference.state"},
  {"ignore-degenerate-elements", "Assume degen elements don't exist", "no"},
  {"exact-residual-update", "Recompute residual from scratch in CG iteration?", "no"},
  {"compute-exact-solution", "Use NL-CG as a starting solution for Newton Raphson process.\n", "no"},
  {"plane-stress", "Use plane stress iso. plane strain", "yes"},
  {0,0,0},
};

void
init_deformation_module ()
{
  init_settings (init_alist, init_str_alist);
}

void
init_after_settings ()
{
  init_lame_parameters ();
}

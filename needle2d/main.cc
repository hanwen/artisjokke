#include <stdio.h>
#include <getopt.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <locale.h>

#include "proto.hh"
#include "mesh-connectivity.hh"
#include "test-mesh.hh"
#include "artisjokke-drawer2.hh"
#include "glut-interface.hh"
#include "debug.hh"
#include "maubach-tree.hh"
#include "setting.hh"
#include "deformation-state.hh"
#include "bottom-fix-hook.hh"
#include "needle-inserter.hh"


static
Number_option_list_entry init_alist[] = {
  {"auto-insert-speed", "Insertion speed (per solution step) of the needle (measured in h-ref)",
   0.1},
  {"auto-insert-y", "Y position of the needle handle", 0.05},
  {"auto-insert-x", "X position of the needle handle", -0.01},
  {"auto-insert-angle", "angle of insert (in degrees)", 0.0},  
  {"auto-insert-depth", "maximum movement of needle", 0.07},
  {"dynamic-friction-factor", "How much lower are slipping forces?", 0.5},
  {"init-force", "add force to initial configuration", 0.0},
  {"initial-h" ,  "Make mesh this fine initially", 0.1},
  {"refinement-h", "Maximum length during refinement meshing ", 2.5e-3},
  {0, 0, 0},

#if 0 
  /*
    These are related to old code...
   */
  {"window-position", "The initial position of the OpenGL window. Negative: don't set" , -1},
  {"window-size", "The initial size of the OpenGL window.", 400},
  {"cube-size", "Size of init cube" , 11},
  {"gravity", "Gravity acceleration", 0.0},


  {"minimal-component-volume", "Delete tetrahedral components with relative volume less than this.", 1e-3},
  {"collapse-degeneracy-passes", "Maximum number of passes for trying to collapse degenerate tetrahedra", 100.0},
  {"degeneracy-stats", "If set, show extensive dumps of degenerate nodes.", 0.0},
   {"improvement-aspect-ratio-threshold", "Aspect ratio threshold below which to subdivide tetras after reading a tet model",0.0},
  {"dihedral-angle-threshold", "Consider dihedral angles of a tetra to be big when larger than this", 0.99},  
  {"eye-distance","eye-distance for stereo viewing.", 0.0},
  {"scalpel-angle" , "angle of default scalpel  configuration (in PI radians)", 0.0},
  {"scalpel-length", "length of scalpel (controls depth of penetration)", 3.0},
  {"scalpel-still-x-angle",  "Rotation for the still point of the scalpel (in PI radians)" ,  0.0},
  {"scalpel-still-y-angle",  "Rotation for the still point of the scalpel (in PI radians)" , 0.0},
  {"scalpel-velocity", "velocity of scalpel", 0.1},
  {"scalpel-x-start",  "X-offset of scalpel still point", 0.0},
  {"tetra-aspect-ratio-threshold", "Consider tetrahedra  with aspect ratio lower to be degenerate", 0.01},
  {"triangle-aspect-ratio-threshold", "Consider triangles  with aspect ratio lower to be degenerate", 0.01},
  {"update-factor","change relaxation rate compared to frame rate",1.0},

  /****************************************************************/
  
  {"do-reposition", "Perform node repositioning?", "yes"},
  {"force-apply-control", "Set initial force from where? {init,gui}", "init"},
  {"0-volume-read", "Try to import 0-volume tetrahedra against all odds", "no"},

#endif

};

static
String_option_list_entry init_str_alist[] = {
  {"normalize-input", "Center and scale model after reading", "no"},
  {"normal-orientation", "orientation of outward normal. {cw,ccw} (arbitrary!)", "ccw"},
  {"fix-planes", "Which planes of the object to fix. Multiple allowed", "down"},
  {"init-write-file", "Write initial configuration to which file? (none)", ""},
  {"relocate-nodes", "Relocate needle nodes?", "no"},
  {"print-mesh", "Print mesh when scenario finishes successfully?", "no"},
  {0,0,0},  
};





void
fpe_handler (int)
{
  fprintf (stderr,"FP exception\n");
  exit (2);
}

void
help ()
{
  fprintf (stderr,
	   "Usage: needle [options]\n"\
	   "  -h      help\n"\
	   "  -oKEY=VALUE \n"
	   "          set option (use -ohelp for option help)\n"
	   "  -SSCEN  run test scenario SCEN\n"
	   "  -d      print gobs of debugging info.\n"\
	   "\n"\
	   );
}

void notice ()
{
  fprintf (stdout, "\nNeedle: 2d needle insertion. Version %s"
	   "\nThis software is copyright (c) Utrecht University 2002.\n"
	   "All rights reserved.\n"
	   "Report bugs to Han-Wen Nienhuys <hanwen@cs.uu.nl>\n",
	   __DATE__ );
}


extern void test_stuff ();


int argument_count;
char * *argument_value ;
char const * scenario = NULL;

bool debug_out = false;

void
my_main (int ac, char *av[])
{
  /*
    Be a little paranoid as Locale settings may mess with the numeric stuff.
   */
  setlocale (LC_NUMERIC, "POSIX");
  
  argument_value = av;
  argument_count = ac;
  signal (SIGFPE, fpe_handler);

  init_deformation_module ();
  init_settings (init_alist, init_str_alist);
  set_number_setting ("young", 34e3* 0.010);

  
  int c;
  while ((c = getopt (ac, av, "S:o:hs:t:dw:")) != -1)
    {
      switch (c)
	{
	case 'd':
	  debug_out = true;
	  break;
	  
	case 'h':
	  help ();
          notice ();
	  exit (0);
	  break;
	  
	case 'S':
	  scenario = strdup (optarg);
	  break;

	case 'o':
	  parse_setting (optarg);
	  break;
	  
	default:		// fixme
	  log_message ("unknown option `%ca\'\n", c);
	  exit (2);
	  break;
	}
    }

  notice ();

  init_after_settings();
  test_stuff ();
}

int
main (int argc, char*argv[])
{
  my_main (argc, argv);
  Maubach_tree* fem = new Maubach_tree ();
  Real size = 0.1;
  
  Deformation_state * def = make_new_deformation_state (fem);
  fem->set_geometry (def, size);

  fprintf (stderr, "Refining ...");
  fflush (stderr);
  refine_uniformly2 (fem, def, 0.1 /get_number_setting ("initial-h" ));
  fprintf (stderr, " done...\n");
  
  Needle_inserter * ins  = new Needle_inserter (dynamic_cast<Maubach_tree*>(fem),
						def);
  Bottom_fix_watcher * bf = new Bottom_fix_watcher (def);
  if (scenario)
    {
      test_scenario (scenario, fem, def, ins);
#ifndef OPENGL
      exit (0);
#endif
    }
  
#ifdef OPENGL
  Artisjokke_drawer2 * md = new Artisjokke_drawer2 (fem, def);
  
  md->needle_inserter_ = ins;
  glut_setup (md, argc, argv);
#else
  while (1)
    {
      def->simulation_body ();
    }
#endif
  //  apply_gravity (def, X_AXIS, -1);
}



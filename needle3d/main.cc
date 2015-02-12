#include <stdio.h>
#include <getopt.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <locale.h>

#include "defo-proto.hh"
#include "mesh-connectivity.hh"
#include "artisjokke-drawer3.hh"
#include "glut-interface.hh"
#include "debug.hh"
#include "setting.hh"
#include "deformation-state.hh"
#include "bottom-fix-hook.hh"
#include "misc3d.hh"
#include "maubach-tree3.hh"
#include "needle-inserter3.hh"
#include "scenario.hh"

static
Number_option_list_entry init_alist[] = {
  {"auto-insert-speed", "Insertion speed (per solution step) of the needle (measured in h-ref)",
   0.1},
  {"auto-insert-x", "X position of the needle handle", -0.01},
  {"auto-insert-y", "Y position of the needle handle", 0.07001},
  {"auto-insert-z", "X position of the needle handle", 0.05001},  
  {"auto-insert-angle", "angle of insert (in degrees)", 0.0},  
  {"auto-insert-depth", "maximum movement of needle", 0.07},
  {"dynamic-friction-factor", "How much lower are slipping forces?", 0.5},
  {"friction-density" , "How much friction along needle (in N/m)", 50}, 
  {"init-force", "add force to initial configuration", 0.0},
  {"initial-level" ,  "Make mesh this fine initially", 0.1},
  {"gravity" , "Gravity acceleration" , 9.8},
  {"refinement-level", "Maximum length during refinement meshing ", 16},
  {"needle-radius", "Radius of needle object", 0.001},
  {"young", "Young modulus", 34e3},
  {0, 0, 0},
};

static
String_option_list_entry init_str_alist[] = {
  {"normalize-input", "Center and scale model after reading", "no"},
  {"normal-orientation", "orientation of outward normal. {cw,ccw} (arbitrary!)", "ccw"},
  {"fix-planes", "Which planes of the object to fix. Multiple allowed", "down"},
  {"init-write-file", "Write initial configuration to which file? (none)", ""},
  {"relocate-nodes", "Relocate needle nodes?", "no"},
  {"print-mesh", "Print mesh when scenario finishes successfully?", "no"},
  {"side-view", "Don't rotate view", "no"},  
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
  fprintf (stdout, "\nNeedle: 3d needle insertion. Version %s"
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
  init_after_settings ();
 }

int
main (int argc, char*argv[])
{
  my_main (argc, argv);

  Maubach_tree3 * fem = new Maubach_tree3 ();
  Real size = 0.1;

  Deformation_state * def = make_new_deformation_state (fem);
  fem->set_geometry (def, size);

  refine_uniformly3 (fem, (int) get_number_setting ("initial-level"),
		     def);
  
  Bottom_fix_watcher * bf = new Bottom_fix_watcher (def);
  Needle_inserter3 * ni = new Needle_inserter3 (fem, def);
  
#if 1
  if (scenario)
    {
      test_scenario (scenario, fem, def, ni);
#ifndef OPENGL
      exit (0);
#endif
    }
#endif
  
#ifdef OPENGL
  Artisjokke_drawer3 * md = new Artisjokke_drawer3 (fem, def, ni );
    
  glut_setup (md, argc, argv);
#else
  while (1)
    {
      def->simulation_body ();
    }
#endif
  //  apply_gravity (def, X_AXIS, -1);
}



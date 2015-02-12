#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "misc2d.hh"
#include "auto-insert.hh"
#include "mesh-connectivity.hh"
#include "mesh-feature.hh"
#include "vector.hh"
#include "node.hh"
#include "debug.hh"
#include "geometry2.hh"
#include "vector-io.hh"
#include "setting.hh"
#include "maubach-tree.hh"
#include "deformation-state.hh"
#include "friction.hh"
#include "big-vector.hh"
#include "convergence-statistics.hh"
#include "needle-inserter.hh"


/*
  Generate nodal forces , and apply to the mesh.
*/
void
node_forces_scenario (Maubach_tree * tree,
		      Deformation_state *def,
		      Needle_inserter * )
{
  Real depth = 0.065 ;
 
  Real h = get_number_setting ("refinement-h")
    <? get_number_setting ("initial-h");

  Real x = 0.0;
  int n = int (depth / h);

  Link_array<Node> nods;  
  for (int i = 0; i < n ; i++)
    {
      Vector2 lookup_loc = Vector2 (x, 0.05);

      Real ref_h =sqrt(2) * h * 1.01;
      refine_around2 (tree, ref_h,  lookup_loc, &reference_location2, def);
      Element * e = tree->locate (lookup_loc, &reference_location2, def);

      Node * nod = 0;
      Real dist = 1e8;
	
      for (int j = 0 ; !nod && j < 3; j++)
	{
	  Vector2 rl = reference_location2 (e->node (j), def);

	  if (euclidian_distance (rl, lookup_loc) < dist
	      && fabs (rl(Y_AXIS) - lookup_loc (Y_AXIS)) < .01 *h)
	    {
	      dist = euclidian_distance (rl, lookup_loc);
	      nod = e->node(j);
	    }
	}

      if (nods.empty ()
	  || nods.top () != nod)
	nods.push (nod);
      x += h;
    }

  Array<Real> pars;
  for (int i = nods.size() ; i--;)
    {
      Node * n = nods[i];
      Real p = depth - reference_location2 (n, def)(X_AXIS);
      pars.push(p);
    }

  Friction_definition  fr;

  Array<Real> forces = compute_friction_forces (pars, &fr);

  for (int i= forces.size(); i--;)
    def->add_force (nods[i], Vector2(forces[i], 0));

  {
    FILE * f = xfopen ("static-needle-node-forces.txt", "w");
    write_invocation_header (f, "# ");
    for (int i= pars.size(); i--;)
      {
	fprintf (f,"%lf %lf\n", pars[i], forces[i]);
      }

    fclose (f);
  }

  
  def->check_for_topology_update();
    
  Real tot_force = big_vector_1_norm (forces.access_array(), forces.size());
  printf ("\nApplying force  %lf to %d nodes\n", tot_force, nods.size());

  def->update_forces ();
  while (!def->good_solution())
    {
      def->simulation_body ();
    }
  write_deformation_state (def, "needle-deformation.state");
      
  Deformation_state * ref = new Deformation_state(2);
  bool succ = read_deformation_state (ref, "static-deformation-reference.state");
  
  if (succ)
    plot_reference_state_distance (tree, def, ref, "error-plot");

  FILE *f  = xfopen ("iteration-statistics.txt", "w");
  write_invocation_header (f, "# ");

  /*
    Twice to ensure proper reads by R.
   */
  fprintf (f, "# ITERATION-COUNT MAX-ITERATION-COUNT\n");

  fprintf (f, "%d %d\n", global_iteration_count,
	   global_max_iteration_count);
  fclose (f);
}


void
angled_node_forces_scenario (Maubach_tree * tree,
			     Deformation_state *def,
			     Needle_inserter * )
{
  Auto_needle_insert * ni = new Auto_needle_insert();
 
  Real h = get_number_setting ("refinement-h")
    <? get_number_setting ("initial-h");
  Real ref_h =sqrt(2) * h * 1.01;
  Real d = 0.0;
  Vector2 tip = ni->handle_ + ni->max_depth_ * ni->dir_;
  bool relocate = get_bool_setting ("relocate-nodes");
  
  while (d < ni->max_depth_)
    {
      d +=  h/2;

      Vector2 x = ni->handle_ + d * ni->dir_;
      refine_around2 (tree, ref_h,  x, &reference_location2, def);
    }

  Array<Element_intersection> ei
    = track_line_through_deformed_mesh (tree, def,
					tip,
					ni->handle_);

  Link_array<Node> nods;
  Array<Real> params;
  for (int i =0 ; i < ei.size (); i++)
    {
      Face  *f = ei[i].exit_face_;
      Real m  = ei[i].exit_param_;

      Vector2 intercept
	= m * reference_location2 (ei[i].exit_face_->node(0), def)
	+ (1-m) * reference_location2 (ei[i].exit_face_->node(1), def);
      
      Node * n = f->node (m > 0.5 ? 0 : 1);
      if (nods.size() && n == nods.top())
	continue;
      
      if (relocate)
	{
	  def->set_reference_location (n, intercept.elts_);
	}
      
      nods.push (n);
      params.push (- ni->dir_ *(deformed_location2 (n,  def) - tip));
    }

  Friction_definition  fr;
  Array<Real>  forces = compute_friction_forces (params, &fr);
  printf ("\nApplying force  %lf to %d nodes\n", 
	  big_vector_1_norm (forces.access_array(), forces.size()),
	  forces.size());


  for (int i = 0; i < nods.size (); i++)
    {
      def->constraints_.change_linear_movement_constraint (nods[i], ni->dir_, 0);
      def->add_force (nods[i], ni->dir_ * forces[i]);
    }

  def->check_for_topology_update();
  def->update_forces ();
  while (!def->good_solution())
    {
      def->simulation_body ();
    }
  write_deformation_state (def, "needle-deformation.state");
      
  Deformation_state * ref = new Deformation_state(2);
  bool succ = read_deformation_state (ref, "angled-deformation-reference.state");
  
  if (succ)
    plot_reference_state_distance (tree, def, ref, "error-plot");

  FILE *f  = xfopen ("iteration-statistics.txt", "w");
  write_invocation_header (f, "# ");

  /*
    Twice to ensure proper reads by R.
   */
  fprintf (f, "# ITERATION-COUNT ITERATION-COUNT\n");

  fprintf (f, "%d %d\n", global_iteration_count,
	   global_iteration_count);
  fclose (f);
}


void
uniform_force (Maubach_tree * tree,
	       Deformation_state *def,
	       Needle_inserter * )
{
  Real traction = 10.0;
  Real total_force =00.0;  
  set<Face*>  fs = tree->boundary ();

  Node * nod = 0;
  for (iterof(i, fs); i != fs.end (); i++)
    {
      Face * f = (*i);
      Vector2 n = simplex_normal2 (f->simplex (), &reference_location2, def);

      if (n(X_AXIS) < 0
	  && fabs (n(Y_AXIS)) < 1e-6)
	{
	  Real l = edge_length2 (f->simplex (), &reference_location2, def);

	  Real force  =l * traction;
	  for (int j =  0; j<2; j++)
	    def->add_force (f->node(j), Vector2 (- force/2, 0));

	  nod = f->node(0) ;
	  total_force += force;
	}
    }
  printf ("\ntotal force %lf\n", total_force);

  def->update_forces();
  while (!def->good_solution())
    {
      def->simulation_body ();
    }

  printf ("Displacement : %lf\n",
	  (deformed_location2( nod, def) - reference_location2 (nod, def))(0));
}


clock_t needle_timer ;

void
finish_auto_insert (Maubach_tree * tree,
		    Deformation_state *def,
		    Needle_inserter * ins)
{
  clock_t dclock =   clock () -needle_timer ;

  Real dt = Real (dclock) / CLOCKS_PER_SEC;
  
  printf  ("\nTook %f seconds for simulating insertion\n", dt );
  printf("%d rearrangements, %f ms per rearrangement.\n", ins->rearrange_count () , 1e3 *dt/ins->rearrange_count ()); 

  FILE*ss =xfopen ("speed-stats.txt", "w");
  write_invocation_header (ss, "# " );
  fprintf (ss, "# TOTAL-TIME REARRANGE-COUNT TOTAL-CG-ITER-COUNT MAX-CG-ITER-COUNT \n");
  fprintf (ss, "%lf %d %d %d\n", dt, ins->rearrange_count(),
	   global_iteration_count,
	   global_max_iteration_count
	   );
  fclose (ss);
  
  ins->dump_needle_forces  ("graph-needle-forces.txt");
  
  write_deformation_state (def, "needle-deformation.state");

  Deformation_state * ref = new Deformation_state(2);
  bool succ = read_deformation_state (ref, "insert-reference.state");

  if (succ)
    plot_reference_state_distance (tree, def, ref, "error-plot");

  mesh2d_print (tree, "deformed-mesh.eps", &deformed_location2, def, ins, 8);
  mesh2d_print (tree, "reference-mesh.eps", &reference_location2, def, ins, 10);
  mesh2d_print_boundary (tree, "deformed-boundary.eps", &deformed_location2, def, ins, 8);
  mesh2d_print_boundary (tree, "reference-boundary.eps", &reference_location2, def, ins, 10);


  Vector2 tiploc = ins->tip_reference_location();
  printf ("Needle tip ended at (%lf,%lf)\n",
	  tiploc(X_AXIS), tiploc(Y_AXIS));
  
}

void
auto_insert_scenario (Maubach_tree *,
		      Deformation_state *def,
		      Needle_inserter *ins)
{
  ins->auto_insert_ = new Auto_needle_insert();
  needle_timer = clock();

  ins->auto_needle_insert();
  
#ifndef OPENGL
  def->update_forces();
  while (!def->good_solution() || ins->auto_insert_)
    {
      def->simulation_body ();
    }
#endif
}



void
read_needle_dump (Maubach_tree * tree,
		  Deformation_state *def,
		  Needle_inserter * )
{
  Real h = get_number_setting ("refinement-h");
  Real ref_h =sqrt(2) * h * 1.01;
  
  char const *fn = "reference-needle-dump.txt";

  FILE *f = xfopen (fn, "r");

  while (1)
    {
      char * line = in_file_read (f);

      if (!line)
	break ;

      Real refx, defx, forcex;
      Real refy, defy, forcey;      
      
      sscanf (line, "%lf %lf %lf %lf %lf %lf\n", &refx, &refy, &defx, &defy, &forcex, &forcey );
      Vector2 loc (refx, refy);
      Vector2 force (forcex, forcey);
      
      refine_around2 (tree, ref_h, loc, &reference_location2, def);
      Element * e = tree->locate (loc, &reference_location2, def);

      Node * nod = 0;
      for (int j = 0 ; e && !nod && j < 3; j++)
        {
          Vector2 rl = reference_location2 (e->node (j), def);

          if (euclidian_distance (rl, loc) <   .01 * h)
            nod = e->node(j); 
        }

      if (!nod)
	{
	  printf ("Huh? Can't find node there. Check refinement settings.\n");  
	  exit (2);
	}
           
      def->add_force (nod, force);
    }
  fclose (f);

#ifndef OPENGL
  def->update_forces();
  while (!def->good_solution())
    {
      def->simulation_body ();
    }
  write_deformation_state (def, "needle-deformation.state");
#endif
}

  
void
file_forces (Maubach_tree * tree,
	     Deformation_state *def,
	     Needle_inserter * )
{
  Real h = get_number_setting ("refinement-h");
  Real ref_h =sqrt(2) * h * 1.01;
  
  char const *fn = "reference-force.txt";
      
  FILE *f = xfopen (fn, "r");

  char * line = in_file_read (f);
  int c;
  sscanf (line, "%d", &c);
  
  for (int i = 0; i < c; i++)
    {
      Vector2 loc = read_vector2 (f);
      Vector2 force = read_vector2 (f);
      Vector2 lookup_loc = loc;
      
      refine_around2 (tree, ref_h, lookup_loc, &reference_location2, def);
      Element * e = tree->locate (lookup_loc, &reference_location2, def);

      if (!e)
        return ; 
      Node * nod = 0;
      for (int j = 0 ; !nod && j < 3; j++)
        {
          Vector2 rl = reference_location2 (e->node (j), def);

          if (euclidian_distance (rl, loc) <   .01 * h)
            nod = e->node(j); 
        }
          
      def->add_force (nod, force);
    }
  fclose (f);
}


void
show_refinement (Maubach_tree * tree,
		 Deformation_state *def,
		 Needle_inserter * )
{
  mesh2d_print (tree, "maubach-no-refinement.eps", reference_location2, def, 0, 10);
  for (int i = 0; i < 10; i++)
    {
      Real h = 0.09 * pow(0.5, (i-1)/2.0);
      
      refine_around2 (tree, h, 
		     Vector2 (0.0702, 0.0351), &reference_location2, def);
      char s [1024];
      sprintf (s, "maubach-refinement-%d.eps", i);
      mesh2d_print (tree, s, reference_location2, def, 0, 10);
    }
  
}


void
test_scenario (char const *ch,
	       Maubach_tree * tree,
	       Deformation_state *def,
	       Needle_inserter * ins)
{
  if (!ch[0])
    return;

  if (0)
    {
      ;
    }
  else if (!strcmp(ch, "synthetic-node-forces"))
    {
      node_forces_scenario (tree, def, ins);
    }
  else if (!strcmp(ch, "angled-forces"))
    {
      angled_node_forces_scenario (tree, def, ins);
    }
  else if (!strcmp(ch, "show-refinement"))
    {
      show_refinement (tree, def, ins);
    }
  else if (!strcmp(ch, "file-node-forces"))
    {
      file_forces (tree, def , ins);
    }
  else if (!strcmp(ch, "read-needle-dump"))
    {
      read_needle_dump (tree, def , ins);
    }  
  else if (!strcmp(ch, "auto-insert"))
    {
      auto_insert_scenario (tree,def,ins);
    }
  else if (!strcmp(ch, "graph-friction"))
    {
      Friction_definition fr ; 
      graph_friction_forces (&fr);
    }
  else if (!strcmp (ch, "uniform-force"))
    {
      uniform_force (tree, def, ins); 
    }
  else
    {
      printf ("Unknown scenario %s, try\n{synthetic-node-forces,angled-forces,auto-insert,graph-friction,uniform-force}\n", ch);
      exit (2);
    }
      
  if (get_bool_setting ("print-mesh"))
    {
      mesh2d_print (tree, "reference-mesh.eps", &reference_location2, def,
		    ins,
		    10);
      mesh2d_print (tree, "deformed-mesh.eps", &deformed_location2, def,
		    ins,
		    10);

      mesh2d_print_boundary (tree, "reference-boundary.eps", &reference_location2, def,
		    ins,
		    10);
      mesh2d_print_boundary (tree, "deformed-boundary.eps", &deformed_location2, def,
		    ins,
		    10);
    }
}



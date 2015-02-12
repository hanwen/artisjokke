#include <stdio.h>
#include <stdlib.h>

#include "vector-io.hh"
#include "deformation-constraints.hh"

void
Deformation_constraints::write_to_file (char const * fn)const
{

  FILE * f = xfopen (fn, "w");

  fprintf (f,"%d\n", position_constraints_.size());
  for (iterof (i, position_constraints_); i != position_constraints_.end (); i++)
    {
        fprintf (f, "%d ", *i);
    }

  fprintf (f,"%d\n", line_constraints_.size());
  for (iterof (i, line_constraints_); i!=line_constraints_.end (); i++)
    {
      Vector2 v = i->second; 
      fprintf (f, "%d %lf %lf %lf\n", i->first, v(0), v(1), v(2));
    }

  fclose (f);
}

bool
Deformation_constraints::read_from_file (char const * fn)
{
  FILE * f = fopen (fn, "r");
  if (!f)
    {
      fprintf (stderr, "State file `%s' not found. Not reading deformation constraints.\n", fn);
      return false;
    }
  else
    {
      fprintf (stderr, "Reading `%s' ...\n", fn);
    }

  int dim = 0;
  fscanf (f,"%d\n", &dim);
  
  int constraints;
  fscanf (f,"%d\n", &constraints);

  position_constraints_.clear ();
  for (int i=0; i < constraints; i++)
    {
      int k;
  
      fscanf (f,"%d ", &k);
      position_constraints_.insert (k);
    }

  line_constraints_.clear ();
  for (int i=0; i < constraints; i++)
    {
      Real *rv = new Real [dim];

      char * line = in_file_read (f);
      int k = strtol (line, &line, 10);
      for (int i = 0; dim; i++)
	{
	  rv [i]= strtod (line , &line);
	}
      
      line_constraints_[k] = rv;
    }

  
  fclose (f);
  return true;
}

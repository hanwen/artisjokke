#include <ctype.h>
#include <stdio.h>


#include "deformation-state.hh"
#include "big-vector.hh"
#include "vector-io.hh"

void
write_big_vector (FILE * f,  Real const *a , int n)
{
  while (n--)
    {
      fprintf (f, "%a ", *a++);
    }
  fputc ('\n', f);
}

void
read_big_vector (FILE* f, Real *a, int n)
{
  int len = n*40;
  char buffer [ len];

  in_file_read_into (f, buffer, len); 
  char *s =buffer;

  char * last_s = buffer;
  while (n--)
    {
      *a++ = strtod (last_s, &s);
      
      assert (s != last_s);
      last_s  = s;
    }
  fscanf (f,"\n");
}



void
write_deformation_state (Deformation_state const*def, char const *fn)
{
  fprintf (stderr, "Writing `%s' ... \n", fn);
  
  FILE * f = fopen (fn, "w");
  if (!f)
    {
      fprintf(stderr, "Can't write file..\n");
      return;
    }
  int n =def->dimension ();

  write_invocation_header (f, "# ");
  fprintf (f, "%d\n",   n);

  write_big_vector (f, def->reference_locations_.access_array (), n);  
  write_big_vector (f, def->displacements_.access_array (), n);
  write_big_vector (f, def->external_force_.access_array (), n);

  fclose (f);

  char s [1024];
  strcpy (s,fn);
  strcat (s,".constraints");
  
  def->constraints_.write_to_file (s);
}




/*
  Try to read deformation state into DEF.

  RETURN

  whether success


  SIDE EFFECT

  def is modified. Message printed if failed. 
 */

bool
read_deformation_state (Deformation_state * def, char const * fn)
{
  FILE * f = fopen (fn, "r");
  if (!f)
    {
      fprintf (stderr, "State file `%s' not found. Not reading deformation state.\n", fn);
      return false;
    }
  else
    {
      fprintf (stderr, "Reading `%s' ...\n", fn);
    }


  char *str = in_file_read (f);
  char * end = str ;
  int n = strtol  (str, &end, 10);
  assert (end != str);

  def->completize_nodal_arrays (n/def->spatial_dimension ());
  
  read_big_vector (f, def->reference_locations_.unsafe_access_array (), n);  
  read_big_vector (f, def->displacements_.unsafe_access_array (),n);
  read_big_vector (f, def->external_force_.unsafe_access_array (),n);

  fclose (f);

  char s [1024];
  strcpy (s,fn);
  strcat (s,".constraints");
  
  def->constraints_.read_from_file (s);
  return true;
}


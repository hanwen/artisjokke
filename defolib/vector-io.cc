#include <ctype.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "vector-io.hh" 


char *
in_file_read (FILE *f)
{
  const int len  = 1024;
  static char s[1024];

  char * p =0; 
  while ((p = fgets (s,len, f)))
    {
      while (*p && isspace (*p))
	p++;

      if (*p && *p != '#')
	break;
    }

  return p;
}

void
in_file_read_into (FILE *f, char *s, int len)
{
  char * p =0; 
  while ((p = fgets (s,len, f)))
    {
      while (*p && isspace (*p))
	p++;

      if (*p && *p != '#')
	break;
    }
}

Vector2
read_vector2 (FILE *f)
{
  char *s1 = in_file_read (f);
  char *s2 = s1;
  Vector2 v;
  
  v(0) = strtod (s1, &s2);
  assert (s1!= s2);
  v(1) = strtod (s2, &s1);
  assert (s1!= s2);

  return v;
}

void
write_vector2 (FILE *f, Vector2 v)
{
  fprintf (f, "%a %a\n", v(0),v(1));
}

void
write_invocation_header (FILE*f, const char * comment )
{
  extern char ** argument_value;
  time_t t (time (0));

  fprintf (f, "%s file created on %s",  comment, ctime(&t));
  fprintf (f, "%s Invoked as: \n", comment); 
  for (char **a =  argument_value;
       *a; a ++)
    {
      fprintf (f, "%s %s \n", comment, *a);
    }
}

FILE *
xfopen (char const *fn, char const *mode)
{
  FILE * f= fopen (fn, mode);

  char what[100];
  strcpy (what , (mode[0] == 'w') ? "writing " : "reading ");

  if (!f)
    {
      fprintf (stderr, "Can't open `%s' for %s.\n", fn, what);
      exit (2);
    }

  what[0] = toupper (what[0]);
  fprintf (stderr, "%s `%s'..\n", what, fn);
  return f;
}

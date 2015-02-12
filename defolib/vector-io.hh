#ifndef VECTOR_IO
#define VECTOR_IO

#include <stdio.h>

#include "vector.hh"
#include "defo-proto.hh" 

Vector2 read_vector2 (FILE *f);
void write_vector2 (FILE *f, Vector2 v);
char * in_file_read (FILE *f);
void in_file_read_into (FILE *f, char *s, int len);
void write_invocation_header (FILE*f, char const*);
FILE * xfopen (char const *, char const *); 
#endif

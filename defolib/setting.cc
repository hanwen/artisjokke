#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "array.hh"
#include "setting.hh"

int
Number_option_list_entry::compare ( Number_option_list_entry const &a,
					       Number_option_list_entry const &b)
{
  return strcmp(a.name ,b.name);
}

void
Number_option_list_entry::print (FILE *out) const
{
  fprintf(out, "  %s -- %s. (default: %lf)\n",
	  name, desc ? desc : "" , value);
}

int 
String_option_list_entry::compare (String_option_list_entry const &a,
				   String_option_list_entry const &b)
{
  return strcmp(a.name ,b.name);
}

void 
String_option_list_entry::print (FILE*out) const
{
  fprintf(out, "  %s -- %s. (default: `%s')\n",
	  name, desc ?desc:"", value?value:"");
}



static Array<Number_option_list_entry> alist;
static Array<String_option_list_entry> salist;

static int
number_setting_index (char const * ch)
{
  for (int i= alist.size ();i--;)
    if (!strcmp (alist[i].name, ch))
      return i;
  return -1;
}

static int
string_setting_index (char const * ch)
{
  for (int i= salist.size ();i--;)
    if (!strcmp (salist[i].name, ch))
      return i;
  return -1;
}

Real
get_number_setting (char const * ch)
{
  int i = number_setting_index (ch);

  if (i <0)
    {
      fprintf (stderr, "No  such setting: %s", ch);
      return 0.0;
    }
  else
    return alist[i].value;
}

const char*
get_string_setting (char const * ch)
{
  int i = string_setting_index (ch);

  if (i <0)
    return "";
  else
    return salist[i].value;
}

bool
get_bool_setting (const char * s)
{
  return !strcmp (get_string_setting (s), "yes");
}


void
set_bool_setting (const char * s, bool b)
{
  set_string_setting (s, b ? "yes" : "no");
}


void
set_number_setting (char const * k, Real r)
{
  int i = number_setting_index (k);
  if (i < 0)
    {
      printf( "Warning: no description found for %s. Perhaps you made a typo?\n", k);
      exit (2);
    }
  else
    alist[i].value = r;
}


void
set_string_setting (char const * k, char const *val)
{
  int i = string_setting_index (k);
  if (i < 0)
    {
      printf( "Warning: no description found for %s. Perhaps you made a typo?\n", k);
      exit (2);

    }
  else
    salist[i].value = strdup (val);
}





void
init_settings (Number_option_list_entry * init_alist,
	       String_option_list_entry * init_str_alist)
{
  for (Number_option_list_entry * p = init_alist; p->name; p++)
    alist.push  (*p);
  for (String_option_list_entry *p = init_str_alist; p->name ; p++)
    salist.push (*p);

  salist.sort (String_option_list_entry::compare);
  alist.sort (Number_option_list_entry::compare);
}

void
print_settings ()
{
  for (int i = 0; i < salist.size (); i++)
    salist[i].print (stdout);

  printf("\n");
  for (int i = 0; i < alist.size (); i++)
    alist[i].print (stdout);
}





void
parse_setting (char  * k)
{
  char * p = strchr (k, '=');
  if(p)
    {
      char key[1024];
      strcpy (key, k);
      key[p -k] =0; 

      p++;
      char *end ; 

      Real r = strtod (p, &end);
      
      //      int succ = sscanf (p, "%lf", &r);
      bool succ = (p != end);
      
      if (succ == 1)
	{
	  set_number_setting (key, r);
	  return ;
	}
      else
	{
	  set_string_setting (key, p);
	  return; 
	}
      
    }
  if (!strcmp (k, "help"))
    {
      print_settings ();
      exit (0);
    }
    
  fprintf (stderr,"Couldn't parse `%s'. Use -ohelp for all options.\n ",  k);
  exit (2);
}


static Real update_factor;
Real
get_update_factor ()
{
  if (!update_factor)
    {
      update_factor = get_number_setting ("update-factor");
      if (!update_factor)
	update_factor = 1.0;
    }

  return update_factor;
}


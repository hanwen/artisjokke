
#ifndef SETTING_HH
#define SETTING_HH

#include <stdio.h>

#include "defo-proto.hh"

extern bool global_test_flag;

struct Number_option_list_entry
{
  char * name;
  char * desc;
  Real value;

  static int compare ( Number_option_list_entry const &a,
		       Number_option_list_entry const &b);
  void print (FILE *out) const;
};

struct String_option_list_entry
{
  char * name;
  char * desc;
  char * value;

  static int compare (String_option_list_entry const &a,
		      String_option_list_entry const &b);
  void print (FILE*out) const;
};



Real get_number_setting (char const * ch);
char const * get_string_setting (char const *ch);
bool get_bool_setting (char const * ch);
void set_bool_setting (char const * ch, bool);
void set_number_setting (char const * k, Real r);
void set_string_setting (char const * k, const char* );
void parse_setting (char  * k);
Real get_update_factor ();
void init_settings (Number_option_list_entry * init_alist,
	       String_option_list_entry * init_str_alist );

void print_settings ();

#endif

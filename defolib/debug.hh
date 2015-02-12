#ifndef DEBUG_HH
#define DEBUG_HH

extern bool debug_out;

#define log_message(x, ...) \
do { if (debug_out) \
  fprintf(stderr, x, ## __VA_ARGS__); \
  fflush (stderr); \
} while (0)

#endif


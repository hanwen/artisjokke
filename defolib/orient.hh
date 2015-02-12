
#ifndef PREDICATE_H
#define PREDICATE_H

/* On some machines, the exact arithmetic routines might be defeated by the  */
/*   use of internal extended precision floating-point registers.  Sometimes */
/*   this problem can be fixed by defining certain values to be volatile,    */
/*   thus forcing them to be stored to memory and rounded off.  This isn't   */
/*   a great solution, though, as it slows the arithmetic down.              */
/*                                                                           */
/* To try this out, write "#define INEXACT volatile" below.  Normally,       */
/*   however, INEXACT should be defined to be nothing.  ("#define INEXACT".) */

#define INEXACT                          /* Nothing */
/*#define INEXACT volatile */

#include <fpu_control.h>


inline void
set_extended_precision (int extended_on)
{
  fpu_control_t curr;
  _FPU_GETCW(curr);

  curr = curr & ~_FPU_EXTENDED;
  if (extended_on)
    {
      curr = curr | _FPU_EXTENDED;
    }
  else
    {
#ifdef PREDICATES_USE_FLOAT
      curr = curr | _FPU_SINGLE;
#else
      curr = curr | _FPU_DOUBLE;
#endif
    }

  _FPU_SETCW(curr);
}


#ifdef PREDICATES_USE_FLOAT

#ifndef REAL_DEFINED
typedef float Real;                      /* float or float */
#endif

#define REALPRINT floatprint
#define REALRAND floatrand
#define NARROWRAND narrowfloatrand
#define UNIFORMRAND uniformfloatrand

#else

#ifndef REAL_DEFINED
typedef double Real;
#endif

#define REALPRINT doubleprint
#define REALRAND doublerand
#define NARROWRAND narrowdoublerand
#define UNIFORMRAND uniformdoublerand

#endif


/* Which of the following two methods of finding the absolute values is      */
/*   fastest is compiler-dependent.  A few compilers can inline and optimize */
/*   the fabs() call; but most will incur the overhead of a function call,   */
/*   which is disastrously slow.  A faster way on IEEE machines might be to  */
/*   mask the appropriate bit, but that's difficult to do in C.              */

#define Absolute(a)  ((a) >= 0.0 ? (a) : -(a))

/* #define Absolute(a)  fabs(a) */



Real robust_insphere(Real const *pa, Real const *pb, Real const *pc, Real const *pd, Real const *pe);
Real robust_orient3d(Real const *pa, Real const *pb, Real const *pc, Real const *pd);
     
#endif
     


#include <fpu_control.h>

#include "defo-proto.hh"

#define INEXACT


/* Which of the following two methods of finding the absolute values is      */
/*   fastest is compiler-dependent.  A few compilers can inline and optimize */
/*   the fabs() call; but most will incur the overhead of a function call,   */
/*   which is disastrously slow.  A faster way on IEEE machines might be to  */
/*   mask the appropriate bit, but that's difficult to do in C.              */

#define Absolute(a)  ((a) >= 0.0 ? (a) : -(a))

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


/* Many of the operations are broken up into two pieces, a main part that    */
/*   performs an approximate operation, and a "tail" that computes the       */
/*   roundoff error of that operation.                                       */
/*                                                                           */
/* The operations Fast_Two_Sum(), Fast_Two_Diff(), Two_Sum(), Two_Diff(),    */
/*   Split(), and Two_Product() are all implemented as described in the      */
/*   reference.  Each of these macros requires certain variables to be       */
/*   defined in the calling routine.  The variables `bvirt', `c', `abig',    */
/*   `_i', `_j', `_k', `_l', `_m', and `_n' are declared `INEXACT' because   */
/*   they store the result of an operation that may incur roundoff error.    */
/*   The input parameter `x' (or the highest numbered `x_' parameter) must   */
/*   also be declared `INEXACT'.                                             */

#define Fast_Two_Sum_Tail(a, b, x, y) \
  bvirt = x - a; \
  y = b - bvirt

#define Fast_Two_Sum(a, b, x, y) \
  x = (Real) (a + b); \
  Fast_Two_Sum_Tail(a, b, x, y)

#define Fast_Two_Diff_Tail(a, b, x, y) \
  bvirt = a - x; \
  y = bvirt - b

#define Fast_Two_Diff(a, b, x, y) \
  x = (Real) (a - b); \
  Fast_Two_Diff_Tail(a, b, x, y)

#define Two_Sum_Tail(a, b, x, y) \
  bvirt = (Real) (x - a); \
  avirt = x - bvirt; \
  bround = b - bvirt; \
  around = a - avirt; \
  y = around + bround

#define Two_Sum(a, b, x, y) \
  x = (Real) (a + b); \
  Two_Sum_Tail(a, b, x, y)

#define Two_Diff_Tail(a, b, x, y) \
  bvirt = (Real) (a - x); \
  avirt = x + bvirt; \
  bround = bvirt - b; \
  around = a - avirt; \
  y = around + bround

#define Two_Diff(a, b, x, y) \
  x = (Real) (a - b); \
  Two_Diff_Tail(a, b, x, y)

#define Split(a, ahi, alo) \
  c = (Real) (splitter * a); \
  abig = (Real) (c - a); \
  ahi = c - abig; \
  alo = a - ahi

#define Two_Product_Tail(a, b, x, y) \
  Split(a, ahi, alo); \
  Split(b, bhi, blo); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

#define Two_Product(a, b, x, y) \
  x = (Real) (a * b); \
  Two_Product_Tail(a, b, x, y)

/* Two_Product_Presplit() is Two_Product() where one of the inputs has       */
/*   already been split.  Avoids redundant splitting.                        */

#define Two_Product_Presplit(a, b, bhi, blo, x, y) \
  x = (Real) (a * b); \
  Split(a, ahi, alo); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

/* Two_Product_2Presplit() is Two_Product() where both of the inputs have    */
/*   already been split.  Avoids redundant splitting.                        */

#define Two_Product_2Presplit(a, ahi, alo, b, bhi, blo, x, y) \
  x = (Real) (a * b); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

/* Square() can be done more quickly than Two_Product().                     */

#define Square_Tail(a, x, y) \
  Split(a, ahi, alo); \
  err1 = x - (ahi * ahi); \
  err3 = err1 - ((ahi + ahi) * alo); \
  y = (alo * alo) - err3

#define Square(a, x, y) \
  x = (Real) (a * a); \
  Square_Tail(a, x, y)

/* Macros for summing expansions of various fixed lengths.  These are all    */
/*   unrolled versions of Expansion_Sum().                                   */

#define Two_One_Sum(a1, a0, b, x2, x1, x0) \
  Two_Sum(a0, b , _i, x0); \
  Two_Sum(a1, _i, x2, x1)

#define Two_One_Diff(a1, a0, b, x2, x1, x0) \
  Two_Diff(a0, b , _i, x0); \
  Two_Sum( a1, _i, x2, x1)

#define Two_Two_Sum(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Sum(a1, a0, b0, _j, _0, x0); \
  Two_One_Sum(_j, _0, b1, x3, x2, x1)

#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Diff(a1, a0, b0, _j, _0, x0); \
  Two_One_Diff(_j, _0, b1, x3, x2, x1)

#define Four_One_Sum(a3, a2, a1, a0, b, x4, x3, x2, x1, x0) \
  Two_One_Sum(a1, a0, b , _j, x1, x0); \
  Two_One_Sum(a3, a2, _j, x4, x3, x2)

#define Four_Two_Sum(a3, a2, a1, a0, b1, b0, x5, x4, x3, x2, x1, x0) \
  Four_One_Sum(a3, a2, a1, a0, b0, _k, _2, _1, _0, x0); \
  Four_One_Sum(_k, _2, _1, _0, b1, x5, x4, x3, x2, x1)

#define Four_Four_Sum(a3, a2, a1, a0, b4, b3, b1, b0, x7, x6, x5, x4, x3, x2, \
                      x1, x0) \
  Four_Two_Sum(a3, a2, a1, a0, b1, b0, _l, _2, _1, _0, x1, x0); \
  Four_Two_Sum(_l, _2, _1, _0, b4, b3, x7, x6, x5, x4, x3, x2)

#define Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b, x8, x7, x6, x5, x4, \
                      x3, x2, x1, x0) \
  Four_One_Sum(a3, a2, a1, a0, b , _j, x3, x2, x1, x0); \
  Four_One_Sum(a7, a6, a5, a4, _j, x8, x7, x6, x5, x4)

#define Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, x9, x8, x7, \
                      x6, x5, x4, x3, x2, x1, x0) \
  Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b0, _k, _6, _5, _4, _3, _2, \
                _1, _0, x0); \
  Eight_One_Sum(_k, _6, _5, _4, _3, _2, _1, _0, b1, x9, x8, x7, x6, x5, x4, \
                x3, x2, x1)

#define Eight_Four_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b4, b3, b1, b0, x11, \
                       x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0) \
  Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, _l, _6, _5, _4, _3, \
                _2, _1, _0, x1, x0); \
  Eight_Two_Sum(_l, _6, _5, _4, _3, _2, _1, _0, b4, b3, x11, x10, x9, x8, \
                x7, x6, x5, x4, x3, x2)

/* Macros for multiplying expansions of various fixed lengths.               */

#define Two_One_Product(a1, a0, b, x3, x2, x1, x0) \
  Split(b, bhi, blo); \
  Two_Product_Presplit(a0, b, bhi, blo, _i, x0); \
  Two_Product_Presplit(a1, b, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, x1); \
  Fast_Two_Sum(_j, _k, x3, x2)

#define Four_One_Product(a3, a2, a1, a0, b, x7, x6, x5, x4, x3, x2, x1, x0) \
  Split(b, bhi, blo); \
  Two_Product_Presplit(a0, b, bhi, blo, _i, x0); \
  Two_Product_Presplit(a1, b, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, x1); \
  Fast_Two_Sum(_j, _k, _i, x2); \
  Two_Product_Presplit(a2, b, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, x3); \
  Fast_Two_Sum(_j, _k, _i, x4); \
  Two_Product_Presplit(a3, b, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, x5); \
  Fast_Two_Sum(_j, _k, x7, x6)

#define Two_Two_Product(a1, a0, b1, b0, x7, x6, x5, x4, x3, x2, x1, x0) \
  Split(a0, a0hi, a0lo); \
  Split(b0, bhi, blo); \
  Two_Product_2Presplit(a0, a0hi, a0lo, b0, bhi, blo, _i, x0); \
  Split(a1, a1hi, a1lo); \
  Two_Product_2Presplit(a1, a1hi, a1lo, b0, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, _1); \
  Fast_Two_Sum(_j, _k, _l, _2); \
  Split(b1, bhi, blo); \
  Two_Product_2Presplit(a0, a0hi, a0lo, b1, bhi, blo, _i, _0); \
  Two_Sum(_1, _0, _k, x1); \
  Two_Sum(_2, _k, _j, _1); \
  Two_Sum(_l, _j, _m, _2); \
  Two_Product_2Presplit(a1, a1hi, a1lo, b1, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _n, _0); \
  Two_Sum(_1, _0, _i, x2); \
  Two_Sum(_2, _i, _k, _1); \
  Two_Sum(_m, _k, _l, _2); \
  Two_Sum(_j, _n, _k, _0); \
  Two_Sum(_1, _0, _j, x3); \
  Two_Sum(_2, _j, _i, _1); \
  Two_Sum(_l, _i, _m, _2); \
  Two_Sum(_1, _k, _i, x4); \
  Two_Sum(_2, _i, _k, x5); \
  Two_Sum(_m, _k, x7, x6)

/* An expansion of length two can be squared more quickly than finding the   */
/*   product of two different expansions of length two, and the result is    */
/*   guaranteed to have no more than six (rather than eight) components.     */

#define Two_Square(a1, a0, x5, x4, x3, x2, x1, x0) \
  Square(a0, _j, x0); \
  _0 = a0 + a0; \
  Two_Product(a1, _0, _k, _1); \
  Two_One_Sum(_k, _1, _j, _l, _2, x1); \
  Square(a1, _j, _1); \
  Two_Two_Sum(_j, _1, _l, _2, x5, x4, x3, x2)

Real splitter;     /* = 2^ceiling(p / 2) + 1.  Used to split floats in half. */
Real epsilon;                /* = 2^(-p).  Used to estimate roundoff errors. */
/* A set of coefficients used to calculate maximum roundoff errors.          */
Real resulterrbound;
Real ccwerrboundA, ccwerrboundB, ccwerrboundC;
Real o3derrboundA, o3derrboundB, o3derrboundC;
Real iccerrboundA, iccerrboundB, iccerrboundC;
Real isperrboundA, isperrboundB, isperrboundC;


/*****************************************************************************/
/*                                                                           */
/*  fast_expansion_sum_zeroelim()   Sum two expansions, eliminating zero     */
/*                                  components from the output expansion.    */
/*                                                                           */
/*  Sets h = e + f.  See the long version of my paper for details.           */
/*                                                                           */
/*  If round-to-even is used (as with IEEE 754), maintains the strongly      */
/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   */
/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      */
/*  properties.                                                              */
/*                                                                           */
/*****************************************************************************/

static int
fast_expansion_sum_zeroelim(int elen, Real *e, int flen, Real *f, Real *h)  /* h cannot be e or f. */
         
        
         
        
        
{
  Real Q;
  INEXACT Real Qnew;
  INEXACT Real hh;
  INEXACT Real bvirt;
  Real avirt, bround, around;
  int eindex, findex, hindex;
  Real enow, fnow;

  enow = e[0];
  fnow = f[0];
  eindex = findex = 0;
  if ((fnow > enow) == (fnow > -enow)) {
    Q = enow;
    enow = e[++eindex];
  } else {
    Q = fnow;
    fnow = f[++findex];
  }
  hindex = 0;
  if ((eindex < elen) && (findex < flen)) {
    if ((fnow > enow) == (fnow > -enow)) {
      Fast_Two_Sum(enow, Q, Qnew, hh);
      enow = e[++eindex];
    } else {
      Fast_Two_Sum(fnow, Q, Qnew, hh);
      fnow = f[++findex];
    }
    Q = Qnew;
    if (hh != 0.0) {
      h[hindex++] = hh;
    }
    while ((eindex < elen) && (findex < flen)) {
      if ((fnow > enow) == (fnow > -enow)) {
        Two_Sum(Q, enow, Qnew, hh);
        enow = e[++eindex];
      } else {
        Two_Sum(Q, fnow, Qnew, hh);
        fnow = f[++findex];
      }
      Q = Qnew;
      if (hh != 0.0) {
        h[hindex++] = hh;
      }
    }
  }
  while (eindex < elen) {
    Two_Sum(Q, enow, Qnew, hh);
    enow = e[++eindex];
    Q = Qnew;
    if (hh != 0.0) {
      h[hindex++] = hh;
    }
  }
  while (findex < flen) {
    Two_Sum(Q, fnow, Qnew, hh);
    fnow = f[++findex];
    Q = Qnew;
    if (hh != 0.0) {
      h[hindex++] = hh;
    }
  }
  if ((Q != 0.0) || (hindex == 0)) {
    h[hindex++] = Q;
  }
  return hindex;
}


/*****************************************************************************/
/*                                                                           */
/*  scale_expansion()   Multiply an expansion by a scalar.                   */
/*                                                                           */
/*  Sets h = be.  See either version of my paper for details.                */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */
/*  properties as well.  (That is, if e has one of these properties, so      */
/*  will h.)                                                                 */
/*                                                                           */
/*****************************************************************************/

static int
scale_expansion(int elen, Real *e, Real b, Real *h)            /* e and h cannot be the same. */
{
  INEXACT Real Q;
  INEXACT Real sum;
  INEXACT Real product1;
  Real product0;
  int eindex, hindex;
  Real enow;
  INEXACT Real bvirt;
  Real avirt, bround, around;
  INEXACT Real c;
  INEXACT Real abig;
  Real ahi, alo, bhi, blo;
  Real err1, err2, err3;

  Split(b, bhi, blo);
  Two_Product_Presplit(e[0], b, bhi, blo, Q, h[0]);
  hindex = 1;
  for (eindex = 1; eindex < elen; eindex++) {
    enow = e[eindex];
    Two_Product_Presplit(enow, b, bhi, blo, product1, product0);
    Two_Sum(Q, product0, sum, h[hindex]);
    hindex++;
    Two_Sum(product1, sum, Q, h[hindex]);
    hindex++;
  }
  h[hindex] = Q;
  return elen + elen;
}

/*****************************************************************************/
/*                                                                           */
/*  scale_expansion_zeroelim()   Multiply an expansion by a scalar,          */
/*                               eliminating zero components from the        */
/*                               output expansion.                           */
/*                                                                           */
/*  Sets h = be.  See either version of my paper for details.                */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */
/*  properties as well.  (That is, if e has one of these properties, so      */
/*  will h.)                                                                 */
/*                                                                           */
/*****************************************************************************/

static int
scale_expansion_zeroelim(int elen, Real *e, Real b, Real *h)   /* e and h cannot be the same. */
{
  INEXACT Real Q, sum;
  Real hh;
  INEXACT Real product1;
  Real product0;
  int eindex, hindex;
  Real enow;
  INEXACT Real bvirt;
  Real avirt, bround, around;
  INEXACT Real c;
  INEXACT Real abig;
  Real ahi, alo, bhi, blo;
  Real err1, err2, err3;

  Split(b, bhi, blo);
  Two_Product_Presplit(e[0], b, bhi, blo, Q, hh);
  hindex = 0;
  if (hh != 0) {
    h[hindex++] = hh;
  }
  for (eindex = 1; eindex < elen; eindex++) {
    enow = e[eindex];
    Two_Product_Presplit(enow, b, bhi, blo, product1, product0);
    Two_Sum(Q, product0, sum, hh);
    if (hh != 0) {
      h[hindex++] = hh;
    }
    Fast_Two_Sum(product1, sum, Q, hh);
    if (hh != 0) {
      h[hindex++] = hh;
    }
  }
  if ((Q != 0.0) || (hindex == 0)) {
    h[hindex++] = Q;
  }
  return hindex;
}

/*****************************************************************************/
/*                                                                           */
/*  compress()   Compress an expansion.                                      */
/*                                                                           */
/*  See the long version of my paper for details.                            */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), then any nonoverlapping expansion is converted to a      */
/*  nonadjacent expansion.                                                   */
/*                                                                           */
/*****************************************************************************/

static int
compress(int elen, Real *e, Real *h)                         /* e and h may be the same. */
{
  Real Q, q;
  INEXACT Real Qnew;
  int eindex, hindex;
  INEXACT Real bvirt;
  Real enow, hnow;
  int top, bottom;

  bottom = elen - 1;
  Q = e[bottom];
  for (eindex = elen - 2; eindex >= 0; eindex--) {
    enow = e[eindex];
    Fast_Two_Sum(Q, enow, Qnew, q);
    if (q != 0) {
      h[bottom--] = Qnew;
      Q = q;
    } else {
      Q = Qnew;
    }
  }
  top = 0;
  for (hindex = bottom + 1; hindex < elen; hindex++) {
    hnow = h[hindex];
    Fast_Two_Sum(hnow, Q, Qnew, q);
    if (q != 0) {
      h[top++] = q;
    }
    Q = Qnew;
  }
  h[top] = Q;
  return top + 1;
}

/*****************************************************************************/
/*                                                                           */
/*  estimate()   Produce a one-word estimate of an expansion's value.        */
/*                                                                           */
/*  See either version of my paper for details.                              */
/*                                                                           */
/*****************************************************************************/

Real
estimate(int elen, Real *e)
{
  Real Q;
  int eindex;

  Q = e[0];
  for (eindex = 1; eindex < elen; eindex++) {
    Q += e[eindex];
  }
  return Q;
}


/*****************************************************************************/
/*                                                                           */
/*  orient2dfast()   Approximate 2D orientation test.  Nonrobust.            */
/*  orient2dexact()   Exact 2D orientation test.  Robust.                    */
/*  orient2dslow()   Another exact 2D orientation test.  Robust.             */
/*  orient2d()   Adaptive exact 2D orientation test.  Robust.                */
/*                                                                           */
/*               Return a positive value if the points pa, pb, and pc occur  */
/*               in counterclockwise order; a negative value if they occur   */
/*               in clockwise order; and zero if they are collinear.  The    */
/*               result is also a rough approximation of twice the signed    */
/*               area of the triangle defined by the three points.           */
/*                                                                           */
/*  Only the first and last routine should be used; the middle two are for   */
/*  timings.                                                                 */
/*                                                                           */
/*  The last three use exact arithmetic to ensure a correct answer.  The     */
/*  result returned is the determinant of a matrix.  In orient2d() only,     */
/*  this determinant is computed adaptively, in the sense that exact         */
/*  arithmetic is used only to the degree it is needed to ensure that the    */
/*  returned value has the correct sign.  Hence, orient2d() is usually quite */
/*  fast, but will run more slowly when the input points are collinear or    */
/*  nearly so.                                                               */
/*                                                                           */
/*****************************************************************************/

Real
orient2dfast(Real *pa, Real *pb, Real *pc)
{
  Real acx, bcx, acy, bcy;

  acx = pa[0] - pc[0];
  bcx = pb[0] - pc[0];
  acy = pa[1] - pc[1];
  bcy = pb[1] - pc[1];
  return acx * bcy - acy * bcx;
}

static Real
orient2dexact(Real *pa, Real *pb, Real *pc)
{
  INEXACT Real axby1, axcy1, bxcy1, bxay1, cxay1, cxby1;
  Real axby0, axcy0, bxcy0, bxay0, cxay0, cxby0;
  Real aterms[4], bterms[4], cterms[4];
  INEXACT Real aterms3, bterms3, cterms3;
  Real v[8], w[12];
  int vlength, wlength;

  INEXACT Real bvirt;
  Real avirt, bround, around;
  INEXACT Real c;
  INEXACT Real abig;
  Real ahi, alo, bhi, blo;
  Real err1, err2, err3;
  INEXACT Real _i, _j;
  Real _0;

  Two_Product(pa[0], pb[1], axby1, axby0);
  Two_Product(pa[0], pc[1], axcy1, axcy0);
  Two_Two_Diff(axby1, axby0, axcy1, axcy0,
               aterms3, aterms[2], aterms[1], aterms[0]);
  aterms[3] = aterms3;

  Two_Product(pb[0], pc[1], bxcy1, bxcy0);
  Two_Product(pb[0], pa[1], bxay1, bxay0);
  Two_Two_Diff(bxcy1, bxcy0, bxay1, bxay0,
               bterms3, bterms[2], bterms[1], bterms[0]);
  bterms[3] = bterms3;

  Two_Product(pc[0], pa[1], cxay1, cxay0);
  Two_Product(pc[0], pb[1], cxby1, cxby0);
  Two_Two_Diff(cxay1, cxay0, cxby1, cxby0,
               cterms3, cterms[2], cterms[1], cterms[0]);
  cterms[3] = cterms3;

  vlength = fast_expansion_sum_zeroelim(4, aterms, 4, bterms, v);
  wlength = fast_expansion_sum_zeroelim(vlength, v, 4, cterms, w);

  return w[wlength - 1];
}

static Real
orient2dslow(Real *pa, Real *pb, Real *pc)
{
  INEXACT Real acx, acy, bcx, bcy;
  Real acxtail, acytail;
  Real bcxtail, bcytail;
  Real negate, negatetail;
  Real axby[8], bxay[8];
  INEXACT Real axby7, bxay7;
  Real deter[16];
  int deterlen;

  INEXACT Real bvirt;
  Real avirt, bround, around;
  INEXACT Real c;
  INEXACT Real abig;
  Real a0hi, a0lo, a1hi, a1lo, bhi, blo;
  Real err1, err2, err3;
  INEXACT Real _i, _j, _k, _l, _m, _n;
  Real _0, _1, _2;

  Two_Diff(pa[0], pc[0], acx, acxtail);
  Two_Diff(pa[1], pc[1], acy, acytail);
  Two_Diff(pb[0], pc[0], bcx, bcxtail);
  Two_Diff(pb[1], pc[1], bcy, bcytail);

  Two_Two_Product(acx, acxtail, bcy, bcytail,
                  axby7, axby[6], axby[5], axby[4],
                  axby[3], axby[2], axby[1], axby[0]);
  axby[7] = axby7;
  negate = -acy;
  negatetail = -acytail;
  Two_Two_Product(bcx, bcxtail, negate, negatetail,
                  bxay7, bxay[6], bxay[5], bxay[4],
                  bxay[3], bxay[2], bxay[1], bxay[0]);
  bxay[7] = bxay7;

  deterlen = fast_expansion_sum_zeroelim(8, axby, 8, bxay, deter);

  return deter[deterlen - 1];
}

static Real
orient2dadapt(Real const *pa, Real const *pb, Real const *pc, Real detsum)
{
  INEXACT Real acx, acy, bcx, bcy;
  Real acxtail, acytail, bcxtail, bcytail;
  INEXACT Real detleft, detright;
  Real detlefttail, detrighttail;
  Real det, errbound;
  Real B[4], C1[8], C2[12], D[16];
  INEXACT Real B3;
  int C1length, C2length, Dlength;
  Real u[4];
  INEXACT Real u3;
  INEXACT Real s1, t1;
  Real s0, t0;

  INEXACT Real bvirt;
  Real avirt, bround, around;
  INEXACT Real c;
  INEXACT Real abig;
  Real ahi, alo, bhi, blo;
  Real err1, err2, err3;
  INEXACT Real _i, _j;
  Real _0;

  acx = (Real) (pa[0] - pc[0]);
  bcx = (Real) (pb[0] - pc[0]);
  acy = (Real) (pa[1] - pc[1]);
  bcy = (Real) (pb[1] - pc[1]);

  Two_Product(acx, bcy, detleft, detlefttail);
  Two_Product(acy, bcx, detright, detrighttail);

  Two_Two_Diff(detleft, detlefttail, detright, detrighttail,
               B3, B[2], B[1], B[0]);
  B[3] = B3;

  det = estimate(4, B);
  errbound = ccwerrboundB * detsum;
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  Two_Diff_Tail(pa[0], pc[0], acx, acxtail);
  Two_Diff_Tail(pb[0], pc[0], bcx, bcxtail);
  Two_Diff_Tail(pa[1], pc[1], acy, acytail);
  Two_Diff_Tail(pb[1], pc[1], bcy, bcytail);

  if ((acxtail == 0.0) && (acytail == 0.0)
      && (bcxtail == 0.0) && (bcytail == 0.0)) {
    return det;
  }

  errbound = ccwerrboundC * detsum + resulterrbound * Absolute(det);
  det += (acx * bcytail + bcy * acxtail)
       - (acy * bcxtail + bcx * acytail);
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  Two_Product(acxtail, bcy, s1, s0);
  Two_Product(acytail, bcx, t1, t0);
  Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
  u[3] = u3;
  C1length = fast_expansion_sum_zeroelim(4, B, 4, u, C1);

  Two_Product(acx, bcytail, s1, s0);
  Two_Product(acy, bcxtail, t1, t0);
  Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
  u[3] = u3;
  C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, C2);

  Two_Product(acxtail, bcytail, s1, s0);
  Two_Product(acytail, bcxtail, t1, t0);
  Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
  u[3] = u3;
  Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, D);

  return(D[Dlength - 1]);
}

Real
robust_orient2d(Real const *pa, Real const *pb, Real const *pc)
{
  Real detleft, detright, det;
  Real detsum, errbound;

  set_extended_precision (0);

  
  detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
  detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
  det = detleft - detright;

  if (detleft > 0.0) {
    if (detright <= 0.0) {
      return det;
    } else {
      detsum = detleft + detright;
    }
  } else if (detleft < 0.0) {
    if (detright >= 0.0) {
      return det;
    } else {
      detsum = -detleft - detright;
    }
  } else {
    return det;
  }

  errbound = ccwerrboundA * detsum;
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  {

    Real adapt = orient2dadapt(pa, pb, pc, detsum);
    set_extended_precision (1);
    

    return adapt;
  }
    
}


#include "defo-proto.hh"
/*****************************************************************************/
/*                                                                           */
/*  orient3dfast()   Approximate 3D orientation test.  Nonrobust.            */
/*  orient3dexact()   Exact 3D orientation test.  Robust.                    */
/*  orient3dslow()   Another exact 3D orientation test.  Robust.             */
/*  orient3d()   Adaptive exact 3D orientation test.  Robust.                */
/*                                                                           */
/*               Return a positive value if the point pd lies below the      */
/*               plane passing through pa, pb, and pc; "below" is defined so */
/*               that pa, pb, and pc appear in counterclockwise order when   */
/*               viewed from above the plane.  Returns a negative value if   */
/*               pd lies above the plane.  Returns zero if the points are    */
/*               coplanar.  The result is also a rough approximation of six  */
/*               times the signed volume of the tetrahedron defined by the   */
/*               four points.                                                */
/*                                                                           */
/*  Only the first and last routine should be used; the middle two are for   */
/*  timings.                                                                 */
/*                                                                           */
/*  The last three use exact arithmetic to ensure a correct answer.  The     */
/*  result returned is the determinant of a matrix.  In orient3d() only,     */
/*  this determinant is computed adaptively, in the sense that exact         */
/*  arithmetic is used only to the degree it is needed to ensure that the    */
/*  returned value has the correct sign.  Hence, orient3d() is usually quite */
/*  fast, but will run more slowly when the input points are coplanar or     */
/*  nearly so.                                                               */
/*                                                                           */
/*****************************************************************************/

Real
orient3dfast(Real const *pa, Real const *pb, Real const *pc, Real const *pd)
{
  Real adx, bdx, cdx;
  Real ady, bdy, cdy;
  Real adz, bdz, cdz;

  adx = pa[0] - pd[0];
  bdx = pb[0] - pd[0];
  cdx = pc[0] - pd[0];
  ady = pa[1] - pd[1];
  bdy = pb[1] - pd[1];
  cdy = pc[1] - pd[1];
  adz = pa[2] - pd[2];
  bdz = pb[2] - pd[2];
  cdz = pc[2] - pd[2];

  return adx * (bdy * cdz - bdz * cdy)
       + bdx * (cdy * adz - cdz * ady)
       + cdx * (ady * bdz - adz * bdy);
}

#if 0
static Real
orient3dexact(Real const *pa, Real const *pb, Real const *pc, Real const *pd)
{
  INEXACT Real axby1, bxcy1, cxdy1, dxay1, axcy1, bxdy1;
  INEXACT Real bxay1, cxby1, dxcy1, axdy1, cxay1, dxby1;
  Real axby0, bxcy0, cxdy0, dxay0, axcy0, bxdy0;
  Real bxay0, cxby0, dxcy0, axdy0, cxay0, dxby0;
  Real ab[4], bc[4], cd[4], da[4], ac[4], bd[4];
  Real temp8[8];
  int templen;
  Real abc[12], bcd[12], cda[12], dab[12];
  int abclen, bcdlen, cdalen, dablen;
  Real adet[24], bdet[24], cdet[24], ddet[24];
  int alen, blen, clen, dlen;
  Real abdet[48], cddet[48];
  int ablen, cdlen;
  Real deter[96];
  int deterlen;
  int i;

  INEXACT Real bvirt;
  Real avirt, bround, around;
  INEXACT Real c;
  INEXACT Real abig;
  Real ahi, alo, bhi, blo;
  Real err1, err2, err3;
  INEXACT Real _i, _j;
  Real _0;

  Two_Product(pa[0], pb[1], axby1, axby0);
  Two_Product(pb[0], pa[1], bxay1, bxay0);
  Two_Two_Diff(axby1, axby0, bxay1, bxay0, ab[3], ab[2], ab[1], ab[0]);

  Two_Product(pb[0], pc[1], bxcy1, bxcy0);
  Two_Product(pc[0], pb[1], cxby1, cxby0);
  Two_Two_Diff(bxcy1, bxcy0, cxby1, cxby0, bc[3], bc[2], bc[1], bc[0]);

  Two_Product(pc[0], pd[1], cxdy1, cxdy0);
  Two_Product(pd[0], pc[1], dxcy1, dxcy0);
  Two_Two_Diff(cxdy1, cxdy0, dxcy1, dxcy0, cd[3], cd[2], cd[1], cd[0]);

  Two_Product(pd[0], pa[1], dxay1, dxay0);
  Two_Product(pa[0], pd[1], axdy1, axdy0);
  Two_Two_Diff(dxay1, dxay0, axdy1, axdy0, da[3], da[2], da[1], da[0]);

  Two_Product(pa[0], pc[1], axcy1, axcy0);
  Two_Product(pc[0], pa[1], cxay1, cxay0);
  Two_Two_Diff(axcy1, axcy0, cxay1, cxay0, ac[3], ac[2], ac[1], ac[0]);

  Two_Product(pb[0], pd[1], bxdy1, bxdy0);
  Two_Product(pd[0], pb[1], dxby1, dxby0);
  Two_Two_Diff(bxdy1, bxdy0, dxby1, dxby0, bd[3], bd[2], bd[1], bd[0]);

  templen = fast_expansion_sum_zeroelim(4, cd, 4, da, temp8);
  cdalen = fast_expansion_sum_zeroelim(templen, temp8, 4, ac, cda);
  templen = fast_expansion_sum_zeroelim(4, da, 4, ab, temp8);
  dablen = fast_expansion_sum_zeroelim(templen, temp8, 4, bd, dab);
  for (i = 0; i < 4; i++) {
    bd[i] = -bd[i];
    ac[i] = -ac[i];
  }
  templen = fast_expansion_sum_zeroelim(4, ab, 4, bc, temp8);
  abclen = fast_expansion_sum_zeroelim(templen, temp8, 4, ac, abc);
  templen = fast_expansion_sum_zeroelim(4, bc, 4, cd, temp8);
  bcdlen = fast_expansion_sum_zeroelim(templen, temp8, 4, bd, bcd);

  alen = scale_expansion_zeroelim(bcdlen, bcd, pa[2], adet);
  blen = scale_expansion_zeroelim(cdalen, cda, -pb[2], bdet);
  clen = scale_expansion_zeroelim(dablen, dab, pc[2], cdet);
  dlen = scale_expansion_zeroelim(abclen, abc, -pd[2], ddet);

  ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
  cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
  deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, deter);

  return deter[deterlen - 1];
}

static Real
orient3dslow(Real const *pa, Real const *pb, Real const *pc, Real const *pd)
{
  INEXACT Real adx, ady, adz, bdx, bdy, bdz, cdx, cdy, cdz;
  Real adxtail, adytail, adztail;
  Real bdxtail, bdytail, bdztail;
  Real cdxtail, cdytail, cdztail;
  Real negate, negatetail;
  INEXACT Real axby7, bxcy7, axcy7, bxay7, cxby7, cxay7;
  Real axby[8], bxcy[8], axcy[8], bxay[8], cxby[8], cxay[8];
  Real temp16[16], temp32[32], temp32t[32];
  int temp16len, temp32len, temp32tlen;
  Real adet[64], bdet[64], cdet[64];
  int alen, blen, clen;
  Real abdet[128];
  int ablen;
  Real deter[192];
  int deterlen;

  INEXACT Real bvirt;
  Real avirt, bround, around;
  INEXACT Real c;
  INEXACT Real abig;
  Real a0hi, a0lo, a1hi, a1lo, bhi, blo;
  Real err1, err2, err3;
  INEXACT Real _i, _j, _k, _l, _m, _n;
  Real _0, _1, _2;

  Two_Diff(pa[0], pd[0], adx, adxtail);
  Two_Diff(pa[1], pd[1], ady, adytail);
  Two_Diff(pa[2], pd[2], adz, adztail);
  Two_Diff(pb[0], pd[0], bdx, bdxtail);
  Two_Diff(pb[1], pd[1], bdy, bdytail);
  Two_Diff(pb[2], pd[2], bdz, bdztail);
  Two_Diff(pc[0], pd[0], cdx, cdxtail);
  Two_Diff(pc[1], pd[1], cdy, cdytail);
  Two_Diff(pc[2], pd[2], cdz, cdztail);

  Two_Two_Product(adx, adxtail, bdy, bdytail,
                  axby7, axby[6], axby[5], axby[4],
                  axby[3], axby[2], axby[1], axby[0]);
  axby[7] = axby7;
  negate = -ady;
  negatetail = -adytail;
  Two_Two_Product(bdx, bdxtail, negate, negatetail,
                  bxay7, bxay[6], bxay[5], bxay[4],
                  bxay[3], bxay[2], bxay[1], bxay[0]);
  bxay[7] = bxay7;
  Two_Two_Product(bdx, bdxtail, cdy, cdytail,
                  bxcy7, bxcy[6], bxcy[5], bxcy[4],
                  bxcy[3], bxcy[2], bxcy[1], bxcy[0]);
  bxcy[7] = bxcy7;
  negate = -bdy;
  negatetail = -bdytail;
  Two_Two_Product(cdx, cdxtail, negate, negatetail,
                  cxby7, cxby[6], cxby[5], cxby[4],
                  cxby[3], cxby[2], cxby[1], cxby[0]);
  cxby[7] = cxby7;
  Two_Two_Product(cdx, cdxtail, ady, adytail,
                  cxay7, cxay[6], cxay[5], cxay[4],
                  cxay[3], cxay[2], cxay[1], cxay[0]);
  cxay[7] = cxay7;
  negate = -cdy;
  negatetail = -cdytail;
  Two_Two_Product(adx, adxtail, negate, negatetail,
                  axcy7, axcy[6], axcy[5], axcy[4],
                  axcy[3], axcy[2], axcy[1], axcy[0]);
  axcy[7] = axcy7;

  temp16len = fast_expansion_sum_zeroelim(8, bxcy, 8, cxby, temp16);
  temp32len = scale_expansion_zeroelim(temp16len, temp16, adz, temp32);
  temp32tlen = scale_expansion_zeroelim(temp16len, temp16, adztail, temp32t);
  alen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t,
                                     adet);

  temp16len = fast_expansion_sum_zeroelim(8, cxay, 8, axcy, temp16);
  temp32len = scale_expansion_zeroelim(temp16len, temp16, bdz, temp32);
  temp32tlen = scale_expansion_zeroelim(temp16len, temp16, bdztail, temp32t);
  blen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t,
                                     bdet);

  temp16len = fast_expansion_sum_zeroelim(8, axby, 8, bxay, temp16);
  temp32len = scale_expansion_zeroelim(temp16len, temp16, cdz, temp32);
  temp32tlen = scale_expansion_zeroelim(temp16len, temp16, cdztail, temp32t);
  clen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t,
                                     cdet);

  ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
  deterlen = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, deter);

  return deter[deterlen - 1];
}
#endif

static Real
orient3dadapt(Real const *pa, Real const *pb, Real const *pc, Real const *pd, Real permanent)
{
  INEXACT Real adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz;
  Real det, errbound;

  INEXACT Real bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
  Real bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
  Real bc[4], ca[4], ab[4];
  INEXACT Real bc3, ca3, ab3;
  Real adet[8], bdet[8], cdet[8];
  int alen, blen, clen;
  Real abdet[16];
  int ablen;
  Real *finnow, *finother, *finswap;
  Real fin1[192], fin2[192];
  int finlength;

  Real adxtail, bdxtail, cdxtail;
  Real adytail, bdytail, cdytail;
  Real adztail, bdztail, cdztail;
  INEXACT Real at_blarge, at_clarge;
  INEXACT Real bt_clarge, bt_alarge;
  INEXACT Real ct_alarge, ct_blarge;
  Real at_b[4], at_c[4], bt_c[4], bt_a[4], ct_a[4], ct_b[4];
  int at_blen, at_clen, bt_clen, bt_alen, ct_alen, ct_blen;
  INEXACT Real bdxt_cdy1, cdxt_bdy1, cdxt_ady1;
  INEXACT Real adxt_cdy1, adxt_bdy1, bdxt_ady1;
  Real bdxt_cdy0, cdxt_bdy0, cdxt_ady0;
  Real adxt_cdy0, adxt_bdy0, bdxt_ady0;
  INEXACT Real bdyt_cdx1, cdyt_bdx1, cdyt_adx1;
  INEXACT Real adyt_cdx1, adyt_bdx1, bdyt_adx1;
  Real bdyt_cdx0, cdyt_bdx0, cdyt_adx0;
  Real adyt_cdx0, adyt_bdx0, bdyt_adx0;
  Real bct[8], cat[8], abt[8];
  int bctlen, catlen, abtlen;
  INEXACT Real bdxt_cdyt1, cdxt_bdyt1, cdxt_adyt1;
  INEXACT Real adxt_cdyt1, adxt_bdyt1, bdxt_adyt1;
  Real bdxt_cdyt0, cdxt_bdyt0, cdxt_adyt0;
  Real adxt_cdyt0, adxt_bdyt0, bdxt_adyt0;
  Real u[4], v[12], w[16];
  INEXACT Real u3;
  int vlength, wlength;
  Real negate;

  INEXACT Real bvirt;
  Real avirt, bround, around;
  INEXACT Real c;
  INEXACT Real abig;
  Real ahi, alo, bhi, blo;
  Real err1, err2, err3;
  INEXACT Real _i, _j, _k;
  Real _0;

  adx = (Real) (pa[0] - pd[0]);
  bdx = (Real) (pb[0] - pd[0]);
  cdx = (Real) (pc[0] - pd[0]);
  ady = (Real) (pa[1] - pd[1]);
  bdy = (Real) (pb[1] - pd[1]);
  cdy = (Real) (pc[1] - pd[1]);
  adz = (Real) (pa[2] - pd[2]);
  bdz = (Real) (pb[2] - pd[2]);
  cdz = (Real) (pc[2] - pd[2]);

  Two_Product(bdx, cdy, bdxcdy1, bdxcdy0);
  Two_Product(cdx, bdy, cdxbdy1, cdxbdy0);
  Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
  bc[3] = bc3;
  alen = scale_expansion_zeroelim(4, bc, adz, adet);

  Two_Product(cdx, ady, cdxady1, cdxady0);
  Two_Product(adx, cdy, adxcdy1, adxcdy0);
  Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
  ca[3] = ca3;
  blen = scale_expansion_zeroelim(4, ca, bdz, bdet);

  Two_Product(adx, bdy, adxbdy1, adxbdy0);
  Two_Product(bdx, ady, bdxady1, bdxady0);
  Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
  ab[3] = ab3;
  clen = scale_expansion_zeroelim(4, ab, cdz, cdet);

  ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
  finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

  det = estimate(finlength, fin1);
  errbound = o3derrboundB * permanent;
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  Two_Diff_Tail(pa[0], pd[0], adx, adxtail);
  Two_Diff_Tail(pb[0], pd[0], bdx, bdxtail);
  Two_Diff_Tail(pc[0], pd[0], cdx, cdxtail);
  Two_Diff_Tail(pa[1], pd[1], ady, adytail);
  Two_Diff_Tail(pb[1], pd[1], bdy, bdytail);
  Two_Diff_Tail(pc[1], pd[1], cdy, cdytail);
  Two_Diff_Tail(pa[2], pd[2], adz, adztail);
  Two_Diff_Tail(pb[2], pd[2], bdz, bdztail);
  Two_Diff_Tail(pc[2], pd[2], cdz, cdztail);

  if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
      && (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)
      && (adztail == 0.0) && (bdztail == 0.0) && (cdztail == 0.0)) {
    return det;
  }

  errbound = o3derrboundC * permanent + resulterrbound * Absolute(det);
  det += (adz * ((bdx * cdytail + cdy * bdxtail)
                 - (bdy * cdxtail + cdx * bdytail))
          + adztail * (bdx * cdy - bdy * cdx))
       + (bdz * ((cdx * adytail + ady * cdxtail)
                 - (cdy * adxtail + adx * cdytail))
          + bdztail * (cdx * ady - cdy * adx))
       + (cdz * ((adx * bdytail + bdy * adxtail)
                 - (ady * bdxtail + bdx * adytail))
          + cdztail * (adx * bdy - ady * bdx));
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  finnow = fin1;
  finother = fin2;

  if (adxtail == 0.0) {
    if (adytail == 0.0) {
      at_b[0] = 0.0;
      at_blen = 1;
      at_c[0] = 0.0;
      at_clen = 1;
    } else {
      negate = -adytail;
      Two_Product(negate, bdx, at_blarge, at_b[0]);
      at_b[1] = at_blarge;
      at_blen = 2;
      Two_Product(adytail, cdx, at_clarge, at_c[0]);
      at_c[1] = at_clarge;
      at_clen = 2;
    }
  } else {
    if (adytail == 0.0) {
      Two_Product(adxtail, bdy, at_blarge, at_b[0]);
      at_b[1] = at_blarge;
      at_blen = 2;
      negate = -adxtail;
      Two_Product(negate, cdy, at_clarge, at_c[0]);
      at_c[1] = at_clarge;
      at_clen = 2;
    } else {
      Two_Product(adxtail, bdy, adxt_bdy1, adxt_bdy0);
      Two_Product(adytail, bdx, adyt_bdx1, adyt_bdx0);
      Two_Two_Diff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0,
                   at_blarge, at_b[2], at_b[1], at_b[0]);
      at_b[3] = at_blarge;
      at_blen = 4;
      Two_Product(adytail, cdx, adyt_cdx1, adyt_cdx0);
      Two_Product(adxtail, cdy, adxt_cdy1, adxt_cdy0);
      Two_Two_Diff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0,
                   at_clarge, at_c[2], at_c[1], at_c[0]);
      at_c[3] = at_clarge;
      at_clen = 4;
    }
  }
  if (bdxtail == 0.0) {
    if (bdytail == 0.0) {
      bt_c[0] = 0.0;
      bt_clen = 1;
      bt_a[0] = 0.0;
      bt_alen = 1;
    } else {
      negate = -bdytail;
      Two_Product(negate, cdx, bt_clarge, bt_c[0]);
      bt_c[1] = bt_clarge;
      bt_clen = 2;
      Two_Product(bdytail, adx, bt_alarge, bt_a[0]);
      bt_a[1] = bt_alarge;
      bt_alen = 2;
    }
  } else {
    if (bdytail == 0.0) {
      Two_Product(bdxtail, cdy, bt_clarge, bt_c[0]);
      bt_c[1] = bt_clarge;
      bt_clen = 2;
      negate = -bdxtail;
      Two_Product(negate, ady, bt_alarge, bt_a[0]);
      bt_a[1] = bt_alarge;
      bt_alen = 2;
    } else {
      Two_Product(bdxtail, cdy, bdxt_cdy1, bdxt_cdy0);
      Two_Product(bdytail, cdx, bdyt_cdx1, bdyt_cdx0);
      Two_Two_Diff(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0,
                   bt_clarge, bt_c[2], bt_c[1], bt_c[0]);
      bt_c[3] = bt_clarge;
      bt_clen = 4;
      Two_Product(bdytail, adx, bdyt_adx1, bdyt_adx0);
      Two_Product(bdxtail, ady, bdxt_ady1, bdxt_ady0);
      Two_Two_Diff(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0,
                  bt_alarge, bt_a[2], bt_a[1], bt_a[0]);
      bt_a[3] = bt_alarge;
      bt_alen = 4;
    }
  }
  if (cdxtail == 0.0) {
    if (cdytail == 0.0) {
      ct_a[0] = 0.0;
      ct_alen = 1;
      ct_b[0] = 0.0;
      ct_blen = 1;
    } else {
      negate = -cdytail;
      Two_Product(negate, adx, ct_alarge, ct_a[0]);
      ct_a[1] = ct_alarge;
      ct_alen = 2;
      Two_Product(cdytail, bdx, ct_blarge, ct_b[0]);
      ct_b[1] = ct_blarge;
      ct_blen = 2;
    }
  } else {
    if (cdytail == 0.0) {
      Two_Product(cdxtail, ady, ct_alarge, ct_a[0]);
      ct_a[1] = ct_alarge;
      ct_alen = 2;
      negate = -cdxtail;
      Two_Product(negate, bdy, ct_blarge, ct_b[0]);
      ct_b[1] = ct_blarge;
      ct_blen = 2;
    } else {
      Two_Product(cdxtail, ady, cdxt_ady1, cdxt_ady0);
      Two_Product(cdytail, adx, cdyt_adx1, cdyt_adx0);
      Two_Two_Diff(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0,
                   ct_alarge, ct_a[2], ct_a[1], ct_a[0]);
      ct_a[3] = ct_alarge;
      ct_alen = 4;
      Two_Product(cdytail, bdx, cdyt_bdx1, cdyt_bdx0);
      Two_Product(cdxtail, bdy, cdxt_bdy1, cdxt_bdy0);
      Two_Two_Diff(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0,
                   ct_blarge, ct_b[2], ct_b[1], ct_b[0]);
      ct_b[3] = ct_blarge;
      ct_blen = 4;
    }
  }

  bctlen = fast_expansion_sum_zeroelim(bt_clen, bt_c, ct_blen, ct_b, bct);
  wlength = scale_expansion_zeroelim(bctlen, bct, adz, w);
  finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                          finother);
  finswap = finnow; finnow = finother; finother = finswap;

  catlen = fast_expansion_sum_zeroelim(ct_alen, ct_a, at_clen, at_c, cat);
  wlength = scale_expansion_zeroelim(catlen, cat, bdz, w);
  finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                          finother);
  finswap = finnow; finnow = finother; finother = finswap;

  abtlen = fast_expansion_sum_zeroelim(at_blen, at_b, bt_alen, bt_a, abt);
  wlength = scale_expansion_zeroelim(abtlen, abt, cdz, w);
  finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                          finother);
  finswap = finnow; finnow = finother; finother = finswap;

  if (adztail != 0.0) {
    vlength = scale_expansion_zeroelim(4, bc, adztail, v);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (bdztail != 0.0) {
    vlength = scale_expansion_zeroelim(4, ca, bdztail, v);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (cdztail != 0.0) {
    vlength = scale_expansion_zeroelim(4, ab, cdztail, v);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }

  if (adxtail != 0.0) {
    if (bdytail != 0.0) {
      Two_Product(adxtail, bdytail, adxt_bdyt1, adxt_bdyt0);
      Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (cdztail != 0.0) {
        Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
    if (cdytail != 0.0) {
      negate = -adxtail;
      Two_Product(negate, cdytail, adxt_cdyt1, adxt_cdyt0);
      Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (bdztail != 0.0) {
        Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
  }
  if (bdxtail != 0.0) {
    if (cdytail != 0.0) {
      Two_Product(bdxtail, cdytail, bdxt_cdyt1, bdxt_cdyt0);
      Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (adztail != 0.0) {
        Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
    if (adytail != 0.0) {
      negate = -bdxtail;
      Two_Product(negate, adytail, bdxt_adyt1, bdxt_adyt0);
      Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (cdztail != 0.0) {
        Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
  }
  if (cdxtail != 0.0) {
    if (adytail != 0.0) {
      Two_Product(cdxtail, adytail, cdxt_adyt1, cdxt_adyt0);
      Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (bdztail != 0.0) {
        Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
    if (bdytail != 0.0) {
      negate = -cdxtail;
      Two_Product(negate, bdytail, cdxt_bdyt1, cdxt_bdyt0);
      Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (adztail != 0.0) {
        Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
  }

  if (adztail != 0.0) {
    wlength = scale_expansion_zeroelim(bctlen, bct, adztail, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (bdztail != 0.0) {
    wlength = scale_expansion_zeroelim(catlen, cat, bdztail, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (cdztail != 0.0) {
    wlength = scale_expansion_zeroelim(abtlen, abt, cdztail, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }

  return finnow[finlength - 1];
}

Real
robust_orient3d(Real const *pa, Real const *pb, Real const *pc, Real const *pd)
{
  Real adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz;
  Real bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
  Real det;
  Real permanent, errbound;

  set_extended_precision (0);
  
  adx = pa[0] - pd[0];
  bdx = pb[0] - pd[0];
  cdx = pc[0] - pd[0];
  ady = pa[1] - pd[1];
  bdy = pb[1] - pd[1];
  cdy = pc[1] - pd[1];
  adz = pa[2] - pd[2];
  bdz = pb[2] - pd[2];
  cdz = pc[2] - pd[2];

  bdxcdy = bdx * cdy;
  cdxbdy = cdx * bdy;

  cdxady = cdx * ady;
  adxcdy = adx * cdy;

  adxbdy = adx * bdy;
  bdxady = bdx * ady;

  det = adz * (bdxcdy - cdxbdy) 
      + bdz * (cdxady - adxcdy)
      + cdz * (adxbdy - bdxady);

  permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * Absolute(adz)
            + (Absolute(cdxady) + Absolute(adxcdy)) * Absolute(bdz)
            + (Absolute(adxbdy) + Absolute(bdxady)) * Absolute(cdz);
  errbound = o3derrboundA * permanent;
  if ((det > errbound) || (-det > errbound)) {
    set_extended_precision (1);    
    return det;
  }


  {
    Real adapt = orient3dadapt(pa, pb, pc, pd, permanent);
    set_extended_precision (1);
    return adapt;
  }
}


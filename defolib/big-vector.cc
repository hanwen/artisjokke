#include <stdio.h> 
#include <string.h>
#include <math.h>
#include <gsl/gsl_cblas.h>

#include "big-vector.hh"

long long vector_flop_count = 0LL ;

void
completize_big_vector (Array<Real> *a, int n)
{
  int prevsize = a->size();

  if (prevsize > n)
    return ;

  a->set_size (n);
  Real * aptr = a->unsafe_access_array();
  memset (aptr + prevsize,  0x0, (n - prevsize)* sizeof(Real));
}

extern "C" {

void
big_vector_axpy (Real * dest,
		  Real a, Real const *x, Real const *y, int n)
{
  for (; n--; )
    {
      *dest++ = a * (*x++  ) + (*y++);
    }
  vector_flop_count += 2*n;
}

Real
big_vector_ip (Real const * p1, Real const * p2, int n)
{
  vector_flop_count += 2*n;    

  Real  const * e1 = p1 + n;

  Real i =0.0;
  while (p1 < e1)
    {
      i += (*p1++) * (*p2++);
    }
 
  return i;
}


void
big_vector_print (Real const *ar, int n)
{
  for (int i = 0 ;  i < n ; i++)
    printf ("%s %lf",  ((i && (i % 3)==0)? " | " : ""),  ar[i]);
  
  printf ("\n");
}




void
big_vector_nullify (Real *ar, int n)
{
  memset (ar, 0, n*sizeof (Real));
}


void
big_vector_add (Real  * dest, Real const * p1, Real const * p2, int n)
{
  for (; n--; )
    {
      *dest ++ = *p1++ + *p2++;
    }
  vector_flop_count += n;    
}

void
big_vector_subtract (Real  * dest, Real const * p1, Real const * p2, int n)
{
  for (; n--; )
    {
      *dest ++ = *p1++ - *p2++;
    }
  vector_flop_count += n;    
}

void
big_vector_pointwise_multiply (Real  * dest, Real const * p1, Real const * p2, int n)
{
 for (; n--; )
    {
      *dest ++ = *p1++ * *p2++;
    }
  vector_flop_count += n;    
}



void
big_vector_negate (Real *dest, Real const*src, int n)
{
 for (; n--; )
   *dest++ = - *src ++;
}

void
big_vector_copy (Real * d , Real const *s , int n)
{
  memcpy  (d, s, n*sizeof (Real));
}

void 
big_vector_pointwise_inverse (Real *d, Real const * s, int n)
{
  for (; n--; )
    {
      *d++ = 1 / (*s++);
    }
  vector_flop_count += n;    
}

void
big_vector_fill (Real * d, Real f , int n)
{
  for (;n--;)
    *d++ = f;
}

void
big_vector_scale (Real * d, Real s, Real *src, int n)
{
  for (; n--; )
    {
      *d++ = s * (*src++);
    }
  vector_flop_count += n;    
}

/*
  Copy SRC to DEST, leaving alone those entries that have ONOFF set to true
 */
void
big_vector_partial_copy (Real *dest, bool const*onoff, Real const *src,int n)
{
  for (; n--; )
    {
      if (! *onoff++)
	{
	  *dest = *src;
	}

      dest++;
      src ++;
    } 
}

/*
  Copy SRC to DEST, setting to 0.0 those entries that have ONOFF set to true.
 */
void
big_vector_partial_nullify (Real * dest, Real const *src, bool const * onoff,  int n)
{
  for (; n--; )
    {
      *dest++ = *onoff++ ? 0.0 : (*src );
      src ++;
    }
}

/*
  Set to 0.0 all entries of SRC for which ONOFF == false.
 */
void
big_vector_partial_inv_nullify (Real * dest, Real const *src, bool const * onoff,  int n)
{
  for (; n--; )
    {
      *dest++ = *onoff++ ?  (*src ) :  0.0 ;
      src ++;
    }
}

Real
big_vector_length_squared (Real const * p1, int n)
{
  Real  const * e1 = p1 + n;

  Real i =0.0;
  while (p1 < e1)
    {
      i += (*p1) * (*p1);

      p1 ++;
    }

  vector_flop_count += 2*n;
  return i;
}


Real
big_vector_length (Real const * s, int n)
{
  return sqrt (big_vector_length_squared (s,n) * n / 3);
}

Real
big_vector_distance (Real const * a, Real  const *b, int n)
{
  Real s =0.0;
  while (n--)
    {
      Real d =  (*a++ - *b++);
      s+= d*d;
    }


  vector_flop_count += 3*n;  
  return sqrt (s);
}

Real
big_vector_max_distance (Real const * a, Real  const *b, int n)
{
  Real m =0.0;
  while (n--)
    {
      Real d = fabs (*a++ - *b++);

      m =  m  >? d;
    }

  vector_flop_count += n;  		//  ? 
  return m;
}

Real
big_vector_1_norm (Real const * a , int n)
{
  Real norm = 0.0;
  while (n--)
    norm += fabs (*a++);
  
  vector_flop_count += n;    
  return norm;
}

bool
big_vector_sane (Real const * a , int n)
{
  bool b = true;
  while (b && n--)
    {
      b = b && !isinf (*a) && !isnan (*a);
      a++;
    }
  return b;
}


bool
big_vector_is_null (Real const * a , int n)
{
  while (n--)
    {
      if (*a++)
	return false;
    }
  return true;
}


} // extern "C"

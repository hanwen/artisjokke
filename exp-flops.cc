/*

Benchmark computation speed of a platform, and derive how expensive an
exp() is.


  compile with

  g++ -O2 -fno-prefetch-loop-arrays -march=pentium3 -o exp-flops exp-flops.cc

  g++ -finline-functions -funroll-loops  -O3 -mcpu=pentium3 -mfpmath=sse,387 -msse  -fno-prefetch-loop-arrays -march=pentium3 -o exp-flops exp-flops.cc   

*/

#include <math.h>
#include <time.h>
#include <stdio.h>

 typedef float Real;
//   typedef double  Real;
int
main()
{
  int O = 40960;
  int O2 = O/5;

  int NO_ARRAYS = 1; 
  /*
    1024 * 8 = 8092 bytes.
  */
  int I = 1024 * 10;
  
  Real f[NO_ARRAYS][I];
  for (int j = 0; j < NO_ARRAYS; j++)
    for (int i = 0; i < I; i++)
      f[j][i] = i * (j*1.0/ (I*O));
  
  clock_t t1 = clock();

  for (int j = 0; j < O; j++)
    {
      int k = j % NO_ARRAYS;
#if 1
      Real  *fp = f[k];
      Real  * end = fp + I;
      
      for (; fp < end; fp++)
	*fp = (1.0 - *fp);
#else
      for  (int i = 0; i< I; i++)
	f[k][i] = (1.0 + f[k][i]);
      
#endif
    }
  clock_t t2 = clock ();

  Real secs = (t2 -t1)*1.0/ CLOCKS_PER_SEC;
  
  printf ("%d ops %lf secs, %5.2lf mflops, throughput %5.1lf mb/s\n", I *O,  secs, 1e-6 *I*O/secs,
	  1e-6 * sizeof (Real) * I * O / secs);
  
  t1 = clock();
  for (int j = 0; j < O2; j++)
    {
      int k = j % NO_ARRAYS;
      for (int i = 0; i < I; i++)
	{
	  //	  f[k][i] = log (fabs (f[k][i]));
	  f[k][i] = exp (- (f[k][i])); 	  
	}
    }
  t2 = clock();

  Real secs2 = (t2 -t1)*1.0/ CLOCKS_PER_SEC;
  
  printf ("%d ops %lf secs, %5.2lf mflops \n", I*O2,  secs2, 1e-6 *I*O2/secs2 );
  printf ("Cost of an exp: %lf flops\n",  I * O * secs2 / (I * O2* secs) );
}

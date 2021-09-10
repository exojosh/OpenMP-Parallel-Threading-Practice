//OpenMP version.  Edit and submit only this file.
/* Enter your details below
 * Name : Josh Limon
 * UCLA ID : 804-984-257
 * Email : jjlemonheads@gmail.com
 */

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "utils.h"

double work_it_par(long *old, long *new, long *super, long *simple, long *fibonacci) {
  int i, j, k;
  int u, v, w;
  int ton = 0;
  long compute_it, moving_average;
  double pi, pi2, x , y, sum, step = 0.0;
  long dot_product=0;
  long nCirc=0;
  long aggregate=1.0;
  double r=1.0;
  int was_smart = 16;

  double tmp, tmp2;
  #pragma omp parallel for 
  for(i=0; i<DIM-1;i++)
  {
    super[i] += simple[i];
  }
  #pragma omp parallel private(tmp)
 {tmp = 0.0;
  #pragma omp for lastprivate(i, moving_average)
  for(i=0; i<DIM-1;i++)
  {
    tmp += super[i]*simple[i];
    //#prgama omp parallel lastprivate(moving_average)    
    moving_average = 0;
    for(ton=i;ton<DIM-1-WINDOW_SIZE;ton++)
    {
      moving_average += simple[ton];
    }
  }
   #pragma omp critical
	dot_product += tmp;
}

  int a_secret = 5;
  fibonacci[0] = 1;
  fibonacci[1] = 1;
  #pragma omp parallel for
  for(i=2; i<DIM-1;i++)
  {
    fibonacci[i]=fibonacci[i-1]+fibonacci[i-2];
    if(i==3)
    {
      printf("\n A secret is: %d",obfuscate_obfuscate_obfuscate(a_secret));
    }
  }

  step = 1.0 / NUM_STEPS;
  #pragma omp parallel private(tmp)
  { tmp = 0.0;
  #pragma omp for private (x)
  for (i=0;i<NUM_STEPS; i++)
  {
    x = (i+0.5)*step;
    tmp += 4.0/(1.0+x*x);
  }
  #pragma omp critical
	sum += tmp;
  }
  pi = step * sum;
  printf("\n %d trials, Riemann flavored pi is %f \n",NUM_STEPS, pi); 

  for(i = 0;i<NUM_TRIALS; i++)
  {
    x = (random()%10000000)/10000000.0;
    y = (random()%10000000)/10000000.0;
    if (( x*x + y*y) <= r*r) {
      nCirc++;
    }
  }
  pi2 = 4.0 * ((double)nCirc/(double)NUM_TRIALS);
  printf("\n %d trials, Monte-Carlo flavored pi is %f \n",NUM_TRIALS, pi2);


  //long func = we_need_the_func() / gimmie_the_func();
  #pragma omp parallel private(tmp)
{ tmp = 0.0;
  #pragma omp for private(j, k)
  for (i=1; i<DIM-1; i++) {
    for (j=1; j<DIM-1; j++) {
      for (k=1; k<DIM-1; k++) {
        compute_it = old[i*DIM*DIM+j*DIM+k] * we_need_the_func();
        tmp+= compute_it / gimmie_the_func();
      }
    }
  }
  #pragma omp critical
  aggregate += tmp;
}
  printf("AGGR:%ld\n",aggregate);

#pragma omp parallel private(tmp2)
{tmp2 = 0;
#pragma omp for private(j, k)
  for (i=1; i<DIM-1; i++) {
    for (j=1; j<DIM-1; j++) {
      for (k=1; k<DIM-1; k++) {
        //new[i*DIM*DIM+j*DIM+k]=0;
        tmp2+=old[(i-1)*DIM*DIM+(j-1)*DIM+(k-1)];
tmp2+=old[(i-1)*DIM*DIM+(j-1)*DIM+(k)];
tmp2+=old[(i-1)*DIM*DIM+(j-1)*DIM+(k+1)];
tmp2+=old[(i-1)*DIM*DIM+(j)*DIM+(k-1)];
tmp2+=old[(i-1)*DIM*DIM+(j)*DIM+(k)];
tmp2+=old[(i-1)*DIM*DIM+(j)*DIM+(k+1)];
tmp2+=old[(i-1)*DIM*DIM+(j+1)*DIM+(k-1)];
tmp2+=old[(i-1)*DIM*DIM+(j+1)*DIM+(k)];
tmp2+=old[(i-1)*DIM*DIM+(j+1)*DIM+(k+1)];
tmp2+=old[(i)*DIM*DIM+(j-1)*DIM+(k-1)];
tmp2+=old[(i)*DIM*DIM+(j-1)*DIM+(k)];
tmp2+=old[(i)*DIM*DIM+(j-1)*DIM+(k+1)];
tmp2+=old[(i)*DIM*DIM+(j)*DIM+(k-1)];
tmp2+=old[(i)*DIM*DIM+(j)*DIM+(k)];
tmp2+=old[(i)*DIM*DIM+(j)*DIM+(k+1)];
tmp2+=old[(i)*DIM*DIM+(j+1)*DIM+(k-1)];
tmp2+=old[(i)*DIM*DIM+(j+1)*DIM+(k)];
tmp2+=old[(i)*DIM*DIM+(j+1)*DIM+(k+1)];
tmp2+=old[(i+1)*DIM*DIM+(j-1)*DIM+(k-1)];
tmp2+=old[(i+1)*DIM*DIM+(j-1)*DIM+(k)];
tmp2+=old[(i+1)*DIM*DIM+(j-1)*DIM+(k+1)];
tmp2+=old[(i+1)*DIM*DIM+(j)*DIM+(k-1)];
tmp2+=old[(i+1)*DIM*DIM+(j)*DIM+(k)];
tmp2+=old[(i+1)*DIM*DIM+(j)*DIM+(k+1)];
tmp2+=old[(i+1)*DIM*DIM+(j+1)*DIM+(k-1)];
tmp2+=old[(i+1)*DIM*DIM+(j+1)*DIM+(k)];
tmp2+=old[(i+1)*DIM*DIM+(j+1)*DIM+(k+1)];         
         tmp2 /=27;
        new[i*DIM*DIM+j*DIM+k] = tmp2;
	tmp2 = 0.0;
      }
    }
  }
}

#pragma omp parallel private(tmp, tmp2)
{tmp = 0.0;
tmp2 = 0.0;
#pragma omp for private(j, k, u)
  for (i=1; i<DIM-1; i++) {
    for (j=1; j<DIM-1; j++) {
      for (k=1; k<DIM-1; k++) {
        u=(new[i*DIM*DIM+j*DIM+k]/100);
        if (u<=0) tmp++;
        if (u>=9) tmp2++;
      }
    }
  }
#pragma omp critical
{
 histogrammy[0]+= tmp;
 histogrammy[9]+= tmp2;
}}
  return (double) (dot_product+moving_average+pi+pi2);


}

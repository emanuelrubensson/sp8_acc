#include <iostream>
#include <iomanip>
#include <malloc.h>

#include <cassert>
#include "sp8_eval.h"
int main(int argc, char* const  argv[]){
  int n=2;
  double *A=(double *)malloc(n*n*sizeof(double));
  double *mem1=(double *)malloc(n*n*sizeof(double));
  double *mem2=(double *)malloc(n*n*sizeof(double));
  double *mem3=(double *)malloc(n*n*sizeof(double));
  double *mem4=(double *)malloc(n*n*sizeof(double));

  double *b=(double *)malloc(9*sizeof(double));
  // Matrix: [0.1 0.2 ; 0.3 0.4];
  A[0]=0.1;
  A[1]=0.2;
  A[2]=0.3;
  A[3]=0.4;


  // polynomial: 1/0!+x*1/1!+x^2*1/2!
  b[0]=1.0;
  for (unsigned int ind = 1; ind<9; ind++)
    b[ind]=b[ind-1]/ind;

  // evaluate it
  sp8_eval(b,n,A, mem1, mem2, mem3);

  printf("%f %f \n %f %f\n",A[0],A[1],A[2],A[3]);
  // Should be:
  //   1.14209   0.260351
  //   0.390526  1.53262

}

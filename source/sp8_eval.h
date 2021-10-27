#include <vector>
#include <cmath>
#include <complex>

void sp8_eval(double *b,int n,
	      double *A, double *mem1, double *mem2, double *mem3){
  double b0=b[0];
  double b1=b[1];
  double b2=b[2];
  double b3=b[3];
  double b4=b[4];
  double b5=b[5];
  double b6=b[6];
  double b7=b[7];
  double b8=b[8];

  if (b8<0){
    // Avoid complex numbers by switching evaluating with minus sign and
    // subsequently reversing
    for (unsigned int ind = 0; ind<9; ind++)
      b[ind]=-b[ind];
    sp8_eval(b,n,A,mem1,mem2,mem3);
    for (unsigned int ind = 0; ind<n*n; ind++)
      A[ind]=-A[ind];
    return;
  }


  // plus minus: Two options possible.
  // "Correct"  option leads to a e2_num_sqrt > 0 in code below?
  double c4 = sqrt(b8);
  double c3 = b7 / (2.0 * c4);
  double d2_plus_e2 = (b6 - c3*c3) / c4;
  double d1 = (b5 - c3 * d2_plus_e2) / c4;

  double e2_num_sqrt =  std::pow(d1 - (c3 / c4) * d2_plus_e2, 2) +
    4.0 * (c3 / c4) * (b3 + (std::pow(c3,2) / c4) * d1 - (c3 / c4) * b4);

  // If e2_num_sqrt<0 we may want to try other sign in c4 ?
  // Here are two options on the square roor: plus minus
  double e2_num = (c3 / c4) * d2_plus_e2 - d1 + sqrt(e2_num_sqrt);
  double e2 = e2_num / (2 * c3 / c4);
  double d2 = d2_plus_e2 - e2;
  double f2 = b2;
  double f1 = b1;
  double f0 = b0;

  double e0 = (b3 - d1 * e2) / c3; // Not explicitly documented by sastre?


  // Stupid debug output
  printf("%f\n",e0);

  // Actual evaluation
  //  A2=A^2;
  //  y02=A2*(c4*A2+c3*A);
  //  y12=(y02+d2*A2+d1*A)*(y02+e2*A2)+e0*y02+f2*A2+f1*A+f0*eye(size(A));

  // Degopt form
  // [0 1  0  0  | 0 1   0 0
  //  0 0  1  0  | 0 x1 x2 0
  //  0 x3 x4 1  | 0 0  x5 1]
  //
  // [x6 x7 x8 x9 1]
  //x1=c3;
  //x2=c4;
  //x3=d1;
  //x4=d2;
  //x5=e2;
  //x6=f0;
  //x7=f1;
  //x8=f2;
  //x9=e0;






  // The matrix operations with memory slot usage
  // mem1 = A^2;
  // mem2 = c4*A2+c3*A;
  // mem3 = mem1 * mem2;
  // mem2 = d1*A+d2*mem1+mem3;
  // A  = f0*I+f1*A+f2*mem1+e0*mem3;
  // mem1 = e2*mem1+mem3;
  // mem3 = mem1*mem2;
  // A = A + mem3;


  return;

}

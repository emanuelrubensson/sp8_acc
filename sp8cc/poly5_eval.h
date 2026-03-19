#ifndef POLY5_EVAL_HEADER
#define POLY5_EVAL_HEADER
#include <vector>
#include <cmath>
namespace sp8 {

template<typename T_matrix>
  void poly_5_eval(std::vector<typename T_matrix::value_type> const & mc,
		   T_matrix & M1,
		   T_matrix & M2,
		   T_matrix & M3) {
  typedef typename T_matrix::value_type scalar_type;
  const scalar_type b0=mc[5], b1=mc[4], b2=mc[3];
  const scalar_type b3=mc[2], b4=mc[1], b5=mc[0];
  const scalar_type d5  = b5;
  const scalar_type b4s = b4/d5;
  const scalar_type b3s = b3/d5;
  const scalar_type b2s = b2/d5;  
  const scalar_type a2  = b4s/2.0 + 1.0;
  const scalar_type c1  = b4s/2.0 - 1.0;
  const scalar_type a1  = -0.5*( b2s+(a2*( b4s*b4s/4.0 - 1.0 -b3s )) );
  const scalar_type c0  = 1.0+b3s-a1-b4s*b4s/4.0;
  const scalar_type e1  = b1 - d5*a1*c0;
  const scalar_type e0  = b0;
  M2.scale_and_add(1.0, a2, M1);      // M2 = M2 + a2*M1       = a2*X + X2
  M3.multiply(M2, M1);                // M3 = M2*M1            = a2*X2 + X3
  M3.scale_and_add(1.0, a1, M1);      // M3 = M3 + a1*M1       = a1*X + a2*X2 + X3
  M2.scale_and_add(1.0, (c1-a2), M1); // M2 = M2 + (c1-a2)*M1  = c1*X + X2
  M2.add_scaled_identity(c0);         // M2 = M2 + c0*I        = c0*I + c1*X + X2
  M1.multiply(M3, M2, d5, e1);        // M1 = e1*M1 + d5*M3*M2 = e1*X +
                                      //  d5*(a1*X + a2*X2 + X3)*(c0*I + c1*X + X2)
  M1.add_scaled_identity(e0);
 }


template<typename T_matrix>
  void poly_5_eval(std::vector<typename T_matrix::value_type> const & mc,
		   T_matrix & X,
		   T_matrix & X2) {
  T_matrix T1;
  poly_5_eval(mc, X, X2, T1);
 }
 
} // end namespace sp8
#endif

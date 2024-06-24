#ifndef SASTRE_POLY8_EVAL_HEADER
#define SASTRE_POLY8_EVAL_HEADER
#include <vector>
#include <cmath>

namespace sp8 {
  /** If the input matrix to sastre_poly_8_eval is symmetric, then all
      intermediate matrices and the output matrix are also symmetric. In
      this case symmetric matrix storage and computation may be used
      throughout.
  */
  template<typename T_matrix>
  void sastre_poly_8_eval(std::vector<typename T_matrix::value_type> const & mc,
			  T_matrix & A,
			  T_matrix & A2,
			  T_matrix & M2,
			  T_matrix & M3) {
    typedef typename T_matrix::value_type scalar_type;
    const scalar_type b0=mc[8], b1=mc[7], b2=mc[6];
    const scalar_type b3=mc[5], b4=mc[4], b5=mc[3];
    const scalar_type b6=mc[2], b7=mc[1], b8=mc[0];
    if (b8<0) {
      // Avoid complex numbers by switching evaluating with minus sign and
      // subsequently reversing
      std::vector<scalar_type> mc_minus(9);
      for (unsigned int ind = 0; ind<9; ind++)      
	mc_minus[ind]=-mc[ind];
      sastre_poly_8_eval(mc_minus, A, A2, M2, M3);
      A.negate();
      return;
    }
  
    // plus minus: Two options possible.
    // "Correct"  option leads to a e2_num_sqrt > 0 in code below?
    const scalar_type c4 = std::sqrt(b8);
    const scalar_type c3 = b7 / (2.0 * c4);
    const scalar_type d2_plus_e2 = (b6 - c3*c3) / c4;
    const scalar_type d1 = (b5 - c3 * d2_plus_e2) / c4;
  
    const scalar_type e2_num_sqrt =  std::pow(d1 - (c3 / c4) * d2_plus_e2, 2) +
      4.0 * (c3 / c4) * (b3 + (std::pow(c3,2) / c4) * d1 - (c3 / c4) * b4);

    // If e2_num_sqrt<0 we may want to try other sign in c4 ?
    // Here are two options on the square roor: plus minus
    const scalar_type e2_num = (c3 / c4) * d2_plus_e2 - d1 + std::sqrt(e2_num_sqrt);
    const scalar_type e2 = e2_num / (2 * c3 / c4);
    const scalar_type d2 = d2_plus_e2 - e2;
    const scalar_type f2 = b2;
    const scalar_type f1 = b1;
    const scalar_type f0 = b0;
  
    const scalar_type e0 = (b3 - d1 * e2) / c3; // Not explicitly
    // documented by sastre?
    const scalar_type one = 1.0;
    M2 = A2;                                    // M2 = A2
    M2.scale_and_add(c4, c3, A);                // M2 = c4*M2 + c3*A
    M3.multiply(A2,M2);                         // M3 = A2*M2
    M2 = M3;                                    // M2 = M3
    M2.scale_and_add(one, d2, A2);              // M2 = M2 + d2*A2
    M2.scale_and_add(one, d1, A);               // M2 = M2 + d1*A
    A.scale_and_add(f1, f2, A2);                // A = f1*A + f2*A2
    A.scale_and_add(one, e0, M3);               // A = A + e0*M3
    A.add_scaled_identity(f0);                  // A = A+f0*I
    A2.scale_and_add(e2, one, M3);              // A2 = e2*A2 + M3
    M3.multiply(A2,M2);                         // M3 = A2*M2
    A.scale_and_add(one, one, M3);              // A = A + M3
  }

  template<typename T_matrix>
  void sastre_poly_8_eval(std::vector<typename T_matrix::value_type> const & mc,
			  T_matrix & A,
			  T_matrix & A2) {
    T_matrix M2;
    T_matrix M3;
    sastre_poly_8_eval(mc, A, A2, M2, M3);
  }


} // end namespace sp8
#endif

#ifndef POLY8_EVAL_HEADER
#define POLY8_EVAL_HEADER
#include <vector>
#include <cmath>
namespace sp8 {
  
  /**
     Modified version of Sastres scheme. 
     - Works also if b7 is zero. 
     - Does not lead to complex coefficients (and complex multiplication)
       for any polynomial coefficients.
     - b8 still has to be nonzero. This is the only restriction on the
       coefficients.
   */
  template<typename T_matrix>
    void poly_8_eval(std::vector<typename T_matrix::value_type> const & mc,
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
      poly_8_eval(mc_minus, A, A2, M2, M3);
      A.negate();
      return;
    }
    // plus minus std::sqrt(b8): Two options possible, both should work fine!
    const scalar_type c4  = std::sqrt(b8);
    const scalar_type c3  = b7 / (2.0 * c4);
    const scalar_type de2 = (b6 - c3*c3) / c4;
    const scalar_type de1 = (b5 - c3 * de2) / c4;
    // Extra degree of freedom d0 added to avoid complex coefficients
    // d0 chosen to make [expression under sqrt] = 1 (this should always
    // be possible if b8 != 0, which we have assumed anyway). Other
    // positive choices may also be fine. However, [expression under
    // sqrt] = 0 makes scheme break if de2 is near 0.
    const scalar_type d0 = (1-de2*de2+4*b4-4*c3*de1)/(4*c4);
    // plus minus : Two options possible, choosing + --> e2-d2 = 1
    // which is good due to later division by e2-d2.
    const scalar_type e2 = 0.5*(de2+1);
    const scalar_type d2 = 0.5*(de2-1);
    const scalar_type e1 = (c3*d0+de1*e2-b3); // /(e2-d2) (which is = 1)
    const scalar_type d1 = de1-e1;
    const scalar_type f2 = b2-(d0*e2+d1*e1);
    const scalar_type f1 = b1-d0*e1;
    const scalar_type f0 = b0;
    const scalar_type one = 1.0;
    M2 = A2;                                    // M2 = A2
    M2.scale_and_add(c4, c3, A);                // M2 = c4*M2 + c3*A
    M3.multiply(A2,M2);                         // M3 = A2*M2
    M2 = M3;                                    // M2 = M3
    M2.scale_and_add(one, d2, A2);              // M2 = M2 + d2*A2
    M2.scale_and_add(one, d1, A);               // M2 = M2 + d1*A    
    M2.add_scaled_identity(d0);                 // M2 = M2 + d0*I
    M3.scale_and_add(one, e1, A);               // M3 = M3 + e1*A
    A.scale_and_add(f1, f2, A2);                // A = f1*A + f2*A2
    A.add_scaled_identity(f0);                  // A = A+f0*I
    A2.scale_and_add(e2, one, M3);              // A2 = e2*A2 + M3
    M3.multiply(M2,A2);                         // M3 = M2*A2
    A.scale_and_add(one, one, M3);              // A = A + M3       
  }    

  template<typename T_matrix>
    void poly_8_eval(std::vector<typename T_matrix::value_type> const & mc,
		     T_matrix & A,
		     T_matrix & A2) {
    T_matrix M2;
    T_matrix M3;
    poly_8_eval(mc, A, A2, M2, M3);
  }
} // end namespace sp8
#endif

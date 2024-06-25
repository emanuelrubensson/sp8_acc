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
		     T_matrix & X,
		     T_matrix & X2,
		     T_matrix & T1,
		     T_matrix & T2) {
    typedef typename T_matrix::value_type scalar_type;
    const scalar_type b0=mc[8], b1=mc[7], b2=mc[6];
    const scalar_type b3=mc[5], b4=mc[4], b5=mc[3];
    const scalar_type b6=mc[2], b7=mc[1], b8=mc[0];
    const scalar_type f4 = b8;
    const scalar_type c1 = b7 / (2.0 * f4);
    const scalar_type t2 = b6/f4 - c1*c1;
    const scalar_type t1 = b5/f4 - c1*t2;
    // Extra degree of freedom d0 added to avoid complex coefficients
    // d0 chosen to make [expression under sqrt] = 1 (this should
    // always be possible if b8 != 0, which we have assumed
    // anyway). Other positive choices may also be fine, but our
    // choice leads to e2-d2=1.
    const scalar_type d0 = (1.0-t2*t2+4*b4/f4-4*c1*t1)/4.0;
    // plus minus : Two options possible for e2, choosing +
    // --> e2-d2=1 which is good due to later division by e2-d2.
    const scalar_type e2 = 0.5*(t2+1);
    const scalar_type d2 = 0.5*(t2-1);
    const scalar_type e1 = c1*d0+t1*e2-b3/f4; // /(e2-d2) (which is = 1)
    const scalar_type d1 = t1-e1;
    const scalar_type f2 = b2-f4*(d0*e2+d1*e1);
    const scalar_type f1 = b1-f4*d0*e1;
    const scalar_type f0 = b0;
    //    const scalar_type one = 1.0;
    T1 = X2;                                    // T1 = X2
    T1.scale_and_add(1.0, c1, X);               // T1 = T1 + c1*X
    T2.multiply(X2,T1);                         // T2 = X2*T1
    T1 = T2;                                    // T1 = T2
    T1.scale_and_add(1.0, d2, X2);              // T1 = T1 + d2*X2
    T1.scale_and_add(1.0, d1, X);               // T1 = T1 + d1*X
    T1.add_scaled_identity(d0);                 // T1 = T1 + d0*I
    T2.scale_and_add(1.0, e1, X);               // T2 = T2 + e1*X
    X.scale_and_add(f1, f2, X2);                // X = f1*X + f2*X2
    X.add_scaled_identity(f0);                  // X = X+f0*I
    X2.scale_and_add(e2, 1.0, T2);              // X2 = e2*X2 + T2
    T2.multiply(T1,X2);                         // T2 = T1*X2
    X.scale_and_add(1.0, f4, T2);               // X = X + f4*T2
  }    

  template<typename T_matrix>
    void poly_8_eval(std::vector<typename T_matrix::value_type> const & mc,
		     T_matrix & X,
		     T_matrix & X2) {
    T_matrix T1;
    T_matrix T2;
    poly_8_eval(mc, X, X2, T1, T2);
  }
} // end namespace sp8
#endif

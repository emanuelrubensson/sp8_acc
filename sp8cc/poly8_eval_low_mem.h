#ifndef POLY8_EVAL_LOW_MEM_HEADER
#define POLY8_EVAL_LOW_MEM_HEADER
#include <vector>
#include <cmath>
namespace sp8 {

template<typename T_matrix>
    void poly_8_eval_low_mem(std::vector<typename T_matrix::value_type> const & mc,
		     T_matrix & M1,
		     T_matrix & M2,
		     T_matrix & M3) {
    //------------------------------------------------------------------------------
    // Precompute c,d,e,f
    //------------------------------------------------------------------------------
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
    //------------------------------------------------------------------------------
    // Precompute r
    //------------------------------------------------------------------------------
    const scalar_type r1 = d1 - 0.5*c1*(d2-0.25*c1*c1)
    const scalar_type r2 = d2 -0.25*c1*c1
    const scalar_type r3 = e1 - d1 - 0.5*c1
    const scalar_type r4 = f_1 - f_2*(e1-d1)
    //------------------------------------------------------------------------------
    // Computation Stage
    //------------------------------------------------------------------------------
    // M2 already contains X^2        // line 3
    M2.scale_and_add(1.0,0.5*c1,M1);  // line 4
    M3.multiply(M2,M2)                // line 5
    M3.scale_and_add(1.0, r1, M1)     // line 6
    M3.scale_and_add(1.0, r2, M2)     // line 7
    M2.scale_and_add(1.0, r3, M1)     // line 8
    M1.scale_and_add(r4, f2, M2);     // line 9 (1/2)
    M1.add_scaled_identity(f0)        // line 9 (2/2)
    M2.scale_and_add(1.0,1.0,M3)      // line 10
    M3.add_scaled_identity(d0)        // line 11
    M1.multiply(M2,M3,f4,1.0)         // line 12
  }   


  template<typename T_matrix>
    void poly_8_eval_low_mem(std::vector<typename T_matrix::value_type> const & mc,
		     T_matrix & X,
		     T_matrix & X2) {
    T_matrix T1;
    poly_8_eval_low_mem(mc, X, X2, T1);
  }
} // end namespace sp8
#endif

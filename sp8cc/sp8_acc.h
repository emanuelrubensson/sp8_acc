#ifndef SP8_ACC_HEADER
#define SP8_ACC_HEADER
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include "Complex_step_solver.h"

namespace sp8 {
  /** Evaluate 8-order polynomial 
   *
   *  The polynomial is parameterized by 7 extremal points
   *  (r1,r2,...,r7) and two scale and shift factors (c1,c2)
   *  such that 
   *  sp8'(x) = c1*(x-r1)*(x-r2)*(x-r3)*(x-r4)*(x-r5)*(x-r6)*(x-r7)
   *  and sp8(0) = c2
   * 
   *  Input: parameters in vector v = [c1, c2, r1, r2, r3, r4, r5, r6, r7]
   *         point of evaluation  x
   */
  template<typename T>
  T sp8_eval(std::vector<T> const & v, T const & x) {
    const T c1 = v[0], c2 = v[1], r1 = v[2];
    const T r2 = v[3], r3 = v[4], r4 = v[5];
    const T r5 = v[6], r6 = v[7], r7 = v[8];
    T x0 = r1*r2;
    T x1 = r7*x0;
    T x2 = r3*r4;
    T x3 = r5*r6;
    T x4 = x2*x3;
    T x5 = x0*x2;
    T x6 = r5*r7;
    T x7 = r6*r7;
    T x8 = r3*r5;
    T x9 = x0*x8;
    T x10 = r4*r5;
    T x11 = x0*x10;
    T x12 = r1*r3;
    T x13 = x10*x12;
    T x14 = r2*r3;
    T x15 = x10*x14;
    T x16 = r1*r4;
    T x17 = r1*r5;
    T x18 = r1*r6;
    T x19 = r2*r4;
    T x20 = r2*r5;
    T x21 = r2*r6;
    T x22 = r3*r6;
    T x23 = r3*r7;
    T x24 = r4*r6;
    T x25 = r4*r7;
    T x26 = r5*x0;
    T x27 = r6*x0;
    T x28 = r6*x12;
    T x29 = r7*x12;
    T x30 = r7*x16;
    T x31 = r6*x14;
    T x32 = r7*x14;
    T x33 = r7*x19;
    T x34 = r7*x2;
    T res = (1.0/8.0)*c1*std::pow(x, 8) - 1.0/7.0*c1*std::pow(x, 7)*(r1 + r2 + r3 + r4 + r5 + r6 + r7) + (1.0/6.0)*c1*std::pow(x, 6)*(r1*r7 + r2*r7 + x0 + x10 + x12 + x14 + x16 + x17 + x18 + x19 + x2 + x20 + x21 + x22 + x23 + x24 + x25 + x3 + x6 + x7 + x8) - 1.0/5.0*c1*std::pow(x, 5)*(r3*x0 + r4*x0 + r4*x12 + r4*x14 + r5*x12 + r5*x14 + r5*x16 + r5*x19 + r5*x2 + r6*x10 + r6*x16 + r6*x17 + r6*x19 + r6*x2 + r6*x20 + r6*x8 + r7*x10 + r7*x17 + r7*x18 + r7*x20 + r7*x21 + r7*x22 + r7*x24 + r7*x3 + r7*x8 + x1 + x26 + x27 + x28 + x29 + x30 + x31 + x32 + x33 + x34) + (1.0/4.0)*c1*std::pow(x, 4)*(x0*x22 + x0*x23 + x0*x24 + x0*x25 + x0*x3 + x0*x6 + x0*x7 + x10*x7 + x11 + x12*x24 + x12*x25 + x12*x3 + x12*x6 + x12*x7 + x13 + x14*x24 + x14*x25 + x14*x3 + x14*x6 + x14*x7 + x15 + x16*x3 + x16*x6 + x16*x7 + x17*x7 + x19*x3 + x19*x6 + x19*x7 + x2*x6 + x2*x7 + x20*x7 + x4 + x5 + x7*x8 + x9) - 1.0/3.0*c1*std::pow(x, 3)*(x1*x10 + x1*x2 + x1*x22 + x1*x24 + x1*x3 + x1*x8 + x10*x27 + x10*x28 + x10*x29 + x10*x31 + x10*x32 + x2*x26 + x2*x27 + x24*x29 + x24*x32 + x27*x8 + x29*x3 + x3*x30 + x3*x32 + x3*x33 + x3*x34) + (1.0/2.0)*c1*std::pow(x, 2)*(x11*x7 + x13*x7 + x15*x7 + x3*x5 + x5*x6 + x5*x7 + x7*x9) - c1*x*x1*x4 + c2;
    return res;
  }

  template<typename T>
  T sp8_prim(std::vector<T> const & v, T const & x) {
    const T c1 = v[0], c2 = v[1], r1 = v[2];
    const T r2 = v[3], r3 = v[4], r4 = v[5];
    const T r5 = v[6], r6 = v[7], r7 = v[8];
    return c1*(x-r1)*(x-r2)*(x-r3)*(x-r4)*(x-r5)*(x-r6)*(x-r7);
  }

  /** Get polynomial corresponding to flipped ([0,1]) interval.
   *
   *  p_flipped(x) = 1-p(1-x)
   *
   * Parameters given by
   *   r_flipped(i) = 1-r(8-i), i = 1,...,7
   *   p_flipped'(1-x) = p'(x) -->
   *     c1_flipped*prod_i(1-r_flipped(i)-x) = p'(x)
   *   p_flipped(0) = 1-p(1) -->
   *     c2_flipped = 1-p(1)
   */
  template<typename T>
  void get_flipped_polynomial(std::vector<T> const & v, std::vector<T> & v_flipped) {
    const T c1 = v[0], c2 = v[1], r1 = v[2];
    const T r2 = v[3], r3 = v[4], r4 = v[5];
    const T r5 = v[6], r6 = v[7], r7 = v[8];
    v_flipped.resize(9);
    // Compute values r_flipped(i) = 1-r(8-i), i = 1,...,7
    for(unsigned int ind = 0; ind < 7; ind++)
      v_flipped[ind+2] = 1-v[8-ind];
    // Find x value in [0,1] as far away from r_i, i=1,...,7 as possible
    // Candidates are 0, (r1+r2)/2, (r2+r3)/2, ..., (r6+r7)/2, 1.0
    T max_dist = r1;
    T x = 0;
    for(unsigned int ind = 0; ind < 6; ind++) {
      T dist = (v[ind+3]-v[ind+2])/2;
      if (dist > max_dist) {
	max_dist = dist;
	x = v[ind+2]+dist;
      }
    }
    if (1-r7 > max_dist) {
      max_dist = 1-r7;
      x = 1.0;
    }
    // Compute product
    T r_product_flipped = 1.0;
    for(unsigned int ind = 0; ind < 7; ind++)
      r_product_flipped = r_product_flipped*(1-v_flipped[ind+2]-x);
    // p_flipped'(1-x) = p'(x)
    T c1_flipped = sp8_prim(v,x)/r_product_flipped;
    T c2_flipped = 1 - sp8_eval(v,(T)1.0);
    v_flipped[0] = c1_flipped;
    v_flipped[1] = c2_flipped;
  }


  template<typename T>
  void get_sp8_monomial_coefficients(std::vector<T> const & v,
				     std::vector<T> & mc) {
    const T c1 = v[0], c2 = v[1], r1 = v[2];
    const T r2 = v[3], r3 = v[4], r4 = v[5];
    const T r5 = v[6], r6 = v[7], r7 = v[8];
    mc.assign(9,0.0); // Fill with 9 zeros
    mc[0] = 1.0; // x8 coeff
    for(int i = 0; i < 7; i++) {
      mc[1] = mc[1] - v[i+2]; // x7 coeff
      for(int j = i+1; j < 7; j++) {
	mc[2] = mc[2] + v[i+2]*v[j+2]; // x6 coeff
	for(int k = j+1; k < 7; k++) {
	  mc[3] = mc[3] - v[i+2]*v[j+2]*v[k+2]; // x5 coeff
	  for(int l = k+1; l < 7; l++) {
	    mc[4] = mc[4] + v[i+2]*v[j+2]*v[k+2]*v[l+2]; // x4 coeff
	    for(int m = l+1; m < 7; m++) {
	      mc[5] = mc[5] - v[i+2]*v[j+2]*v[k+2]*v[l+2]*v[m+2]; // x3 coeff
	      for(int n = m+1; n < 7; n++)
		mc[6] = mc[6] + v[i+2]*v[j+2]*v[k+2]*v[l+2]*v[m+2]*v[n+2]; // x2 coeff
	    }
	  }
	}
      }
    }
    mc[7] = -r1*r2*r3*r4*r5*r6*r7; // x coeff
    for(int ind = 0; ind < 8; ind++)
      mc[ind] = mc[ind]*c1*(1.0/(8.0-(T)ind));
    mc[8] = c2;
  }

  template<typename value_type>
  struct Homotopy_solver_base {
    virtual int solve(value_type const L_target,
		      value_type const H_target,
		      std::vector<value_type> & v) = 0;
  };

  template<class Objective_function>
  struct Homotopy_solver: public Homotopy_solver_base<typename Objective_function::value_type> {
    typedef typename Objective_function::value_type value_type;
    int error_code_incorrect_order(bool ok_left, bool ok_right) {
      assert(!ok_left || !ok_right);
      if (!ok_left && !ok_right) return 3;
      if (!ok_left)  return -2;
      if (!ok_right) return  2;
      return 10; // To suppress compiler warning, this point should not be reached
    }
    /*
      Return codes:
      0: success
      1: failure in step_newton
      -2: incorrect order of roots to the left
      +2: incorrect order of roots to the right
      3: incorrect order of roots to the left and right
      4: failed to reach relative error < sqrt(eps) within maxiter
      5: reached maxiter in homotopy part
    */
    int solve(value_type L_target,
	      value_type H_target,
	      std::vector<value_type> & v) {
      assert(Objective_function::left != 0); // can be handled by flipping the interval
      if (Objective_function::right == 0)
	H_target = 1; // Not used in objective function in this case
      bool ok_left, ok_right = true;
      value_type v_maxabs;
      value_type dv_maxabs;
      v = Objective_function::v_start;
      value_type L = Objective_function::L_start;
      value_type H = Objective_function::H_start;
      Objective_function objfun(L, H);
      Complex_step_solver<Objective_function> solver(objfun);
      value_type gap_target = H_target - L_target;
      std::vector<value_type> tmp;
      int count = 0;
      while (L != L_target || H != H_target) {
	count++;
	// closest root to the left of L
	value_type root_near_L = v[Objective_function::left+1];
	// closest root to the right of H
	// or 1 if no roots to the right
	value_type root_near_H = Objective_function::right == 0 ?
	  1 : v[9-Objective_function::right];
	// max step length half distance to closest root
	value_type max_step_L = (L - root_near_L)/2;
	value_type max_step_H = (root_near_H - H)/2;
	value_type L_save = L;
	if (L_target <= L)
	  L = std::max(L_target, L - max_step_L);
	else
	  L = std::min({L_target,
	      std::max(L, H-0.9*gap_target), // keep distance to H
	      L + max_step_L});
	if (H_target >= H)
	  H = std::min(H_target, H + max_step_H);
	else
	  H = std::max({H_target,
	      std::min(H, L_save+0.9*gap_target), // keep distance to L
	      H - max_step_H});
	objfun.L = L;
	objfun.H = H;
	int info = solver.step_newton(v,tmp,v_maxabs,dv_maxabs);
	if (info != 0)
	  return 1;
	v.swap(tmp);
	if ( !objfun.correct_order_of_roots(v, ok_left, ok_right) )
	  return error_code_incorrect_order(ok_left, ok_right);
	if (count > 10000)
	  // Count usually between 0 and 300
	  return 5;
      }
      value_type dv_relative = dv_maxabs/v_maxabs;
      value_type dv_relative_prev;
      int max_iter = 10;
      for (int ind = 0;ind < max_iter;ind++) {
	dv_relative_prev = dv_relative;
	int info = solver.step_newton(v,tmp,v_maxabs,dv_maxabs);
	if (info != 0)
	  return 1;
	v.swap(tmp);
	if ( !objfun.correct_order_of_roots(v, ok_left, ok_right) )
	  return error_code_incorrect_order(ok_left, ok_right);
	dv_relative = dv_maxabs/v_maxabs;
	if (dv_relative < std::sqrt(std::numeric_limits<value_type>::epsilon())) {
	  // Ok, we are at least "half way" to machine epsilon
	  // Break if error is not decreasing substantially
	  if (dv_relative > dv_relative_prev/2)
	    break;
	}
      }
      if (dv_relative >= std::sqrt(std::numeric_limits<value_type>::epsilon()))
	return 4;
      // Add ones in the end (for 4-0, 5-0, and 6-0 cases).
      v.resize(9, 1.0);
      return 0;
    }
  };

  template<typename T2>
  struct Objective_Fun43 {
    typedef T2 value_type;
    static const int left  = 4;
    static const int right = 3;
    static std::vector<value_type> const v_start;
    static const value_type constexpr L_start = 0.55;
    static const value_type constexpr H_start = 0.65;
    value_type L;
    value_type H;
    Objective_Fun43(value_type const & L, value_type const & H)
      :L(L),H(H) {}
    template<typename T1>
    void fun(std::vector<T1> const & v, std::vector<T1> & result) const {
      const T1 c1 = v[0], c2 = v[1], r1 = v[2];
      const T1 r2 = v[3], r3 = v[4], r4 = v[5];
      const T1 r5 = v[6], r6 = v[7], r7 = v[8];
      auto p = [&](T1 x){return sp8_eval(v,x);};
      result.resize(9);
      result[0] = p(0.0);
      result[1] = p(r2);
      result[2] = p(r4);
      result[3] = p(r1) - p(r3);
      result[4] = p(r1) - p(L);
      result[5] = p(H)  - p(1);
      result[6] = p(r6) - p(1);
      result[7] = p(r5) - 1.0;
      result[8] = p(r7) - 1.0;
    }
    /** Check order of roots
     *
     *  For 4-3, the correct order is
     *  0<r1<r2<r3<r4<L<H<r5<r6<r7<1
     */
    bool correct_order_of_roots(std::vector<T2> const & v,
				bool & ok_left, bool & ok_right) const {
      const T2 c1 = v[0], c2 = v[1], r1 = v[2];
      const T2 r2 = v[3], r3 = v[4], r4 = v[5];
      const T2 r5 = v[6], r6 = v[7], r7 = v[8];
      ok_left  = (0<r1) && (r1<r2) && (r2<r3) && (r3<r4) && (r4<L);
      ok_right = (H<r5) && (r5<r6) && (r6<r7) && (r7<1.0);
      return ok_left && ok_right && L<H;
    }
  };
  template<typename T2>
  std::vector<T2> const Objective_Fun43<T2>::v_start =
    {-3.828330354372377e+04, // c1
     0.000000000000000e+00,  // c2
     3.630457073236520e-02,  // r1
     1.391234727632784e-01,  // ...
     2.902225497081609e-01,
     4.558182467167558e-01,
     7.341809614537667e-01,
     8.676371783744240e-01,
     9.650762995000088e-01}; // r7

  template<typename T2>
  struct Objective_Fun52 {
    typedef T2 value_type;
    static const int left  = 5;
    static const int right = 2;
    static std::vector<value_type> const v_start;
    static const value_type constexpr L_start = 0.75;
    static const value_type constexpr H_start = 0.85;
    value_type L;
    value_type H;
    Objective_Fun52(value_type const & L, value_type const & H)
      :L(L),H(H) {}
    template<typename T1>
    void fun(std::vector<T1> const & v, std::vector<T1> & result) const {
      const T1 c1 = v[0], c2 = v[1], r1 = v[2];
      const T1 r2 = v[3], r3 = v[4], r4 = v[5];
      const T1 r5 = v[6], r6 = v[7], r7 = v[8];
      auto p = [&](T1 x){return sp8_eval(v,x);};
      result.resize(9);
      result[0] = p(0.0) - p(L);
      result[1] = p(r2)  - p(L);
      result[2] = p(r4)  - p(L);
      result[3] = p(r1);
      result[4] = p(r3);
      result[5] = p(r5);
      result[6] = p(r6)  - 1.0;
      result[7] = p(1.0) - 1.0;
      result[8] = p(r7)  - p(H);
    }
    /** Check order of roots
     *
     *  For 5-2, the correct order is
     *  0<r1<r2<r3<r4<r5<L<H<r6<r7<1
     */
    bool correct_order_of_roots(std::vector<T2> const & v,
				bool & ok_left, bool & ok_right) const {
      const T2 c1 = v[0], c2 = v[1], r1 = v[2];
      const T2 r2 = v[3], r3 = v[4], r4 = v[5];
      const T2 r5 = v[6], r6 = v[7], r7 = v[8];
      ok_left  = (0<r1) && (r1<r2) && (r2<r3) && (r3<r4) && (r4<r5) && (r5<L);
      ok_right = (H<r6) && (r6<r7) && (r7<1.0);
      return ok_left && ok_right && L<H;
    }
  };
  template<typename T2>
  std::vector<T2> const Objective_Fun52<T2>::v_start =
    {3.616269094358620e+04,  // c1
     2.515360660658102e-01,  // c2
     3.734847296159297e-02,  // r1
     1.435436349152357e-01,  // ...
     3.017655595515090e-01,
     4.859232582684261e-01,
     6.598969344170539e-01,
     8.998903128223554e-01,
     9.714353568216227e-01}; // r7

  template<typename T2>
  struct Objective_Fun61 {
    typedef T2 value_type;
    static const int left  = 6;
    static const int right = 1;
    static std::vector<value_type> const v_start;
    static const value_type constexpr L_start = 0.88;
    static const value_type constexpr H_start = 0.98;
    value_type L;
    value_type H;
    Objective_Fun61(value_type const & L, value_type const & H)
      :L(L),H(H) {}
    template<typename T1>
    void fun(std::vector<T1> const & v, std::vector<T1> & result) const {
      const T1 c1 = v[0], c2 = v[1], r1 = v[2];
      const T1 r2 = v[3], r3 = v[4], r4 = v[5];
      const T1 r5 = v[6], r6 = v[7], r7 = v[8];
      auto p = [&](T1 x){return sp8_eval(v,x);};
      result.resize(9);
      result[0] = p(0.0);
      result[1] = p(r2);
      result[2] = p(r4);
      result[3] = p(r6);
      result[4] = p(r1)  - p(L);
      result[5] = p(r3)  - p(L);
      result[6] = p(r5)  - p(L);
      result[7] = p(1.0) - p(H);
      result[8] = p(r7)  - 1.0;
    }
    /** Check order of roots
     *
     *  For 6-1, the correct order is
     *  0<r1<r2<r3<r4<r5<r6<L<H<r7<1
     */
    bool correct_order_of_roots(std::vector<T2> const & v,
				bool & ok_left, bool & ok_right) const {
      const T2 c1 = v[0], c2 = v[1], r1 = v[2];
      const T2 r2 = v[3], r3 = v[4], r4 = v[5];
      const T2 r5 = v[6], r6 = v[7], r7 = v[8];
      ok_left  =
	(0<r1) && (r1<r2) && (r2<r3) && (r3<r4) && (r4<r5) && (r5<r6) && (r6<L);
      ok_right = (H<r7) && (r7<1.0);
      return ok_left && ok_right && L<H;
    }
  };
  template<typename T2>
  std::vector<T2> const Objective_Fun61<T2>::v_start =
    {-2.703277318672865e+04, // c1
     0.000000000000000e+00,  // c2
     3.704945885349221e-02,  // r1
     1.424978475605410e-01,  // r2
     3.000736741629149e-01,  // ...
     4.852354821229843e-01,
     6.682580170649393e-01,
     8.156966140526589e-01,
     9.904826516137375e-01}; // r7

  template<typename T2>
  struct Objective_Fun70 {
    typedef T2 value_type;
    static const int left  = 7;
    static const int right = 0;
    static std::vector<value_type> const v_start;
    static const value_type constexpr L_start = 0.94;
    static const value_type constexpr H_start = 1.0;
    value_type L;
    value_type H;
    Objective_Fun70(value_type const & L, value_type const & H)
      :L(L),H(H) {}
    template<typename T1>
    void fun(std::vector<T1> const & v, std::vector<T1> & result) const {
      const T1 c1 = v[0], c2 = v[1], r1 = v[2];
      const T1 r2 = v[3], r3 = v[4], r4 = v[5];
      const T1 r5 = v[6], r6 = v[7], r7 = v[8];
      auto p = [&](T1 x){return sp8_eval(v,x);};
      result.resize(9);
      result[0] = p(0.0) - p(L);
      result[1] = p(r2)  - p(L);
      result[2] = p(r4)  - p(L);
      result[3] = p(r6)  - p(L);
      result[4] = p(r1);
      result[5] = p(r3);
      result[6] = p(r5);
      result[7] = p(r7);
      result[8] = p(1.0) - 1.0;
    }
    /** Check order of roots
     *
     *  For 7-0, the correct order is
     *  0<r1<r2<r3<r4<r5<r6<r7<L<H<1
     */
    bool correct_order_of_roots(std::vector<T2> const & v,
				bool & ok_left, bool & ok_right) const {
      const T2 c1 = v[0], c2 = v[1], r1 = v[2];
      const T2 r2 = v[3], r3 = v[4], r4 = v[5];
      const T2 r5 = v[6], r6 = v[7], r7 = v[8];
      ok_left  = (0<r1) && (r1<r2) && (r2<r3) && (r3<r4) && (r4<r5)
	&& (r5<r6) && (r6<r7) && (r7<L);
      ok_right = true;
      return ok_left && ok_right && L<H;
    }
  };
  template<typename T2>
  std::vector<T2> const Objective_Fun70<T2>::v_start =
    {1.518407434056799e+04, // c1
     7.061569273738105e-02, // c2
     3.577661971967767e-02, // r1
     1.376598128422552e-01, // ...
     2.901387867882566e-01,
     4.699999999997491e-01,
     6.498612132111929e-01,
     8.023401871574832e-01,
     9.042233802799536e-01};// r7

  template<typename T2>
  struct Objective_Fun10 {
    typedef T2 value_type;
    static const int left  = 1;
    static const int right = 0;
    static std::vector<value_type> const v_start;
    static const value_type constexpr L_start = 0.3;
    static const value_type constexpr H_start = 1.0;
    value_type L;
    value_type H;
    Objective_Fun10(value_type const & L, value_type const & H)
      :L(L),H(H) {}
    template<typename T1>
    void fun(std::vector<T1> const & v, std::vector<T1> & result) const {
      const T1 c1 = v[0], c2 = v[1], r1 = v[2];
      std::vector<T1> v_extended = {c1,c2,r1,1.0,1.0,1.0,1.0,1.0,1.0};
      auto p = [&](T1 x){return sp8_eval(v_extended,x);};
      result.resize(3);
      result[0] = p(0.0) - p(L);
      result[1] = p(r1);
      result[2] = p(1.0) - 1.0;
    }
    /** Check order of roots
     *
     *  For 1-0, the correct order is
     *  0<r1<L<H<r2=r3=r4=r5=r6=r7=1
     */
    bool correct_order_of_roots(std::vector<T2> const & v,
				bool & ok_left, bool & ok_right) const {
      const T2 c1 = v[0], c2 = v[1], r1 = v[2];
      // Note that v has only 3 elements in the 1-0 case (r2=r3=r4=r5=r6=r7=1)
      ok_left  = (0<r1) && (r1<L);
      ok_right = true;
      return ok_left && ok_right && L<H;
    }
  };
  template<typename T2>
  std::vector<T2> const Objective_Fun10<T2>::v_start =
    {131.7706407785533429,  // c1
     0.5565331390208555,    // c2
     0.1014418841062514};   // r1

  template<typename T2>
  struct Objective_Fun20 {
    typedef T2 value_type;
    static const int left  = 2;
    static const int right = 0;
    static std::vector<value_type> const v_start;
    static const value_type constexpr L_start = 0.4;
    static const value_type constexpr H_start = 1.0;
    value_type L;
    value_type H;
    Objective_Fun20(value_type const & L, value_type const & H)
      :L(L),H(H) {}
    template<typename T1>
    void fun(std::vector<T1> const & v, std::vector<T1> & result) const {
      const T1 c1 = v[0], c2 = v[1], r1 = v[2];
      const T1 r2 = v[3];
      std::vector<T1> v_extended = {c1,c2,r1,r2,1.0,1.0,1.0,1.0,1.0};
      auto p = [&](T1 x){return sp8_eval(v_extended,x);};
      result.resize(4);
      result[0] = p(0.0);
      result[1] = p(r2);
      result[2] = p(r1) - p(L);
      result[3] = p(1)  - 1.0;
    }
    /** Check order of roots
     *
     *  For 2-0, the correct order is
     *  0<r1<r2<L<H<r3=r4=r5=r6=r7=1
     */
    bool correct_order_of_roots(std::vector<T2> const & v,
				bool & ok_left, bool & ok_right) const {
      const T2 c1 = v[0], c2 = v[1], r1 = v[2];
      const T2 r2 = v[3];
      // Note that v has only 4 elements in the 2-0 case (r3=r4=r5=r6=r7=1)
      ok_left  = (0<r1) && (r1<r2) && (r2<L);
      ok_right = true;
      return ok_left && ok_right && L<H;
    }
  };
  template<typename T2>
  std::vector<T2> const Objective_Fun20<T2>::v_start =
    {-796.3368359148582840, // c1
     0.0,                   // c2
     0.0645147150177223,    // r1
     0.2420578253664362};   // r2

  template<typename T2>
  struct Objective_Fun30 {
    typedef T2 value_type;
    static const int left  = 3;
    static const int right = 0;
    static std::vector<value_type> const v_start;
    static const value_type constexpr L_start = 0.5;
    static const value_type constexpr H_start = 1.0;
    value_type L;
    value_type H;
    Objective_Fun30(value_type const & L, value_type const & H)
      :L(L),H(H) {}
    template<typename T1>
    void fun(std::vector<T1> const & v, std::vector<T1> & result) const {
      const T1 c1 = v[0], c2 = v[1], r1 = v[2];
      const T1 r2 = v[3], r3 = v[4]; 
      std::vector<T1> v_extended = {c1,c2,r1,r2,r3,1.0,1.0,1.0,1.0};
      auto p = [&](T1 x){return sp8_eval(v_extended,x);};
      result.resize(5);
      result[0] = p(0.0) - p(L);
      result[1] = p(r2)  - p(L);
      result[2] = p(r1);
      result[3] = p(r3);
      result[4] = p(1.0) - 1.0;
    }
    /** Check order of roots
     *
     *  For 3-0, the correct order is
     *  0<r1<r2<r3<L<H<r4=r5=r6=r7=1
     */
    bool correct_order_of_roots(std::vector<T2> const & v,
				bool & ok_left, bool & ok_right) const {
      const T2 c1 = v[0], c2 = v[1], r1 = v[2];
      const T2 r2 = v[3], r3 = v[4];
      // Note that v has only 5 elements in the 3-0 case (r4=r5=r6=r7=1)
      ok_left  = (0<r1) && (r1<r2) && (r2<r3) && (r3<L);
      ok_right = true;
      return ok_left && ok_right && L<H;
    }
  };
  template<typename T2>
  std::vector<T2> const Objective_Fun30<T2>::v_start =
    {2632.4596551462145726,  // c1
     0.1848064209714597,     // c2
     0.0491352255195878,     // r1
     0.1869223450229846,     // r2
     0.3779197335216983};    // r3

  template<typename T2>
  struct Objective_Fun40 {
    typedef T2 value_type;
    static const int left  = 4;
    static const int right = 0;
    static std::vector<value_type> const v_start;
    static const value_type constexpr L_start = 0.5;
    static const value_type constexpr H_start = 1.0;
    value_type L;
    value_type H;
    Objective_Fun40(value_type const & L, value_type const & H)
      :L(L),H(H) {}
    template<typename T1>
    void fun(std::vector<T1> const & v, std::vector<T1> & result) const {
      const T1 c1 = v[0], c2 = v[1], r1 = v[2];
      const T1 r2 = v[3], r3 = v[4], r4 = v[5];
      std::vector<T1> v_extended = {c1,c2,r1,r2,r3,r4,1.0,1.0,1.0};
      auto p = [&](T1 x){return sp8_eval(v_extended,x);};
      result.resize(6);
      result[0] = p(0.0);
      result[1] = p(r2);
      result[2] = p(r4);
      result[3] = p(r1) - p(r3);
      result[4] = p(r1) - p(L);
      result[5] = p(1)  - 1.0;
    }
    /** Check order of roots
     *
     *  For 4-0, the correct order is
     *  0<r1<r2<r3<r4<L<H<r5=r6=r7=1
     */
    bool correct_order_of_roots(std::vector<T2> const & v,
				bool & ok_left, bool & ok_right) const {
      const T2 c1 = v[0], c2 = v[1], r1 = v[2];
      const T2 r2 = v[3], r3 = v[4], r4 = v[5];
      // Note that v has only 6 elements in the 4-0 case (r5=r6=r7=1)
      ok_left  = (0<r1) && (r1<r2) && (r2<r3) && (r3<r4) && (r4<L);
      ok_right = true;
      return ok_left && ok_right && L<H;
    }
  };
  template<typename T2>
  std::vector<T2> const Objective_Fun40<T2>::v_start =
    {-3082.84741694587,     // c1
     -2.52435489670693e-29, // c2
     0.0366703565717832,	  // r1
     0.139392693385197,	  // r2
     0.285519637000744,	  // r3
     0.431541718562603};	  // r4

  template<typename T2>
  struct Objective_Fun50 {
    typedef T2 value_type;
    static const int left  = 5;
    static const int right = 0;
    static std::vector<value_type> const v_start;
    static const value_type constexpr L_start = 0.5;
    static const value_type constexpr H_start = 1.0;
    value_type L;
    value_type H;
    Objective_Fun50(value_type const & L, value_type const & H)
      :L(L),H(H) {}
    template<typename T1>
    void fun(std::vector<T1> const & v, std::vector<T1> & result) const {
      const T1 c1 = v[0], c2 = v[1], r1 = v[2];
      const T1 r2 = v[3], r3 = v[4], r4 = v[5];
      const T1 r5 = v[6];
      std::vector<T1> v_extended = {c1,c2,r1,r2,r3,r4,r5,1.0,1.0};
      auto p = [&](T1 x){return sp8_eval(v_extended,x);};
      result.resize(7);
      result[0] = p(0.0) - p(L);
      result[1] = p(r2)  - p(L);
      result[2] = p(r4)  - p(L);
      result[3] = p(r1);
      result[4] = p(r3);
      result[5] = p(r5);
      result[6] = p(1.0) - 1.0;
    }
    /** Check order of roots
     *
     *  For 5-0, the correct order is
     *  0<r1<r2<r3<r4<r5<L<H<r6=r7=1
     */
    bool correct_order_of_roots(std::vector<T2> const & v,
				bool & ok_left, bool & ok_right) const {
      const T2 c1 = v[0], c2 = v[1], r1 = v[2];
      const T2 r2 = v[3], r3 = v[4], r4 = v[5];
      const T2 r5 = v[6];
      // Note that v has only 7 elements in the 5-0 case (r6=r7=1)
      ok_left  = (0<r1) && (r1<r2) && (r2<r3) && (r3<r4) && (r4<r5) && (r5<L);
      ok_right = true;
      return ok_left && ok_right && L<H;
    }
  };
  template<typename T2>
  std::vector<T2> const Objective_Fun50<T2>::v_start =
    {1974.00119732597,    // c1
     0.00267129139635331, // c2
     0.028633198688273,   // r1
     0.109223055107629,   // ...
     0.226156115705252,
     0.354106988313488,
     0.458452413536081};  // r5

  template<typename T2>
  struct Objective_Fun60 {
    typedef T2 value_type;
    static const int left  = 6;
    static const int right = 0;
    static std::vector<value_type> const v_start;
    static const value_type constexpr L_start = 0.88;
    static const value_type constexpr H_start = 1.0;
    value_type L;
    value_type H;
    Objective_Fun60(value_type const & L, value_type const & H)
      :L(L),H(H) {}
    template<typename T1>
    void fun(std::vector<T1> const & v, std::vector<T1> & result) const {
      const T1 c1 = v[0], c2 = v[1], r1 = v[2];
      const T1 r2 = v[3], r3 = v[4], r4 = v[5];
      const T1 r5 = v[6], r6 = v[7];
      std::vector<T1> v_extended = {c1,c2,r1,r2,r3,r4,r5,r6,1.0};
      auto p = [&](T1 x){return sp8_eval(v_extended,x);};
      result.resize(8);
      result[0] = p(0.0);
      result[1] = p(r2);
      result[2] = p(r4);
      result[3] = p(r6);
      result[4] = p(r1)  - p(L);
      result[5] = p(r3)  - p(L);
      result[6] = p(r5)  - p(L);
      result[7] = p(1.0) - 1.0;
    }
    /** Check order of roots
     *
     *  For 6-0, the correct order is
     *  0<r1<r2<r3<r4<r5<r6<L<H<r7=1
     */
    bool correct_order_of_roots(std::vector<T2> const & v,
				bool & ok_left, bool & ok_right) const {
      const T2 c1 = v[0], c2 = v[1], r1 = v[2];
      const T2 r2 = v[3], r3 = v[4], r4 = v[5];
      const T2 r5 = v[6], r6 = v[7];
      ok_left  =
	(0<r1) && (r1<r2) && (r2<r3) && (r3<r4) && (r4<r5) && (r5<r6) && (r6<L);
      ok_right = true;
      return ok_left && ok_right && L<H;
    }
  };
  template<typename T2>
  std::vector<T2> const Objective_Fun60<T2>::v_start =
    {-22303.4615521158,     // c1
     -1.00974195868148e-28, // c2
     0.03716257529248,      // r1
     0.14292553399174,      // r2
     0.300942293257609,	  // ..
     0.486538327322085,
     0.669757555005894,
     0.816704819712282};	  // r6


  template<typename value_type>
  Homotopy_solver_base<value_type>* homotopy_solver_factory(unsigned int roots_left, unsigned int roots_right) {
    assert(roots_left  <= 7);
    assert(roots_right <= 7);
    Homotopy_solver_base<value_type>* solver;
    auto switch_pair = [](unsigned int x, unsigned int y){return (x<<3)+y;};
    switch(switch_pair(roots_left, roots_right)) {
    case switch_pair(1,0):
      solver = new Homotopy_solver<Objective_Fun10<value_type> >;
      break;
    case switch_pair(2,0):
      solver = new Homotopy_solver<Objective_Fun20<value_type> >;
      break;
    case switch_pair(3,0):
      solver = new Homotopy_solver<Objective_Fun30<value_type> >;
      break;
    case switch_pair(4,0):
      solver = new Homotopy_solver<Objective_Fun40<value_type> >;
      break;
    case switch_pair(5,0):
      solver = new Homotopy_solver<Objective_Fun50<value_type> >;
      break;
    case switch_pair(6,0):
      solver = new Homotopy_solver<Objective_Fun60<value_type> >;
      break;
    case switch_pair(7,0):
      solver = new Homotopy_solver<Objective_Fun70<value_type> >;
      break;
    case switch_pair(4,3):
      solver = new Homotopy_solver<Objective_Fun43<value_type> >;
      break;
    case switch_pair(5,2):
      solver = new Homotopy_solver<Objective_Fun52<value_type> >;
      break;
    case switch_pair(6,1):
      solver = new Homotopy_solver<Objective_Fun61<value_type> >;
      break;
    default:
      std::cout << "homotopy_solver_factory not implemented for <left-right> = <" << roots_left << "-" << roots_right << ">" << std::endl;
      std::exit(1);
    }
    return solver;
  }

  struct SP8_spec {
    int left;       /**< number of stationary points to the left of L  */
    int right;      /**< number of stationary points to the right of H */
    bool acc_left;  /**< use acceleration to the left?                 */
    bool acc_right; /**< use acceleration to the right?                */
    SP8_spec(int left, int right, bool acc_left, bool acc_right)
      : left(left), right(right), acc_left(acc_left), acc_right(acc_right) {
      // Acceleration to the left not possible with left == 0.
      assert( !(acc_left  && left  == 0) );
      // Likewise to the right.
      assert( !(acc_right && right == 0) );
    }
    SP8_spec reversed() const {
      return SP8_spec(right, left, acc_right, acc_left);
    }
  };

  template<typename T>
  void get_sp8_params_no_acc(const int left, const int right,
                             std::vector<T> & v,
                             T & Csp8, T & qsp8) {
    auto switch_pair = [](unsigned int x, unsigned int y){return (x<<3)+y;};
    switch(switch_pair(left, right)) {
    case switch_pair(4,3):
      v.assign({-280.0, 0, 0, 0, 0, 0, 1, 1, 1});
      Csp8 = 82; qsp8 = 4;
      break;
    case switch_pair(3,4):
      v.assign({280.0, 0, 0, 0, 0, 1, 1, 1, 1});
      Csp8 = 82; qsp8 = 4;
      break;
    case switch_pair(5,2):
      v.assign({168.0, 0, 0, 0, 0, 0, 0, 1, 1});
      Csp8 = 56; qsp8 = 3;
      break;
    case switch_pair(2,5):
      v.assign({-168.0, 0, 0, 0, 1, 1, 1, 1, 1});
      Csp8 = 56; qsp8 = 3;
      break;
    case switch_pair(6,1):
      Csp8 = 28; qsp8 = 2;
      v.assign({-56.0, 0, 0, 0, 0, 0, 0, 0, 1});
      break;
    case switch_pair(1,6):
      v.assign({56.0, 0, 0, 1, 1, 1, 1, 1, 1});
      Csp8 = 28; qsp8 = 2;
      break;
    case switch_pair(7,0):
      v.assign({8.0, 0, 0, 0, 0, 0, 0, 0, 0});
      Csp8 = std::numeric_limits<double>::infinity(); qsp8 = 1;
      break;
    case switch_pair(0,7):
      v.assign({-8.0, 0, 1, 1, 1, 1, 1, 1, 1});
      Csp8 = std::numeric_limits<double>::infinity(); qsp8 = 1;
      break;
    default:
      std::cout << "sp8 without acceleration does not exist with <left-right> = <" << left << "-" << right << ">" << std::endl;
      std::exit(1);
    }
  }
  template<typename T>
  void get_sp8_params_no_acc(const int left, const int right,
			     std::vector<T> & v) {
    T Csp8 = 0; T qsp8 = 1; // dummy variables
    get_sp8_params_no_acc(left, right, v, Csp8, qsp8);
  }

  /** Get coefficients of eigth degree polynomial. User decides location
      of extremal points and whether to use acceleration to the left
      and/or right via the sp8_spec parameter.
      @param L        lumo eigenvalue within [0,1] interval
      @param H        homo eigenvalue within [0,1] interval (L<H)
      @param sp8_spec see definition of SP8_spec
      @param v        output vector with 7 extremal points (r1,..,r7) and
      two scale and shift factors, see sp8.
  */
  template<typename T>
  int get_sp8_params(const T L, const T H,
		     const SP8_spec sp8_spec,
		     std::vector<T> & v,
		     T & Csp8, T & qsp8) {
    assert(0 <= L);
    assert(L <  H);
    assert(H <= 1);
    // By default params for stopping criteria set to not stop
    Csp8 = std::numeric_limits<double>::infinity();
    qsp8 = 1;
    if (!sp8_spec.acc_left && !sp8_spec.acc_right) {
      get_sp8_params_no_acc(sp8_spec.left, sp8_spec.right, v, Csp8, qsp8);
      // Params for stopping criteria updated in this case
      return 0;
    }
    int left  = sp8_spec.left;
    int right = sp8_spec.right;
    if (!sp8_spec.acc_left)
      left = 0;
    if (!sp8_spec.acc_right)
      right = 0;
    if (left < right) {
      std::vector<T> v_tmp;
      int info = get_sp8_params(1-H, 1-L, sp8_spec.reversed(),  v_tmp);
      if (info != 0)
	return -info; // reversed
      get_flipped_polynomial(v_tmp, v); // v is flipped
      return 0;
    }
    Homotopy_solver_base<double>* solver = homotopy_solver_factory<double>(left, right);
    int info = solver->solve(L, H, v);
    delete solver;
    return info;
  }
  template<typename T>
  int get_sp8_params(const T L, const T H,
		     const SP8_spec sp8_spec,
		     std::vector<T> & v) {
    T Csp8 = 0; T qsp8 = 1; // dummy variables
    return get_sp8_params(L, H, sp8_spec, v, Csp8, qsp8);
  }

  /** Choose number of extremal points to the left and right based on
   *  which one gives the largest slope at mu = (H+L)/2.  The break
   *  points (0.692, 0.855, 0.963) come from numerical experiments.
   *  Experiments: small gap and acceleration based on exact H and L
   *  values.
   */
  template<typename T>
  void get_no_of_extremal_points_left_right(const T mu,
					    int & left, int & right) {
    if (mu < 0.5) {
      get_no_of_extremal_points_left_right(1-mu, right, left);
      return;
    }
    // Ok, only need to handle mu in [0.5, 1]
    if (mu < 0.692) { // Use <4-3>
      left  = 4;
      right = 3;
      return;
    }
    if (mu < 0.855) { // Use <5-2>
      left  = 5;
      right = 2;
      return;
    }
    if (mu < 0.963) { // Use <6-1>
      left  = 6;
      right = 1;
      return;
    }
    // Use <7-0>
    left  = 7;
    right = 0;
  }

  /** Get coefficients of eigth degree polynomial. This version chooses
      location of extremal points based on location of mu = (H+L)/2. The
      polynomial with the largest slope at mu is chosen.
      @param L      lumo eigenvalue within [0,1] interval
      @param H      homo eigenvalue within [0,1] interval (L<H)
      @param v      output vector with 7 extremal points (r1,..,r7) and
      two scale and shift factors, see sp8.
  */
  template<typename T>
  void get_sp8_params_max_slope(const T L, const T H, std::vector<T> & v) {
    assert(0 <= L);
    assert(L <  H);
    assert(H <= 1);
    const T limit_value = 0.01;
    int left, right = 0;
    get_no_of_extremal_points_left_right((L+H)/2, left, right);
    // Use acceleration when L is away from 0 and H away from 1
    // but not when there is no stationary points to the left and
    // right respectively
    bool acc_left  = L > limit_value && left > 0;
    bool acc_right = H < 1-limit_value && right > 0;
    SP8_spec sp8_spec = {left, right, acc_left, acc_right};
    int info = get_sp8_params(L, H, sp8_spec, v);
    if (info != 0) {
      // get_sp8_params is not expected to fail in this case...
      std::cerr << "get_sp8_params failed with info = " << info
		<< " with params L = " << L
		<< ", H = " << H
		<< "<" << sp8_spec.left
		<< "," << sp8_spec.right
		<< ">" << std::endl;
      std::exit(1);
    }
  }

  template<typename T>
  void get_sp8_params_max_gap(const T L_outer, const T L_inner,
			      const T H_inner, const T H_outer,
			      const T kappa,
			      std::vector<T> & v,
			      T & Csp8,
			      T & qsp8) {
    T Csp8_tmp = std::numeric_limits<double>::infinity();
    T qsp8_tmp = 1.0;
    int left_start = 0;
    int left_end   = 7;
    if ((L_inner < kappa) && (H_inner > 1-kappa)) {
      // Only <3,4> and <4,3> are candidate polynomials close to
      // convergence.
      left_start = 3;
      left_end   = 4;
    }
    T gap_old = H_inner-L_inner;
    assert(gap_old > 0);
    T gap = 0;
    std::vector<T> v_tmp;
    for(int left = left_start; left <= left_end; left++) {
      int right = 7 - left;
      // Use acceleration to the left when L_outer is away from 0
      // but not when there are no stationary points to the left
      bool acc_left  = L_outer > kappa && left > 0;
      // and similarly to the right
      bool acc_right = H_outer < 1-kappa && right > 0;
      SP8_spec sp8_spec = {left, right, acc_left, acc_right};
      // Outer bounds used for acceleration, a precautious approach and
      // necessary if we want homo and lumo to be isolated for eigenvalue
      // estimation
      int info = get_sp8_params(L_outer, H_outer, sp8_spec, v_tmp, Csp8_tmp, qsp8_tmp);
      if (info != 0) {
	// Attempt to turn off acceleration
	if (info !=  2)  // if problem not to the right (only)
	  sp8_spec.acc_left  = false;
	if (info != -2) // if problem not to the left (only)
	  sp8_spec.acc_right = false;
	std::cout << "Rerun without acceleration triggered! info = " << info << std::endl;
	info = get_sp8_params(L_outer, H_outer, sp8_spec, v_tmp, Csp8_tmp, qsp8_tmp);
      }
      if (info != 0) {
	std::cerr << "Warning: get_sp8_params_max_gap: get_sp8_params failed with info = " << info
		  << " with params L_outer = " << L_outer
		  << ", L_inner = " << L_inner
		  << ", H_inner = " << H_inner
		  << ", H_outer = " << H_outer
		  << ", <" << sp8_spec.left
		  << "," << sp8_spec.right
		  << "> "
		  << acc_left << ", "
		  << acc_right
		  << std::endl;
	continue;
      }
      // Ok now we can compute the updated gap_new = H_inner_new-L_inner_new
      T gap_tmp = sp8_eval(v_tmp,H_inner) - sp8_eval(v_tmp,L_inner);
      if (gap_tmp > gap) {
	gap = gap_tmp;
	v = v_tmp;
	Csp8 = Csp8_tmp;
	qsp8 = qsp8_tmp;
      }
    }
    if (gap == 0) {
      std::cerr << "get_sp8_params_max_gap: get_sp8_params failed for all possible configurations, exiting." << std::endl;
      std::exit(1);
    }
    if (gap_old < 1-std::sqrt(std::numeric_limits<T>::epsilon()) && gap <= gap_old) {
      std::cerr << "get_sp8_params_max_gap: "
		<< " with params L_outer = " << L_outer
		<< ", L_inner = " << L_inner
		<< ", H_inner = " << H_inner
		<< ", H_outer = " << H_outer
		<< ", gap_old = " << gap_old
		<< ", gap = " << gap
		<< std::endl;
      std::cerr << "get_sp8_params_max_gap: get_sp8_params failed to produce polynomial giving increasing gap, exiting." << std::endl;
      std::exit(1);
    }
  }

} // end namespace
#endif

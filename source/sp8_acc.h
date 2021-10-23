#ifndef SP8_ACC_HEADER
#define SP8_ACC_HEADER
#include <vector>
#include <cmath>
#include <complex>
#include "Complex_step_solver.h"
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
T sp8(std::vector<T> const & v, T const & x) {
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

template<class Objective_function>
struct Homotopy_solver {
  typedef typename Objective_function::value_type value_type;
  void operator()(value_type const L_target,
		  value_type const H_target,
		  std::vector<value_type> & v,
		  value_type & v_maxabs,
		  value_type & dv_maxabs) {
    v = Objective_function::v_start;
    value_type L = Objective_function::L_start;
    value_type H = Objective_function::H_start;
    int nsteps = 100;
    value_type L_step = (L_target-L)/nsteps;
    value_type H_step = (H_target-H)/nsteps;
    std::vector<value_type> tmp;
    for (int ind = 0; ind<nsteps; ind++) {
      L = L+L_step;
      H = H+H_step;
      Objective_function objfun(L,H);
      Complex_step_solver<Objective_function> solver(objfun);
      solver.step_newton(v,tmp,v_maxabs,dv_maxabs);
      v.swap(tmp);
    }
    Objective_function objfun(L,H);
    Complex_step_solver<Objective_function> solver(objfun);
    for (int ind = 0; ind<4; ind++) {
      solver.step_newton(v,tmp,v_maxabs,dv_maxabs);
      v.swap(tmp);
    }
  }
};

template<typename T2>
struct Objective_Fun43 {
  typedef T2 value_type;
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
    auto p = [&](T1 x){return sp8(v,x);};
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
#endif

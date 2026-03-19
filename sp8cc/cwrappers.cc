#include <array>
#include "sp8_acc.h"
#include "Matrix_proxy.h"
#include "sastre_poly8_eval.h"
#include "poly8_eval.h"
#include "poly8_eval_low_mem.h"
#include "poly5_eval.h"

template<typename T_scalar>
void matrix_multiply(T_scalar const * ap, T_scalar const * bp, T_scalar * cp, int n) {
  sp8::Matrix_proxy<T_scalar const> const A(ap,n);
  sp8::Matrix_proxy<T_scalar const> const B(bp,n);
  sp8::Matrix_proxy<T_scalar> C(cp,n);
  C.multiply(A,B);
}

template<typename T, std::size_t N>
bool none_is_equal(std::array<T*,N> const & ptrs) {
  for(unsigned int ind=0;ind<N;ind++)
    for(unsigned int ind2=ind+1;ind2<N;ind2++)
      if ( ptrs[ind] == ptrs[ind2] )
	return false;
  return true;
}

template<typename T_scalar>
void sastre_poly_8_eval(T_scalar const * mc, T_scalar * ap, T_scalar * a2p,
			T_scalar * m2p, T_scalar * m3p, int n) {
  sp8::Matrix_proxy<T_scalar> A(ap,n);
  sp8::Matrix_proxy<T_scalar> A2(a2p,n);
  sp8::Matrix_proxy<T_scalar> M2(m2p,n);
  sp8::Matrix_proxy<T_scalar> M3(m3p,n);
  std::vector<T_scalar> mc_vec;
  mc_vec.resize(9);
  std::copy_n(mc, 9, mc_vec.begin());
  sp8::sastre_poly_8_eval(mc_vec, A, A2, M2, M3);
}

template<typename T_scalar>
void poly_8_eval(T_scalar const * mc, T_scalar * ap, T_scalar * a2p,
		 T_scalar * m2p, T_scalar * m3p, int n) {
  assert( none_is_equal(std::array<T_scalar*,4>({ap,a2p,m2p,m3p})) );
  sp8::Matrix_proxy<T_scalar> A(ap,n);
  sp8::Matrix_proxy<T_scalar> A2(a2p,n);
  sp8::Matrix_proxy<T_scalar> M2(m2p,n);
  sp8::Matrix_proxy<T_scalar> M3(m3p,n);
  std::vector<T_scalar> mc_vec;
  mc_vec.resize(9);
  std::copy_n(mc, 9, mc_vec.begin());
  sp8::poly_8_eval(mc_vec, A, A2, M2, M3);
}

template<typename T_scalar>
void poly_8_eval_low_mem(T_scalar const * mc, T_scalar * ap, T_scalar * a2p,
		 T_scalar * t1p, int n) {
  assert( none_is_equal(std::array<T_scalar*,3>({ap,a2p,t1p})) );
  sp8::Matrix_proxy<T_scalar> A(ap,n);
  sp8::Matrix_proxy<T_scalar> A2(a2p,n);
  sp8::Matrix_proxy<T_scalar> T1(t1p,n);
  std::vector<T_scalar> mc_vec;
  mc_vec.resize(9);
  std::copy_n(mc, 9, mc_vec.begin());
  sp8::poly_8_eval_low_mem(mc_vec, A, A2, T1);
}

template<typename T_scalar>
void poly_5_eval(T_scalar const * mc, T_scalar * ap, T_scalar * a2p,
		 T_scalar * m2p, int n) {
  assert( none_is_equal(std::array<T_scalar*,3>({ap,a2p,m2p})) );
  sp8::Matrix_proxy<T_scalar> A(ap,n);
  sp8::Matrix_proxy<T_scalar> A2(a2p,n);
  sp8::Matrix_proxy<T_scalar> M2(m2p,n);
  std::vector<T_scalar> mc_vec;
  mc_vec.resize(6);
  std::copy_n(mc, 6, mc_vec.begin());
  sp8::poly_5_eval(mc_vec, A, A2, M2);
}

template<typename T_scalar>
T_scalar trace_XmX2(T_scalar const * ap, int n) {
  sp8::Matrix_proxy<T_scalar const> const A(ap,n);
  return A.trace_XmX2();
}

extern "C" {
  /** v_output should point to array with 9 elements */
  void get_sp8_params_no_acc(const int left, const int right,
			     double* v_output,
			     double* C_q_sp8) {
    std::vector<double> v;
    sp8::get_sp8_params_no_acc(left, right, v, C_q_sp8[0], C_q_sp8[1]);
    assert(v.size() == 9);
    std::copy_n(v.begin(), 9, v_output);
  }

  /** v_output should point to array with 9 elements */
  int get_sp8_params(const double L, const double H,
                     const int left, const int right,
                     const bool acc_left, const bool acc_right,
                     double* v_output) {
    sp8::SP8_spec sp8_spec = {left, right, acc_left, acc_right};
    std::vector<double> v;
    int info = sp8::get_sp8_params(L, H, sp8_spec, v);
    if (info != 0)
      return info;
    assert(v.size() == 9);
    std::copy_n(v.begin(), 9, v_output);
    return 0;
  }

  /** v_output should point to array with 9 elements */
  void get_sp8_params_max_slope(const double L, const double H,
				double* v_output) {
    std::vector<double> v;
    sp8::get_sp8_params_max_slope(L, H, v);
    assert(v.size() == 9);
    std::copy_n(v.begin(), 9, v_output);
  }

  /** v_output should point to array with 9 elements */
  void get_sp8_params_max_gap(const double L_outer, const double L_inner,
                              const double H_inner, const double H_outer,
			      const double kappa,
                              double* v_output,
			      double* C_q_sp8) {
    std::vector<double> v;
    sp8::get_sp8_params_max_gap(L_outer, L_inner, H_inner, H_outer, kappa, v, C_q_sp8[0], C_q_sp8[1]);
    assert(v.size() == 9);
    std::copy_n(v.begin(), 9, v_output);
  }

  /** v_input should point to array with 9 elements.
   *  mc      should point to array with 9 elements.
   */
  void get_sp8_monomial_coefficients(double const * v_input, double * mc) {
    std::vector<double> v;
    std::vector<double> mc_vec;
    v.resize(9);
    std::copy_n(v_input, 9, v.begin());
    sp8::get_sp8_monomial_coefficients(v, mc_vec);
    std::copy_n(mc_vec.begin(), 9, mc);
  }

  double sp8_eval(double const * v_input, double const x) {
    std::vector<double> v;
    v.resize(9);
    std::copy_n(v_input, 9, v.begin());
    return sp8::sp8_eval(v, x);
  }

  double sp8_prim(double const * v_input, double const x) {
    std::vector<double> v;
    v.resize(9);
    std::copy_n(v_input, 9, v.begin());
    return sp8::sp8_prim(v, x);
  }

  void matmul_single(float const * ap, float const * bp, float * cp, int n) {
    matrix_multiply(ap, bp, cp, n);
  }
  void matmul_double(double const * ap, double const * bp, double * cp, int n) {
    matrix_multiply(ap, bp, cp, n);
  }

  void sastre_poly_8_eval_single(float const * mc, float * ap, float * a2p,
				 float * m2p, float * m3p, int n) {
    sastre_poly_8_eval(mc, ap, a2p, m2p, m3p, n);
  }
  void sastre_poly_8_eval_double(double const * mc, double * ap, double * a2p,
				 double * m2p, double * m3p, int n) {
    sastre_poly_8_eval(mc, ap, a2p, m2p, m3p, n);
  }

  void poly_8_eval_single(float const * mc, float * ap, float * a2p,
			  float * m2p, float * m3p, int n) {
    poly_8_eval(mc, ap, a2p, m2p, m3p, n);
  }
  void poly_8_eval_double(double const * mc, double * ap, double * a2p,
			  double * m2p, double * m3p, int n) {
    poly_8_eval(mc, ap, a2p, m2p, m3p, n);
  }

  void poly_8_eval_low_mem_single(float const * mc, float * ap, float * a2p,
			  float * m2p, int n) {
    poly_8_eval_low_mem(mc, ap, a2p, m2p, n);
  }
  void poly_8_eval_low_mem_double(double const * mc, double * ap, double * a2p,
			  double * m2p, int n) {
    poly_8_eval_low_mem(mc, ap, a2p, m2p, n);
  }

  void poly_5_eval_single(float const * mc, float * ap, float * a2p,
			  float * m2p, int n) {
    poly_5_eval(mc, ap, a2p, m2p, n);
  }
  void poly_5_eval_double(double const * mc, double * ap, double * a2p,
			  double * m2p, int n) {
    poly_5_eval(mc, ap, a2p, m2p, n);
  }
  
  float trace_XmX2_single(float const * ap, int n) {
    return trace_XmX2(ap, n);
  }
  double trace_XmX2_double(double const * ap, int n) {
    return trace_XmX2(ap, n);
  }

} // extern "C"

#ifndef COMPLEX_STEP_SOLVER_HEADER
#define COMPLEX_STEP_SOLVER_HEADER
#include <cassert>
#include <vector>
#include <complex>
#include <limits>
#include <algorithm>

extern "C" void dgesv_(const int *n, const int *nrhs, double *a, const int *lda, int *ipiv, double *b, const int *ldb, int *info);
extern "C" void sgesv_(const int *n, const int *nrhs, float *a, const int *lda, int *ipiv, float *b, const int *ldb, int *info);
void gesv(const int *n, const int *nrhs, double *a, const int *lda,
	  int *ipiv, double *b, const int *ldb, int *info) {
  dgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}
void gesv(const int *n, const int *nrhs, float *a, const int *lda,
	  int *ipiv, float *b, const int *ldb, int *info) {
  sgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}

template<class Objective_function>
struct Complex_step_solver {
  typedef typename Objective_function::value_type value_type;
  Objective_function const & obj_fun;
  Complex_step_solver(Objective_function const & obj_fun) : obj_fun(obj_fun) {}
  void partial_derivative(std::vector<value_type> const & v,
			  int const & i,
			  std::vector<value_type> & result) {
    std::vector<std::complex<value_type> > v_complex_step;
    v_complex_step.reserve(v.size());
    v_complex_step.assign(v.begin(), v.end());
    value_type h = std::numeric_limits<value_type>::epsilon();
    v_complex_step[i] += std::complex<value_type>(0,h);
    std::vector<std::complex<value_type> > tmp;
    obj_fun.fun(v_complex_step, tmp);
    result.resize( tmp.size() );
    for (unsigned int ind = 0; ind<tmp.size(); ind++) 
      result[ind] = tmp[ind].imag()/h;
  }

  void get_jacobian(std::vector<value_type> const & v,
		    std::vector<value_type> & jacobian_colwise) {
    jacobian_colwise.resize(v.size()*v.size());
    for (unsigned int ind = 0; ind<v.size(); ind++) {
      std::vector<value_type> partial(v.size());
      partial_derivative(v, ind, partial);
      for (unsigned int ind2 = 0; ind2<partial.size(); ind2++)
	jacobian_colwise[ind*v.size()+ind2] = partial[ind2];    
    }
  }

  void step_newton(std::vector<value_type> const & v,
		   std::vector<value_type>  & vnew,
		   value_type & v_maxabs,
		   value_type & dv_maxabs) {
    std::vector<value_type> jacobian_colwise;
    get_jacobian(v, jacobian_colwise);
    // Get rhs vector in vnew
    obj_fun.fun(v, vnew);
    int n = v.size();
    int nrhs = 1;
    value_type* a = &jacobian_colwise[0];
    int lda = n;
    std::vector<int> ipiv_vec(n);
    int* ipiv = &ipiv_vec[0];
    value_type* b = &vnew[0];
    int ldb = n;
    int info;
    gesv(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    assert(info == 0);
    // Ok we have vnew = inv(J)*fun(v)
    auto compare_abs = [](value_type a, value_type b) { return std::abs(a) < std::abs(b); };
    // Compute max absolute element of correction
    dv_maxabs = std::abs( *std::max_element(vnew.begin(), vnew.end(), compare_abs) );
    // Now take the Newton step
    for (unsigned int ind = 0; ind<v.size(); ind++) 
      vnew[ind] = v[ind]-vnew[ind];
    // Compute max absolute element of approximate solution
    v_maxabs = std::abs( *std::max_element(vnew.begin(), vnew.end(), compare_abs) );
  }
};
#endif

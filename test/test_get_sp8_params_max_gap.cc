#include <iostream>
#include <iomanip>

#include <cassert>
#include "sp8_acc.h"

template<typename T>
static void test_get_sp8_params_max_gap(const T L_outer, const T L_inner,
					const T H_inner, const T H_outer) {
  std::cout << "Running test_get_sp8_params_max_gap with Lout, Lin, Hin, Hout : "
	    << std::setw(6) << L_outer << ", "
	    << std::setw(6) << L_inner << ", "
	    << std::setw(6) << H_inner << ", "
	    << std::setw(6) << H_outer
	    << std::endl;
  int n_values = 100;
  std::vector<T> v;
  T kappa = 0.01;
  T Csp8;
  T qsp8;
  sp8::get_sp8_params_max_gap(L_outer, L_inner, H_inner, H_outer, kappa,
			      v, Csp8, qsp8);
  // Check if values are in [0,1]
  for (int ind = 0; ind <= n_values; ind++) {
    T x = ind/T(n_values);
    T y = sp8::sp8_eval(v,x);
    assert( y > -std::sqrt(std::numeric_limits<T>::epsilon()) );
    assert( y < 1+std::sqrt(std::numeric_limits<T>::epsilon()) );
  }
  // Check no mixing of occ-vir parts and that gap is increasing
  T vir_max = 0;
  T occ_min = 1;
  for (int ind = 0; ind <= n_values; ind++) {
    T x = L_inner*ind/T(n_values);
    T y = sp8::sp8_eval(v,x);
    vir_max = std::max(vir_max,y);
  }
  for (int ind = 0; ind <= n_values; ind++) {
    T x = H_inner+(1-H_inner)*ind/T(n_values);
    T y = sp8::sp8_eval(v,x);
    occ_min = std::min(occ_min,y);
  }
  assert(occ_min > vir_max);                 // No mixing
  assert(occ_min-vir_max > H_inner-L_inner); // Increasing gap  
}

int main(int argc, char* const  argv[]){
  // See that the program does not fail for certain input some of
  // which may be challenging for the solver.
  test_get_sp8_params_max_gap<double>(0.5,   0.55, 0.66, 0.77);
  test_get_sp8_params_max_gap<double>(0.011,  0.9, 0.91, 0.92);
  test_get_sp8_params_max_gap<double>(0.0,   0.01, 0.02, 0.03);
  test_get_sp8_params_max_gap<double>(0.0,   0.0001, 0.0002, 0.0003);
  test_get_sp8_params_max_gap<double>(0.0,   0.0001, 0.0002, 1.0);
  test_get_sp8_params_max_gap<double>(0.8,   0.95, 0.96, 0.99);
  test_get_sp8_params_max_gap<double>(0.011,  0.5, 0.51, 0.989);
  test_get_sp8_params_max_gap<double>(0.011,  0.011, 0.988, 0.989);
  test_get_sp8_params_max_gap<double>(0.011,  0.011, 0.012, 0.989);
  test_get_sp8_params_max_gap<double>(0.25,   0.6, 0.65, 0.995);
  test_get_sp8_params_max_gap<double>(0.005,  0.6, 0.65, 0.9);
  test_get_sp8_params_max_gap<double>(0.0,   0.49, 0.51, 1.0);
  test_get_sp8_params_max_gap<double>(0.48,  0.49, 0.51, 0.52);
  std::cout << argv[0] << " finished OK!" << std::endl;
  return 0;
}

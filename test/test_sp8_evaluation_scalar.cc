#include <iostream>
#include <iomanip>

#include <cassert>
#include "sp8_acc.h"
#include "sastre_poly8_eval.h"

/** Test checking that three different evaluations of an sp8 polynomial
    gives the same value. Polynomial evaluated directly from
    original form, via calculation of the canonical form as sum of
    monomials, and using Sastres method.
    @param v coefficients of sp8 polynomial in original representation
    @param x point of evaluation
*/
template<typename T>
static void test_sp8_evaluation(std::vector<T> const &  v, T const x) {
  T sp8_poly_value = sp8::sp8_eval(v,x);
  std::vector<T> mc;
  sp8::get_sp8_monomial_coefficients(v, mc);
  assert(mc.size() == 9);
  T sp8_poly_value_monomial = 0;
  T x_pow = 1;
  for(int ind = 8; ind >= 0; ind--) {
    sp8_poly_value_monomial += mc[ind]*x_pow;
    x_pow = x_pow*x;
  }
  T diff_sp8_mono = std::abs(sp8_poly_value - sp8_poly_value_monomial);
  T sp8_poly_value_sastre = x;
  T x2 = x*x;
  sp8::sastre_poly_8_eval(mc, sp8_poly_value_sastre, x2);
  T diff_sp8_sastre = std::abs(sp8_poly_value - sp8_poly_value_sastre);
  //  std::cout << std::setprecision(15);
  //  std::cout << "sastre: " << sp8_poly_value_sastre << std::endl;
  //  std::cout << "monom : " << sp8_poly_value_monomial << std::endl;
  //  std::cout << "sp8   : " << sp8_poly_value << std::endl;
  assert(diff_sp8_mono   < std::sqrt(std::numeric_limits<T>::epsilon()));
  assert(diff_sp8_sastre < std::sqrt(std::numeric_limits<T>::epsilon()));
}

int main(int argc, char* const  argv[]){
  typedef double real;
  real x = 0.5;
  // <4,3> L,H : 0.5, 0.6
  std::vector<real> v43 = {
    -36328.6859565887,
    -8.50474862353282e-31,
    0.0341080318412086,
    0.130423556398877,
    0.27072434777572,
    0.420393998539823,
    0.703968652953286,
    0.855051714130911,
    0.961997546899386
  };
  // 3.616e+04, 2.515e-01, 0.0373, 0.143, 0.301, 0.485, 0.659, 0.899, 0.971};
  test_sp8_evaluation(v43,x);
  test_sp8_evaluation(v43,(real)0.);
  test_sp8_evaluation(v43,(real)1.);

  // <5,2> L,H : 0.5, 0.6
  std::vector<real> v52 = {
    13838.0776158176,
    0.0122046264023502,
    0.0277339412746366,
    0.106051726505751,
    0.22062904850672,
    0.348187204925736,
    0.455602321566929,
    0.793044473864489,
    0.947094433741181
  };
  test_sp8_evaluation(v52,x);
  test_sp8_evaluation(v52,(real)0.);
  test_sp8_evaluation(v52,(real)1.);

  std::cout << argv[0] << " finished OK!" << std::endl;
  return 0;
}

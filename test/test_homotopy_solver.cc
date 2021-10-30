#include <iostream>
#include <iomanip>

#include <cassert>
#include "sp8_acc.h"

template<typename Objective_fun>
static void test_homotopy_solver(typename Objective_fun::value_type L_target,
				 typename Objective_fun::value_type H_target) {
  typedef typename Objective_fun::value_type value_type;
  std::vector<value_type> v;
  Homotopy_solver<Objective_fun> solver;
  solver(L_target, H_target, v);
  assert(v.size() == 9);
  // Check that all residual elements abs values are below some tolerance
  value_type residual_tolerance = std::sqrt(std::numeric_limits<value_type>::epsilon()); // this is ad hoc
  Objective_fun objfun(L_target,H_target);
  std::vector<value_type> tmp;
  objfun.fun(v,tmp);
  assert(std::all_of(tmp.begin(), tmp.end(), [residual_tolerance](value_type el){return std::abs(el)<residual_tolerance;}));
  // Check that order of roots is correct
  assert( objfun.correct_order_of_roots(v) );
  //  for (int ind = 0;ind < v.size();ind++)    
  //    std::cout << std::setprecision(15) << v[ind] << std::endl;
}

template<typename T>
static void test_flip() {
    std::vector<T> v = {3.616269094358620e+04, 2.515360660658102e-01, 3.734847296159297e-02, 1.435436349152357e-01, 3.017655595515090e-01, 4.859232582684261e-01, 6.598969344170539e-01, 8.998903128223554e-01, 9.714353568216227e-01};
    std::vector<T> v_flipped;
    get_flipped_polynomial(v, v_flipped);
    T x = 0.384;
    T ref = 1-sp8(v,1-x);
    T value = sp8(v_flipped,x);
    assert( std::abs(ref-value) < std::sqrt(std::numeric_limits<T>::epsilon()));
}

int main(int argc, char* const  argv[]){
  double L_target = 0.5;
  double H_target = 0.6;
  test_homotopy_solver<Objective_Fun43<double> >(L_target, H_target);
  test_homotopy_solver<Objective_Fun52<double> >(L_target, H_target);
  test_homotopy_solver<Objective_Fun40<double> >(L_target, H_target);
  test_homotopy_solver<Objective_Fun50<double> >(L_target, H_target);
  L_target = 0.94;
  H_target = 0.95;
  test_homotopy_solver<Objective_Fun61<double> >(L_target, H_target);
  test_homotopy_solver<Objective_Fun60<double> >(L_target, H_target);
  L_target = 0.97;
  H_target = 0.98;
  test_homotopy_solver<Objective_Fun70<double> >(L_target, H_target);

  test_flip<double>();

  std::cout << argv[0] << " finished OK!" << std::endl;
}

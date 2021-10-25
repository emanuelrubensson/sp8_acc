#include <iostream>
#include <cassert>
#include "sp8_acc.h"

template<typename Objective_fun>
static void test_homotopy_solver(typename Objective_fun::value_type L_target,
				 typename Objective_fun::value_type H_target,
				 typename Objective_fun::value_type dv_tolerance_relative) {
  typedef typename Objective_fun::value_type value_type;
  std::vector<value_type> v;
  value_type v_maxabs;
  value_type dv_maxabs;
  Homotopy_solver<Objective_fun> solver;
  solver(L_target, H_target, dv_tolerance_relative, v, v_maxabs,dv_maxabs);
  // Check solution
  // Check that step size satisfies tolerance
  assert(dv_maxabs/v_maxabs < dv_tolerance_relative);
  // Check that all residual elements abs values are below some tolerance
  value_type residual_tolerance = dv_tolerance_relative; // this is ad hoc
  Objective_fun objfun(L_target,H_target);
  std::vector<value_type> tmp;
  objfun.fun(v,tmp);
  assert(std::all_of(tmp.begin(), tmp.end(), [residual_tolerance](value_type el){return std::abs(el)<residual_tolerance;}));
  //  for (int ind = 0;ind < v.size();ind++)    
  //    std::cout << v[ind] << std::endl;
}

int main(int argc, char* const  argv[]){
  double L_target = 0.5;
  double H_target = 0.6;
  double dv_tolerance_relative = 1e-10;
  test_homotopy_solver<Objective_Fun40<double> >(L_target,
						 H_target,
						 dv_tolerance_relative);
  test_homotopy_solver<Objective_Fun43<double> >(L_target,
						 H_target,
						 dv_tolerance_relative);
  test_homotopy_solver<Objective_Fun52<double> >(L_target,
						 H_target,
						 dv_tolerance_relative);
  std::cout << argv[0] << " finished OK!" << std::endl;
}

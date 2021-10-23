#include <iostream>
#include <cassert>
#include "sp8_acc.h"

int main(int argc, char* const  argv[]){
  double L_target = 0.5;
  double H_target = 0.6;
  std::vector<double> v;
  double v_maxabs;
  double dv_maxabs;
  Homotopy_solver<Objective_Fun43<double> > solver;
  solver(L_target,H_target,v,v_maxabs,dv_maxabs);
  // Check solution
  // Check step size
  assert(dv_maxabs/v_maxabs < 1e-10);
  // Check that all residual elements abs values are below 1e-10
  Objective_Fun43<double> objfun(L_target,H_target);
  std::vector<double> tmp;
  objfun.fun(v,tmp);
  assert(std::all_of(tmp.begin(), tmp.end(), [](double el){return std::abs(el)<1e-10;}));
  
  std::cout << argv[0] << " finished OK!" << std::endl;
  return 0;
}

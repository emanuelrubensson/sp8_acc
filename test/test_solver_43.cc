#include <iostream>
#include <cassert>
#include "sp8_acc.h"
#include "Complex_step_solver.h"

int main(int argc, char* const  argv[]){
  double L_target = 0.5;
  double H_target = 0.6;
  std::vector<double> v = Objective_Fun43<double>::v_start;
  double L_guess        = Objective_Fun43<double>::L_start;
  double H_guess        = Objective_Fun43<double>::H_start;
  int nsteps = 100;
  double L_step = (L_target-L_guess)/nsteps;
  double H_step = (H_target-H_guess)/nsteps;
  double L = L_guess;
  double H = H_guess;
  std::vector<double> tmp;
  double v_maxabs;
  double dv_maxabs;
  for (int ind = 0; ind<nsteps; ind++) {
    L = L+L_step;
    H = H+H_step;
    Objective_Fun43<double> objfun(L,H);
    Complex_step_solver<Objective_Fun43<double> > solver(objfun);
    solver.step_newton(v,tmp,v_maxabs,dv_maxabs);
    v.swap(tmp);
  }
  Objective_Fun43<double> objfun(L,H);
  Complex_step_solver<Objective_Fun43<double> > solver(objfun);
  for (int ind = 0; ind<4; ind++) {
    solver.step_newton(v,tmp,v_maxabs,dv_maxabs);
    v.swap(tmp);
    //    std::cout << "max(abs(v)), max(abs(dv)) = " << v_maxabs << ", " << dv_maxabs << std::endl;
  }
  // Check step size
  assert(dv_maxabs/v_maxabs < 1e-10);
  // Check that all residual elements abs values are below 1e-10
  objfun.fun(v,tmp);
  assert(std::all_of(tmp.begin(), tmp.end(), [](double el){return std::abs(el)<1e-10;}));
  
  std::cout << argv[0] << " finished OK!" << std::endl;
  return 0;
}

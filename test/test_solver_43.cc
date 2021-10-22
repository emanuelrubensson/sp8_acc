#include <iostream>
#include "sp8_acc.h"
#include "Complex_step_solver.h"

int main() {
  double L_target = 0.5;
  double H_target = 0.6;
  std::vector<double> v_guess_43 = {-3.828330354372377e+04, // c1
                                    0.000000000000000e+00,  // c2
                                    3.630457073236520e-02,
                                    1.391234727632784e-01,
                                    2.902225497081609e-01,
                                    4.558182467167558e-01,
                                    7.341809614537667e-01,
                                    8.676371783744240e-01,
                                    9.650762995000088e-01}; // for L,H=0.55, 0.65
  double L_guess = 0.55;
  double H_guess = 0.65;
  int nsteps = 100;
  double L_step = (L_target-L_guess)/nsteps;
  double H_step = (H_target-H_guess)/nsteps;
  double L = L_guess;
  double H = H_guess;
  std::vector<double> v = v_guess_43;
  std::vector<double> tmp;
  for (int ind = 0; ind<nsteps; ind++) {
    L = L+L_step;
    H = H+H_step;
    Objective_Fun43<double> objfun(L,H);
    Complex_step_solver<Objective_Fun43<double> > solver(objfun);
    solver.step_newton(v,tmp);
    v.swap(tmp);
  }
  Objective_Fun43<double> objfun(L,H);
  Complex_step_solver<Objective_Fun43<double> > solver(objfun);
  for (int ind = 0; ind<4; ind++) {
    solver.step_newton(v,tmp);
    v.swap(tmp);
  }
  
  for (unsigned int ind = 0; ind<v.size(); ind++) 
    std::cout << v[ind] << std::endl;
  return 0;
}

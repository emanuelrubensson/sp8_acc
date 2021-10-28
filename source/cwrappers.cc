#include "sp8_acc.h"

extern "C" {
  /** v_output should point to array with 9 elements */
  void get_sp8_params(const double L, const double H,
		      const int roots_left, const int roots_right,
		      double* v_output) {
    std::vector<double> v;
    Homotopy_solver_base<double>* solver = homotopy_solver_factory<double>(roots_left, roots_right);
    (*solver)(L, H, v);
    delete solver;
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
    get_sp8_monomial_coefficients(v, mc_vec);
    std::copy_n(mc_vec.begin(), 9, mc);
  }
}

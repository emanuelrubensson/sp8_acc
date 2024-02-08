#ifndef SASTRE_POLY8_EVAL_HEADER
#define SASTRE_POLY8_EVAL_HEADER
#include <vector>
#include <cmath>

template<typename T_matrix_scalar>
struct Scalar_proxy {
  T_matrix_scalar x;
  Scalar_proxy(T_matrix_scalar const x) :x(x) {}
  inline Scalar_proxy operator=(Scalar_proxy const & other) {
    x = other.x;
    return *this;
  }
  inline static Scalar_proxy<T_matrix_scalar> multiply_symm_product(Scalar_proxy<T_matrix_scalar> const & A,
								    Scalar_proxy<T_matrix_scalar> const & B) {
    return Scalar_proxy(A.x*B.x);
  }
  inline void scale_and_add(T_matrix_scalar const a,
			    T_matrix_scalar const b,
			    Scalar_proxy<T_matrix_scalar> const & M1) {
    x = a*x + b*M1.x;
  }
  inline void add_scaled_identity(T_matrix_scalar const a) {
    x = x+a;
  }
  inline void negate() {
    x = -x;
  }
};


template<typename T_scalar, typename T_matrix_scalar>
void sastre_poly_8_eval(std::vector<T_scalar> const & mc,
			T_matrix_scalar & x,
			std::true_type) {
  Scalar_proxy X(x);
  sastre_poly_8_eval(mc, X);
  x = X.x;
}

template<typename T_scalar, typename T_matrix>
void sastre_poly_8_eval(std::vector<T_scalar> const & mc,
			T_matrix & A,
			std::false_type) {
  const T_scalar b0=mc[8], b1=mc[7], b2=mc[6];
  const T_scalar b3=mc[5], b4=mc[4], b5=mc[3];
  const T_scalar b6=mc[2], b7=mc[1], b8=mc[0];
  if (b8<0) {
    // Avoid complex numbers by switching evaluating with minus sign and
    // subsequently reversing
    std::vector<T_scalar> mc_minus(9);
    for (unsigned int ind = 0; ind<9; ind++)      
      mc_minus[ind]=-mc[ind];
    sastre_poly_8_eval(mc_minus, A);
    A.negate();
    return;
  }
  
  // plus minus: Two options possible.
  // "Correct"  option leads to a e2_num_sqrt > 0 in code below?
  const T_scalar c4 = std::sqrt(b8);
  const T_scalar c3 = b7 / (2.0 * c4);
  const T_scalar d2_plus_e2 = (b6 - c3*c3) / c4;
  const T_scalar d1 = (b5 - c3 * d2_plus_e2) / c4;
  
  const T_scalar e2_num_sqrt =  std::pow(d1 - (c3 / c4) * d2_plus_e2, 2) +
    4.0 * (c3 / c4) * (b3 + (std::pow(c3,2) / c4) * d1 - (c3 / c4) * b4);

  // If e2_num_sqrt<0 we may want to try other sign in c4 ?
  // Here are two options on the square roor: plus minus
  const T_scalar e2_num = (c3 / c4) * d2_plus_e2 - d1 + std::sqrt(e2_num_sqrt);
  const T_scalar e2 = e2_num / (2 * c3 / c4);
  const T_scalar d2 = d2_plus_e2 - e2;
  const T_scalar f2 = b2;
  const T_scalar f1 = b1;
  const T_scalar f0 = b0;
  
  const T_scalar e0 = (b3 - d1 * e2) / c3; // Not explicitly
					   // documented by sastre?
  const T_scalar one = 1.0;
  T_matrix M1 = T_matrix::multiply_symm_product(A, A); // M1 = A*A
  T_matrix M2 = M1;                           // M2 = M1
  M2.scale_and_add(c4, c3, A);                // M2 = c4*M2 + c3*A
  T_matrix M3 = T_matrix::multiply_symm_product(M1,M2); // M3 = M1*M2
  M2 = M3;                                    // M2 = M3
  M2.scale_and_add(one, d2, M1);              // M2 = M2 + d2*M1
  M2.scale_and_add(one, d1, A);               // M2 = M2 + d1*A
  A.scale_and_add(f1, f2, M1);                // A = f1*A + f2*M1
  A.scale_and_add(one, e0, M3);               // A = A + e0*M3
  A.add_scaled_identity(f0);                  // A = A+f0*I
  M1.scale_and_add(e2, one, M3);              // M1 = e2*M1 + M3
  M3 = T_matrix::multiply_symm_product(M1,M2); // M3 = M1*M2
  A.scale_and_add(one, one, M3);              // A = A + M3
}


template<typename T_scalar, typename T_matrix>
void sastre_poly_8_eval(std::vector<T_scalar> const & mc,
			T_matrix & x) {
  sastre_poly_8_eval(mc, x, std::is_floating_point<T_matrix>());
}

#endif

#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "sastre_poly8_eval.h"
#include "poly8_eval.h"
#include "poly8_eval_low_mem.h"
extern "C" void dgemm_(const char *ta,const char *tb,
		       const int *n, const int *k, const int *l,
		       const double *alpha,const double *A,const int *lda,
		       const double *B, const int *ldb,
		       const double *beta, double *C, const int *ldc);
extern "C" void sgemm_(const char *ta,const char *tb,
		       const int *n, const int *k, const int *l,
		       const float *alpha,const float *A,const int *lda,
		       const float *B, const int *ldb,
		       const float *beta, float *C, const int *ldc);
void gemm(const char *ta,const char *tb,
	  const int *n, const int *k, const int *l,
	  const double *alpha,const double *A,const int *lda,
	  const double *B, const int *ldb,
	  const double *beta, double *C, const int *ldc) {
  dgemm_(ta, tb, n, k, l, alpha, A, lda, B, ldb, beta, C, ldc);
}
void gemm(const char *ta,const char *tb,
	  const int *n, const int *k, const int *l,
	  const float *alpha,const float *A,const int *lda,
	  const float *B, const int *ldb,
	  const float *beta, float *C, const int *ldc) {
  sgemm_(ta, tb, n, k, l, alpha, A, lda, B, ldb, beta, C, ldc);
}
  
template<typename T_scalar>
struct Matrix {
  typedef T_scalar value_type;
  unsigned int n;
  std::vector<T_scalar> elements; // n*n
  Matrix(unsigned int n = 0) :n(n), elements(n*n) {}
  Matrix const & operator=(Matrix const & other) {
    this->n = other.n;
    this->elements.resize(n*n);
    std::copy_n(other.elements.begin(), n*n, this->elements.begin());
    return *this;
  }  
  inline void print_elements() const {
    std::cout << std::setprecision(15);
    for (int ind = 0; ind < n*n;ind++)
      std::cout << elements[ind] << "  ";
    std::cout << std::endl;
  }
  void multiply(Matrix<T_scalar> const & A,
		Matrix<T_scalar> const & B,
		value_type const alpha = 1.0,
		value_type const beta = 0.0) {
    int n = A.n;
    this->n = n;
    this->elements.resize(n*n);
    gemm("N", "N", &n, &n, &n, &alpha, &A.elements[0], &n,
	 &B.elements[0], &n, &beta, &this->elements[0], &n);
  }
  void scale_and_add(T_scalar const a,
		     T_scalar const b,
		     Matrix<T_scalar> const & A) {
    assert(this->n == A.n);
    for(unsigned int ind = 0; ind < this->n*this->n;ind++)
      this->elements[ind] = a*this->elements[ind] + b*A.elements[ind];
  }
  void add_scaled_identity(T_scalar const a) {
    for(unsigned int ind = 0; ind < this->n*this->n; ind+=this->n+1)
      this->elements[ind] += a;
  }
  void negate() {
    for(unsigned int ind = 0; ind < this->n*this->n;ind++)
      this->elements[ind] = -this->elements[ind];
  }
};

template<typename T>
static T maxabs_diff(std::vector<T> const & x, std::vector<T> const & y) {
  T maxdiff = 0;
  assert( x.size() == y.size() );
  for(unsigned int ind = 0; ind < x.size(); ind++) {
    T diff = std::abs(x[ind] - y[ind]);
    if (diff > maxdiff)
      maxdiff = diff;
  }
  return maxdiff;
}

template<typename T>
static void test_poly8_evaluation(std::vector<T> const & mc, Matrix<T> const & A, Matrix<T> const & sol_ref) {
  bool verbose_output = false;
  {
    Matrix<T> poly_value_sastre = A;
    Matrix<T> A2;
    A2.multiply(A, A);                          // A2 = A*A
    sp8::sastre_poly_8_eval(mc, poly_value_sastre, A2);
    if (verbose_output)
      poly_value_sastre.print_elements();
    assert( maxabs_diff(poly_value_sastre.elements, sol_ref.elements) < std::sqrt(std::numeric_limits<T>::epsilon()) );
  }
  {
    Matrix<T> poly_value = A;
    Matrix<T> A2;
    A2.multiply(A, A);                          // A2 = A*A
    sp8::poly_8_eval(mc, poly_value, A2);
    if (verbose_output)
      poly_value.print_elements();
    assert( maxabs_diff(poly_value.elements, sol_ref.elements) < std::sqrt(std::numeric_limits<T>::epsilon()) );
  }
  {
    Matrix<T> poly_value_low_mem = A;
    Matrix<T> A2;
    A2.multiply(A, A);                          // A2 = A*A
    sp8::poly_8_eval_low_mem(mc, poly_value_low_mem, A2);
    if (verbose_output)
      poly_value_low_mem.print_elements();
    assert( maxabs_diff(poly_value_low_mem.elements, sol_ref.elements) < std::sqrt(std::numeric_limits<T>::epsilon()) );
  }
}

int main(int argc, char* const  argv[]){
  {
    Matrix<double> A(2);
    A.elements[0] = 0.1;
    A.elements[1] = 0.2;
    A.elements[2] = 0.3;
    A.elements[3] = 0.4;
    std::vector<double> mc = {2.48016e-05, 0.000198413, 0.00138889, 0.00833333, 0.0416667, 0.166667, 0.5, 1, 1};
    Matrix<double> sol_ref(2);
    sol_ref.elements[0] = 1.1420935592955;
    sol_ref.elements[1] = 0.260350762687821;
    sol_ref.elements[2] = 0.390526144031731;
    sol_ref.elements[3] = 1.53261970332723;
    test_poly8_evaluation(mc, A, sol_ref);
  }
  {
    // Symmetric matrix without symmetric storage
    Matrix<double> A(2);
    A.elements[0] = 0.1;
    A.elements[1] = 0.3;
    A.elements[2] = 0.3;
    A.elements[3] = 0.4;
    std::vector<double> mc = {2.48016e-05, 0.000198413, 0.00138889, 0.00833333, 0.0416667, 0.166667, 0.5, 1, 1};
    Matrix<double> sol_ref(2);
    sol_ref.elements[0] = 1.160696024142013;
    sol_ref.elements[1] = 0.392471025439357;
    sol_ref.elements[2] = 0.392471025439357;
    sol_ref.elements[3] = 1.553167049581370;
    test_poly8_evaluation(mc, A, sol_ref);
  }
  std::cout << argv[0] << " finished OK!" << std::endl;
  return 0;
}



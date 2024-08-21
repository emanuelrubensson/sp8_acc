#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "sastre_poly8_eval.h"
#include "poly8_eval.h"
#include "poly8_eval_low_mem.h"
#include "test_utils.h"

template<typename T_scalar>
struct Matrix_symmetric {
  typedef T_scalar value_type;
  unsigned int n;
  std::vector<T_scalar> elements; // n*n, but only upper triangle is referenced 
  Matrix_symmetric(unsigned int n = 0) :n(n), elements(n*n) {}
  Matrix_symmetric const & operator=(Matrix_symmetric const & other) {
    this->n = other.n;
    this->elements.resize(n*n);
    std::copy_n(other.elements.begin(), other.elements.size(), this->elements.begin());
    return *this;
  }
  T_scalar const & operator()(const int i, const int j) const {
    if (i <= j)
      // Element in upper triangle
      return this->elements[i+j*this->n];
    // Element in lower triangle
    return this->elements[j+i*this->n];
  }
  T_scalar & operator()(const int i, const int j) {
    if (i <= j)
      return this->elements[i+j*this->n];
    return this->elements[j+i*this->n];
  }
  void multiply(Matrix_symmetric<T_scalar> const & A,
		Matrix_symmetric<T_scalar> const & B,
		T_scalar const alpha = 1.0,
		T_scalar const beta = 0.0) {
    int n = A.n;
    if (beta == 0.0) {
      this->n = n;
      this->elements.resize(n*n);
      std::fill(this->elements.begin(), this->elements.end(), 0.0);
    }
    else
      assert(this->n == n);
    for (int colC = 0; colC<n; colC++)
      for (int rowC = 0; rowC<=colC; rowC++) {
	(*this)(rowC,colC) *= beta;
	for (int k = 0; k<n; k++)
	  (*this)(rowC,colC) += alpha*A(rowC,k)*B(k,colC);
      }
  }
  void scale_and_add(T_scalar const a,
		     T_scalar const b,
		     Matrix_symmetric<T_scalar> const & A) {
    assert(this->n == A.n);
    for (int col = 0; col<n; col++)
      for (int row = 0; row<=col; row++)
	(*this)(row,col) = a*(*this)(row,col) + b*A(row,col);
  }
  void add_scaled_identity(T_scalar const a) {
    for(unsigned int ind = 0; ind < this->n*this->n; ind+=this->n+1)
      this->elements[ind] += a;
  }
  void negate() {
    for (int col = 0; col<n; col++)
      for (int row = 0; row<=col; row++)
	(*this)(row,col) = -(*this)(row,col);
  }
  void print_elements() const {
    std::cout << std::setprecision(15);
    for(unsigned int ind = 0; ind < this->n*this->n;ind++) 
      std::cout << std::setprecision(15) << this->elements[ind] << "  ";
    std::cout << std::endl;
  }
};

template<typename T>
static void test_poly8_evaluation(std::vector<T> const & mc, Matrix_symmetric<T> const & A, Matrix_symmetric<T> const & sol_ref) {
  bool verbose_output = true;
  {
    Matrix_symmetric<T> sp8_poly_value_sastre = A;
    Matrix_symmetric<T> A2;
    A2.multiply(A, A);                          // A2 = A*A
    sp8::sastre_poly_8_eval(mc, sp8_poly_value_sastre, A2);
    if (verbose_output)
      sp8_poly_value_sastre.print_elements();
    assert( maxabs_diff(sp8_poly_value_sastre.elements, sol_ref.elements) < std::sqrt(std::numeric_limits<T>::epsilon()) );
  }
  {
    Matrix_symmetric<T> sp8_poly_value = A;
    Matrix_symmetric<T> A2;
    A2.multiply(A, A);                          // A2 = A*A
    sp8::poly_8_eval(mc, sp8_poly_value, A2);
    if (verbose_output)
      sp8_poly_value.print_elements();
    assert( maxabs_diff(sp8_poly_value.elements, sol_ref.elements) < std::sqrt(std::numeric_limits<T>::epsilon()) );
  }
  {
    Matrix_symmetric<T> sp8_poly_value_low_mem = A;
    Matrix_symmetric<T> A2;
    A2.multiply(A, A);                          // A2 = A*A
    sp8::poly_8_eval_low_mem(mc, sp8_poly_value_low_mem, A2);
    if (verbose_output)
      sp8_poly_value_low_mem.print_elements();
    assert( maxabs_diff(sp8_poly_value_low_mem.elements, sol_ref.elements) < std::sqrt(std::numeric_limits<T>::epsilon()) );
  }
}

int main(int argc, char* const  argv[]){
  Matrix_symmetric<double> A(2);
  A.elements[0] = 0.1;
  // A.elements[1] = 0; only upper triangle referenced
  A.elements[2] = 0.3;
  A.elements[3] = 0.4;
  std::vector<double> mc = {2.48016e-05, 0.000198413, 0.00138889, 0.00833333, 0.0416667, 0.166667, 0.5, 1, 1};
  Matrix_symmetric<double> sol_ref(2);
  sol_ref.elements[0] = 1.160696024142013;
  // sol_ref.elements[1] = 0
  sol_ref.elements[2] = 0.392471025439357;
  sol_ref.elements[3] = 1.553167049581370;
  test_poly8_evaluation(mc, A, sol_ref);
  std::cout << argv[0] << " finished OK!" << std::endl;
  return 0;
}

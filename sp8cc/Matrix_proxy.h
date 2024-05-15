#ifndef MATRIX_PROXY_HEADER
#define MATRIX_PROXY_HEADER
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
namespace sp8 {
  template<typename T_scalar>
  struct Matrix_proxy {
    T_scalar* elements;
    int n;
    Matrix_proxy(T_scalar* elements, int n)
      :elements(elements), n(n)
    {}
    Matrix_proxy(const Matrix_proxy&) = delete;
    inline Matrix_proxy const & operator=(Matrix_proxy const & other) {
      assert(this->n == other.n);
      std::copy_n(other.elements, n*n, this->elements);
      return *this;
    }
    inline void multiply(Matrix_proxy<T_scalar> const & A,
			 Matrix_proxy<T_scalar> const & B) {
      const T_scalar one  = 1.0;
      const T_scalar zero = 0.0;
      int n = this->n;
      assert(n   == A.n);
      assert(B.n == A.n);
      gemm("N", "N", &n, &n, &n, &one, A.elements, &n,
	   B.elements, &n, &zero, this->elements, &n);
    }
    inline void scale_and_add(T_scalar const a,
			      T_scalar const b,
			      Matrix_proxy<T_scalar> const & A) {
      assert(this->n == A.n);
      for(unsigned int ind = 0; ind < this->n*this->n;ind++)
	this->elements[ind] = a*this->elements[ind] + b*A.elements[ind];
    }
    inline void add_scaled_identity(T_scalar const a) {
      for(unsigned int ind = 0; ind < this->n*this->n; ind+=this->n+1)
	this->elements[ind] += a;
    }
    inline void negate() {
      for(unsigned int ind = 0; ind < this->n*this->n; ind++)
	this->elements[ind] = -this->elements[ind];
    }
  };

} // end namespace sp8
#endif

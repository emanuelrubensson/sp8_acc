#ifndef SCALAR_PROXY_HEADER
#define SCALAR_PROXY_HEADER

/** Proxy template for handling the special case where the 'matrix' is
    a single floating point number. This template also provides a
    template for the implementation of Matrix classes to be used with
    the polynomial evaluation schemes.    
*/
template<typename T_scalar>
struct Scalar_proxy {
  typedef T_scalar value_type;
  T_scalar x;
Scalar_proxy(T_scalar const x = T_scalar()) :x(x) {}
  inline Scalar_proxy operator=(Scalar_proxy const & other) {
    x = other.x;
    return *this;
  }
    inline void multiply(Scalar_proxy<T_scalar> const & A,
			 Scalar_proxy<T_scalar> const & B,
			 value_type const alpha = 1.0,
			 value_type const beta = 0.0) {
      this->x = beta*this->x + alpha*A.x*B.x;
    }
  inline void scale_and_add(T_scalar const a,
			    T_scalar const b,
			    Scalar_proxy<T_scalar> const & M1) {
    x = a*x + b*M1.x;
  }
  inline void add_scaled_identity(T_scalar const a) {
    x = x+a;
  }
  inline void negate() {
    x = -x;
  }
};

#endif

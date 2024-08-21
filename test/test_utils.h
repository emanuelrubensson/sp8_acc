#ifndef TESTUTILS_HEADER
#define TESTUTILS_HEADER
#include <vector>
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

#endif

# SP8-ACC Library

The SP8-ACC library implements methods for computing matrix step
functions by recursive polynomial expansion. It includes the SP8-ACC
algorithm, algorithms for efficient matrix polynomial evaluation, and
reference implementations of related methods.

## Features

- **Generic header-only C++ template** for evaluating arbitrary
  degree-eight matrix polynomials using only three non-scalar matrix
  multiplications while requiring storage for only three matrices. The
  template is parameterized by a matrix class implementing scaling,
  addition, and multiplication operations.

- **Python implementation** of the SP8-ACC algorithm for recursive
  polynomial expansion of the matrix step function.

- **C++ tools** for generating the coefficients of the component
  polynomials used by the SP8-ACC algorithm.

- **Reference implementations** of several alternative matrix
  polynomial evaluation schemes and recursive polynomial expansion
  methods for comparison.

## Repository contents

The repository contains

- `sp8cc/`: Header-only C++ library for matrix polynomial evaluation,
generation of coefficients for the SP8-ACC component polynomials, and
related utilities. Also contains C wrappers and a lightweight proxy
matrix class for matrices stored in contiguous column-major format,
primarily intended for the C interface and as a reference
implementation for integrating other matrix libraries.

- `sp8py/`: Python implementation of the recursive expansion
  algorithms built on the C wrappers in sp8cc.

- `test/` and `testpy/`: Unit tests for both the C++ and Python implementations.

- `plots/`: Python code to generate figures.

## Reference

The algorithms implemented in this repository are described in the paper

> Emanuel H. Rubensson, Elias Jarlebring, and Gustaf Lorentzon,
> *Recursive expansion of the matrix step function using polynomials of
> degree eight*, 2026.
> https://doi.org/10.48550/arXiv.2606.24701

## Dependencies

The core C++ library is header-only and has no external dependencies.
Building the remaining components, including the wrappers, Python
bindings, and tests, requires the LAPACK linear algebra library.

## Getting started

Build the library and run the tests

```bash
make
make check
```

Generate the figures from the paper

```bash
cd plots
make
```

/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2016 Imperial College London
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#define _SCL_SECURE_NO_WARNINGS

#include "mirtk/Transformation.h"
#include "mirtk/MultiLevelTransformation.h"
#include "mirtk/Matrix.h"

#include "boost/numeric/ublas/lu.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/vector.hpp"

#include <cstdint>
#include <tuple>
#include <algorithm>


// =============================================================================
// Newton-Raphson solver
// =============================================================================

namespace boost { namespace math { namespace tools {

namespace detail {

template<class T>
bool invert_matrix(const boost::numeric::ublas::matrix<T> & operand,
  boost::numeric::ublas::matrix<T> & result)
{
  namespace ublas = boost::numeric::ublas;
  typedef ublas::matrix<T> matrix_type;
  typedef typename matrix_type::size_type size_type;
  typedef typename matrix_type::value_type value_type;
  typedef ublas::permutation_matrix<size_type> permutation_matrix_type;
  typedef ublas::identity_matrix<value_type> identity_matrix_type;

  // Create temporary copy of operand.
  matrix_type A(operand);
  // Create a permutation matrix for the LU-factorization.
  permutation_matrix_type P(A.size1());
  // Perform LU-factorization.
  size_type singular = lu_factorize(A, P);
  // Set inverse to identity.
  result.assign(identity_matrix_type(A.size1()));
  // Substitute back to compute the inverse.
  lu_substitute(A, P, result);

  return (singular == 0);
}

} // namespace detail

template<class F>
bool multivariate_newton_raphson_iterate(F f,
  typename F::vector_type & x, int digits,
  uintmax_t max_iter=100)
{
  typedef typename F::vector_type vector_type;
  typedef typename F::matrix_type matrix_type;
  typedef typename vector_type::value_type value_type;
  typedef typename vector_type::size_type size_type;

  value_type factor = static_cast<value_type>(ldexp(1.0, 1-digits));
  size_type n = x.size();
  vector_type d_x(n), e_x(n), f0(n);
  matrix_type f1(n, n), f1_inv(n, n);
  uintmax_t iter=1;

  do
  {
    std::tie(f0, f1) = f(x);
    if(!detail::invert_matrix<value_type>(f1, f1_inv))
      return false;
    d_x = boost::numeric::ublas::prod(f1_inv, f0);
    x -= d_x;
    for(size_type ni = 0; ni < n; ni++)
      // FIXME: use zip iterator.
      e_x(ni) = fabs(d_x(ni)) - fabs(x(ni) * factor);
  }
  while((iter++ < max_iter) &&
         std::any_of(e_x.begin(), e_x.end(), [](value_type val){return val > 0;}));

  return true;
}


} } } // namespace boost::math::tools


namespace mirtk {

// =============================================================================
// Inverse global transformation
// =============================================================================

struct GlobalInverseTransformFunctor
{
  typedef boost::numeric::ublas::vector<double> vector_type;
  typedef boost::numeric::ublas::matrix<double> matrix_type;
  typedef std::tuple<vector_type, matrix_type>  return_type;

  GlobalInverseTransformFunctor(const Transformation * transformation,
    double x, double y, double z, double t, double t0)
    : transformation_(transformation), x_(x), y_(y), z_(z), t_(t), t0_(t0)
  {
  }

  return_type operator()(const vector_type & x)
  {
    typedef typename vector_type::value_type value_type;

    value_type x1 = x(0), x2 = x(1), x3 = x(2);
    transformation_->GlobalTransform(x1, x2, x3, t_, t0_);
    vector_type f0(3);
    f0(0) = x1 - x_;
    f0(1) = x2 - y_;
    f0(2) = x3 - z_;

    Matrix J(3, 3);
    transformation_->GlobalJacobian(J, x(0), x(1), x(2), t_, t0_);
    matrix_type f1(3, 3);
    for (int icol = 0; icol < 3; ++icol)
    for (int irow = 0; irow < 3; ++irow) {
      f1(irow, icol) = J(irow, icol);
    }

    return std::make_tuple(f0, f1);
  }

  const Transformation * transformation_;
  double x_;
  double y_;
  double z_;
  double t_;
  double t0_;
};

// -----------------------------------------------------------------------------
bool EvaluateGlobalInverse(const Transformation *T,
                           double &x, double &y, double &z, double t, double t0)
{
  using namespace boost::math::tools;
  typedef GlobalInverseTransformFunctor functor_type;
  typedef typename functor_type::vector_type vector_type;

  functor_type fun(T, x, y, z, t, t0);
  vector_type sol(3); sol(0) = x;  sol(1) = y;  sol(2) = z;
  bool success = multivariate_newton_raphson_iterate<functor_type>(fun, sol, 6);
  x = sol(0); y = sol(1); z = sol(2);
  return success;
}

// =============================================================================
// Inverse local transformation
// =============================================================================

struct LocalInverseTransformFunctor
{
  typedef boost::numeric::ublas::vector<double> vector_type;
  typedef boost::numeric::ublas::matrix<double> matrix_type;
  typedef std::tuple<vector_type, matrix_type>  return_type;

  LocalInverseTransformFunctor(const Transformation * transformation,
    double x, double y, double z, double t, double t0)
    : transformation_(transformation), x_(x), y_(y), z_(z), t_(t), t0_(t0)
  {
  }

  return_type operator()(vector_type x)
  {
    typedef typename vector_type::value_type value_type;

    value_type x1 = x(0), x2 = x(1), x3 = x(2);
    transformation_->LocalTransform(x1, x2, x3, t_, t0_);
    vector_type f0(3);
    f0(0) = x1 - x_;
    f0(1) = x2 - y_;
    f0(2) = x3 - z_;

    Matrix J(3, 3);
    transformation_->LocalJacobian(J, x(0), x(1), x(2), t_, t0_);
    matrix_type f1(3, 3);
    for (int icol = 0; icol < 3; ++icol)
    for (int irow = 0; irow < 3; ++irow) {
      f1(irow, icol) = J(irow, icol);
    }

    return std::make_tuple(f0, f1);
  }

  const Transformation * transformation_;
  double x_;
  double y_;
  double z_;
  double t_;
  double t0_;
};

// -----------------------------------------------------------------------------
bool EvaluateLocalInverse(const Transformation *T,
                          double &x, double &y, double &z, double t, double t0)
{
  using namespace boost::math::tools;
  typedef LocalInverseTransformFunctor functor_type;
  typedef typename functor_type::vector_type vector_type;

  functor_type fun(T, x, y, z, t, t0);
  vector_type sol(3); sol(0) = x;  sol(1) = y;  sol(2) = z;
  bool success = multivariate_newton_raphson_iterate<functor_type>(fun, sol, 6);
  x = sol(0); y = sol(1); z = sol(2);
  return success;
}

// =============================================================================
// Inverse transformation
// =============================================================================

struct InverseTransformFunctor
{
  typedef boost::numeric::ublas::vector<double> vector_type;
  typedef boost::numeric::ublas::matrix<double> matrix_type;
  typedef std::tuple<vector_type, matrix_type>  return_type;

  InverseTransformFunctor(const Transformation * transformation,
    double x, double y, double z, double t, double t0)
    : transformation_(transformation), x_(x), y_(y), z_(z), t_(t), t0_(t0)
  {
  }

  return_type operator()(vector_type x)
  {
    typedef typename vector_type::value_type value_type;

    value_type x1 = x(0), x2 = x(1), x3 = x(2);
    transformation_->Transform(x1, x2, x3, t_, t0_);
    vector_type f0(3);
    f0(0) = x1 - x_;
    f0(1) = x2 - y_;
    f0(2) = x3 - z_;

    Matrix J(3, 3);
    transformation_->Jacobian(J, x(0), x(1), x(2), t_, t0_);
    matrix_type f1(3, 3);
    for (int icol = 0; icol < 3; ++icol)
    for (int irow = 0; irow < 3; ++irow) {
      f1(irow, icol) = J(irow, icol);
    }

    return std::make_tuple(f0, f1);
  }

  const Transformation * transformation_;
  double x_;
  double y_;
  double z_;
  double t_;
  double t0_;
};

// -----------------------------------------------------------------------------
bool EvaluateInverse(const Transformation *T,
                     double &x, double &y, double &z, double t, double t0)
{
  using namespace boost::math::tools;
  typedef InverseTransformFunctor functor_type;
  typedef typename functor_type::vector_type vector_type;

  functor_type fun(T, x, y, z, t, t0);
  vector_type sol(3); sol(0) = x;  sol(1) = y;  sol(2) = z;
  bool success = multivariate_newton_raphson_iterate<functor_type>(fun, sol, 6);
  x = sol(0); y = sol(1); z = sol(2);
  return success;
}

// =============================================================================
// Inverse multi-level transformation
// =============================================================================

struct InverseMultiLevelTransformFunctor
{
  typedef boost::numeric::ublas::vector<double> vector_type;
  typedef boost::numeric::ublas::matrix<double> matrix_type;
  typedef std::tuple<vector_type, matrix_type>  return_type;

  InverseMultiLevelTransformFunctor(const MultiLevelTransformation * transformation,
    double x, double y, double z, int m, int n, double t, double t0)
    : transformation_(transformation), x_(x), y_(y), z_(z), m_(m), n_(n), t_(t), t0_(t0)
  {
  }

  return_type operator()(vector_type x)
  {
    typedef typename vector_type::value_type value_type;

    value_type x1 = x(0), x2 = x(1), x3 = x(2);
    transformation_->Transform(m_, n_, x1, x2, x3, t_, t0_);
    vector_type f0(3);
    f0(0) = x1 - x_;
    f0(1) = x2 - y_;
    f0(2) = x3 - z_;

    Matrix J(3, 3);
    transformation_->Jacobian(m_, n_, J, x(0), x(1), x(2), t_, t0_);
    matrix_type f1(3, 3);
    for (int icol = 0; icol < 3; ++icol)
    for (int irow = 0; irow < 3; ++irow) {
      f1(irow, icol) = J(irow, icol);
    }

    return std::make_tuple(f0, f1);
  }

  const MultiLevelTransformation * transformation_;
  double x_;
  double y_;
  double z_;
  int m_;
  int n_;
  double t_;
  double t0_;
};

// -----------------------------------------------------------------------------
bool EvaluateInverse(const MultiLevelTransformation *T, int m, int n,
                     double &x, double &y, double &z, double t, double t0)
{
  using namespace boost::math::tools;
  typedef InverseMultiLevelTransformFunctor functor_type;
  typedef typename functor_type::vector_type vector_type;

  functor_type fun(T, x, y, z, m, n, t, t0);
  vector_type sol(3); sol(0) = x;  sol(1) = y;  sol(2) = z;
  bool success = multivariate_newton_raphson_iterate<functor_type>(fun, sol, 6);
  x = sol(0); y = sol(1); z = sol(2);
  return success;
}


} // namespace mirtk

/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#include "NumericsTest.h"

#include "mirtk/Polynomial.h"
using namespace mirtk;

// =============================================================================
// Construction
// =============================================================================

// -----------------------------------------------------------------------------
TEST(Polynomial, Construction)
{
  // p(x) = a1 x^2 + a2 x + a3
  {
    Polynomial polynomial(1, 2);
    #ifndef NDEBUG
      polynomial.Print(1);
    #endif
    ASSERT_EQ(polynomial.Dimension(), 1);
    ASSERT_EQ(polynomial.Order(), 2);
    ASSERT_EQ(polynomial.NumberOfTerms(), 3);
    ASSERT_EQ(polynomial.NumberOfTerms(), polynomial.Coefficients().Rows());
    ASSERT_EQ(polynomial.NumberOfTerms(), polynomial.ModelTerms().Rows());
    ASSERT_EQ(polynomial.Dimension(),     polynomial.ModelTerms().Cols());
    ASSERT_EQ(polynomial.Exponent(0), 2);
    ASSERT_EQ(polynomial.Exponent(1), 1);
    ASSERT_EQ(polynomial.Exponent(2), 0);
    for (int i = 0; i < polynomial.NumberOfTerms(); ++i) {
      ASSERT_EQ(polynomial.Coefficient(i), .0);
    }
  }

  // p(x) = a1 x^2 + a2 x y + a3 x + a4 y^2 + a5 y + a6
  {
    Polynomial polynomial(2, 2);
    #ifndef NDEBUG
      polynomial.Print(1);
    #endif
    ASSERT_EQ(polynomial.Dimension(), 2);
    ASSERT_EQ(polynomial.Order(), 2);
    ASSERT_EQ(polynomial.NumberOfTerms(), 6);
    ASSERT_EQ(polynomial.NumberOfTerms(), polynomial.Coefficients().Rows());
    ASSERT_EQ(polynomial.NumberOfTerms(), polynomial.ModelTerms().Rows());
    ASSERT_EQ(polynomial.Dimension(),     polynomial.ModelTerms().Cols());
    ASSERT_EQ(polynomial.Exponent(0, 0), 2);
    ASSERT_EQ(polynomial.Exponent(0, 1), 0);
    ASSERT_EQ(polynomial.Exponent(1, 0), 1);
    ASSERT_EQ(polynomial.Exponent(1, 1), 1);
    ASSERT_EQ(polynomial.Exponent(2, 0), 1);
    ASSERT_EQ(polynomial.Exponent(2, 1), 0);
    ASSERT_EQ(polynomial.Exponent(3, 0), 0);
    ASSERT_EQ(polynomial.Exponent(3, 1), 2);
    ASSERT_EQ(polynomial.Exponent(4, 0), 0);
    ASSERT_EQ(polynomial.Exponent(4, 1), 1);
    ASSERT_EQ(polynomial.Exponent(5, 0), 0);
    ASSERT_EQ(polynomial.Exponent(5, 1), 0);
    for (int i = 0; i < polynomial.NumberOfTerms(); ++i) {
      ASSERT_EQ(polynomial.Coefficient(i), .0);
    }
  }
}

// -----------------------------------------------------------------------------
TEST(Polynomial, Exponent)
{
  Polynomial polynomial(4, 3);
  #ifndef NDEBUG
    polynomial.Print(1);
  #endif
  for (int i = 0; i < polynomial.NumberOfTerms(); ++i)
  for (int j = 0; j < polynomial.Dimension(); ++j) {
    ASSERT_EQ(polynomial.Exponent(i, j), polynomial.ModelTerms()(i, j));
  }
}

// -----------------------------------------------------------------------------
TEST(Polynomial, IsConstant)
{
  Polynomial polynomial(4, 2);
  #ifndef NDEBUG
    polynomial.Print(1);
  #endif
  for (int i = 0; i < polynomial.NumberOfTerms() - 1; ++i) {
    ASSERT_FALSE(polynomial.IsConstant(i));
  }
  ASSERT_TRUE(polynomial.IsConstant(polynomial.NumberOfTerms() - 1));
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
TEST(Polynomial, Evaluate)
{
  // p(x) = a1 x^2 + a2 x + a3
  {
    Polynomial polynomial(1, 2);
    polynomial.Coefficient(0,  1.0);
    polynomial.Coefficient(1,  0.1);
    polynomial.Coefficient(2, -4.0);
    ASSERT_EQ(polynomial.Evaluate(0.0), -4.0);
    ASSERT_EQ(polynomial.Evaluate(1.0), -2.9);
  }
}

// =============================================================================
// Polynomial regression
// =============================================================================

// -----------------------------------------------------------------------------
TEST(Polynomial, Fit)
{
  // p(x) = a1 x^2 + a2 x + a3
  {
    Polynomial polynomial(1, 2);
    // Use Polynomial::Fit(const Matrix &, const Vector &)
    Matrix x(3, 1);
    Vector y(3);
    x(0, 0) = 0.0, y(0) = -4.0;
    x(1, 0) = 1.0, y(1) = -2.9;
    x(2, 0) = 5.0, y(2) =  8.5;
    const double rms = polynomial.Fit(x, y);
    #ifndef NDEBUG
      cout << "  RMS error = " << rms << endl;
    #endif
    EXPECT_DOUBLE_EQ(rms, .0);
    for (int i = 0; i < x.Rows(); ++i) {
      EXPECT_DOUBLE_EQ(polynomial.Evaluate(x(i, 0)), y(i));
    }
    EXPECT_DOUBLE_EQ(polynomial.Coefficient(2), y(0));
    // Use Polynomial::Fit(const Vector &, const Vector &) instead
    Vector v(3);
    for (int i = 0; i < v.Rows(); ++i) v(i) = x(i, 0);
    Polynomial polynomial2(polynomial.Order());
    EXPECT_DOUBLE_EQ(polynomial2.Fit(v, y), rms);
    for (int i = 0; i < polynomial.NumberOfTerms(); ++i) {
      EXPECT_DOUBLE_EQ(polynomial2.Coefficient(i), polynomial.Coefficient(i));
    }
  }
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

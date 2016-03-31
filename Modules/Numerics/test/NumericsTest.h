/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
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

#ifndef MIRTK_NumericsTest_H
#define MIRTK_NumericsTest_H

#include <cmath>
#include "gtest/gtest.h"


namespace testing { namespace mirtk {


// Helper template function for comparing floating-points.
//
// Template parameter:
//
//   RawType: the raw floating-point type (either float or double)
template <typename RawType>
AssertionResult CmpHelperFloatingPointEQ(const char* expected_expression,
                                         const char* actual_expression,
                                         RawType expected,
                                         RawType actual) {
  if (::std::abs(actual - expected) < static_cast<RawType>(1e-3)) {
    return AssertionSuccess();
  }

  ::std::stringstream expected_ss;
  expected_ss << std::setprecision(std::numeric_limits<RawType>::digits10 + 2)
  << expected;

  ::std::stringstream actual_ss;
  actual_ss << std::setprecision(std::numeric_limits<RawType>::digits10 + 2)
  << actual;

  return testing::internal::EqFailure(expected_expression,
                   actual_expression,
                   testing::internal::StringStreamToString(&expected_ss),
                   testing::internal::StringStreamToString(&actual_ss),
                   false);
}


} } // namespace testing::mirtk


#undef EXPECT_FLOAT_EQ
#define EXPECT_FLOAT_EQ(expected, actual) \
  EXPECT_PRED_FORMAT2(::testing::mirtk::CmpHelperFloatingPointEQ<double>, \
                      expected, actual)

#undef EXPECT_DOUBLE_EQ
#define EXPECT_DOUBLE_EQ(expected, actual) \
  EXPECT_PRED_FORMAT2(::testing::mirtk::CmpHelperFloatingPointEQ<double>, \
                      expected, actual)

#define EXPECT_MATRIX_EQ(expected, actual) \
  EXPECT_EQ(expected.Rows(), actual.Rows()); \
  EXPECT_EQ(expected.Cols(), actual.Cols()); \
  for (int _c = 0; _c < expected.Cols(); ++_c) \
  for (int _r = 0; _r < expected.Rows(); ++_r) \
    EXPECT_DOUBLE_EQ(expected(_r, _c), actual(_r, _c))


#endif // MIRTK_NumericsTest_H

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

#include "FuzzyCorrespondenceUtils.h"

#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"


namespace mirtk { namespace FuzzyCorrespondenceUtils {


// =============================================================================
// Sinkhorn-Knopp algorithm
// =============================================================================

// -----------------------------------------------------------------------------
struct NormalizeRows
{
  WeightMatrix *_Weight;

  void operator()(const blocked_range<int> &rows) const
  {
    double wsum;
    for (int r = rows.begin(); r != rows.end(); ++r) {
      wsum = _Weight->RowSum(r);
      _Weight->ScaleRow(r, 1.0 / wsum);
    }
  }

  static void Run(WeightMatrix &weight, int m = -1)
  {
    if (m < 0) m = weight.Rows();
    NormalizeRows body;
    body._Weight = &weight;
    MIRTK_START_TIMING();
    parallel_for(blocked_range<int>(0, m), body);
    MIRTK_DEBUG_TIMING(7, "row normalization");
  }
};

// -----------------------------------------------------------------------------
struct NormalizeCols
{
  WeightMatrix *_Weight;

  void operator()(const blocked_range<int> &cols) const
  {
    double wsum;
    for (int c = cols.begin(); c != cols.end(); ++c) {
      wsum = _Weight->ColumnSum(c);
      _Weight->ScaleColumn(c, 1.0 / wsum);
    }
  }

  static void Run(WeightMatrix &weight, int n = -1)
  {
    if (n < 0) n = weight.Cols();
    NormalizeCols body;
    body._Weight = &weight;
    MIRTK_START_TIMING();
    parallel_for(blocked_range<int>(0, n), body);
    MIRTK_DEBUG_TIMING(7, "col normalization");
  }
};

// -----------------------------------------------------------------------------
class RowNormalizationError
{
  WeightMatrix *_Weight;
  double        _Error;

public:

  RowNormalizationError(WeightMatrix *weight = NULL)
  :
    _Weight(weight), _Error(.0)
  {}

  RowNormalizationError(const RowNormalizationError &lhs, split)
  :
    _Weight(lhs._Weight), _Error(.0)
  {}

  void join(const RowNormalizationError &rhs)
  {
    _Error += rhs._Error;
  }

  void operator()(const blocked_range<int> &re)
  {
    for (int i = re.begin(); i != re.end(); ++i) {
      _Error += pow(1.0 - _Weight->RowSum(i), 2);
    }
  }

  static double Run(WeightMatrix &weight, int m = -1)
  {
    if (m < 0) m = weight.Rows();
    RowNormalizationError body(&weight);
    MIRTK_START_TIMING();
    parallel_reduce(blocked_range<int>(0, m), body);
    MIRTK_DEBUG_TIMING(7, "evaluating row normalization error");
    return body._Error / m;
  }
};

// -----------------------------------------------------------------------------
class ColumnNormalizationError
{
  WeightMatrix *_Weight;
  double        _Error;

public:

  ColumnNormalizationError(WeightMatrix *weight = NULL)
  :
    _Weight(weight), _Error(.0)
  {}

  ColumnNormalizationError(const ColumnNormalizationError &lhs, split)
  :
    _Weight(lhs._Weight), _Error(.0)
  {}

  void join(const ColumnNormalizationError &rhs)
  {
    _Error += rhs._Error;
  }

  void operator()(const blocked_range<int> &re)
  {
    for (int i = re.begin(); i != re.end(); ++i) {
      _Error += pow(1.0 - _Weight->ColumnSum(i), 2);
    }
  }

  static double Run(WeightMatrix &weight, int n = -1)
  {
    if (n < 0) n = weight.Cols();
    ColumnNormalizationError body(&weight);
    MIRTK_START_TIMING();
    parallel_reduce(blocked_range<int>(0, n), body);
    MIRTK_DEBUG_TIMING(7, "evaluating col normalization error");
    return body._Error / n;
  }
};

// -----------------------------------------------------------------------------
void NormalizeWeights(WeightMatrix &weight, int m, int n, int maxit, double maxerr)
{
  // Precompute column/row indices
  weight.Index();
  // Iterative row and column normalization
  if (maxerr <= .0) {
    for (int i = 0; i < maxit; ++i) {
      NormalizeCols::Run(weight, n);
      NormalizeRows::Run(weight, m);
    }
  } else {
    if (weight.Layout() == WeightMatrix::CRS) {
      for (int i = 0; i < maxit; ++i) {
        NormalizeRows::Run(weight, m);
        NormalizeCols::Run(weight, n);
        if (RowNormalizationError::Run(weight, m) < maxerr) break;
      }
    } else {
      for (int i = 0; i < maxit; ++i) {
        NormalizeCols::Run(weight, n);
        NormalizeRows::Run(weight, m);
        if (ColumnNormalizationError::Run(weight, n) < maxerr) break;
      }
    }
  }
}


} } // namespace mirtk::FuzzyCorrespondenceUtils

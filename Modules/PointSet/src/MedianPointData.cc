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

#include "mirtk/MedianPointData.h"

#include "mirtk/Algorithm.h"
#include "mirtk/Parallel.h"


namespace mirtk {


// =============================================================================
// Auxiliaries
// =============================================================================

namespace MedianPointDataUtils {


// -----------------------------------------------------------------------------
/// Median filter given point data components
struct MedianFilter
{
  vtkDataArray     *_Input;
  vtkDataArray     *_Output;
  EdgeConnectivity *_Neighbors;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    Array<double> values;
    int median, nbrPts;
    const int  *nbrIds;

    for (auto ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      _Neighbors->GetConnectedPoints(ptId, nbrPts, nbrIds);
      if (nbrPts > 0) {
        median = nbrPts / 2;
        values.resize(nbrPts + 1);
        for (int j = 0; j < _Input->GetNumberOfComponents(); ++j) {
          values[nbrPts] = _Input->GetComponent(ptId, j);
          for (int i = 0; i < nbrPts; ++i) {
            values[i] = _Input->GetComponent(nbrIds[i], j);
          }
          _Output->SetComponent(ptId, j, NthElement(values, median));
        }
      } else {
        for (int j = 0; j < _Input->GetNumberOfComponents(); ++j) {
          _Output->SetComponent(ptId, j, _Input->GetComponent(ptId, j));
        }
      }
    }
  }
};


} // namespace MedianPointDataUtils
using namespace MedianPointDataUtils;

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void MedianPointData::CopyAttributes(const MedianPointData &)
{
}

// -----------------------------------------------------------------------------
MedianPointData::MedianPointData()
{
}

// -----------------------------------------------------------------------------
MedianPointData::MedianPointData(const MedianPointData &other)
:
  PointDataFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
MedianPointData &MedianPointData::operator =(const MedianPointData &other)
{
  if (this != &other) {
    PointDataFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
MedianPointData::~MedianPointData()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void MedianPointData::Execute()
{
  MedianFilter filter;
  filter._Input     = _InputData;
  filter._Output    = _OutputData;
  filter._Neighbors = _Neighbors.get();
  parallel_for(blocked_range<int>(0, static_cast<int>(_Input->GetNumberOfPoints())), filter);
}


} // namespace mirtk

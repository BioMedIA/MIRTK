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

#include "mirtk/DilatePointData.h"

#include "mirtk/Parallel.h"


namespace mirtk {


// =============================================================================
// Auxiliaries
// =============================================================================

namespace DilatePointDataUtils {


// -----------------------------------------------------------------------------
struct DilateScalars
{
  vtkDataArray           *_Input;
  vtkDataArray           *_Output;
  const EdgeConnectivity *_Neighbors;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    int        nbrPts;
    const int *nbrIds;
    double     value;

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      for (int j = 0; j < _Input->GetNumberOfComponents(); ++j) {
        value = _Input->GetComponent(ptId, j);
        _Neighbors->GetConnectedPoints(ptId, nbrPts, nbrIds);
        for (int i = 0; i < nbrPts; ++i) {
          value = max(value, _Input->GetComponent(nbrIds[i], j));
        }
        _Output->SetComponent(ptId, j, value);
      }
    }
  }
};


} // namespace DilatePointDataUtils
using namespace DilatePointDataUtils;

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void DilatePointData::CopyAttributes(const DilatePointData &other)
{
  _Iterations = other._Iterations;
}

// -----------------------------------------------------------------------------
DilatePointData::DilatePointData()
:
  _Iterations(1)
{
}

// -----------------------------------------------------------------------------
DilatePointData::DilatePointData(const DilatePointData &other)
:
  PointDataFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
DilatePointData &DilatePointData::operator =(const DilatePointData &other)
{
  if (this != &other) {
    PointDataFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
DilatePointData::~DilatePointData()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void DilatePointData::Execute()
{
  vtkSmartPointer<vtkDataArray> arr = _InputData;
  vtkSmartPointer<vtkDataArray> res = _OutputData;
  DilateScalars body;
  body._Input     = arr;
  body._Output    = res;
  body._Neighbors = _Neighbors.get();
  for (int iter = 0; iter < _Iterations; ++iter) {
    if (iter == 1) {
      arr.TakeReference(res->NewInstance());
      arr->SetNumberOfComponents(res->GetNumberOfComponents());
      arr->SetNumberOfTuples(res->GetNumberOfTuples());
      arr->SetName(res->GetName());
      for (int j = 0; j < res->GetNumberOfComponents(); ++j) {
        arr->SetComponentName(j, res->GetComponentName(j));
        arr->CopyComponent(j, res, j);
      }
      body._Input = arr;
    } else if (iter > 1) {
      swap(body._Input, body._Output);
    }
    parallel_for(blocked_range<int>(0, static_cast<int>(_Input->GetNumberOfPoints())), body);
  }
  if (body._Output != _OutputData) {
    for (int j = 0; j < body._Output->GetNumberOfComponents(); ++j) {
      _OutputData->CopyComponent(j, body._Output, j);
    }
  }
}


} // namespace mirtk

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

#include "mirtk/DilateCellData.h"

#include "mirtk/Parallel.h"


namespace mirtk {


// =============================================================================
// Auxiliaries
// =============================================================================

namespace DilateCellDataUtils {


// -----------------------------------------------------------------------------
struct DilateScalars
{
  CellDataFilter *_Filter;
  vtkDataArray   *_Input;
  vtkDataArray   *_Output;

  void operator ()(const blocked_range<int> &cellIds) const
  {
    UnorderedSet<int> nbrIds;
    double            value;

    for (int cellId = cellIds.begin(); cellId != cellIds.end(); ++cellId) {
      for (int j = 0; j < _Input->GetNumberOfComponents(); ++j) {
        value = _Input->GetComponent(cellId, j);
        _Filter->GetNodeNeighbors(cellId, nbrIds);
        for (auto nbrId : nbrIds) {
          value = max(value, _Input->GetComponent(nbrId, j));
        }
        _Output->SetComponent(cellId, j, value);
      }
    }
  }
};


} // namespace DilateCellDataUtils
using namespace DilateCellDataUtils;

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void DilateCellData::CopyAttributes(const DilateCellData &other)
{
  _Iterations = other._Iterations;
}

// -----------------------------------------------------------------------------
DilateCellData::DilateCellData()
:
  _Iterations(1)
{
}

// -----------------------------------------------------------------------------
DilateCellData::DilateCellData(const DilateCellData &other)
:
  CellDataFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
DilateCellData &DilateCellData::operator =(const DilateCellData &other)
{
  if (this != &other) {
    CellDataFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
DilateCellData::~DilateCellData()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void DilateCellData::Execute()
{
  vtkSmartPointer<vtkDataArray> arr = _InputData;
  vtkSmartPointer<vtkDataArray> res = _OutputData;
  DilateScalars body;
  body._Filter = this;
  body._Input  = arr;
  body._Output = res;
  for (int iter = 0; iter < _Iterations; ++iter) {
    if (iter == 1) {
      arr.TakeReference(res->NewInstance());
      arr->DeepCopy(res);
      arr->SetName(res->GetName());
      arr->CopyComponentNames(res);
      body._Input = arr;
    } else if (iter > 1) {
      swap(body._Input, body._Output);
    }
    parallel_for(blocked_range<int>(0, static_cast<int>(_Input->GetNumberOfCells())), body);
  }
  if (body._Output != _OutputData) {
    _OutputData->DeepCopy(body._Output);
  }
}


} // namespace mirtk

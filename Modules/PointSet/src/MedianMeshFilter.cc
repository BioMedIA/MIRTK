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

#include "mirtk/MedianMeshFilter.h"

#include "mirtk/Algorithm.h"
#include "mirtk/EdgeConnectivity.h"
#include "mirtk/Parallel.h"

#include "vtkPointData.h"
#include "vtkDataArray.h"


namespace mirtk {


// =============================================================================
// Auxiliaries
// =============================================================================

namespace MedianMeshFilterUtils {


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


} // namespace MedianMeshFilterUtils
using namespace MedianMeshFilterUtils;

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void MedianMeshFilter::CopyAttributes(const MedianMeshFilter &other)
{
  _Connectivity  = other._Connectivity;
  _DataArray     = other._DataArray;
  _DataArrayName = other._DataArrayName;
}

// -----------------------------------------------------------------------------
MedianMeshFilter::MedianMeshFilter()
:
  _Connectivity(1)
{
}

// -----------------------------------------------------------------------------
MedianMeshFilter::MedianMeshFilter(const MedianMeshFilter &other)
:
  MeshFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
MedianMeshFilter &MedianMeshFilter::operator =(const MedianMeshFilter &other)
{
  if (this != &other) {
    MeshFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
MedianMeshFilter::~MedianMeshFilter()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void MedianMeshFilter::Initialize()
{
  // Initialize base class
  MeshFilter::Initialize();

  // Get input point data array
  if (_DataArray) {
    if (_DataArrayName.empty()) {
      if (_DataArray->GetName() == nullptr) {
        Throw(ERR_RuntimeError, __FUNCTION__, "Input data array must have a name or specify output data array name");
      }
      _DataArrayName = _DataArray->GetName();
    }
  } else {
    if (_DataArrayName.empty()) {
      Throw(ERR_RuntimeError, __FUNCTION__, "No point data array named or set to process");
    }
    _DataArray = _Input ->GetPointData()->GetArray(_DataArrayName.c_str());
    if (_DataArray == nullptr) {
      Throw(ERR_RuntimeError, __FUNCTION__, "Input mesh has no point data array named ", _DataArrayName);
    }
  }

  // Add new output point data array
  vtkPointData * const outputPD = _Output->GetPointData();
  vtkSmartPointer<vtkDataArray> result;
  result = NewPointArray(_DataArrayName.c_str(), _DataArray->GetNumberOfComponents(), _DataArray->GetDataType());
  outputPD->RemoveArray(_DataArrayName.c_str());
  outputPD->AddArray(result);
}

// -----------------------------------------------------------------------------
void MedianMeshFilter::Execute()
{
  MedianFilter filter;
  EdgeConnectivity neighbors;
  neighbors.Initialize(_Input, _Connectivity, _EdgeTable.get());
  filter._Input     = _DataArray;
  filter._Output    = _Output->GetPointData()->GetArray(_DataArrayName.c_str());
  filter._Neighbors = &neighbors;
  parallel_for(blocked_range<int>(0, static_cast<int>(_Input->GetNumberOfPoints())), filter);
}


} // namespace mirtk

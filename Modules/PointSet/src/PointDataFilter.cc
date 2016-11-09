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

#include "mirtk/PointDataFilter.h"

#include "mirtk/PointSetUtils.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void PointDataFilter::CopyAttributes(const PointDataFilter &other)
{
  _Connectivity = other._Connectivity;
  _Radius       = other._Radius;
  _Neighbors    = other._Neighbors;
  _DataName     = other._DataName;
  _InputData    = other._InputData;
  _OutputData   = _Output->GetPointData()->GetArray(_DataName.c_str());
}

// -----------------------------------------------------------------------------
PointDataFilter::PointDataFilter()
:
  _Connectivity(1),
  _Radius(0.)
{
}

// -----------------------------------------------------------------------------
PointDataFilter::PointDataFilter(const PointDataFilter &other)
:
  MeshFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
PointDataFilter &PointDataFilter::operator =(const PointDataFilter &other)
{
  if (this != &other) {
    MeshFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
PointDataFilter::~PointDataFilter()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void PointDataFilter::Initialize()
{
  // Initialize base class
  MeshFilter::Initialize();

  // Get input point data array
  if (_InputData) {
    if (_DataName.empty()) {
      if (_InputData->GetName() == nullptr) {
        Throw(ERR_RuntimeError, __FUNCTION__, "Input data array must have a name or specify output data array name");
      }
      _DataName = _InputData->GetName();
    }
  } else if (_DataName.empty()) {
    _InputData = _Input->GetPointData()->GetScalars();
    Throw(ERR_RuntimeError, __FUNCTION__, "Input mesh has no active scalars point data array, data array name required");
  } else {
    _InputData = _Input->GetPointData()->GetArray(_DataName.c_str());
    if (_InputData == nullptr) {
      _InputData = GetArrayByCaseInsensitiveName(_Input->GetPointData(), _DataName.c_str());
      if (_InputData == nullptr) {
        Throw(ERR_RuntimeError, __FUNCTION__, "Input mesh has no point data array named ", _DataName);
      }
    }
  }

  // Add new output point data array
  int attr = -1;
  vtkPointData * const outputPD = _Output->GetPointData();
  for (int i = 0; i < outputPD->GetNumberOfArrays(); ++i) {
    if (outputPD->GetArray(i) == _InputData) {
      attr = outputPD->IsArrayAnAttribute(i);
      break;
    }
  }
  _OutputData = NewPointArray(_DataName.c_str(), _InputData->GetNumberOfComponents(), _InputData->GetDataType());
  for (int j = 0; j < _InputData->GetNumberOfComponents(); ++j) {
    _OutputData->SetComponentName(j, _InputData->GetComponentName(j));
  }
  outputPD->RemoveArray(_DataName.c_str());
  const int idx = outputPD->AddArray(_OutputData);
  if (attr >= 0) outputPD->SetActiveAttribute(idx, attr);

  // Pre-compute sets of neighboring nodes
  if (!_Neighbors) {
    InitializeEdgeTable();
    if (_Radius > 0.) {
      _Neighbors = NewShared<EdgeConnectivity>(_Input, _Radius, _EdgeTable.get());
    } else {
      _Neighbors = NewShared<EdgeConnectivity>(_Input, _Connectivity, _EdgeTable.get());
    }
  }
}


} // namespace mirtk

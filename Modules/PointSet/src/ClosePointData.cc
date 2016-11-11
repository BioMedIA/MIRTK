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

#include "mirtk/ClosePointData.h"

#include "mirtk/DilatePointData.h"
#include "mirtk/ErodePointData.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void ClosePointData::CopyAttributes(const ClosePointData &other)
{
  _Iterations = other._Iterations;
}

// -----------------------------------------------------------------------------
ClosePointData::ClosePointData()
:
  _Iterations(1)
{
}

// -----------------------------------------------------------------------------
ClosePointData::ClosePointData(const ClosePointData &other)
:
  PointDataFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
ClosePointData &ClosePointData::operator =(const ClosePointData &other)
{
  if (this != &other) {
    PointDataFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
ClosePointData::~ClosePointData()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void ClosePointData::Initialize()
{
  // Nothing to be done
}

// -----------------------------------------------------------------------------
void ClosePointData::Execute()
{
  DilatePointData dilate;
  dilate.Input(_Input);
  dilate.InputData(_InputData);
  dilate.DataName(_DataName);
  dilate.Connectivity(_Connectivity);
  dilate.Radius(_Radius);
  dilate.EdgeTable(_EdgeTable);
  dilate.Neighbors(_Neighbors);
  dilate.Iterations(_Iterations);
  dilate.Run();

  _EdgeTable = dilate.EdgeTable();
  _Neighbors = dilate.Neighbors();

  ErodePointData erode;
  erode.Input(dilate.Output());
  erode.InputData(dilate.OutputData());
  erode.DataName(_DataName);
  erode.EdgeTable(_EdgeTable);
  erode.Neighbors(_Neighbors);
  erode.Iterations(_Iterations);
  erode.Run();

  _Output     = erode.Output();
  _OutputData = erode.OutputData();
}


} // namespace mirtk

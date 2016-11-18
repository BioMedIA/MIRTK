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

#include "mirtk/OpenPointData.h"

#include "mirtk/DilatePointData.h"
#include "mirtk/ErodePointData.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void OpenPointData::CopyAttributes(const OpenPointData &other)
{
  _Iterations = other._Iterations;
}

// -----------------------------------------------------------------------------
OpenPointData::OpenPointData()
:
  _Iterations(1)
{
}

// -----------------------------------------------------------------------------
OpenPointData::OpenPointData(const OpenPointData &other)
:
  PointDataFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
OpenPointData &OpenPointData::operator =(const OpenPointData &other)
{
  if (this != &other) {
    PointDataFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
OpenPointData::~OpenPointData()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void OpenPointData::Initialize()
{
  // Nothing to be done
}

// -----------------------------------------------------------------------------
void OpenPointData::Execute()
{
  ErodePointData erode;
  erode.Input(_Input);
  erode.InputData(_InputData);
  erode.DataName(_DataName);
  erode.Connectivity(_Connectivity);
  erode.Radius(_Radius);
  erode.EdgeTable(_EdgeTable);
  erode.Neighbors(_Neighbors);
  erode.Iterations(_Iterations);
  erode.Run();

  _EdgeTable = erode.EdgeTable();
  _Neighbors = erode.Neighbors();

  DilatePointData dilate;
  dilate.Input(erode.Output());
  dilate.InputData(erode.OutputData());
  dilate.EdgeTable(_EdgeTable);
  dilate.Neighbors(_Neighbors);
  dilate.Iterations(_Iterations);
  dilate.Run();

  _Output     = dilate.Output();
  _OutputData = dilate.OutputData();
}


} // namespace mirtk

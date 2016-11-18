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

#include "mirtk/OpenCellData.h"

#include "mirtk/DilateCellData.h"
#include "mirtk/ErodeCellData.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void OpenCellData::CopyAttributes(const OpenCellData &other)
{
  _Iterations = other._Iterations;
}

// -----------------------------------------------------------------------------
OpenCellData::OpenCellData()
:
  _Iterations(1)
{
}

// -----------------------------------------------------------------------------
OpenCellData::OpenCellData(const OpenCellData &other)
:
  CellDataFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
OpenCellData &OpenCellData::operator =(const OpenCellData &other)
{
  if (this != &other) {
    CellDataFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
OpenCellData::~OpenCellData()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void OpenCellData::Initialize()
{
  // Nothing to be done
}

// -----------------------------------------------------------------------------
void OpenCellData::Execute()
{
  ErodeCellData erode;
  erode.Input(_Input);
  erode.InputData(_InputData);
  erode.DataName(_DataName);
  erode.Iterations(_Iterations);
  erode.Run();

  DilateCellData dilate;
  dilate.Input(erode.Output());
  dilate.InputData(erode.OutputData());
  dilate.Iterations(_Iterations);
  dilate.Run();

  _Output     = dilate.Output();
  _OutputData = dilate.OutputData();
}


} // namespace mirtk

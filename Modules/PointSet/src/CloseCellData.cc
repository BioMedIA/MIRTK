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

#include "mirtk/CloseCellData.h"

#include "mirtk/DilateCellData.h"
#include "mirtk/ErodeCellData.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void CloseCellData::CopyAttributes(const CloseCellData &other)
{
  _Iterations = other._Iterations;
}

// -----------------------------------------------------------------------------
CloseCellData::CloseCellData()
:
  _Iterations(1)
{
}

// -----------------------------------------------------------------------------
CloseCellData::CloseCellData(const CloseCellData &other)
:
  CellDataFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
CloseCellData &CloseCellData::operator =(const CloseCellData &other)
{
  if (this != &other) {
    CellDataFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
CloseCellData::~CloseCellData()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void CloseCellData::Initialize()
{
  // Nothing to be done
}

// -----------------------------------------------------------------------------
void CloseCellData::Execute()
{
  DilateCellData dilate;
  dilate.Input(_Input);
  dilate.InputData(_InputData);
  dilate.DataName(_DataName);
  dilate.Iterations(_Iterations);
  dilate.Run();

  ErodeCellData erode;
  erode.Input(dilate.Output());
  erode.InputData(dilate.OutputData());
  erode.Iterations(_Iterations);
  erode.Run();

  _Output     = erode.Output();
  _OutputData = erode.OutputData();
}


} // namespace mirtk

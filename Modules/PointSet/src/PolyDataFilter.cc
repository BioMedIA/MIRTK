/*
 * Medical Image Registration ToolKit (MMIRTK)
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

#include "mirtk/PolyDataFilter.h"

#include "mirtk/Profiling.h"
#include "mirtk/EdgeTable.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
PolyDataFilter::PolyDataFilter()
:
  _EdgeTable(NULL),
  _EdgeTableOwner(false),
  _DoublePrecision(false)
{
}

// -----------------------------------------------------------------------------
void PolyDataFilter::CopyAttributes(const PolyDataFilter &other)
{
  _Input = other._Input;

  if (_EdgeTableOwner) delete _EdgeTable;
  _EdgeTable      = NULL;
  _EdgeTableOwner = false;
  if (other._EdgeTable) {
    _EdgeTableOwner = other._EdgeTableOwner;
    _EdgeTable      = (_EdgeTableOwner ? new class EdgeTable(*other._EdgeTable) : other._EdgeTable);
  }

  if (other._Output) {
    _Output = vtkSmartPointer<vtkPolyData>::New();
    _Output->DeepCopy(other._Output);
  } else {
    _Output = NULL;
  }

  _DoublePrecision = other._DoublePrecision;
}

// -----------------------------------------------------------------------------
PolyDataFilter::PolyDataFilter(const PolyDataFilter &other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
PolyDataFilter &PolyDataFilter::operator =(const PolyDataFilter &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
PolyDataFilter::~PolyDataFilter()
{
  if (_EdgeTableOwner) Delete(_EdgeTable);
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void PolyDataFilter::Run()
{
  MIRTK_START_TIMING();
  {
    MIRTK_START_TIMING();
    this->Initialize();
    MIRTK_DEBUG_TIMING(2, this->NameOfClass() << "::Initialize");
  }
  {
    MIRTK_START_TIMING();
    this->Execute();
    MIRTK_DEBUG_TIMING(2, this->NameOfClass() << "::Execute");
  }
  {
    MIRTK_START_TIMING();
    this->Finalize();
    MIRTK_DEBUG_TIMING(2, this->NameOfClass() << "::Finalize");
  }
  MIRTK_DEBUG_TIMING(1, this->NameOfClass());
}

// -----------------------------------------------------------------------------
void PolyDataFilter::InitializeEdgeTable()
{
  if (_EdgeTable == NULL) {
    _EdgeTable      = new class EdgeTable(_Input);
    _EdgeTableOwner = true;
  }
}

// -----------------------------------------------------------------------------
void PolyDataFilter::Initialize()
{
  // Check input
  if (!_Input) {
    cerr << this->NameOfClass() << "::Initialize: Input surface mesh not set!" << endl;
    exit(1);
  }

  // By default, set output to be shallow copy of input
  _Output = vtkSmartPointer<vtkPolyData>::New();
  _Output->ShallowCopy(_Input);
}

// -----------------------------------------------------------------------------
void PolyDataFilter::Finalize()
{
  // Destroy edge table
  if (_EdgeTableOwner) {
    delete _EdgeTable;
    _EdgeTable      = NULL;
    _EdgeTableOwner = false;
  }
}


} // namespace mirtk

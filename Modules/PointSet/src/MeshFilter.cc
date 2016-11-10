/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#include "mirtk/MeshFilter.h"

#include "mirtk/Config.h" // MIRTK_USE_FLOAT_BY_DEFAULT
#include "mirtk/Profiling.h"

#include "mirtk/Vtk.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void MeshFilter::CopyAttributes(const MeshFilter &other)
{
  _Input           = other._Input;
  _EdgeTable       = other._EdgeTable;
  _DoublePrecision = other._DoublePrecision;

  if (other._Output) {
    _Output.TakeReference(other._Output->NewInstance());
    _Output->ShallowCopy(other._Output);
  } else {
    _Output = nullptr;
  }
}

// -----------------------------------------------------------------------------
MeshFilter::MeshFilter()
:
#if MIRTK_USE_FLOAT_BY_DEFAULT
  _DoublePrecision(false)
#else
  _DoublePrecision(true)
#endif
{
}

// -----------------------------------------------------------------------------
MeshFilter::MeshFilter(const MeshFilter &other)
:
  Object(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
MeshFilter &MeshFilter::operator =(const MeshFilter &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
MeshFilter::~MeshFilter()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void MeshFilter::Run()
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
void MeshFilter::Initialize()
{
  // Check input
  if (_Input == nullptr) {
    Throw(ERR_LogicError, __FUNCTION__, "Input mesh not set!");
  }

  // Build mesh links
  _Input->BuildLinks();

  // By default, set output to be shallow copy of input
  _Output.TakeReference(_Input->NewInstance());
  _Output->ShallowCopy(_Input);

  // Ensure we don't modify shallow copy of input cell types when
  // marking cells as deleted... also need to build cells/links
  _Output->DeleteCells();
  _Output->BuildLinks();
}

// -----------------------------------------------------------------------------
void MeshFilter::Finalize()
{
}

// =============================================================================
// Auxiliaries
// =============================================================================

// ------------------------------------------------------------------------------
vtkSmartPointer<vtkDataArray>
MeshFilter::NewArray(const char *name, vtkIdType n, int c, int type) const
{
  if (type == VTK_VOID) type = (_DoublePrecision ? VTK_DOUBLE : VTK_FLOAT);
  vtkSmartPointer<vtkDataArray> array = NewVTKDataArray(type);
  array->SetName(name);
  array->SetNumberOfComponents(c);
  array->SetNumberOfTuples(n);
  return array;
}

// ------------------------------------------------------------------------------
vtkSmartPointer<vtkDataArray>
MeshFilter::NewArray(const char *name, int c, int type) const
{
  return NewArray(name, _Input->GetNumberOfPoints(), c, type);
}

// ------------------------------------------------------------------------------
vtkSmartPointer<vtkDataArray>
MeshFilter::NewPointArray(const char *name, int c, int type) const
{
  return NewArray(name, _Input->GetNumberOfPoints(), c, type);
}

// ------------------------------------------------------------------------------
vtkSmartPointer<vtkDataArray>
MeshFilter::NewCellArray(const char *name, int c, int type) const
{
  return NewArray(name, _Input->GetNumberOfCells(), c, type);
}

// -----------------------------------------------------------------------------
void MeshFilter::InitializeEdgeTable()
{
  if (_EdgeTable == nullptr || static_cast<vtkIdType>(_EdgeTable->NumberOfPoints()) != _Input->GetNumberOfPoints()) {
    _EdgeTable = NewShared<class EdgeTable>(_Input);
  }
}


} // namespace mirtk

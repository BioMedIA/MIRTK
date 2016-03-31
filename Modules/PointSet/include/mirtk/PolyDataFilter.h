/*
 * Medical Image Registration ToolKit (MIRTK)
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

#ifndef MIRTK_PolyDataFilter_H
#define MIRTK_PolyDataFilter_H

#include <mirtkObject.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>


namespace mirtk {


class EdgeTable;


/**
 * Base class for filters which process polygonal surface meshes
 */
class PolyDataFilter : public Object
{
  mirtkAbstractMacro(PolyDataFilter);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Input surface mesh
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Input);

  /// Precomputed edge table of input surface mesh
  mirtkPublicAggregateMacro(const class EdgeTable, EdgeTable);

  /// Whether edge table is to be destroyed by this filter
  mirtkPublicAttributeMacro(bool, EdgeTableOwner);

  /// Output surface mesh (NULL if filter does not produce polygonal output)
  mirtkReadOnlyAttributeMacro(vtkSmartPointer<vtkPolyData>, Output);

  /// Whether to output floating point results in double precision
  /// If \c true, output data arrays are of type vtkDoubleArray.
  /// If \c false, output data arrays are of type vtkFloatArray.
  mirtkPublicAttributeMacro(bool, DoublePrecision);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Copy attributes of this filter from another instance
  void CopyAttributes(const PolyDataFilter &);

protected:

  /// Default constructor
  PolyDataFilter();

  /// Copy constructor
  PolyDataFilter(const PolyDataFilter &);

  /// Assignment operator
  PolyDataFilter &operator =(const PolyDataFilter &);

  /// Destructor
  virtual ~PolyDataFilter();

  // ---------------------------------------------------------------------------
  // Execution

public:

  /// Run filter
  ///
  /// \note Most filters assume that irtkPolyData::BuildLinks of the input
  ///       surface mesh has been invoked before execution of the filter.
  ///       See the documentation/implementation of each individual filter.
  ///       Generally assume that this is required unless documented otherwise.
  virtual void Run();

protected:

  /// To be called by subclass in Initialize when _EdgeTable is needed
  void InitializeEdgeTable();

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Execute filter
  virtual void Execute() = 0;

  /// Finalize filter execution
  virtual void Finalize();

  // ---------------------------------------------------------------------------
  // Alternative VTK-like API

public:

  /// Enable/disable double precision output
  mirtkOnOffMacro(DoublePrecision);

  /// Set input surface mesh
  void SetInputData(vtkPolyData *);

  /// Set input surface mesh
  void SetInput(vtkPolyData *);

  /// Run filter
  void Update();

  /// Get output surface mesh
  vtkPolyData *GetOutput();

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Alternative VTK-like API
// =============================================================================

// -----------------------------------------------------------------------------
inline void PolyDataFilter::SetInputData(vtkPolyData *poly)
{
  this->Input(poly);
}

// -----------------------------------------------------------------------------
inline void PolyDataFilter::SetInput(vtkPolyData *poly)
{
  this->Input(poly);
}

// -----------------------------------------------------------------------------
inline void PolyDataFilter::Update()
{
  this->Run();
}

// -----------------------------------------------------------------------------
inline vtkPolyData *PolyDataFilter::GetOutput()
{
  return this->Output();
}


} // namespace mirtk

#endif // MIRTK_PolyDataFilter

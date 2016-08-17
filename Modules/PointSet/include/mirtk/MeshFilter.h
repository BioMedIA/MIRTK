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

#ifndef MIRTK_MeshFilter_H
#define MIRTK_MeshFilter_H

#include "mirtk/Object.h"

#include "mirtk/EdgeTable.h"
#include "mirtk/Memory.h"

#include "vtkType.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"


namespace mirtk {


/**
 * Base class for filters which process discrete surface and/or volumetric meshes
 */
class MeshFilter : public Object
{
  mirtkAbstractMacro(MeshFilter);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Input mesh
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Input);

  /// Precomputed edge table of input surface mesh
  mirtkPublicAttributeMacro(SharedPtr<const class EdgeTable>, EdgeTable);

  /// Output surface mesh (NULL if filter does not produce polygonal output)
  mirtkReadOnlyAttributeMacro(vtkSmartPointer<vtkPolyData>, Output);

  /// Whether to output floating point results in double precision
  /// If \c true,  output data arrays are of type vtkDoubleArray.
  /// If \c false, output data arrays are of type vtkFloatArray.
  mirtkPublicAttributeMacro(bool, DoublePrecision);

  /// Copy attributes of this filter from another instance
  void CopyAttributes(const MeshFilter &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  MeshFilter();

  /// Copy constructor
  MeshFilter(const MeshFilter &);

  /// Assignment operator
  MeshFilter &operator =(const MeshFilter &);

  /// Destructor
  virtual ~MeshFilter();

  // ---------------------------------------------------------------------------
  // Execution

public:

  /// Run filter
  virtual void Run();

protected:

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

  // ---------------------------------------------------------------------------
  // Auxiliaries

protected:

  /// To be called by subclass in Initialize when _EdgeTable is needed
  virtual void InitializeEdgeTable();

  /// Allocate new data set attributes array
  ///
  /// \param[in] name Name of data array.
  /// \param[in] n    Number of tuples.
  /// \param[in] c    Number of components.
  /// \param[in] type Type of VTK array. When VTK_VOID, a floating point array
  ///                 is allocated with either single or double precision depending
  ///                 on the DoublePrecision flag of this mesh filter.
  ///
  /// \returns New floating point array
  vtkSmartPointer<vtkDataArray> NewArray(const char *name, vtkIdType n, int c,
                                         int type = VTK_VOID) const;

  /// Allocate new point data array
  ///
  /// \param[in] name Name of data array.
  /// \param[in] c    Number of components.
  /// \param[in] type Type of VTK array. When VTK_VOID, a floating point array
  ///                 is allocated with either single or double precision depending
  ///                 on the DoublePrecision flag of this mesh filter.
  ///
  /// \returns New floating point array
  ///
  /// \deprecated Use NewPointArray instead.
  vtkSmartPointer<vtkDataArray> NewArray(const char *name, int c = 1,
                                         int type = VTK_VOID) const;

  /// Allocate new point data array
  ///
  /// \param[in] name Name of data array.
  /// \param[in] c    Number of components.
  /// \param[in] type Type of VTK array. When VTK_VOID, a floating point array
  ///                 is allocated with either single or double precision depending
  ///                 on the DoublePrecision flag of this mesh filter.
  ///
  /// \returns New floating point array
  vtkSmartPointer<vtkDataArray> NewPointArray(const char *name, int c = 1,
                                              int type = VTK_VOID) const;

  /// Allocate new cell data array
  ///
  /// \param[in] name Name of data array.
  /// \param[in] c    Number of components.
  /// \param[in] type Type of VTK array. When VTK_VOID, a floating point array
  ///                 is allocated with either single or double precision depending
  ///                 on the DoublePrecision flag of this mesh filter.
  ///
  /// \returns New floating point array
  vtkSmartPointer<vtkDataArray> NewCellArray(const char *name, int c = 1,
                                             int type = VTK_VOID) const;
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Alternative VTK-like API
// =============================================================================

// -----------------------------------------------------------------------------
inline void MeshFilter::SetInputData(vtkPolyData *poly)
{
  this->Input(poly);
}

// -----------------------------------------------------------------------------
inline void MeshFilter::SetInput(vtkPolyData *poly)
{
  this->Input(poly);
}

// -----------------------------------------------------------------------------
inline void MeshFilter::Update()
{
  this->Run();
}

// -----------------------------------------------------------------------------
inline vtkPolyData *MeshFilter::GetOutput()
{
  return this->Output();
}


} // namespace mirtk

#endif // MIRTK_MeshFilter_H

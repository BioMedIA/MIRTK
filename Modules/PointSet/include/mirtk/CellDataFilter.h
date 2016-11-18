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

#ifndef MIRTK_CellDataFilter_H
#define MIRTK_CellDataFilter_H

#include "mirtk/MeshFilter.h"

#include "mirtk/Memory.h"
#include "mirtk/SparseMatrix.h"

#include "vtkSmartPointer.h"
#include "vtkDataArray.h"
#include "vtkCellData.h"


namespace mirtk {


/**
 * Filter for mesh cell data
 */
class CellDataFilter : public MeshFilter
{
  mirtkAbstractMacro(CellDataFilter);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Name of (input and) output cell data array
  ///
  /// When an input _DataArray is given, this name is assigned to the respective
  /// output cell data array. Otherwise, the name of the input data array must
  /// be set and is used also for the output data array. When no input _DataArray
  /// is given, the name of the cell data array of the input mesh must be set.
  mirtkPublicAttributeMacro(string, DataName);

  /// Input cell data array
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, InputData);

  /// Output cell data array
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, OutputData);

  /// Copy attributes of this filter from another instance
  void CopyAttributes(const CellDataFilter &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  CellDataFilter();

  /// Copy constructor
  CellDataFilter(const CellDataFilter &);

  /// Assignment operator
  CellDataFilter &operator =(const CellDataFilter &);

  /// Destructor
  virtual ~CellDataFilter();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

public:

  /// Get cell neighbors sharing a node with the specified cell
  void GetNodeNeighbors(int, UnorderedSet<int> &) const;

  /// Get cell neighbors sharing an edge with the specified cell
  void GetEdgeNeighbors(int, UnorderedSet<int> &) const;

};


} // namespace mirtk

#endif // MIRTK_CellDataFilter_H

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

#ifndef MIRTK_PointDataFilter_H
#define MIRTK_PointDataFilter_H

#include "mirtk/MeshFilter.h"

#include "mirtk/Memory.h"
#include "mirtk/EdgeConnectivity.h"

#include "vtkSmartPointer.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"


namespace mirtk {


/**
 * Filter for mesh node data
 */
class PointDataFilter : public MeshFilter
{
  mirtkAbstractMacro(PointDataFilter);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Maximum edge-connectivity of neighboring nodes
  ///
  /// Used instead of _Radius attribute when radius is non-positive value.
  /// Ignored when pre-computed _Neighbors edge-connectivity set.
  mirtkPublicAttributeMacro(int, Connectivity);

  /// Maximum point distance of neighboring points
  ///
  /// Used instead of _Connectivity attribute when set to positive value.
  /// Ignored when pre-computed _Neighbors edge-connectivity set.
  mirtkPublicAttributeMacro(double, Radius);

  /// Set of considered neighboring nodes for each mesh node
  ///
  /// Computed based on either _Connectivity or _Radius when no pre-computed
  /// edge-connectivity is specified before the filter is executed.
  mirtkPublicAttributeMacro(SharedPtr<EdgeConnectivity>, Neighbors);

  /// Name of (input and) output point data array
  ///
  /// When an input _DataArray is given, this name is assigned to the respective
  /// output point data array. Otherwise, the name of the input data array must
  /// be set and is used also for the output data array. When no input _DataArray
  /// is given, the name of the point data array of the input mesh must be set.
  mirtkPublicAttributeMacro(string, DataName);

  /// Input point data array
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, InputData);

  /// Output point data array
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, OutputData);

  /// Copy attributes of this filter from another instance
  void CopyAttributes(const PointDataFilter &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  PointDataFilter();

  /// Copy constructor
  PointDataFilter(const PointDataFilter &);

  /// Assignment operator
  PointDataFilter &operator =(const PointDataFilter &);

  /// Destructor
  virtual ~PointDataFilter();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

};


} // namespace mirtk

#endif // MIRTK_PointDataFilter_H

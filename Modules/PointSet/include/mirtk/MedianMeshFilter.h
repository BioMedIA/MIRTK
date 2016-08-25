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

#ifndef MIRTK_MedianMeshFilter_H
#define MIRTK_MedianMeshFilter_H

#include "mirtk/MeshFilter.h"

#include "vtkSmartPointer.h"
#include "vtkDataArray.h"


namespace mirtk {


/**
 * Component-wise median filter for mesh node data
 */
class MedianMeshFilter : public MeshFilter
{
  mirtkObjectMacro(MedianMeshFilter);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Maximum edge-connectivity of neighboring nodes
  mirtkPublicAttributeMacro(int, Connectivity);

  /// Input point data array
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, DataArray);

  /// Name of (input and) output point data array
  ///
  /// When an input _DataArray is given, this name is assigned to the respective
  /// output point data array. Otherwise, the name of the input data array must
  /// be set and is used also for the output data array. When no input _DataArray
  /// is given, the name of the point data array of the input mesh must be set.
  mirtkPublicAttributeMacro(string, DataArrayName);

  /// Copy attributes of this filter from another instance
  void CopyAttributes(const MedianMeshFilter &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  MedianMeshFilter();

  /// Copy constructor
  MedianMeshFilter(const MedianMeshFilter &);

  /// Assignment operator
  MedianMeshFilter &operator =(const MedianMeshFilter &);

  /// Destructor
  virtual ~MedianMeshFilter();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Execute filter
  virtual void Execute();

};


} // namespace mirtk

#endif // MIRTK_MedianMeshFilter_H

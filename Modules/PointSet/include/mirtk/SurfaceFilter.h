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

#ifndef MIRTK_SurfaceFilter_H
#define MIRTK_SurfaceFilter_H

#include "mirtk/MeshFilter.h"


namespace mirtk {


/**
 * Base class for filters which process polygonal surface meshes
 */
class SurfaceFilter : public MeshFilter
{
  mirtkAbstractMacro(SurfaceFilter);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this filter from another instance
  void CopyAttributes(const SurfaceFilter &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  SurfaceFilter();

  /// Copy constructor
  SurfaceFilter(const SurfaceFilter &);

  /// Assignment operator
  SurfaceFilter &operator =(const SurfaceFilter &);

  /// Destructor
  virtual ~SurfaceFilter();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

};


} // namespace mirtk

#endif // MIRTK_SurfaceFilter_H

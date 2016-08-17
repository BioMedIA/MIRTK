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

#ifndef MIRTK_SurfacePatches_H
#define MIRTK_SurfacePatches_H

#include "mirtk/SurfaceFilter.h"

#include "mirtk/Array.h"


namespace mirtk {


/**
 * Labels connected surface patches
 *
 * Each surface patch consists of interior cells with edges that have only
 * one cell edge neigbhor. Only boundary cells of a surface patch have
 * an edge with no or multiple cell edge neigbhors.
 *
 * Other definitions of patch borders based on surface attributes may be
 * considered by specializations (subclasses) of this filter.
 */
class SurfacePatches : public SurfaceFilter
{
  mirtkObjectMacro(SurfacePatches);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Enumeration of possible orderings of surface patches
  enum Ordering
  {
    NoOrdering,   ///< No sorting of output patches
    LargestFirst, ///< Sort by decreasing patch size
    SmallestFirst ///< Sort by increasing patch size
  };

  // ---------------------------------------------------------------------------
  // Attributes
private:

  /// Ordering of surface patch labels
  mirtkPublicAttributeMacro(enum Ordering, Ordering);

  /// Number of surface patches
  mirtkReadOnlyAttributeMacro(int, NumberOfPatches);

  /// Sizes of each surface patch in number of cells
  mirtkReadOnlyAttributeMacro(Array<int>, PatchSize);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const SurfacePatches &);

  // ---------------------------------------------------------------------------
  // Construction/destruction
public:

  /// Constructor
  SurfacePatches();

  /// Copy constructor
  SurfacePatches(const SurfacePatches &);

  /// Assignment operator
  SurfacePatches &operator =(const SurfacePatches &);

  /// Destructor
  virtual ~SurfacePatches();

  // ---------------------------------------------------------------------------
  // Execution
protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Execute filter
  virtual void Execute();

  /// Finalize filter output
  virtual void Finalize();

  // ---------------------------------------------------------------------------
  // Output
public:

  /// Get cell data array storing surface patch labels
  vtkDataArray *GetLabelsArray() const;
};


} // namespace mirtk

#endif // MIRKT_SurfacePatches_H

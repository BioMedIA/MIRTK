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

#ifndef MIRTK_ClosestPointLabel_H
#define MIRTK_ClosestPointLabel_H

#include "mirtk/ClosestPoint.h"

#include "mirtk/UnorderedMap.h"

#include "vtkSmartPointer.h"
#include "vtkIdTypeArray.h"


namespace mirtk {


/**
 * Closest point correspondence map
 *
 * This class establishes point correspondence based on a known parcellation
 * of the surfaces, i.e., the discrete "Labels" point data array
 * (cf. MIRTK polydataassignlabels tool). A source point is said to correspond
 * to a target point if it is the closest point on the source surface which
 * has the same label as the target point. If no such point is found,
 * the target point is marked as outlier.
 */
class ClosestPointLabel : public ClosestPoint
{
  mirtkObjectMacro(ClosestPointLabel);

  typedef UnorderedMap<int, int> LabelToComponentIndexMap;

  // ---------------------------------------------------------------------------
  // Attributes

  /// IDs of points on the target surface closest to each respective surface label
  mirtkAttributeMacro(vtkSmartPointer<vtkIdTypeArray>, TargetIds);

  /// Map of target surface label to target ID array component index
  mirtkAttributeMacro(LabelToComponentIndexMap, TargetComponent);

  /// IDs of points on the source surface closest to each respective surface label
  mirtkAttributeMacro(vtkSmartPointer<vtkIdTypeArray>, SourceIds);

  /// Map of target surface label to target ID array component index
  mirtkAttributeMacro(LabelToComponentIndexMap, SourceComponent);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  ClosestPointLabel();

  /// Copy constructor
  ClosestPointLabel(const ClosestPointLabel &);

  /// Copy construct a new instance
  virtual PointCorrespondence *NewInstance() const;

  /// Destructor
  virtual ~ClosestPointLabel();

  /// Type enumeration value
  virtual TypeId Type() const;

  // ---------------------------------------------------------------------------
  // Correspondences

protected:

  /// Common (re-)initialization steps of this class
  /// \note Must be a non-virtual function!
  void Init();

public:

  /// Initialize correspondence map after input and parameters are set
  virtual void Initialize();

  /// Reinitialize correspondence map after change of input topology
  virtual void Reinitialize();

  /// Update correspondence map
  virtual void Update();

};


} // namespace mirtk

#endif // MIRTK_ClosestPointLabel_H

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

#ifndef MIRTK_ClosestPoint_H
#define MIRTK_ClosestPoint_H

#include "mirtk/PointCorrespondence.h"

#include "mirtk/Array.h"


namespace mirtk {


/**
 * Closest point correspondence map
 */
class ClosestPoint : public PointCorrespondence
{
  mirtkObjectMacro(ClosestPoint);

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// Mark correspondences with distance greater than the mean distance plus
  /// _Sigma times the standard deviation of the current point distances as outliers
  mirtkPublicAttributeMacro(double, Sigma);

  /// Mark correspondences with distance greater than this value as outliers
  /// Unused if _Sigma is set to a non-negative value
  mirtkPublicAttributeMacro(double, MaxDistance);

  /// Maximum squared distance
  mirtkAttributeMacro(double, MaxSquaredDistance);

  /// Indices of target points corresponding to the respective source points
  mirtkAttributeMacro(Array<int>, TargetIndex);

  /// Distance of closest target points from source samples
  mirtkAttributeMacro(Array<double>, TargetDistance);

  /// Indices of source points corresponding to the respective target points
  mirtkAttributeMacro(Array<int>, SourceIndex);

  /// Distance of closest source points from target samples
  mirtkAttributeMacro(Array<double>, SourceDistance);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  ClosestPoint();

  /// Copy constructor
  ClosestPoint(const ClosestPoint &);

  /// Copy construct a new instance
  virtual PointCorrespondence *NewInstance() const;

  /// Destructor
  virtual ~ClosestPoint();

  /// Type enumeration value
  virtual TypeId Type() const;

  // ---------------------------------------------------------------------------
  // Parameters

  // Import other overloads
  using PointCorrespondence::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameter key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Correspondences

  /// Initialize correspondence map
  virtual void Initialize();

  /// Update correspondence map
  virtual void Update();

  /// Update correspondence map after convergence
  virtual bool Upgrade();

  /// Get untransformed target point corresponding to i-th source (sample) point
  virtual bool GetInputTargetPoint(int, Point &) const;

  /// Get untransformed source point corresponding to i-th target (sample) point
  virtual bool GetInputSourcePoint(int, Point &) const;

  /// Get (transformed) target point corresponding to i-th source (sample) point
  virtual bool GetTargetPoint(int, Point &) const;

  /// Get (transformed) source point corresponding to i-th target (sample) point
  virtual bool GetSourcePoint(int, Point &) const;

  /// Get index of target point corresponding to i-th source (sample) point
  ///
  /// \returns Index of corresponding target point and -1 if point is an outlier
  ///          or its corresponding point is not a target vertex position.
  virtual int GetTargetIndex(int) const;

  /// Get index of source point corresponding to i-th target (sample) point
  ///
  /// \returns Index of corresponding source point and -1 if point is an outlier
  ///          or its corresponding point is not a source vertex position.
  virtual int GetSourceIndex(int) const;

};


} // namespace mirtk

#endif // MIRTK_ClosestPoint_H

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

#ifndef MIRTK_SpectralMatch_H
#define MIRTK_SpectralMatch_H

#include "mirtk/PointCorrespondence.h"

#include "mirtk/PointSet.h"


namespace mirtk {


class PointLocator;


/**
 * Correspondence map based on diffeomorphic spectral matching
 */
class SpectralMatch : public PointCorrespondence
{
  mirtkObjectMacro(SpectralMatch);

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// Number of spectral nearest neighbors to consider
  mirtkPublicAttributeMacro(int, NumberOfNeighbors);

  /// Standard deviation of Gaussian interpolation kernel
  mirtkPublicAttributeMacro(double, Sigma);

  /// Target points corresponding to source points
  mirtkAttributeMacro(PointSet, TargetPoints);

  /// Source points corresponding to target points
  mirtkAttributeMacro(PointSet, SourcePoints);

  PointLocator *_TargetLocator;
  PointLocator *_SourceLocator;

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  SpectralMatch();

  /// Copy constructor
  SpectralMatch(const SpectralMatch &);

  /// Copy construct a new instance
  virtual PointCorrespondence *NewInstance() const;

  /// Destructor
  virtual ~SpectralMatch();

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

};


} // namespace mirtk

#endif // MIRTK_SpectralMatch_H

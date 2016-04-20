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

#ifndef MIRTK_FiducialMatch_H
#define MIRTK_FiducialMatch_H

#include "mirtk/PointCorrespondence.h"

#include "mirtk/Array.h"


namespace mirtk {


/**
 * Point correspondence based on the fiducial indices
 */
class FiducialMatch : public PointCorrespondence
{
  mirtkObjectMacro(FiducialMatch);

  /// Optional map of corresponding target point indices
  mirtkPublicAttributeMacro(Array<int>, TargetIndex);

  /// Optional map of corresponding source point indices
  mirtkPublicAttributeMacro(Array<int>, SourceIndex);

  /// Name of text file listing corresponding point indices
  mirtkPublicAttributeMacro(string, CorrespondenceMap);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  FiducialMatch();

  /// Copy constructor
  FiducialMatch(const FiducialMatch &);

  /// Copy construct a new instance
  virtual PointCorrespondence *NewInstance() const;

  /// Destructor
  virtual ~FiducialMatch();

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
protected:

  /// Check if correspondence map has a valid entry for every target point
  void ValidateCorrespondenceMap(const RegisteredPointSet *,
                                 const RegisteredPointSet *,
                                 const Array<int> &,
                                 const char * = NULL) const;

  /// Compute inverse correspondence map
  void InvertCorrespondenceMap(const RegisteredPointSet *,
                               const RegisteredPointSet *,
                               const Array<int> &, Array<int> &) const;

public:

  /// Initialize correspondence map
  virtual void Initialize();

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

#endif // MIRTK_FiducialMatch_H

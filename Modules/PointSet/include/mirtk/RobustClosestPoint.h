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

#ifndef MIRTK_RobustClosestPoint_H
#define MIRTK_RobustClosestPoint_H

#include "mirtk/FuzzyCorrespondence.h"


namespace mirtk {


/**
 * Improved (iterative) closest point (ICP) correspondences
 *
 * The point correspondences returned by an instance of this class implement
 * an improvied iterative closest point matching with outlier rejection, the
 * reference method used by
 *
 *   Chui and Rangarajan, "A new point matching algorithm for non-rigid registration",
 *   Computer Vision and Image Understanding, 89(2-3), pp. 114–141, 2003.
 *
 * for comparison to the proposed robust point matching (RPM) algorithm.
 * This algorithm is realized by RobustPointMatch.
 *
 * \sa FiducialRegistrationError
 */
class RobustClosestPoint : public FuzzyCorrespondence
{
  mirtkObjectMacro(RobustClosestPoint);

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// Factor by which standard deviation of point distances is multiplied
  /// for outlier rejection during (iterative) closest point matching
  mirtkPublicAttributeMacro(double, Sigma);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Default constructor
  RobustClosestPoint();

  /// Construct correspondence map and initialize it
  RobustClosestPoint(const RegisteredPointSet *,
                     const RegisteredPointSet *);

  /// Copy constructor
  RobustClosestPoint(const RobustClosestPoint &);

  /// Copy construct a new instance
  virtual PointCorrespondence *NewInstance() const;

  /// Destructor
  virtual ~RobustClosestPoint();

  /// Type enumeration value
  virtual TypeId Type() const;

  // ---------------------------------------------------------------------------
  // Parameters

  // Import other overloads
  using FuzzyCorrespondence::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameter key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Correspondences
protected:

  /// (Re-)calculate weights of correspondence links
  virtual void CalculateWeights();

};


} // namespace mirtk

#endif // MIRTK_RobustClosestPoint_H 

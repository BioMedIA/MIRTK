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

#ifndef MIRTK_RobustPointMatch_H
#define MIRTK_RobustPointMatch_H

#include "mirtk/FuzzyCorrespondence.h"


namespace mirtk {


/**
 * Robust point match (RPM) correspondences
 *
 * The point correspondences returned by an instance of this class correspond
 * to the cluster centers of the Robust Point Matching (RPM) algorithm of
 *
 *   Chui and Rangarajan, "A new point matching algorithm for non-rigid registration",
 *   Computer Vision and Image Understanding, 89(2-3), pp. 114–141, 2003.
 *
 * where the smoothness regularization term is separate. Therefore, the weight
 * of the smoothness term is not directly coupled to the annealing temperature
 * used by this class, unlike the original RPM proposal.
 *
 * \sa FiducialRegistrationError
 */
class RobustPointMatch : public FuzzyCorrespondence
{
  mirtkObjectMacro(RobustPointMatch);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Initial temperature of the deterministic annealing process
  mirtkPublicAttributeMacro(double, InitialTemperature);

  /// Annealing rate by which the temperature is multiplied or number of NN
  ///
  /// If set to a negative value, the average mean distance to the
  /// corresponding number of nearest neighbors is used as next temperature.
  mirtkPublicAttributeMacro(double, AnnealingRate);

  /// Final temperature which stops the deterministic annealing process
  mirtkPublicAttributeMacro(double, FinalTemperature);

  /// Current temperature of the deterministic annealing process
  mirtkPublicAttributeMacro(double, Temperature);

  /// Variance of extra features (resp. of their differences)
  mirtkPublicAttributeMacro(double, VarianceOfFeatures);

  /// Cluster for outliers in target (i.e., centroid of source points!)
  mirtkReadOnlyAttributeMacro(Point, TargetOutlierCluster);

  /// Cluster for outliers in source (i.e., centroid of target points!)
  mirtkReadOnlyAttributeMacro(Point, SourceOutlierCluster);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Default constructor
  RobustPointMatch();

  /// Construct correspondence map and initialize it
  RobustPointMatch(const RegisteredPointSet *,
                   const RegisteredPointSet *);

  /// Copy constructor
  RobustPointMatch(const RobustPointMatch &);

  /// Copy construct a new instance
  virtual PointCorrespondence *NewInstance() const;

  /// Destructor
  virtual ~RobustPointMatch();

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

  /// Initialize correspondence map
  virtual void Initialize();

  /// Update correspondence map after convergence
  virtual bool Upgrade();

protected:

  /// Initialize annealing process
  virtual void InitializeAnnealing();

  /// (Re-)calculate weights of correspondence links
  virtual void CalculateWeights();

};


} // namespace mirtk

#endif // MIRTK_RobustPointMatch_H

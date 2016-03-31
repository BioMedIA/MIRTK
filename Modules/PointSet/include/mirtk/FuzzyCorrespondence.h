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

#ifndef MIRTK_FuzzyCorrespondence_H
#define MIRTK_FuzzyCorrespondence_H

#include "mirtk/PointCorrespondence.h"

#include "mirtk/Array.h"
#include "mirtk/PointSet.h"
#include "mirtk/SparseMatrix.h"


namespace mirtk {


/**
 * Fuzzy point correspondences
 *
 * The point correspondences realized by this class are inspired by the
 * softassign approach as in particular employed by the Robust Point Matching
 * (RPM) algorithm proposed in
 *
 *   Chui and Rangarajan, "A new point matching algorithm for non-rigid registration",
 *   Computer Vision and Image Understanding, 89(2-3), pp. 114–141, 2003.
 *
 * This class is the common base of the different variants of softassign
 * point correspondences, two of which implementing those presented in the
 * MATLAB Demo by the authors of aforementioned paper.
 *
 * \sa FiducialRegistrationError
 */
class FuzzyCorrespondence : public PointCorrespondence
{
  mirtkAbstractMacro(FuzzyCorrespondence);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Sparse matrix type used for fuzzy correspondence weights
  typedef GenericSparseMatrix<float> WeightMatrix;

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// Minimum correspondence weight (i.e., threshold used to truncate radial weight function)
  mirtkPublicAttributeMacro(double, MinWeight);

  /// Weight matrix of fuzzy point correspondences
  mirtkReadOnlyAttributeMacro(WeightMatrix, Weight);

  /// Untransformed target clusters corresponding to source points
  mirtkAttributeMacro(PointSet, InputTargetClusters);

  /// Untransformed source clusters corresponding to target points
  mirtkAttributeMacro(PointSet, InputSourceClusters);

  /// (Transformed) Target clusters corresponding to source points
  mirtkAttributeMacro(PointSet, TargetClusters);

  /// (Transformed) Source clusters corresponding to target points
  mirtkAttributeMacro(PointSet, SourceClusters);

  /// Whether a given target cluster is an outlier
  mirtkAttributeMacro(Array<bool>, TargetOutlier);

  /// Whether a given source cluster is an outlier
  mirtkAttributeMacro(Array<bool>, SourceOutlier);

  /// Standard deviation of normally distributed noise with zero mean
  mirtkPublicAttributeMacro(double, GaussianNoise);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Default constructor
  FuzzyCorrespondence();

  /// Copy constructor
  FuzzyCorrespondence(const FuzzyCorrespondence &);

public:

  /// Destructor
  virtual ~FuzzyCorrespondence();

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

  /// Reinitialize correspondence map after change of input topology
  virtual void Reinitialize();

protected:

  /// Common (re-)initialization steps of this class
  /// \note Must be a non-virtual function!
  void Init();

public:

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

  // ---------------------------------------------------------------------------
  // Internal
protected:

  /// (Re-)calculate correspondence weights with optional additional outlier row/column
  virtual void CalculateWeights() = 0;

  /// Add normally distributed noise to correspondence weights
  virtual void AddGaussianNoise();

  /// Normalize correspondence matrix using Sinkhorn-Knopp algorithm
  virtual void NormalizeWeights();

  /// (Re-)calculate cluster centers from correspondence weights
  virtual void CalculateClusters();

  // ---------------------------------------------------------------------------
  // Debugging

  /// Write input of data fidelity term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

};


} // namespace mirtk

#endif // MIRTK_FuzzyCorrespondence_H

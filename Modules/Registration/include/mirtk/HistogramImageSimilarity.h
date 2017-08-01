/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2017 Andreas Schuh
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

#ifndef MIRTK_HistogramImageSimilarity_H
#define MIRTK_HistogramImageSimilarity_H

#include "mirtk/ImageSimilarity.h"

#include "mirtk/Math.h"
#include "mirtk/Histogram2D.h"
#include "mirtk/Parallel.h"


namespace mirtk {


/**
 * Base class for probabilistic image similarity measures
 *
 * Subclasses of this intensity-based image similarity measure compute similarity
 * from the joint and marginal probabilities of the intensities in the images.
 * An estimate of the probabilities is obtained using a joint histogram and
 * cubic B-spline Parzen Windows for a continuous representation.
 */
class HistogramImageSimilarity : public ImageSimilarity
{
  mirtkAbstractMacro(HistogramImageSimilarity);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Type of joint histogram bin values
  #if MIRKT_USE_FLOAT_BY_DEFAULT
    typedef float  BinType;
  #else
    typedef double BinType;
  #endif

  /// Type of joint histogram
  typedef Histogram2D<BinType> JointHistogramType;

  /// Get default number of bins for given input image
  ///
  /// \param[in] image         Intensity image.
  /// \param[in] min_intensity Pre-computed minimum intensity value.
  /// \param[in] max_intensity Pre-computed maximum intensity value.
  ///
  /// \returns Number of histogram bins.
  static int DefaultNumberOfBins(const BaseImage *image,
                                 double min_intensity = NaN,
                                 double max_intensity = NaN);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Joint histogram of raw intensity samples (no Parzen window function)
  ///
  /// This class allocates a joint histogram instance and takes care of updating
  /// it upon each Update call. When an external _Samples histogram is provided,
  /// however, the owner of this external histogram is reponsible for updating
  /// and this histogram before the Update of this image similarity.
  mirtkLooseComponentMacro(JointHistogramType, Samples);

  /// Joint histogram of image intensities (cubic B-spline Parzen windows)
  mirtkAttributeMacro(JointHistogramType, Histogram);

  /// Whether to use a Parzen window function
  mirtkPublicAttributeMacro(bool, UseParzenWindow);

  /// Whether to pad histogram when performing Parzen window estimation
  mirtkPublicAttributeMacro(bool, PadHistogram);

  /// Number of histogram bins for target image intensities
  mirtkPublicAttributeMacro(int, NumberOfTargetBins);

  /// Number of histogram bins for source image intensities
  mirtkPublicAttributeMacro(int, NumberOfSourceBins);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const HistogramImageSimilarity &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  HistogramImageSimilarity(const char * = "", double = 1.0);

  /// Copy constructor
  HistogramImageSimilarity(const HistogramImageSimilarity &);

  /// Assignment operator
  HistogramImageSimilarity &operator =(const HistogramImageSimilarity &);

  /// Destructor
  virtual ~HistogramImageSimilarity();

  // ---------------------------------------------------------------------------
  // Parameters
protected:

  /// Set parameter value from string
  virtual bool SetWithPrefix(const char *, const char *);

public:

  // Do not hide other overload
  using ImageSimilarity::Parameter;

  /// Get parameter key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize similarity measure once input and parameters have been set
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Update moving image and internal state of similarity measure
  virtual void Update(bool = true);

  /// Exclude region from similarity evaluation
  ///
  /// Called by ApproximateGradient \b before the registered image region of
  /// the transformed image is updated.
  virtual void Exclude(const blocked_range3d<int> &);

  /// Include region in similarity evaluation
  ///
  /// Called by ApproximateGradient \b after the registered image region of
  /// the transformed image is updated.
  virtual void Include(const blocked_range3d<int> &);

protected:

  /// Update estimated joint histogram after raw samples histogram has changed
  virtual void UpdateHistogram();

  // ---------------------------------------------------------------------------
  // Debugging

public:

  /// Print debug information
  virtual void Print(Indent = 0) const;

  /// Write input of data fidelity term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

};


} // namespace mirtk

#endif // MIRTK_HistogramImageSimilarity_H

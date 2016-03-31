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

#ifndef MIRTK_DummyImageSimilarity_H
#define MIRTK_DummyImageSimilarity_H

#include "mirtk/ImageSimilarity.h"


namespace mirtk {


/**
 * Illustrative dummy implementation of custom image similarity measure
 *
 * An instance of an image (dis-)similarity measure has two instances of
 * RegisteredImage as input. The more similar these images are, the
 * lower the image \b dissimilarity value returned by the Evaluate member
 * function should be. The RegisteredImage instances are owned by
 * an image dissimilarity measure instance and are not shared used by any
 * other instance. Therefore, the image dissimilarity instance may modify
 * these images. The RegisteredImage::Update functions apply the
 * current transformation estimate and resample the original (downsampled)
 * input images on a common image grid shared by both RegisteredImage
 * instances. Note that even though the registered image instances are
 * referred to as _Target and _Source, respectively, either of these
 * two images can be moving or fixed. The image dissimilarity measure
 * can check the RegisteredImage::Transformation() to see whether an
 * image is moving (!= NULL) or fixed (== NULL). Usually, the source image
 * is transformed while the target image remains fixed. In case of an
 * inverse consistent registration, the target image is transformed by the
 * inverse transformation while the source image remains fixed instead.
 * In case of a symmetric registration, both input images are transformed
 * by half the forward and backward (i.e., inverse) transformation, respectively.
 *
 * In order for this image dissimilarity to be instantiated, add an enumeration
 * value to SimilarityMeasure defined in Common/include/mirtkEnums.h.
 * Don't forget to also extend the corresponding FromString and ToString
 * template specializations for SimilarityMeasure. Next, modify the static
 * factory method ImageSimilarity::New at the top of ImageSimilarity.cc.
 *
 * \sa ImageSimilarity, RegisteredImage, SumOfSquaredIntensityDifferences
 */
class DummyImageSimilarity : public ImageSimilarity
{
  mirtkObjectMacro(DummyImageSimilarity);

  // ---------------------------------------------------------------------------
  // Attributes
  //
  // To declare class attributes (i.e., member variables) use one of the following
  // macros defined in mirtkObject.h:
  // - mirtkAttributeMacro:         member variable with protected getter/setter
  // - mirtkPublicAttributeMacro:   member variable with public getter/setter
  // - mirtkReadOnlyAttributeMacro: member variable with public getter
  // - mirtkComponentMacro:         pointer to object owned by this  instance
  // - mirtkAggregateMacro:         pointer to object owned by other instance

  /// A boolean option
  ///
  /// This macro adds a protected member "bool _Option" and public setter
  /// function "void Option(bool)" and getter function "bool Option() const".
  mirtkPublicAttributeMacro(bool, Option);

  /// A floating point parameter
  ///
  /// This macro adds a protected member "double _Parameter" and public setter
  /// function "void Parameter(double)" and getter function "double Parameter() const".
  mirtkPublicAttributeMacro(double, Parameter);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  ///
  /// \param[in] name Name of the image dissimilarity term as specified by the
  ///                 user in the configuration file. For example,
  ///                 "Energy function = SIM[name](I(1), I(2) o T)".
  DummyImageSimilarity(const char *name = "");

  /// Copy constructor
  DummyImageSimilarity(const DummyImageSimilarity &);

  /// Assignment operator
  DummyImageSimilarity &operator =(const DummyImageSimilarity &);

  /// Destructor
  ~DummyImageSimilarity();

  // ---------------------------------------------------------------------------
  // Parameter

  // Import other overloads of Parameter member function to not hide them
  using ImageSimilarity::Parameter;

  /// Set parameter value from string
  ///
  /// \param[in] name  Name of parameter as given in the configuration file.
  /// \param[in] value Value to set as given in the configuration file.
  ///
  /// \returns Whether the given parameter is known by this class and its value
  ///          has been set to the given value successfully.
  virtual bool Set(const char *name, const char *value);

  /// Get parameter name/value pairs
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize dissimilarity measure
  virtual void Initialize();

  /// Update moving input image(s) and internal state of dissimilarity measure
  ///
  /// \param[in] gradient Whether the image dissimilarity gradient will be evaluated
  ///                     next using ImageSimilarity::EvaluateGradient.
  ///                     Otherwise, only the following evaluation of the image
  ///                     dissimilarity must be prepared by this instance.
  virtual void Update(bool gradient = true);

  // ---------------------------------------------------------------------------
  // Evaluation
protected:

  /// Exclude region from dissimilarity evaluation
  ///
  /// Called by ImageSimilarity::ApproximateGradient \b before the registered
  /// image region of the transformed image is updated.
  ///
  /// \param[in] region Image region defining those voxels in the input images
  ///                   which are to be excluded from the dissimilarity measure
  ///                   because the intensities of these voxels may change.
  virtual void Exclude(const blocked_range3d<int> &region);

  /// Include region in dissimilarity evaluation
  ///
  /// Called by ImageSimilarity::ApproximateGradient \b after the registered
  /// image region of the transformed image was updated.
  ///
  /// \param[in] region Image region defining those voxels in the input images
  ///                   which were updated and are to be included in the
  ///                   dissimilarity measure again.
  virtual void Include(const blocked_range3d<int> &region);

  /// Evaluate dissimilarity of images
  virtual double Evaluate();

  /// Evaluate non-parametric/voxel-based dissimilarity gradient w.r.t the given image
  ///
  /// This function is called by ImageSimilarity::EvaluateGradient once for each
  /// moving input image. This image is passed as first argument such that this
  /// function can adjust, for example, the sign of the gradient vectors.
  virtual bool NonParametricGradient(const RegisteredImage *, GradientImageType *);

};


} // namespace mirtk

#endif // MIRTK_DummyImageSimilarity_H

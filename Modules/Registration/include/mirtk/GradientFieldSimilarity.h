/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Stefan Pszczolkowski Parraguez, Andreas Schuh
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

#ifndef MIRTK_GradientFieldSimilarity_H
#define MIRTK_GradientFieldSimilarity_H

#include "mirtk/ImageSimilarity.h"


namespace mirtk {


/**
 * Base class for gradient field similarity measures
 *
 * Subclasses of this image similarity measure evaluate similarity of two images
 * based on their intensity gradient.
 */
class GradientFieldSimilarity : public ImageSimilarity
{
  mirtkObjectMacro(GradientFieldSimilarity);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Noise parameter for target image
  mirtkPublicAttributeMacro(bool, IgnoreJacobianGradientWrtDOFs);

  /// Transformed gradient of the target image
  /// Used only if the target image is being transformed
  mirtkAttributeMacro(GenericImage<RegisteredImage::VoxelType>, TargetTransformedGradient);

  /// Transformed gradient of the source image
  /// Used only if the source image is being transformed
  mirtkAttributeMacro(GenericImage<RegisteredImage::VoxelType>, SourceTransformedGradient);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Constructor
  GradientFieldSimilarity(const char * = "", double = 1.0);

  /// Copy constructor
  GradientFieldSimilarity(const GradientFieldSimilarity &);

  /// Destructor
  virtual ~GradientFieldSimilarity();

  // ---------------------------------------------------------------------------
  // Parameters

public:

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  // Import other overloads
  using ImageSimilarity::Parameter;

  /// Get parameter name/value map
  virtual ParameterList Parameter() const;

protected:

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize similarity measure once input and parameters have been set
  /// \param[in] domain Image domain on which the similarity is evaluated.
  virtual void InitializeInput(const ImageAttributes &domain);

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Convert non-parametric similarity gradient into gradient
  /// w.r.t transformation parameters
  ///
  /// This function calls Transformation::ParametricGradient of the
  /// transformation to apply the chain rule in order to obtain the similarity
  /// gradient w.r.t the transformation parameters. It adds the weighted gradient
  /// to the final registration energy gradient.
  ///
  /// \param[in]    image       Transformed image.
  /// \param[in]    np_gradient Voxel-wise non-parametric gradient.
  /// \param[inout] gradient    Gradient to which the computed parametric gradient
  ///                           is added, after multiplication by the given \p weight.
  /// \param[in]    weight      Weight of image similarity.
  virtual void ParametricGradient(const RegisteredImage *image,
                                  GradientImageType     *np_gradient,
                                  double                *gradient,
                                  double                 weight);

  /// Update moving input image(s) and internal state of similarity measure
  virtual void Update(bool = true);

  /// Reorient transformed image gradient according to dI(y)/dy * dy/dx
  void ReorientGradient(RegisteredImage *, bool = false);

  /// Multiply similarity gradient by 2nd order derivatives of transformed image
  ///
  /// \param[in]     image    Transformed image
  /// \param[in,out] gradient Input must be the gradient of the image similarity
  ///                         w.r.t. the transformed \p image gradient. Output is
  ///                         the voxel-wise gradient of the similarity w.r.t. T(x).
  void MultiplyByImageHessian(const RegisteredImage *image,
                              GradientImageType     *gradient);

};


} // namespace mirtk

#endif // MIRTK_GradientFieldSimilarity_H 

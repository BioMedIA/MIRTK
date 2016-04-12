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

#ifndef MIRTK_ImageSimilarity_H
#define MIRTK_ImageSimilarity_H

#include "mirtk/DataFidelity.h"
#include "mirtk/SimilarityMeasure.h"
#include "mirtk/Parallel.h"
#include "mirtk/FreeFormTransformation.h"
#include "mirtk/RegisteredImage.h"


namespace mirtk {


/**
 * Base class for image similarity measures
 *
 * The lower the value of the similarity measure, the more similar the images are.
 * It may therefore more precisely be referred to as dissimilarity measure, but
 * both terms are for historic reasons used interchangeably in this framework.
 *
 * If the transformed image has a foreground region defined, we assume that this
 * corresponds to the region which needs to be matched with the untransformed
 * image. Therefore we make use of the entire transformed foreground region to
 * evaluate similarity and in particular the gradient forces defined within
 * this region. Note that outside this source foreground region, the gradient
 * is always zero and thus cannot be used to drive the image registration.
 *
 * Otherwise, the region of interest is generally defined by the foreground of
 * the untransformed image. For each such target voxel we want to find a suitable
 * corresponding voxel in the other image. Therefore, we restrict the similarity
 * evaluation to this region. In case of a symmetric transformation of both
 * images, the union of the foreground of both images defines the region for
 * which similarity is evaluated and forces are computed.
 *
 * Subclasses which implement a particular similarity measure should call
 * the IsForeground member function to decide whether or not to consider
 * a given voxel of the grid on which the registered images are defined.
 *
 * \note The similarity measure owns the registered input images and may thus
 *       modify these to improve the runtime of each Update step.
 */
class ImageSimilarity : public DataFidelity
{
  mirtkAbstractMacro(ImageSimilarity);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Voxel type of registered images
  typedef RegisteredImage::VoxelType           VoxelType;

  /// Type of similarity gradient image
  typedef RegisteredImage::GradientImageType   GradientImageType;

  /// Type of similarity gradient components
  typedef GradientImageType::VoxelType         GradientType;

  // ---------------------------------------------------------------------------
  // Attributes

  /// (Transformed) Target image
  mirtkLooseComponentMacro(RegisteredImage, Target);

  /// (Transformed) Source image
  mirtkLooseComponentMacro(RegisteredImage, Source);

  /// Domain on which to evaluate similarity
  mirtkPublicAttributeMacro(ImageAttributes, Domain);

  /// Mask which defines arbitrary domain on which the similarity is evaluated
  mirtkPublicAggregateMacro(BinaryImage,     Mask);

  /// Memory for (non-parametric) similarity gradient w.r.t target transformation
  mirtkComponentMacro(GradientImageType, GradientWrtTarget);

  /// Memory for (non-parametric) similarity gradient w.r.t source transformation
  mirtkComponentMacro(GradientImageType, GradientWrtSource);

  /// Memory for (parametric) similarity gradient
  mirtkComponentMacro(double, Gradient);

  /// Number of voxels per registered image
  mirtkPublicAttributeMacro(int, NumberOfVoxels);

  /// Approximate gradient using finite differences even if the similarity
  /// measure implements the NonParametricGradient function
  mirtkPublicAttributeMacro(bool, UseApproximateGradient);

  /// Voxel-wise gradient preconditioning sigma used to supress noise.
  /// A non-positive value disables the voxel-wise preconditioning all together.
  ///
  /// Zikic, D., Baust, M., Kamen, A., & Navab, N. A General Preconditioning
  /// Scheme for Difference Measures in Deformable Registration. In ICCV 2011.
  mirtkPublicAttributeMacro(double, VoxelWisePreconditioning);

  /// Node-based (M)FFD control point gradient preconditioning sigma used to
  /// supress noise. A non-positive value disables the node-based
  /// preconditioning all together.
  ///
  /// Zikic, D., Baust, M., Kamen, A., & Navab, N. A General Preconditioning
  /// Scheme for Difference Measures in Deformable Registration. In ICCV 2011.
  mirtkPublicAttributeMacro(double, NodeBasedPreconditioning);

  /// Skip initialization of target image
  mirtkPublicAttributeMacro(bool, SkipTargetInitialization);

  /// Skip initialization of source image
  mirtkPublicAttributeMacro(bool, SkipSourceInitialization);

  /// Whether Update has not been called since initialization
  mirtkAttributeMacro(bool, InitialUpdate);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const ImageSimilarity &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  ImageSimilarity(const char * = "", double = 1.0);

  /// Copy constructor
  ImageSimilarity(const ImageSimilarity &);

  /// Assignment operator
  ImageSimilarity &operator =(const ImageSimilarity &);

public:

  /// Instantiate specified similarity measure
  static ImageSimilarity *New(SimilarityMeasure,
                              const char * = "", double = 1.0);

  /// Destructor
  virtual ~ImageSimilarity();

  // ---------------------------------------------------------------------------
  // Initialization

protected:

  /// Initialize similarity measure once input and parameters have been set
  /// \param[in] domain Image domain on which the similarity is evaluated.
  virtual void InitializeInput(const ImageAttributes &domain);

public:

  /// Initialize similarity measure once input and parameters have been set
  virtual void Initialize();

  /// Release input target image
  void ReleaseTarget();

  /// Release input source image
  void ReleaseSource();

  // ---------------------------------------------------------------------------
  // Parameters
protected:

  /// Set parameter value from string
  virtual bool SetWithoutPrefix(const char *, const char *);

public:

  // Import other overloads
  using DataFidelity::Parameter;

  /// Get parameter key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Update moving input image(s) and internal state of similarity measure
  virtual void Update(bool = true);

  /// Whether to evaluate similarity at specified voxel
  bool IsForeground(int) const;

  /// Whether to evaluate similarity at specified voxel
  bool IsForeground(int, int, int) const;

  /// Exclude region from similarity evaluation
  ///
  /// Called by ApproximateGradient \b before the registered image region of
  /// the transformed image is updated.
  ///
  /// Override in subclass for more efficient ApproximateGradient evaluation
  /// If the analytic derivation is possible and thus the NonParametericGradient
  /// function is overriden instead, this function is not used. Otherwise,
  /// if this function is not overriden to update the similarity measure,
  /// the Evaluate function must re-evaluate the similarity for all voxels.
  ///
  /// \sa ApproximateGradient, Include
  virtual void Exclude(const blocked_range3d<int> &);

  /// Include region in similarity evaluation
  ///
  /// Called by ApproximateGradient \b after the registered image region of
  /// the transformed image is updated.
  ///
  /// Override in subclass for more efficient ApproximateGradient evaluation
  /// If the analytic derivation is possible and thus the NonParametericGradient
  /// function is overriden instead, this function is not used. Otherwise,
  /// if this function is not overriden to update the similarity measure,
  /// the Evaluate function must re-evaluate the similarity for all voxels.
  ///
  /// \sa ApproximateGradient, Exclude
  virtual void Include(const blocked_range3d<int> &);

protected:

  /// Multiply voxel-wise similarity gradient by transformed image gradient
  ///
  /// This function is intended for use by subclass implementations to compute
  /// the NonParametericGradient. It applies the chain rule to compute
  /// \f$\frac{dSimilarity}{dy} = \frac{dSimilarity}{dI} * \frac{dI}{dy}\f$, given
  /// \f$\frac{dSimilarity}{dI}\f$ as input, where \f$y = T(x)\f$.
  ///
  /// \param[in]     image    Transformed image.
  /// \param[in,out] gradient Input must be the gradient of the image similarity
  ///                         w.r.t. the transformed \p image in x. Output is the
  ///                         voxel-wise gradient of the similarity w.r.t. T(x).
  void MultiplyByImageGradient(const RegisteredImage *image,
                               GradientImageType     *gradient);

  /// Compute voxel-wise non-parametric similarity gradient w.r.t the given image
  ///
  /// Must be implemented by subclasses to compute the similarity gradient.
  /// The base class implementation can be used to convert the similarity gradient
  /// computed w.r.t transformed image (i.e., dSimilarity/dI) to a voxel-wise
  /// non-parametric gradient (i.e., dSimilarity/dT). Note that the input must
  /// be a scalar field only. The base class copies the x component of the input
  /// gradient image to the y and z components before applying the chain rule.
  ///
  /// \param[in]  image    Transformed image
  /// \param[out] gradient Non-parametric similarity gradient.
  ///
  /// \returns Whether voxel-wise similarity gradient has been computed.
  virtual bool NonParametricGradient(const RegisteredImage *image,
                                     GradientImageType     *gradient);

  /// Normalize voxel-wise non-parametric similarity gradient
  ///
  /// Zikic, D., Baust, M., Kamen, A., & Navab, N. A General Preconditioning
  /// Scheme for Difference Measures in Deformable Registration. In ICCV 2011.
  ///
  /// \param[inout] gradient Non-parametric similarity gradient.
  virtual void NormalizeGradient(GradientImageType *gradient);

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

  /// Approximate similarity gradient using finite differences
  ///
  /// If the image similarity does not provide an implementation of the
  /// NonParametricGradient function, the similarity gradient is
  /// approximated instead using finite differences.
  ///
  /// \param[in]     image    Transformed image.
  /// \param[in]     ffd      Free-form deformation.
  /// \param[in,out] gradient Gradient to which the computed parametric gradient
  ///                         is added, after multiplication by the given \p weight.
  /// \param[in]     step     Step size to use for finite differences.
  /// \param[in]     weight   Weight of image similarity.
  void ApproximateGradient(RegisteredImage        *image,
                           FreeFormTransformation *ffd,
                           double                 *gradient,
                           double                  step,
                           double                  weight);

  /// Approximate similarity gradient using finite differences
  ///
  /// If the image similarity does not provide an implementation of the
  /// NonParametricGradient function, the similarity gradient is
  /// approximated instead using finite differences.
  ///
  /// \param[in]     image    Transformed image.
  /// \param[in,out] gradient Gradient to which the computed parametric gradient
  ///                         is added, after multiplication by the given \p weight.
  /// \param[in]     step     Step size to use for finite differences.
  /// \param[in]     weight   Weight of image similarity.
  void ApproximateGradient(RegisteredImage *image,
                           double          *gradient,
                           double           step,
                           double           weight);

  /// Normalize node-based similarity gradient
  ///
  /// Zikic, D., Baust, M., Kamen, A., & Navab, N. A General Preconditioning
  /// Scheme for Difference Measures in Deformable Registration. In ICCV 2011.
  ///
  /// This function applies the normalization for FFD transformations to the
  /// gradient vectors of the control point coefficients. It does nothing for
  /// non-FFD transformations. In case of a multi-level FFD with more than one
  /// active level, it furthermore normalizes the gradient vectors across levels.
  ///
  /// \param[in]     image    Transformed image.
  /// \param[in,out] gradient Parametric similarity gradient.
  virtual void NormalizeGradient(const RegisteredImage *image, double *gradient);

  /// Evaluate similarity gradient
  ///
  /// This function calls the virtual NonParametricGradient function to be
  /// implemented by subclasses for each transformed input image to obtain
  /// the voxel-wise similarity gradient. It then converts this gradient into
  /// a gradient w.r.t the transformation parameters using the ParametricGradient.
  ///
  /// If both target and source are transformed by different transformations,
  /// the resulting gradient vector contains first the derivative values w.r.t
  /// the parameters of the target transformation followed by those computed
  /// w.r.t the parameters of the source transformation. If both images are
  /// transformed by the same transformation, the sum of the derivative values
  /// is added to the resulting gradient vector. This is in particular the case
  /// for a velocity based transformation model which is applied to deform both
  /// images "mid-way". Otherwise, only one input image is transformed
  /// (usually the source) and the derivative values of only the respective
  /// transformation parameters added to the gradient vector.
  ///
  /// \sa NonParametricGradient, ParametricGradient
  ///
  /// \param[in]     image       Transformed image.
  /// \param[in,out] np_gradient Memory for voxel-wise non-parametric gradient.
  /// \param[in,out] gradient    Gradient to which the computed gradient of the
  ///                            image similarity is added after multiplying by
  ///                            the given similarity \p weight.
  /// \param[in]     step        Step size to use for finite differences.
  /// \param[in]     weight      Weight of image similarity.
  virtual void EvaluateGradient(RegisteredImage     *image,
                                GradientImageType   *&np_gradient,
                                double               *gradient,
                                double step, double  weight);

  /// Evaluate similarity gradient
  ///
  /// This function calls the virtual NonParametricGradient function to be
  /// implemented by subclasses for each transformed input image to obtain
  /// the voxel-wise similarity gradient. It then converts this gradient into
  /// a gradient w.r.t the transformation parameters using the ParametricGradient.
  ///
  /// If both target and source are transformed by different transformations,
  /// the resulting gradient vector contains first the derivative values w.r.t
  /// the parameters of the target transformation followed by those computed
  /// w.r.t the parameters of the source transformation. If both images are
  /// transformed by the same transformation, the sum of the derivative values
  /// is added to the resulting gradient vector. This is in particular the case
  /// for a velocity based transformation model which is applied to deform both
  /// images "mid-way". Otherwise, only one input image is transformed
  /// (usually the source) and the derivative values of only the respective
  /// transformation parameters added to the gradient vector.
  ///
  /// \sa NonParametricGradient, ParametricGradient
  ///
  /// \param[in,out] gradient Gradient to which the computed gradient of the
  ///                         image similarity is added after multiplying by
  ///                         the given similarity \p weight.
  /// \param[in]     step     Step size to use for finite differences.
  /// \param[in]     weight   Weight of image similarity.
  virtual void EvaluateGradient(double *gradient, double step, double weight);

  // ---------------------------------------------------------------------------
  // Debugging

public:

  /// Print debug information
  virtual void Print(Indent = 0) const;

  /// Write input of data fidelity term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

  /// Write gradient of data fidelity term w.r.t each transformed input
  virtual void WriteGradient(const char *, const char *) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline bool ImageSimilarity::IsForeground(int idx) const
{
  // Never evaluate similarity outside explicitly specified domain
  if (_Mask && !_Mask->Get(idx)) return false;
  // If both images are transformed (symmetric registration)...
  if (_Target->Transformation() && _Source->Transformation()) {
    // ... evaluate within union of foreground regions
    return _Target->IsForeground(idx) || _Source->IsForeground(idx);
  // If second image is transformed...
  } else if (_Source->Transformation()) {
    // ... try to match masked region of transformed image if given
    if (_Source->HasMask()) return _Source->IsForeground(idx);
    // ... otherwise evaluate similarity for each untransformed voxel,
    //     but consider that similarity gradient is anyway zero for voxels
    //     that are mapped outside the transformed image foreground
    return _Source->IsForeground(idx) && _Target->IsForeground(idx);
  // If first image is transformed same as above with reversed order though
  } else {
    if (_Target->HasMask()) return _Target->IsForeground(idx);
    return _Target->IsForeground(idx) && _Source->IsForeground(idx);
  }
}

// -----------------------------------------------------------------------------
inline bool ImageSimilarity::IsForeground(int i, int j, int k) const
{
  // Never evaluate similarity outside explicitly specified domain
  if (_Mask && !_Mask->Get(i, j, k)) return false;
  // If both images are transformed (symmetric registration)...
  if (_Target->Transformation() && _Source->Transformation()) {
    // ... evaluate within union of foreground regions
    return _Target->IsForeground(i, j, k) || _Source->IsForeground(i, j, k);
  // If second image is transformed...
  } else if (_Source->Transformation()) {
    // ... try to match masked region of transformed image if given
    if (_Source->HasMask()) return _Source->IsForeground(i, j, k);
    // ... otherwise evaluate similarity for each untransformed voxel,
    //     but consider that similarity gradient is anyway zero for voxels
    //     that are mapped outside the transformed image foreground
    return _Source->IsForeground(i, j, k) && _Target->IsForeground(i, j, k);
  // If first image is transformed same as above with reversed order though
  } else {
    if (_Target->HasMask()) return _Target->IsForeground(i, j, k);
    return _Target->IsForeground(i, j, k) && _Source->IsForeground(i, j, k);
  }
}


} // namespace mirtk

#endif // MIRTK_ImageSimilarity_H

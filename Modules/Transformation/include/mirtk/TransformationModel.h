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

#ifndef MIRTK_TransformationModel_H
#define MIRTK_TransformationModel_H

#include "mirtk/String.h"
#include "mirtk/Array.h"
#include "mirtk/ImageAttributes.h"
#include "mirtk/TransformationType.h"


namespace mirtk {


// -----------------------------------------------------------------------------
/// Enumeration of transformation models
///
/// A transformation model is implemented by one or more transformation
/// classes and futhermore may be used to determine additional registration
/// settings, such as hard and soft transformation constraints. Thus, the
/// transformation model differs in semantics from the transformation type.
///
/// \sa TransformationType, ToTransformationType
enum TransformationModel
{
  TM_Unknown,                       ///< Unknown/invalid transformation model
  // Add new enumeration values below
  TM_Rigid,                         ///< Linear transformation with up to  6 DoFs (rotate, translate)
  TM_Similarity,                    ///< Linear transformation with up to  7 DoFs (rotate, translate, global scale)
  TM_Affine,                        ///< Linear transformation with up to 12 DoFs (rotate, translate, scale, skew)
  TM_LinearFFD,                     ///< Displacement field with linear interpolation
  TM_BSplineFFD,                    ///< Displacement field with B-spline interpolation
  TM_BSplineStatFFD,                ///< Displacement field with B-spline interpolation using a statistical model
  TM_BSplineSVFFD,                  ///< Stationary velocity field with B-spline interpolation
  TM_BSplineTDFFD,                  ///< Non-stationary velocity field with B-spline interpolation
  // Add new enumeration values above
  TM_Last                           ///< Number of available transformation models + 1
};

// -----------------------------------------------------------------------------
template <>
inline string ToString(const TransformationModel &m, int w, char c, bool left)
{
  const char *str;
  switch (m) {
    case TM_Rigid:           str = "Rigid"; break;
    case TM_Similarity:      str = "Similarity"; break;
    case TM_Affine:          str = "Affine"; break;
    case TM_LinearFFD:       str = "LinearFFD"; break;
    case TM_BSplineFFD:      str = "BSplineFFD"; break;
    case TM_BSplineStatFFD:  str = "BSplineStatFFD"; break;
    case TM_BSplineSVFFD:    str = "BSplineSVFFD"; break;
    case TM_BSplineTDFFD:    str = "BSplineTDFFD"; break;
    default:                 str = "Unknown"; break;
  }
  return ToString(str, w, c, left);
}

// -----------------------------------------------------------------------------
inline string ToPrettyString(const TransformationModel &m)
{
  switch (m) {
    case TM_Rigid:           return "rigid transformation";
    case TM_Similarity:      return "similarity transformation";
    case TM_Affine:          return "affine transformation";
    case TM_LinearFFD:       return "non-parametric displacement field";
    case TM_BSplineFFD:      return "free-form deformation";
    case TM_BSplineStatFFD:  return "statistical free-form deformation";
    case TM_BSplineSVFFD:    return "parametric stationary velocity field transformation";
    case TM_BSplineTDFFD:    return "temporal diffeomorphic free-form deformation";
    default:                 return "unknown transformation";
  }
}

// -----------------------------------------------------------------------------
template <>
inline bool FromString(const char *str, TransformationModel &m)
{
  m = TM_Unknown;
  // Alternative names of transformation models
  if      (strcmp(str, "FFD")   == 0) m = TM_BSplineFFD;
  else if (strcmp(str, "SVFFD") == 0) m = TM_BSplineSVFFD;
  else if (strcmp(str, "TDFFD") == 0) m = TM_BSplineTDFFD;
  // Default names of transformation models
  if (m == TM_Unknown) {
    m = static_cast<TransformationModel>(TM_Last - 1);
    while (m != TM_Unknown) {
      if (ToString(m) == str) break;
      m = static_cast<TransformationModel>(m - 1);
    }
  }
  return (m != TM_Unknown);
}

// -----------------------------------------------------------------------------
/// Whether a given transformation model is linear
inline bool IsLinear(TransformationModel model)
{
  return model == TM_Rigid || model == TM_Similarity || model == TM_Affine;
}

// -----------------------------------------------------------------------------
/// Whether any given transformation model is linear
inline bool IsLinear(const Array<TransformationModel> &model)
{
  for (size_t i = 0; i < model.size(); ++i) {
    if (IsLinear(model[i])) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Whether any given transformation model is non-linear
inline bool IsNonLinear(const Array<TransformationModel> &model)
{
  for (size_t i = 0; i < model.size(); ++i) {
    if (!IsLinear(model[i])) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Whether a given transformation model is a FFD with a linear interpolation kernel
inline bool IsLinearFFD(TransformationModel model)
{
  return model == TM_LinearFFD;
}

// -----------------------------------------------------------------------------
/// Whether any given transformation model is a FFD with a linear interpolation kernel
inline bool IsLinearFFD(const Array<TransformationModel> &model)
{
  for (size_t i = 0; i < model.size(); ++i) {
    if (IsLinearFFD(model[i])) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Whether a given transformation model is diffeomorphic (velocity parameterization)
inline bool IsDiffeo(TransformationModel model)
{
  return IsLinear(model) || model == TM_BSplineSVFFD || model == TM_BSplineTDFFD;
}

// -----------------------------------------------------------------------------
/// Whether any given transformation model is a diffeomorphic model
inline bool IsDiffeo(const Array<TransformationModel> &model)
{
  for (size_t i = 0; i < model.size(); ++i) {
    if (IsDiffeo(model[i])) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Whether a given transformation model is 3D+t
inline bool IsSpatioTemporal(TransformationModel model)
{
  return model == TM_BSplineTDFFD || model == TM_BSplineSVFFD;
}

// -----------------------------------------------------------------------------
/// Whether any given transformation model is 3D+t
inline bool IsSpatioTemporal(const Array<TransformationModel> &model)
{
  for (size_t i = 0; i < model.size(); ++i) {
    if (IsSpatioTemporal(model[i])) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Get type of (default) transformation which implements a specific model
///
/// The mapping from transformation model to transformtion type is not
/// one-to-one. More then one transformation can be suitable for a
/// transformation model. This function defines the default type used
/// for each model. The base registration filter implementations make use
/// of it, but a specialized registration filter can choose another
/// transformation for a given model if desired.
///
/// For example, see ImageRegistrationFilter::TransformationType.
inline TransformationType
ToTransformationType(TransformationModel    model,
                     const ImageAttributes &domain)
{
  // 2D/3D
  if (domain._t == 1) {
    switch (model) {
      case TM_LinearFFD:     return TRANSFORMATION_LINEAR_FFD_3D_v3;
      case TM_BSplineFFD:    return TRANSFORMATION_BSPLINE_FFD_3D_v3;
      default: break;
    }
  // 4D
  } else {
    switch (model) {
      case TM_LinearFFD:     return TRANSFORMATION_LINEAR_FFD_4D_v2;
      case TM_BSplineFFD:    return TRANSFORMATION_BSPLINE_FFD_4D_v2;
      default: break;
    }
  }
  // nD
  switch (model) {
    case TM_Rigid:           return TRANSFORMATION_RIGID;
    case TM_Similarity:      return TRANSFORMATION_SIMILARITY;
    case TM_Affine:          return TRANSFORMATION_AFFINE;
    case TM_BSplineStatFFD:  return TRANSFORMATION_BSPLINE_FFD_STATISTICAL;
    case TM_BSplineSVFFD:    return TRANSFORMATION_BSPLINE_FFD_SV_v5;
    case TM_BSplineTDFFD:    return TRANSFORMATION_BSPLINE_FFD_TD_v3;
    default:                 return TRANSFORMATION_UNKNOWN;
  }
}

// -----------------------------------------------------------------------------
/// Enumeration of available multi-level transformation modes
enum MFFDMode
{
  MFFD_Default,  ///< Choose suitable default multi-level transformation model
  MFFD_None,     ///< Use single transformation without additional global or local transformations
  MFFD_Sum,      ///< One transformation for each resolution level with additive composition
  MFFD_Fluid,    ///< One transformation for each resolution level with fluid composition
  MFFD_LogSum    ///< Additive multi-level stationary velocity field
};

// -----------------------------------------------------------------------------
template <>
inline string ToString(const MFFDMode &m, int w, char c, bool left)
{
  const char *str;
  switch (m) {
    case MFFD_Default:  str = "Default"; break;
    case MFFD_None:     str = "None"; break;
    case MFFD_Sum:      str = "Sum"; break;
    case MFFD_LogSum:   str = "LogSum"; break;
    case MFFD_Fluid:    str = "Fluid"; break;
    default:            str = "Unknown"; break;
  }
  return ToString(str, w, c, left);
}

// -----------------------------------------------------------------------------
template <>
inline bool FromString(const char *str, MFFDMode &m)
{
  if      (strcmp(str, "Default")  == 0) m = MFFD_Default;
  else if (strcmp(str, "None")     == 0) m = MFFD_None;
  else if (strcmp(str, "Sum")      == 0) m = MFFD_Sum;
  else if (strcmp(str, "LogSum")   == 0) m = MFFD_LogSum;
  else if (strcmp(str, "Fluid")    == 0) m = MFFD_Fluid;
  else return false;
  return true;
}


} // namespace mirtk

#endif // MIRTK_RegistrationFilter_H

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

#ifndef MIRTK_EnergyMeasure_H
#define MIRTK_EnergyMeasure_H

#include "mirtk/String.h"


namespace mirtk {


// -----------------------------------------------------------------------------
/// Enumeration of all available energy terms
enum EnergyMeasure
{
  EM_Unknown, ///< Unknown/invalid energy term
  // Add new enumeration values below
  // ---------------------------------------------------------------------------
  // Image (dis-)similarity measures (cf. ImageSimilarity subclasses)
  SIM_Begin,

    EM_JE,                      ///< Joint entropy
    EM_CC,                      ///< Cross-correlation
    EM_MI,                      ///< Mutual information
    EM_NMI,                     ///< Normalized mutual information
    EM_SSD,                     ///< Sum of squared differences
    EM_CR_XY,                   ///< Correlation ratio
    EM_CR_YX,                   ///< Correlation ratio
    EM_LC,                      ///< Label consistency
    EM_K,                       ///< Kappa statistic
    EM_ML,                      ///< Maximum likelihood
    EM_NGF_COS,                 ///< Cosine of normalzed gradient field
    EM_NCC,                     ///< Normalized cross-correlation
    EM_LNCC = EM_NCC,           ///< Local normalized cross-correlation
    EM_CoVar,                   ///< Covariance
    EM_PSNR,                    ///< Peak signal-to-noise ratio

  SIM_End,
  // ---------------------------------------------------------------------------
  // Point set distance measures (cf. PointSetDistance subclasses)
  PDM_Begin,

    EM_FRE,                     ///< Fiducial registration error (FRE) measure
    EM_CorrespondenceDistance,  ///< Point correspondence distance measure
    EM_CurrentsDistance,        ///< Distance measure based on currents representation
    EM_VarifoldDistance,        ///< Distance measure based on varifold representation

  PDM_End,
  // ---------------------------------------------------------------------------
  // External point set forces (cf. ExternalForce subclasses)
  EFT_Begin,

    EM_BalloonForce,            ///< Balloon/inflation force
    EM_ImageEdgeForce,          ///< Image edge force
    EM_ImageEdgeDistance,       ///< Image edge distance
    EM_ImplicitSurfaceDistance, ///< Implicit surface distance force

  EFT_End,
  // ---------------------------------------------------------------------------
  // Internal point set forces (cf. InternalForce subclasses)
  IFT_Begin,

    EM_MetricDistortion,        ///< Minimize metric distortion
    EM_Stretching,              ///< Stretching force (rest edge length)
    EM_Curvature,               ///< Minimize curvature of point set surface
    EM_QuadraticCurvature,      ///< Quadratic fit of neighor to tangent plane distance
    EM_GaussCurvature,          ///< Gauss curvature constraint
    EM_MeanCurvature,           ///< Mean curvature constraint
    EM_MaximumCurvature,        ///< Maximum curvature constraint
    EM_NonSelfIntersection,     ///< Repels too close non-neighboring triangles
    EM_RepulsiveForce,          ///< Repels too close non-neighboring nodes
    EM_InflationForce,          ///< Inflate point set surface
    EM_SpringForce,             ///< Spring force
    EM_NormalForce,             ///< Constant force in normal direction

  IFT_End,
  // ---------------------------------------------------------------------------
  // Transformation regularization terms (cf. TransformationConstraint subclasses)
  CM_Begin,

    EM_VolumePreservation,      ///< Volume preservation constraint
    EM_TopologyPreservation,    ///< Topology preservation constraint
    EM_Sparsity,                ///< Default sparsity constraint
    EM_BendingEnergy,           ///< Thin-plate spline bending energy
    EM_LinearElasticity,        ///< Linear elastic energy
    EM_L0Norm,                  ///< Sparsity constraint based on l0-norm
    EM_L1Norm,                  ///< Sparsity constraint based on l1-norm
    EM_L2Norm,                  ///< Sparsity constraint based on l2-norm
    EM_SqLogDetJac,             ///< Squared logarithm of the Jacobian determinant
    EM_NegDetJac,               ///< Penalise negative Jacobian determinant

  CM_End,

  // ---------------------------------------------------------------------------
  // Others
  EM_MeanSquaredDisplacementError,  ///< Mean squared deviation from given deformation

  // ---------------------------------------------------------------------------
  // Add new enumeration values above
  EM_Last ///< Number of enumeration values + 1
};

// -----------------------------------------------------------------------------
/// Convert energy measure enumeration value to string
template <>
inline string ToString(const EnergyMeasure &value, int w, char c, bool left)
{
  const char *str;
  switch (value) {
    // ---------------------------------------------------------------------------
    // Image (dis-)similarity measures
    case EM_JE:      str = "JE"; break;
    case EM_CC:      str = "CC"; break;
    case EM_MI:      str = "MI"; break;
    case EM_NMI:     str = "NMI"; break;
    case EM_SSD:     str = "SSD"; break;
    case EM_CR_XY:   str = "CR_XY"; break;
    case EM_CR_YX:   str = "CR_YX"; break;
    case EM_LC:      str = "LC"; break;
    case EM_K:       str = "K"; break;
    case EM_ML:      str = "ML"; break;
    case EM_NGF_COS: str = "NGF_COS"; break;
    case EM_NCC:     str = "NCC"; break;
    case EM_CoVar:   str = "CoVar"; break;
    case EM_PSNR:    str = "PSNR"; break;

    // ---------------------------------------------------------------------------
    // Point set distance measures
    case EM_FRE:                    str = "FRE"; break;
    case EM_CorrespondenceDistance: str = "PCD"; break;
    case EM_CurrentsDistance:       str = "CurrentsDistance"; break;
    case EM_VarifoldDistance:       str = "VarifoldDistance"; break;

    // ---------------------------------------------------------------------------
    // External point set forces
    case EM_BalloonForce:            str = "BalloonForce"; break;
    case EM_ImageEdgeForce:          str = "ImageEdgeForce"; break;
    case EM_ImageEdgeDistance:       str = "ImageEdgeDistance"; break;
    case EM_ImplicitSurfaceDistance: str = "ImplicitSurfaceDistance"; break;

    // ---------------------------------------------------------------------------
    // Internal point set forces
    case EM_MetricDistortion:    str = "MetricDistortion"; break;
    case EM_Stretching:          str = "Stretching"; break;
    case EM_Curvature:           str = "Curvature"; break;
    case EM_QuadraticCurvature:  str = "QuadraticCurvature"; break;
    case EM_GaussCurvature:      str = "GaussCurvature"; break;
    case EM_MeanCurvature:       str = "MeanCurvature"; break;
    case EM_MaximumCurvature  :  str = "MaximumCurvature"; break;
    case EM_NonSelfIntersection: str = "NSI"; break;
    case EM_RepulsiveForce:      str = "Repulsion"; break;
    case EM_InflationForce:      str = "Inflation"; break;
    case EM_SpringForce:         str = "Spring"; break;
    case EM_NormalForce:         str = "NormalForce"; break;

    // -------------------------------------------------------------------------
    // Transformation constraints
    case EM_BendingEnergy:        str = "BE"; break;
    case EM_LinearElasticity:     str = "LE"; break;
    case EM_VolumePreservation:   str = "VP"; break;
    case EM_TopologyPreservation: str = "TP"; break;
    case EM_Sparsity:             str = "Sparsity"; break;
    case EM_L0Norm:               str = "L0"; break;
    case EM_L1Norm:               str = "L1"; break;
    case EM_L2Norm:               str = "L2"; break;
    case EM_SqLogDetJac:          str = "SqLogDetJac"; break;
    case EM_NegDetJac:            str = "NegDetJac"; break;

    // -------------------------------------------------------------------------
    // Others
    case EM_MeanSquaredDisplacementError: str = "MSDE"; break;

    // ---------------------------------------------------------------------------
    // Unknown/invalid enumeration value
    default: str = "Unknown";
  }
  return ToString(str, w, c, left);
}

// -----------------------------------------------------------------------------
/// Convert energy measure enumeration value to human-friendly descriptive string
inline string ToPrettyString(const EnergyMeasure &value, int w = 0, char c = ' ', bool left = true)
{
  const char *str;
  switch (value) {
    // ---------------------------------------------------------------------------
    // Image (dis-)similarity measures
    case EM_JE:      str = "Joint entropy"; break;
    case EM_CC:      str = "Normalized cross-correlation"; break;
    case EM_MI:      str = "Mutual information"; break;
    case EM_NMI:     str = "Normalized mutual information"; break;
    case EM_SSD:     str = "Sum of squared differences"; break;
    case EM_CR_XY:   str = "Correlation ratio C(X|Y)"; break;
    case EM_CR_YX:   str = "Correlation ratio C(Y|X)"; break;
    case EM_LC:      str = "Label consistency"; break;
    case EM_K:       str = "Kappa statistic"; break;
    case EM_ML:      str = "Maximum likelihood"; break;
    case EM_NGF_COS: str = "Cosine of normalized gradient field"; break;
    case EM_NCC:     str = "(Local) Normalized cross-correlation"; break;
    case EM_CoVar:   str = "Covariance"; break;
    case EM_PSNR:    str = "Peak signal-to-noise ratio"; break;

    // ---------------------------------------------------------------------------
    // Point set distance measures
    case EM_FRE:                    str = "Fiducial registration error"; break;
    case EM_CorrespondenceDistance: str = "Point correspondence distance"; break;
    case EM_CurrentsDistance:       str = "Currents distance"; break;
    case EM_VarifoldDistance:       str = "Varifold distance"; break;

    // ---------------------------------------------------------------------------
    // External point set forces
    case EM_BalloonForce:               str = "Balloon force"; break;
    case EM_ImageEdgeForce:             str = "Image edge force"; break;
    case EM_ImageEdgeDistance:          str = "Image edge distance"; break;
    case EM_ImplicitSurfaceDistance:    str = "Implicit surface distance"; break;

    // ---------------------------------------------------------------------------
    // Internal point set forces
    case EM_MetricDistortion:    str = "Metric distortion"; break;
    case EM_Stretching:          str = "Stretching"; break;
    case EM_Curvature:           str = "Curvature"; break;
    case EM_QuadraticCurvature:  str = "Quadratic curvature"; break;
    case EM_GaussCurvature:      str = "Gauss curvature"; break;
    case EM_MeanCurvature:       str = "Mean curvature"; break;
    case EM_MaximumCurvature:    str = "Maximum curvature"; break;
    case EM_NonSelfIntersection: str = "Non-self intersection"; break;
    case EM_RepulsiveForce:      str = "Repulsion"; break;
    case EM_InflationForce:      str = "Inflation"; break;
    case EM_SpringForce:         str = "Spring"; break;
    case EM_NormalForce:         str = "Normal"; break;

    // -------------------------------------------------------------------------
    // Transformation constraints
    case EM_BendingEnergy:        str = "Bending energy"; break;
    case EM_LinearElasticity:     str = "Linear elasticity"; break;
    case EM_VolumePreservation:   str = "Volume preservation"; break;
    case EM_TopologyPreservation: str = "Topology preservation"; break;
    case EM_Sparsity:             str = "Sparsity constraint"; break;
    case EM_L0Norm:               str = "l0 norm"; break;
    case EM_L1Norm:               str = "l1 norm"; break;
    case EM_L2Norm:               str = "l2 norm"; break;
    case EM_SqLogDetJac:          str = "Squared logarithm of Jacobian determinant"; break;
    case EM_NegDetJac:            str = "Negative Jacobian determinant penalty"; break;

    // -------------------------------------------------------------------------
    // Others
    case EM_MeanSquaredDisplacementError: str = "Mean squared displacement error"; break;

    // ---------------------------------------------------------------------------
    // Unknown/invalid enumeration value
    default: str = "Unknown energy term";
  }
  return ToString(str, w, c, left);
}

// -----------------------------------------------------------------------------
/// Convert energy measure string to enumeration value
template <>
inline bool FromString(const char *str, EnergyMeasure &value)
{
  const string lstr = ToLower(Trim(str));

  value = EM_Unknown;

  // ---------------------------------------------------------------------------
  // Alternative names for image (dis-)similarity measures
  if (value == EM_Unknown) {
    if      (lstr == "mse") value = EM_SSD;
    else if (lstr == "meansquarederror") value = EM_SSD;
    else if (lstr == "mean squared error") value = EM_SSD;
    else if (lstr == "lcc") value = EM_NCC;
    else if (lstr == "lncc") value = EM_NCC;
    else if (lstr == "kappa") value = EM_K;
    else if (lstr == "correlation ratio xy" ||
             lstr == "correlationratioxy") value = EM_CR_XY;
    else if (lstr == "correlation ratio" ||
             lstr == "correlationratio" ||
             lstr == "correlation ratio yx" ||
             lstr == "correlationratioyx") value = EM_CR_YX;
  }

  // ---------------------------------------------------------------------------
  // Alternative names for point set distance measures
  if (value == EM_Unknown) {
    if (lstr == "fiducial error" || lstr == "fiducialerror" ||
        lstr == "landmark error" || lstr == "landmarkerror" ||
        lstr == "landmark registration error" || lstr == "landmarkregistrationerror") {
      value = EM_FRE;
    } else if (lstr == "correspondence distance" || lstr == "correspondencedistance") {
      value = EM_CorrespondenceDistance;
    }
  }

  // ---------------------------------------------------------------------------
  // Alternative names for external point set forces
  if (value == EM_Unknown) {
    if (lstr == "edge force"    || lstr == "edgeforce")    value = EM_ImageEdgeForce;
    if (lstr == "edge distance" || lstr == "edgedistance") value = EM_ImageEdgeDistance;
  }

  // ---------------------------------------------------------------------------
  // Alternative names for internal point set forces
  if (value == EM_Unknown) {
    if      (strcmp(str, "EdgeLength")          == 0) value = EM_Stretching;
    else if (strcmp(str, "MetricDistortion")    == 0) value = EM_MetricDistortion;
    else if (strcmp(str, "Bending")             == 0) value = EM_Curvature;
    else if (strcmp(str, "SurfaceBending")      == 0) value = EM_Curvature;
    else if (strcmp(str, "SurfaceCurvature")    == 0) value = EM_Curvature;
    else if (strcmp(str, "GaussianCurvature")   == 0) value = EM_GaussCurvature;
    else if (strcmp(str, "RepulsiveForce")      == 0) value = EM_RepulsiveForce;
    else if (strcmp(str, "NonSelfIntersection") == 0) value = EM_NonSelfIntersection;
    else if (strcmp(str, "InflationForce")      == 0) value = EM_InflationForce;
    else if (strcmp(str, "SurfaceInflation")    == 0) value = EM_InflationForce;
  }

  // ---------------------------------------------------------------------------
  // Alternative names for transformation regularization terms
  if (value == EM_Unknown) {
    if      (lstr == "jac") value = EM_SqLogDetJac;
    else if (lstr == "logjac") value = EM_SqLogDetJac;
    else if (lstr == "negjac") value = EM_NegDetJac;
    else if (lstr == "elasticity") value = EM_LinearElasticity;
  }

  // ---------------------------------------------------------------------------
  // Convert default names of energy measures
  // (cf. ToString(EnergyMeasure) and ToPrettyString(EnergyMeasure))
  if (value == EM_Unknown) {
    string pretty;
    value = static_cast<EnergyMeasure>(EM_Last - 1);
    while (value != EM_Unknown) {
      if (lstr == ToLower(ToString(value))) break;
      pretty = ToLower(ToPrettyString(value));
      if (lstr == pretty || lstr == TrimAll(pretty, " -")) break;
      value = static_cast<EnergyMeasure>(value - 1);
    }
  }

  return (value != EM_Unknown);
}


} // namespace mirtk


namespace std {

/// Compute hash value from EnergyMeasure enumeration value
template<>
struct hash<mirtk::EnergyMeasure> {
    size_t operator()(const mirtk::EnergyMeasure &enum_value) const {
        return std::hash<int>()(enum_value);
    }
};


} // namespace std


#endif // MIRTK_EnergyMeasure_H

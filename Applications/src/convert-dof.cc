/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/IOConfig.h"
#include "mirtk/GenericImage.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/Transformations.h"

#if MIRTK_IO_WITH_NIfTI
#  include "mirtk/NiftiImageInfo.h"
#  include "mirtk/NiftiImageReader.h"
#endif // MIRTK_IO_WITH_NIfTI


using namespace mirtk;


////////////////////////////////////////////////////////////////////////////////
// Help
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <input> <output> [options]\n";
  cout << "       " << name << " <input> <output> -input-format  image [options] (3D+3 input  image)\n";
  cout << "       " << name << " <input> <output> -output-format image [options] (3D+3 output image)\n";
  cout << "       " << name << " <dx> <dy> <dz> <output> -input-format  image [options] (3x3D input  image)\n";
  cout << "       " << name << " <input> <dx> <dy> <dz>  -output-format image [options] (3x3D output image)\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Converts between different transformation file formats:\n";
  cout << "\n";
  cout << "  =========================   =========================================================================================\n";
  cout << "  unknown                     Unknown, try to guess it from file header/type.\n";
  cout << "  disp_world|disp|image       Dense displacement field image with world space displacement vectors [mm].\n";
  cout << "  disp_voxel                  Dense displacement field image with target space displacement vectors [voxel].\n";
  cout << "  disp_itk|itk_disp           Dense displacement field image in ITK format [mm].\n";
  cout << "  svf_world|svf               Dense velocity field image with world space vectors [mm]. Requires input SVFFD.\n";
  cout << "  svf_voxel                   Dense velocity field image with target space vectors [voxel]. Requires input SVFFD.\n";
  cout << "  svf_itk|itk_svf             Dense velocity field image in ITK format (e.g., LogDemons) [mm]\n";
  cout << "  mirtk                       MIRTK transformation file format.\n";
  cout << "  mirtk_rigid|rigid           Rigid MIRTK transformation file format (6 DoFs).\n";
  cout << "  mirtk_similarity            Similarity MIRTK transformation file format (7 DoFs).\n";
  cout << "  mirtk_affine|affine         Affine MIRTK transformation file format (12 DoFs).\n";
  cout << "  mirtk_linear_ffd            Linear free-form deformation.\n";
  cout << "  mirtk_linear_svffd          Linear free-form deformation parameterized by stationary velocity field.\n";
  cout << "  mirtk_linear_tdffd          Linear free-form deformation parameterized by non-stationary velocity field.\n";
  cout << "  mirtk_bspline_ffd           Cubic B-spline free-form deformation.\n";
  cout << "  mirtk_bspline_svffd         Cubic B-spline free-form deformation parameterized by stationary velocity field.\n";
  cout << "  mirtk_bspline_tdffd         Cubic B-spline free-form deformation parameterized by non-stationary velocity field.\n";
  cout << "  irtk                        IRTK transformation file format.\n";
  cout << "  irtk_rigid                  Rigid IRTK transformation file format (6 DoFs).\n";
  cout << "  irtk_affine                 Affine IRTK transformation file format (12 DoFs).\n";
  cout << "  irtk_linear_ffd             Linear IRTK free-form deformation.\n";
  cout << "  irtk_bspline_ffd            Cubic B-spline IRTK free-form deformation.\n";
  cout << "  mni_xfm|xfm                 Linear FreeSurfer transformation (.xfm file).\n";
  cout << "  fsl                         Guess/choose FSL output file format.\n";
  cout << "  flirt                       FSL FLIRT output file format.\n";
  cout << "  fnirt                       FSL FNIRT output file format.\n";
  cout << "  nreg                        Guess/choose Nifty Reg transformation output file format.\n";
  cout << "  aladin                      Nifty Reg Aladin output file format.\n";
  cout << "  f3d                         Nifty Reg reg_f3d output file format with nifti1.intent_p1 code.\n";
  cout << "  f3d_def_field               Nifty Reg reg_f3d output image deformation  field.\n";
  cout << "  f3d_disp_field              Nifty Reg reg_f3d output image displacement field.\n";
  cout << "  f3d_spline_grid             Nifty Reg reg_f3d output control point displacement field.\n";
  cout << "  f3d_def_vel_field           Nifty Reg reg_f3d output image deformation  field as stationary velocity field.\n";
  cout << "  f3d_disp_vel_field          Nifty Reg reg_f3d output image displacement field as stationary velocity field.\n";
  cout << "  f3d_spline_vel_grid         Nifty Reg reg_f3d output control point velocity field.\n";
  cout << "  elastix|elastix_ffd         Elastix BSplineTransform parameters text file.\n";
  cout << "  dramms                      DRAMMS deformation field.\n";
  cout << "  star_ccm                    Output file suitable for import in STAR CCM+.\n";
  cout << "  star_ccm_table              Point displacements as STAR CCM+ XYZ Table.\n";
  cout << "  star_ccm_table_xyz          Transformed points  as STAR CCM+ XYZ Table.\n";
  cout << "  table|csv|tsv               ASCII table of target point coordinates and displacement vectors.\n";
  cout << "  table_xyz|csv_xyz|tsv_xyz   ASCII table of transformed point coordinates.\n";
  cout << "  =========================   =========================================================================================\n";
#if !MIRTK_IO_WITH_NIfTI
  cout << "\n";
  cout << "  Cannot convert from/to the following formats because the I/O module is missing NIfTI support:\n";
  cout << "\n";
  cout << "    f3d*, fnirt\n";
  cout << "\n";
#endif // !MIRTK_IO_WITH_NIfTI
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Input  transformation file.\n";
  cout << "  output   Output transformation file.\n";
  cout << "\n";
  cout << "Options:\n";
  cout << "  -input-format  <format>    Format of input file. (default: unknown)\n";
  cout << "  -output-format <format>    Format of output file. (default: .nii[.gz] -> disp_world, .dof[.gz] -> mirtk)\n";
  cout << "  -format <format>           Short for :option:`-output-format`.\n";
  cout << "  -dofin <fname>             Affine transformation component in MIRTK format which\n";
  cout << "                             shall be removed from Nifty Reg's FFD such that the\n";
  cout << "                             resulting transformation is a MIRTK MFFD consisting\n";
  cout << "                             of this affine component plus the remaining non-linear\n";
  cout << "                             component of the input FFD. (default: none)\n";
  cout << "  -target <fname>            Target image. Required for from/to FSL format conversion.\n";
  cout << "                             Also required when converting linear transformations to displacment\n";
  cout << "                             fields or (M)FFD transformation types. (default: none)\n";
  cout << "  -source <fname>            Source image. Required for from/to FSL format conversion. (default: none)\n";
  cout << "  -points <fname>            Input point set. Used for CSV/TSV and STAR-CCM+ output of transformed\n";
  cout << "                             points. By default, all target or FFD lattice points are transformed.\n";
  cout << "  -Tt <time>                 Time point of target image. Used by 3D+t, TD, and SV FFDs.\n";
  cout << "  -Ts <time>                 Time point of source image. Used by 3D+t, TD, and SV FFDs.\n";
  cout << "  -ds <value>                Output control point spacing. (default: input spacing)\n";
  cout << "  -dx <value>                Output control point spacing in x dimension. (default: input spacing)\n";
  cout << "  -dy <value>                Output control point spacing in y dimension. (default: input spacing)\n";
  cout << "  -dz <value>                Output control point spacing in z dimension. (default: input spacing)\n";
  cout << "  -dt <value>                Temporal sampling used for CSV/TSV and STAR-CCM+ output. (default: input spacing)\n";
  cout << "  -t1, -tmin <value>         Lower time interval limit for output of CSV/TSV and STAR-CCM+ table. (default: -inf)\n";
  cout << "  -t2, -tmax <value>         Upper time interval limit for output of CSV/TSV and STAR-CCM+ table. (default: +inf)\n";
  cout << "  -scaling-steps <int>       Number of scaling and squaring steps.\n";
  cout << "                             In case of f3d_*_vel_*, use nifti1.intent_p2 by default.\n";
  cout << "  -xyz_units (m|mm|mu)       Spatial units of original target NIfTI header\n";
  cout << "                             if ignored by Nifty Reg's reg_f3d. (default: mm)\n";
  cout << "  -delimiter <string>        Delimiting character sequence to use for CSV/TSV or STAR-CCM+ Table output.\n";
  cout << "  -precision <int>           Number of decimal digits for ASCII output of floating point values. (default: 5)\n";
  PrintStandardOptions(cout);
  cout << endl;
}

////////////////////////////////////////////////////////////////////////////////
// Auxiliaries
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
/// Open file with the specified mode
FILE *OpenFile(const char *fname, const char *mode)
{
  FILE *fp = nullptr;
#ifdef WINDOWS
  if (fopen_s(&fp, fname, mode) != 0) fp = nullptr;
#else
  fp = fopen(fname, mode);
#endif
  return fp;
}

////////////////////////////////////////////////////////////////////////////////
// Types
////////////////////////////////////////////////////////////////////////////////

namespace mirtk {


// -----------------------------------------------------------------------------
/// Enumeration of supported transformation file formats/packages
enum TransformationFileFormat
{
  Format_Unknown,
  // Dense displacement field image
  Format_WorldDisplacement, ///< Displacement vectors in world units (mm)
  Format_VoxelDisplacement, ///< Displacement vectors in voxel units
  Format_WorldSVF,          ///< Stationary velocity field in world units (mm)
  Format_VoxelSVF,          ///< Stationary velocity field in voxel units
  // ITK
  Format_ITKDisplacement,
  Format_ITKSVF,
  // MIRTK
  Format_MIRTK,
  Format_MIRTK_Rigid,
  Format_MIRTK_Similarity,
  Format_MIRTK_Affine,
  Format_MIRTK_LinearFFD,
  Format_MIRTK_LinearSVFFD,
  Format_MIRTK_LinearTDFFD,
  Format_MIRTK_BSplineFFD,
  Format_MIRTK_BSplineSVFFD,
  Format_MIRTK_BSplineTDFFD,
  // IRTK
  Format_IRTK,
  Format_IRTK_Rigid,
  Format_IRTK_Affine,
  Format_IRTK_LinearFFD,
  Format_IRTK_BSplineFFD,
  // FSL
  Format_FSL,
  Format_FSL_FLIRT,
  Format_FSL_FNIRT_Displacement,
  Format_FSL_WarpAbsolute,
  Format_FSL_WarpRelative,
  Format_FSL_DctCoefficients,
  Format_FSL_CubicSpline,
  Format_FSL_QuadraticSpline,
  // Nifty Reg
  Format_NREG,
  Format_Aladin,
  Format_F3D,
  Format_F3D_DEF_FIELD,
  Format_F3D_DISP_FIELD,
  Format_F3D_SPLINE_GRID,
  Format_F3D_DEF_VEL_FIELD,
  Format_F3D_DISP_VEL_FIELD,
  Format_F3D_SPLINE_VEL_GRID,
  // Elastix
  Format_Elastix,
  Format_Elastix_FFD,
  // FreeSurfer
  Format_MNI,
  Format_MNI_XFM,
  Format_MNI_M3Z,
  // DRAMMS
  Format_DRAMMS,
  // STAR-CCM+
  Format_STAR_CCM,
  Format_STAR_CCM_Table,
  Format_STAR_CCM_Table_XYZ,
  // Table
  Format_CSV,
  Format_CSV_XYZ,
  Format_TSV,
  Format_TSV_XYZ,
  // Last enumeration entry
  Format_Last
};

// -----------------------------------------------------------------------------
/// Convert transformation file format enumeration value to string
template <>
inline string ToString(const TransformationFileFormat &format, int w, char c, bool left)
{
  const char *str;
  switch (format) {
    case Format_WorldDisplacement:        str = "disp_world"; break;
    case Format_VoxelDisplacement:        str = "disp_voxel"; break;
    case Format_WorldSVF:                 str = "svf_world"; break;
    case Format_VoxelSVF:                 str = "svf_voxel"; break;
    case Format_ITKDisplacement:          str = "itk_disp"; break;
    case Format_ITKSVF:                   str = "itk_svf"; break;
    case Format_MIRTK:                    str = "mirtk"; break;
    case Format_MIRTK_Rigid:              str = "mirtk_rigid"; break;
    case Format_MIRTK_Similarity:         str = "mirtk_similarity"; break;
    case Format_MIRTK_Affine:             str = "mirtk_affine"; break;
    case Format_MIRTK_LinearFFD:          str = "mirtk_linear_ffd"; break;
    case Format_MIRTK_LinearSVFFD:        str = "mirtk_linear_svffd"; break;
    case Format_MIRTK_LinearTDFFD:        str = "mirtk_linear_tdffd"; break;
    case Format_MIRTK_BSplineFFD:         str = "mirtk_bspline_ffd"; break;
    case Format_MIRTK_BSplineSVFFD:       str = "mirtk_bspline_svffd"; break;
    case Format_MIRTK_BSplineTDFFD:       str = "mirtk_bspline_tdffd"; break;
    case Format_IRTK:                     str = "irtk"; break;
    case Format_IRTK_Rigid:               str = "irtk_rigid"; break;
    case Format_IRTK_Affine:              str = "irtk_affine"; break;
    case Format_IRTK_LinearFFD:           str = "irtk_linear_ffd"; break;
    case Format_IRTK_BSplineFFD:          str = "irtk_bspline_ffd"; break;
    case Format_FSL:                      str = "fsl"; break;
    case Format_FSL_FLIRT:                str = "firt"; break;
    case Format_FSL_FNIRT_Displacement:   str = "fnirt_disp"; break;
    case Format_FSL_WarpAbsolute:         str = "fsl_warp_abs"; break;
    case Format_FSL_WarpRelative:         str = "fsl_warp_rel"; break;
    case Format_FSL_DctCoefficients:      str = "fsl_dct_coeff"; break;
    case Format_FSL_CubicSpline:          str = "fsl_cubic_spline"; break;
    case Format_FSL_QuadraticSpline:      str = "fsl_quadratic_spline"; break;
    case Format_NREG:                     str = "nreg"; break;
    case Format_Aladin:                   str = "aladin"; break;
    case Format_F3D:                      str = "f3d"; break;
    case Format_F3D_DEF_FIELD:            str = "f3d_def_field"; break;
    case Format_F3D_DISP_FIELD:           str = "f3d_disp_field"; break;
    case Format_F3D_SPLINE_GRID:          str = "f3d_spline_grid"; break;
    case Format_F3D_DEF_VEL_FIELD:        str = "f3d_def_vel_field"; break;
    case Format_F3D_DISP_VEL_FIELD:       str = "f3d_disp_vel_field"; break;
    case Format_F3D_SPLINE_VEL_GRID:      str = "f3d_spline_vel_grid"; break;
    case Format_MNI:                      str = "mni"; break;
    case Format_MNI_XFM:                  str = "mni_xfm"; break;
    case Format_MNI_M3Z:                  str = "mni_m3z"; break;
    case Format_Elastix:                  str = "elastix"; break;
    case Format_Elastix_FFD:              str = "elastix_ffd"; break;
    case Format_DRAMMS:                   str = "dramms"; break;
    case Format_STAR_CCM:                 str = "star_ccm"; break;
    case Format_STAR_CCM_Table:           str = "star_ccm_table"; break;
    case Format_STAR_CCM_Table_XYZ:       str = "star_ccm_table_xyz"; break;
    case Format_CSV:                      str = "csv"; break;
    case Format_CSV_XYZ:                  str = "csv_xyz"; break;
    case Format_TSV:                      str = "tsv"; break;
    case Format_TSV_XYZ:                  str = "tsv_xyz"; break;
    default:                              str = "unknown"; break;
  }
  return ToString(str, w, c, left);
}

// -----------------------------------------------------------------------------
/// Convert transformation file format string to enumeration value
template <>
inline bool FromString(const char *str, TransformationFileFormat &format)
{
  string format_name = ToLower(Trim(str));
  format = Format_Unknown;

  // Alternative format names
  if      (format_name == "disp" || format_name == "image") format = Format_WorldDisplacement;
  else if (format_name == "svf") format = Format_WorldSVF;
  else if (format_name == "disp_itk") format = Format_ITKDisplacement;
  else if (format_name == "svf_itk") format = Format_ITKSVF;
  else if (format_name == "rigid") format = Format_MIRTK_Rigid;
  else if (format_name == "similarity") format = Format_MIRTK_Similarity;
  else if (format_name == "affine") format = Format_MIRTK_Affine;
  else if (format_name == "linear_ffd") format = Format_MIRTK_LinearFFD;
  else if (format_name == "linear_svffd") format = Format_MIRTK_LinearSVFFD;
  else if (format_name == "linear_tdffd") format = Format_MIRTK_LinearTDFFD;
  else if (format_name == "bspline_ffd") format = Format_MIRTK_BSplineFFD;
  else if (format_name == "bspline_svffd") format = Format_MIRTK_BSplineSVFFD;
  else if (format_name == "bspline_tdffd") format = Format_MIRTK_BSplineTDFFD;
  else if (format_name == "mirtk_ffd"   || format_name == "ffd")   format = Format_MIRTK_BSplineFFD;
  else if (format_name == "mirtk_svffd" || format_name == "svffd") format = Format_MIRTK_BSplineSVFFD;
  else if (format_name == "mirtk_tdffd" || format_name == "tdffd") format = Format_MIRTK_BSplineTDFFD;
  else if (format_name == "fsl_warp" || format_name == "warp") format = Format_FSL_WarpRelative;
  else if (format_name == "flirt") format = Format_FSL_FLIRT;
  else if (format_name == "fnirt") format = Format_FSL_FNIRT_Displacement;
  else if (format_name == "niftk" || format_name == "niftyreg") format = Format_NREG;
  else if (format_name == "freesurfer") format = Format_MNI;
  else if (format_name == "freesurfer_xfm" || format_name == "xfm") format = Format_MNI_XFM;
  else if (format_name == "freesurfer_m3z" || format_name == "m3z") format = Format_MNI_M3Z;
  else if (format_name.compare(0, 4, "reg_") == 0) format_name = format_name.substr(4);
  else if (format_name == "star-ccm" || format_name == "star-ccm+") format = Format_STAR_CCM;
  else if (format_name == "star-ccm table" || format_name == "star-ccm+ table") format = Format_STAR_CCM_Table;
  else if (format_name == "table") format = Format_TSV;
  else if (format_name == "table_xyz") format = Format_TSV_XYZ;

  // Default format names (cf. ToString(const TransformationFileFormat &))
  if (format == Format_Unknown) {
    format = static_cast<TransformationFileFormat>(Format_Last - 1);
    while (format != Format_Unknown) {
      if (ToString(format) == format_name) break;
      format = static_cast<TransformationFileFormat>(format - 1);
    }
  }

  return (format != Format_Unknown);
}

// -----------------------------------------------------------------------------
TransformationType ToMIRTKTransformationType(TransformationFileFormat format)
{
  switch (format) {
    case Format_MIRTK_Rigid:         return TRANSFORMATION_RIGID;
    case Format_MIRTK_Similarity:    return TRANSFORMATION_SIMILARITY;
    case Format_MIRTK_Affine:        return TRANSFORMATION_AFFINE;
    case Format_MIRTK_LinearFFD:     return TRANSFORMATION_LINEAR_FFD_3D;
    case Format_MIRTK_LinearSVFFD:   return TRANSFORMATION_LINEAR_FFD_SV;
    case Format_MIRTK_LinearTDFFD:   return TRANSFORMATION_LINEAR_FFD_TD;
    case Format_MIRTK_BSplineFFD:    return TRANSFORMATION_BSPLINE_FFD_3D;
    case Format_MIRTK_BSplineSVFFD:  return TRANSFORMATION_BSPLINE_FFD_SV;
    case Format_MIRTK_BSplineTDFFD:  return TRANSFORMATION_BSPLINE_FFD_TD;
    default:                         return TRANSFORMATION_UNKNOWN;
  }
}

// -----------------------------------------------------------------------------
/// NiftyReg transformation types
///
/// \sa _reg_maths.h of NiftyReg source code
enum F3DTransformationType
{
  F3D_TYPE_UNKNOWN    = 42,
  F3D_DEF_FIELD       = 0, ///< Deformation  field
  F3D_DISP_FIELD      = 1, ///< Displacement field
  F3D_SPLINE_GRID     = 2, ///< Deformation  field parameterized by (B-)spline
  F3D_DEF_VEL_FIELD   = 3, ///< Deformation  field parameterized by velocities
  F3D_DISP_VEL_FIELD  = 4, ///< Displacement field parameterized by velocities
  F3D_SPLINE_VEL_GRID = 5  ///< Deformation  field parameterized by (B-)spline velocities
};

// -----------------------------------------------------------------------------
F3DTransformationType ToF3DTransformationType(TransformationFileFormat format)
{
  switch (format) {
    case Format_F3D_DEF_FIELD:       return F3D_DEF_FIELD;
    case Format_F3D_DISP_FIELD:      return F3D_DISP_FIELD;
    case Format_F3D_SPLINE_GRID:     return F3D_SPLINE_GRID;
    case Format_F3D_DEF_VEL_FIELD:   return F3D_DEF_VEL_FIELD;
    case Format_F3D_DISP_VEL_FIELD:  return F3D_DISP_VEL_FIELD;
    case Format_F3D_SPLINE_VEL_GRID: return F3D_SPLINE_VEL_GRID;
    default:                         return F3D_TYPE_UNKNOWN;
  }
}

// -----------------------------------------------------------------------------
/// FSL NIfTI intent codes (see $FSLDIR/src/warpfns/fnirt_file_reader.h)
enum FSLIntentCode
{
  FSL_FNIRT_DISPLACEMENT_FIELD      = 2006,
  FSL_CUBIC_SPLINE_COEFFICIENTS     = 2007,
  FSL_DCT_COEFFICIENTS              = 2008,
  FSL_QUADRATIC_SPLINE_COEFFICIENTS = 2009
};

// -----------------------------------------------------------------------------
/// Integration parameters of output SV/TD FFD
struct FFDIMParams
{
  FFDIM  method;
  double t1;
  double t2;
  int    minsteps;
  int    maxsteps;
  double tol;

  FFDIMParams()
  :
    method(FFDIM_Unknown),
    t1(.0), t2(1.0),
    minsteps(0), maxsteps(0),
    tol(1.0e-3)
  {}

  double T() const
  {
    return (fequal(t1, t2) ? 1.0 : (t2 - t1));
  }

  int MinNumberOfSteps() const
  {
    int n;
    if (minsteps <= maxsteps) {
      n = (minsteps > 0 ? minsteps : (maxsteps > 0 ? maxsteps : 0));
    } else {
      n = (maxsteps > 0 ? maxsteps : (minsteps > 0 ? minsteps : 0));
    }
    return (n > 0 ? n : 10);
  }

  int MaxNumberOfSteps() const
  {
    int n;
    if (minsteps > maxsteps) {
      n = (minsteps > 0 ? minsteps : (maxsteps > 0 ? maxsteps : 0));
    } else {
      n = (maxsteps > 0 ? maxsteps : (minsteps > 0 ? minsteps : 0));
    }
    return (n > 0 ? n : 100);
  }

  double MinTimeStep() const
  {
    return T() / MaxNumberOfSteps();
  }

  double MaxTimeStep() const
  {
    return T() / MinNumberOfSteps();
  }

  double Tolerance() const
  {
    return (tol > .0 ? tol : .0);
  }
};

// -----------------------------------------------------------------------------
/// Parameters of Baker-Campbell-Hausdorff (BCH) method
///
/// These parameters are used for converting a transformation parameterized
/// by displacements into a transformation parameterized by velocities.
struct BCHParams
{
  int  nsteps;
  int  nterms;
  bool smooth;

  BCHParams() : nsteps(8), nterms(3), smooth(false) {}
};


} // namespace mirtk

////////////////////////////////////////////////////////////////////////////////
// Read transformation
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// MIRTK
// =============================================================================

// -----------------------------------------------------------------------------
/// Read transformation from MIRTK transformation file
Transformation *ReadMIRTK(const char *fname)
{
  return Transformation::New(fname);
}

// =============================================================================
// IRTK
// =============================================================================

// -----------------------------------------------------------------------------
/// Read transformation from IRTK transformation file
Transformation *ReadIRTK(const char *fname)
{
  return Transformation::New(fname);
}

// =============================================================================
// Image
// =============================================================================

// -----------------------------------------------------------------------------
/// Read displacement field from component images
GenericImage<double> ConcatenateWorldDisplacement(const char *dx_name, const char *dy_name, const char *dz_name)
{
  GenericImage<double> tmp(dx_name);
  GenericImage<double> disp(tmp.Attributes(), 3);
  for (int k = 0; k < disp.Z(); ++k)
  for (int j = 0; j < disp.Y(); ++j)
  for (int i = 0; i < disp.X(); ++i) {
    disp(i, j, k, 0) = tmp(i, j, k);
  }
  tmp.Read(dy_name);
  for (int k = 0; k < disp.Z(); ++k)
  for (int j = 0; j < disp.Y(); ++j)
  for (int i = 0; i < disp.X(); ++i) {
    disp(i, j, k, 1) = tmp(i, j, k);
  }
  tmp.Read(dz_name);
  for (int k = 0; k < disp.Z(); ++k)
  for (int j = 0; j < disp.Y(); ++j)
  for (int i = 0; i < disp.X(); ++i) {
    disp(i, j, k, 2) = tmp(i, j, k);
  }
  tmp.Clear();
  return disp;
}

// -----------------------------------------------------------------------------
/// Convert dense displacement field to linear FFD
Transformation *ToLinearFFD(const GenericImage<double> &disp,
                            double dx = .0, double dy = .0, double dz = .0)
{
  UniquePtr<LinearFreeFormTransformation> ffd;

  if (dx <= .0) dx = disp.XSize();
  if (dy <= .0) dy = disp.YSize();
  if (dz <= .0) dz = disp.ZSize();

  if (fequal(dx, disp.XSize()) && fequal(dy, disp.YSize()) && fequal(dz, disp.ZSize())) {

    ffd.reset(new LinearFreeFormTransformation3D(disp, true));

  } else {

    double x1 = 0., y1 = 0., z1 = 0.;
    double x2 = disp.X() - 1;
    double y2 = disp.Y() - 1;
    double z2 = disp.Z() - 1;
    disp.ImageToWorld(x1, y1, z1);
    disp.ImageToWorld(x2, y2, z2);
    double ax[3], ay[3], az[3];
    disp.GetOrientation(ax, ay, az);
    ffd.reset(new LinearFreeFormTransformation3D(x1, y1, z1,
                                                 x2, y2, z2,
                                                 dx, dy, dz,
                                                 ax, ay, az));

    double x, y, z, v[3];
    GenericLinearInterpolateImageFunction<GenericImage<double> > d;
    d.Input(&disp);
    d.Initialize();
    for (int k = 0; k < ffd->Z(); ++k)
    for (int j = 0; j < ffd->Y(); ++j)
    for (int i = 0; i < ffd->X(); ++i) {
      x = i, y = j, z = k;
      ffd->LatticeToWorld(x, y, z);
      disp.WorldToImage(x, y, z);
      d.Evaluate(v, x, y, z);
      ffd->Put(i, j, k, v[0], v[1], v[2]);
    }
  }

  UniquePtr<MultiLevelFreeFormTransformation> dof;
  dof.reset(new MultiLevelFreeFormTransformation());
  dof->PushLocalTransformation(ffd.release());
  return dof.release();
}

// -----------------------------------------------------------------------------
/// Convert dense stationary velocity field to cubic B-spline SVFFD
BSplineFreeFormTransformationSV *ToBSplineSVFFD(GenericImage<double> &velo,
                                                double dx = .0, double dy = .0, double dz = .0)
{
  if (dx <= .0) dx = velo.XSize();
  if (dy <= .0) dy = velo.YSize();
  if (dz <= .0) dz = velo.ZSize();
  UniquePtr<BSplineFreeFormTransformationSV> svffd;
  svffd.reset(new BSplineFreeFormTransformationSV(velo.Attributes(), dx, dy, dz));
  svffd->ApproximateVelocitiesAsNew(velo);
  return svffd.release();
}

// -----------------------------------------------------------------------------
/// Read transformation from dense world space displacement field image
Transformation *ReadWorldDisplacement(const char *fname)
{
  GenericImage<double> disp(fname);
  return ToLinearFFD(disp);
}

// -----------------------------------------------------------------------------
/// Read transformation from dense world space displacement field component images
Transformation *ReadWorldDisplacement(const char *dx_name, const char *dy_name, const char *dz_name)
{
  return ToLinearFFD(ConcatenateWorldDisplacement(dx_name, dy_name, dz_name));
}

// -----------------------------------------------------------------------------
/// Convert displacements in voxel space to world space
void ConvertVoxelToWorldDisplacement(GenericImage<double> &disp)
{
  double x0, y0, z0, x1, y1, z1;
  for (int k = 0; k < disp.Z(); ++k)
  for (int j = 0; j < disp.Y(); ++j)
  for (int i = 0; i < disp.X(); ++i) {
    x0 = i, y0 = j, z0 = k;
    x1 = i + disp(i, j, k, 0);
    y1 = j + disp(i, j, k, 1);
    z1 = k + disp(i, j, k, 2);
    disp.ImageToWorld(x0, y0, z0);
    disp.ImageToWorld(x1, y1, z1);
    disp(i, j, k, 0) = x1 - x0;
    disp(i, j, k, 1) = y1 - y0;
    disp(i, j, k, 2) = z1 - z0;
  }
}

// -----------------------------------------------------------------------------
/// Read transformation from dense voxel space displacement field image
Transformation *ReadVoxelDisplacement(const char *fname)
{
  GenericImage<double> disp(fname);
  ConvertVoxelToWorldDisplacement(disp);
  return ToLinearFFD(disp);
}

// -----------------------------------------------------------------------------
/// Read transformation from dense voxel space displacement field component images
Transformation *ReadVoxelDisplacement(const char *dx_name, const char *dy_name, const char *dz_name)
{
  GenericImage<double> disp = ConcatenateWorldDisplacement(dx_name, dy_name, dz_name);
  ConvertVoxelToWorldDisplacement(disp);
  return ToLinearFFD(disp);
}

// -----------------------------------------------------------------------------
/// Read transformation from dense world space stationary velocity field image
Transformation *ReadWorldSVF(const char *fname)
{
  GenericImage<double> velo(fname);
  return ToBSplineSVFFD(velo);
}

// -----------------------------------------------------------------------------
/// Read transformation from dense world space stationary velocity field component images
Transformation *ReadWorldSVF(const char *dx_name, const char *dy_name, const char *dz_name)
{
  auto velo = ConcatenateWorldDisplacement(dx_name, dy_name, dz_name);
  return ToBSplineSVFFD(velo);
}

// -----------------------------------------------------------------------------
/// Read transformation from dense voxel space stationary velocity field image
Transformation *ReadVoxelSVF(const char *fname)
{
  GenericImage<double> velo(fname);
  ConvertVoxelToWorldDisplacement(velo);
  return ToBSplineSVFFD(velo);
}

// -----------------------------------------------------------------------------
/// Read transformation from dense voxel space stationary velocity field component images
Transformation *ReadVoxelSVF(const char *dx_name, const char *dy_name, const char *dz_name)
{
  GenericImage<double> velo = ConcatenateWorldDisplacement(dx_name, dy_name, dz_name);
  ConvertVoxelToWorldDisplacement(velo);
  return ToBSplineSVFFD(velo);
}

// =============================================================================
// ITK
// =============================================================================

// -----------------------------------------------------------------------------
/// Convert dense vector field from ITK to MIRTK format
/// \sa itk::NiftiImageIO::SetNIfTIOrientationFromImageIO
void ITKToMIRTKVectorField(GenericImage<double> &vf)
{
  if (vf.T() != 2 && vf.T() != 3) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Can only convert 2D/3D vector fields");
  }
  for (int l = 0; l < 2; ++l)
  for (int k = 0; k < vf.Z(); ++k)
  for (int j = 0; j < vf.Y(); ++j)
  for (int i = 0; i < vf.X(); ++i) {
    vf(i, j, k, l) *= -1.;
  }
}

// -----------------------------------------------------------------------------
/// Read displacement field from image written by ITK
Transformation *ReadITKDisplacement(const char *fname)
{
  GenericImage<double> disp(fname);
  ITKToMIRTKVectorField(disp);
  return ToLinearFFD(disp);
}

// -----------------------------------------------------------------------------
/// Read displacement field from component images written by ITK
Transformation *ReadITKDisplacement(const char *dx_name, const char *dy_name, const char *dz_name)
{
  auto disp = ConcatenateWorldDisplacement(dx_name, dy_name, dz_name);
  ITKToMIRTKVectorField(disp);
  return ToLinearFFD(disp);
}

// -----------------------------------------------------------------------------
/// Read stationary velocity field from image written by ITK
Transformation *ReadITKSVF(const char *fname)
{
  GenericImage<double> velo(fname);
  ITKToMIRTKVectorField(velo);
  // FIXME: Use LinearFreeFormTransformationSV
  UniquePtr<BSplineFreeFormTransformationSV> svffd(ToBSplineSVFFD(velo));
  svffd->NumberOfSteps(16);
  return svffd.release();
}

// -----------------------------------------------------------------------------
/// Read stationary velocity field from component images written by ITK
Transformation *ReadITKSVF(const char *dx_name, const char *dy_name, const char *dz_name)
{
  auto velo = ConcatenateWorldDisplacement(dx_name, dy_name, dz_name);
  ITKToMIRTKVectorField(velo);
  // FIXME: Use LinearFreeFormTransformationSV
  UniquePtr<BSplineFreeFormTransformationSV> svffd(ToBSplineSVFFD(velo));
  svffd->NumberOfSteps(16);
  return svffd.release();
}

// =============================================================================
// FSL
// =============================================================================

// -----------------------------------------------------------------------------
Matrix ReadFLIRTMatrix(const char *fname)
{
  ifstream ifs(fname);
  if (!ifs) {
    Warning("Failed to open file for reading: " << fname);
    return false;
  }
  Matrix m(4, 4);
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j) {
    ifs >> m(i, j);
  }
  return m;
}

// -----------------------------------------------------------------------------
/// Read transformation from FLIRT output file
Transformation *ReadFLIRT(const char *fname, ImageAttributes target, ImageAttributes source)
{
  Matrix A = ReadFLIRTMatrix(fname);

  if (target.GetWorldToLatticeOrientation().Det() > 0.) {
    if (verbose) cout << "qform determinant of target image is positive, reflecting x axis" << endl;
    target._xaxis[0] = -target._xaxis[0];
    target._xaxis[1] = -target._xaxis[1];
    target._xaxis[2] = -target._xaxis[2];
  } else {
    if (verbose) cout << "qform determinant of target image is negative, do not reflect x axis" << endl;
  }
  if (source.GetWorldToLatticeOrientation().Det() > 0.) {
    if (verbose) cout << "qform determinant of source image is positive, reflecting x axis" << endl;
    source._xaxis[0] = -source._xaxis[0];
    source._xaxis[1] = -source._xaxis[1];
    source._xaxis[2] = -source._xaxis[2];
  } else {
    if (verbose) cout << "qform determinant of source image is negative, do not reflect x axis" << endl;
  }

  Matrix w2t = target.GetWorldToImageMatrix();
  Matrix s2w = source.GetImageToWorldMatrix();

  Matrix S1(4, 4);
  S1(0, 0) = target._dx;
  S1(1, 1) = target._dy;
  S1(2, 2) = target._dz;
  S1(3, 3) = 1.0;

  Matrix S2(4, 4);
  S2(0, 0) = 1.0 / source._dx;
  S2(1, 1) = 1.0 / source._dy;
  S2(2, 2) = 1.0 / source._dz;
  S2(3, 3) = 1.0;

  UniquePtr<AffineTransformation> dof(new AffineTransformation());
  dof->PutMatrix(s2w * S2 * A.Inverse() * S1 * w2t);
  return dof.release();
}

// -----------------------------------------------------------------------------
/// Read transformation from FNIRT output file
Transformation *ReadFNIRTDisplacement(const char *fname)
{
  return nullptr;
}

// =============================================================================
// Nifty Reg
// =============================================================================

// -----------------------------------------------------------------------------
/// Read transformation from Aladin output file
Transformation *ReadAladin(const char *fname)
{
  Matrix m(4, 4);
  FILE *f = OpenFile(fname, "r");
  #ifdef WINDOWS
  #  define fscanf fscanf_s
  #endif // WINDOWS
  if (fscanf(f, "%lf %lf %lf %lf\n", &m(0, 0), &m(0, 1), &m(0, 2), &m(0, 3)) != 4 ||
      fscanf(f, "%lf %lf %lf %lf\n", &m(1, 0), &m(1, 1), &m(1, 2), &m(1, 3)) != 4 ||
      fscanf(f, "%lf %lf %lf %lf\n", &m(2, 0), &m(2, 1), &m(2, 2), &m(2, 3)) != 4) {
    fclose(f);
    FatalError("File does not appear to be a valid Aladin output file: " << fname);
    exit(1);
  }
  #ifdef WINDOWS
  #  undef fscanf
  #endif // WINDOWS
  m(3, 0) = m(3, 1) = m(3, 2) = .0; m(3, 3) = 1.0;
  fclose(f);
  UniquePtr<AffineTransformation> dof(new AffineTransformation);
  dof->PutMatrix(m);
  return dof.release();
}

// -----------------------------------------------------------------------------
/// Read transformation from reg_f3d output file
Transformation *ReadF3D(const char *fname, const char *dofin_name = NULL,
                        int xyz_units = 0, int steps = 0,
                        F3DTransformationType type = F3D_TYPE_UNKNOWN)
{
  #if MIRTK_IO_WITH_NIfTI
    if (!NiftiImageReader::CheckHeader(fname)) {
      FatalError("Input file is not a NIfTI image file, cannot be a NiftyReg F3D output -cpp file!");
    }

    // Read NIfTI header
    NiftiImageInfo hdr(fname);
    if (type == F3D_TYPE_UNKNOWN) {
      if (hdr.intent_code == NIFTI_INTENT_VECTOR) {
        // latest
        if (hdr.intent_name == "NREG_TRANS") {
          type = static_cast<F3DTransformationType>(static_cast<int>(hdr.intent_p1));
        // v1.3.9
        } else if (hdr.intent_name == "NREG_CPP_FILE") {
          if (hdr.descrip == "Control point position from NiftyReg (reg_f3d)") {
            // reg_f3d without -vel option
            type = F3D_SPLINE_GRID;
          } else if (hdr.descrip == "Velocity field grid from NiftyReg (reg_f3d2)") {
            // reg_f3d with -vel and -sym options
            type = F3D_SPLINE_VEL_GRID;
          }
        } else if (hdr.intent_name == "NREG_VEL_STEP") {
          // reg_f3d with -vel option, but without -sym
          if (hdr.descrip == "Velocity field grid from NiftyReg (reg_f3d2)") {
            type = F3D_SPLINE_VEL_GRID;
          }
        }
      }
      if (type == F3D_TYPE_UNKNOWN) {
        FatalError("Cannot determine format of input F3D vector field from NIfTI intent code!");
      }
    }
    if (steps == 0) {
      if (hdr.intent_name == "NREG_VEL_STEP") {
        steps = static_cast<int>(hdr.intent_p1);
      } else if (hdr.intent_name == "NREG_TRANS") {
        steps = static_cast<int>(hdr.intent_p2);
      }
      if (steps == 0) {
        steps = 6;
      }
    }

    bool displacement = (type == F3D_DISP_FIELD || type == F3D_DISP_VEL_FIELD);

    // Read input vector field
    GenericImage<double> cpp(fname);

    // Change/correct spatial units to mm (NiftyReg ignores xyzt_units)
    double scale = 1.0;
    if      (xyz_units == NIFTI_UNITS_METER)  scale = 1.0e+3;
    else if (xyz_units == NIFTI_UNITS_MICRON) scale = 1.0e-3;
    if (scale != 1.0) {
      double dx, dy, dz;
      cpp.PutOrigin(cpp.GetOrigin() * scale);
      cpp.GetPixelSize(dx, dy, dz);
      cpp.PutPixelSize(dx * scale, dy * scale, dz * scale);
      cpp *= scale;
    }

    // Initialize output transformation
    UniquePtr<MultiLevelFreeFormTransformation> mffd(new MultiLevelFreeFormTransformation());

    // Subtract -dofin displacement
    if (dofin_name) {
      AffineTransformation dof;
      dof.Read(dofin_name);
      for (int k = 0; k < cpp.Z(); ++k)
      for (int j = 0; j < cpp.Y(); ++j)
      for (int i = 0; i < cpp.X(); ++i) {
        double dx = i, dy = j, dz = k;
        cpp.ImageToWorld(dx, dy, dz);
        dof.Displacement(dx, dy, dz);
        cpp(i, j, k, 0) -= dx;
        cpp(i, j, k, 1) -= dy;
        if (cpp.T() > 2) cpp(i, j, k, 2) -= dz;
      }
      mffd->GetGlobalTransformation()->PutMatrix(dof.GetMatrix());
    }

    // Construct transformation from vector field
    switch (type) {
      case F3D_DEF_FIELD:
      case F3D_DISP_FIELD: {
        mffd->PushLocalTransformation(new LinearFreeFormTransformation(cpp, displacement));
        break;
      }
      case F3D_SPLINE_GRID: {
        mffd->PushLocalTransformation(new BSplineFreeFormTransformation3D(cpp, displacement));
        break;
      }
      case F3D_SPLINE_VEL_GRID: {
        BSplineFreeFormTransformationSV *ffd = new BSplineFreeFormTransformationSV(cpp, displacement);
        ffd->NumberOfSteps(static_cast<int>(pow(2, steps)));
        ffd->IntegrationMethod(FFDIM_FastSS);
        mffd->PushLocalTransformation(ffd);
        break;
      }
      default: {
        FatalError("Cannot convert given F3D vector field format or memory allocation failed!");
        exit(1);
      }
    }

    return mffd.release();
  #else // MIRTK_IO_WITH_NIfTI
    Warning("Cannot read F3D output file without NIfTI module");
    return NULL;
  #endif // MIRTK_IO_WITH_NIfTI
}

// =============================================================================
// FreeSurfer
// =============================================================================

// -----------------------------------------------------------------------------
/// Read affine transformation from FreeSurfer .xfm file
Transformation *ReadXFM(const char *fname)
{
  Matrix m(4, 4);

  ifstream ifs(fname);
  if (!ifs.is_open()) {
    Warning("Failed to open file for reading: " << fname);
    return nullptr;
  }
  string line;
  if (!getline(ifs, line) || line != "MNI Transform File") {
    Warning("Input file is not a MNI Transform File: " << fname);
    return nullptr;
  }
  while (getline(ifs, line)) {
    if (line == "Linear_Transform = ") break;
  }
  if (ifs.eof()) {
    Warning("Could not find \"Linear_Transform = \" on a single line before the matrix entries");
    return nullptr;
  }
  for (int i = 0; i < 3; ++i) {
    if (!getline(ifs, line)) {
      Warning("Failed to read affine matrix entries from file " << fname);
      return nullptr;
    }
    istringstream is(line);
    is >> m(i, 0) >> m(i, 1) >> m(i, 2) >> m(i, 3);
  }
  m(3, 0) = m(3, 1) = m(3, 2) = .0; m(3, 3) = 1.0;
  ifs.close();
  UniquePtr<AffineTransformation> dof(new AffineTransformation);
  dof->PutMatrix(m);
  return dof.release();
}

// =============================================================================
// elastix
// =============================================================================

// -----------------------------------------------------------------------------
/// Read transformation from elastix output file
Transformation *ReadElastix(const char *fname, const ImageAttributes &target)
{
  int ndims = 0;
  int ndofs = 0;
  int order = 3;
  UniquePtr<double> dofs;
  ImageAttributes grid;
  bool cyclic = false;

  ifstream ifs(fname);
  if (!ifs.is_open()) {
    Warning("Failed to open file for reading: " << fname);
    return nullptr;
  }
  string line;
  Array<string> part;
  if (!getline(ifs, line) || line != "(Transform \"BSplineTransform\")") {
    Warning("Elastix input: Expected first line to specify BSplineTransform type");
    return nullptr;
  }
  while (getline(ifs, line)) {
    line = Trim(line);
    if (line.empty()) continue;
    if (line[0] == '/' && line[1] == '/') continue;
    if (line[0] != '(' || line[line.length()-1] != ')') {
      Warning("Elastix input: Expected lines to be enclosed in parentheses");
      return nullptr;
    }
    line = line.substr(1, line.length() - 2);
    part = Split(line, ' ', 0, true, true);
    if (part.size() < 2) {
      Warning("Elastix input: Expected lines to contain at least two entries, a key and a value");
      return nullptr;
    }
    if (part[0] == "NumberOfParameters") {
      if (part.size() != 2 || !FromString(part[1], ndofs) || ndofs <= 0) {
        Warning("Elastix input: Expected NumberOfParameters to have exactly one positive integer value");
        return nullptr;
      }
      dofs.reset(new double[ndofs]);
    } else if (part[0] == "TransformParameters") {
      if (part.size() != static_cast<size_t>(ndofs + 1)) {
        Warning("Elastix input: Expected " << ndofs << " TransformParameters, got " << part.size() - 1);
        return nullptr;
      }
      for (int i = 0; i < ndofs; ++i) {
        if (!FromString(part[i+1], dofs.get()[i])) {
          Warning("Elastix input: Failed to parse TransformParameters value");
          return nullptr;
        }
      }
    } else if (part[0] == "GridSize") {
      if (ndims == 0) {
        ndims = static_cast<int>(part.size() - 1);
      } else if (static_cast<int>(part.size() - 1) != ndims) {
        Warning("Elastix input: Expected GridSize to have " << ndims << " values, got " << part.size() - 1);
        return nullptr;
      }
      if (ndims < 2 || ndims > 4) {
        Warning("Elastix input: Can only read transformation with 2, 3, or 4 dimensions");
        return nullptr;
      }
      for (int i = 0, n; i < ndims; ++i) {
        if (!FromString(part[i+1], n) || n <= 0) {
          Warning("Elastix input: Failed to parse GridSize value");
          return nullptr;
        }
        if      (i == 0) grid._x = n;
        else if (i == 1) grid._y = n;
        else if (i == 2) grid._z = n;
        else if (i == 3) grid._t = n;
      }
    } else if (part[0] == "GridIndex") {
      if (ndims == 0) {
        ndims = static_cast<int>(part.size() - 1);
      } else if (static_cast<int>(part.size() - 1) != ndims) {
        Warning("Elastix input: Expected GridIndex to have " << ndims << " values, got " << part.size() - 1);
        return nullptr;
      }
      if (ndims < 2 || ndims > 4) {
        Warning("Elastix input: Can only read transformation with 2, 3, or 4 dimensions");
        return nullptr;
      }
      for (int i = 0, n; i < ndims; ++i) {
        if (!FromString(part[i+1], n)) {
          Warning("Elastix input: Failed to parse GridIndex value");
          return nullptr;
        }
        if (n != 0) {
          Warning("Elastix input: Only GridIndex equal to 0 supported");
          return nullptr;
        }
      }
    } else if (part[0] == "GridSpacing") {
      if (ndims == 0) {
        ndims = static_cast<int>(part.size() - 1);
      } else if (static_cast<int>(part.size() - 1) != ndims) {
        Warning("Elastix input: Expected GridSpacing to have " << ndims << " values, got " << part.size() - 1);
        return nullptr;
      }
      if (ndims < 2 || ndims > 4) {
        Warning("Elastix input: Can only read transformation with 2, 3, or 4 dimensions");
        return nullptr;
      }
      double s;
      for (int i = 0; i < ndims; ++i) {
        if (!FromString(part[i+1], s) || s <= 0.) {
          Warning("Elastix input: Failed to parse GridSpacing value");
          return nullptr;
        }
        if      (i == 0) grid._dx = s;
        else if (i == 1) grid._dy = s;
        else if (i == 2) grid._dz = s;
        else if (i == 3) grid._dt = s;
      }
    } else if (part[0] == "GridOrigin") {
      if (ndims == 0) {
        ndims = static_cast<int>(part.size() - 1);
      } else if (static_cast<int>(part.size() - 1) != ndims) {
        Warning("Elastix input: Expected GridOrigin to have " << ndims << " values, got " << part.size() - 1);
        return nullptr;
      }
      if (ndims < 2 || ndims > 4) {
        Warning("Elastix input: Can only read transformation with 2, 3, or 4 dimensions");
        return nullptr;
      }
      double x;
      for (int i = 0; i < ndims; ++i) {
        if (!FromString(part[i+1], x)) {
          Warning("Elastix input: Failed to parse GridOrigin value");
          return nullptr;
        }
        if      (i == 0) grid._xorigin = x;
        else if (i == 1) grid._yorigin = x;
        else if (i == 2) grid._zorigin = x;
        else if (i == 3) grid._torigin = x;
      }
    } else if (part[0] == "GridDirection") {
      if (static_cast<int>(part.size() - 1) != 9) {
        Warning("Elastix input: Expected GridDirection to have 9 values, got " << part.size() - 1);
        return nullptr;
      }
      double v;
      for (int i = 0; i < 9; ++i) {
        if (!FromString(part[i+1], v)) {
          Warning("Elastix input: Failed to parse GridDirection value");
          return nullptr;
        }
        if      (i == 0) grid._xaxis[0] = v;
        else if (i == 1) grid._yaxis[0] = v;
        else if (i == 2) grid._zaxis[0] = v;
        else if (i == 3) grid._xaxis[1] = v;
        else if (i == 4) grid._yaxis[1] = v;
        else if (i == 5) grid._zaxis[1] = v;
        else if (i == 6) grid._xaxis[2] = v;
        else if (i == 7) grid._yaxis[2] = v;
        else if (i == 8) grid._zaxis[2] = v;
      }
    } else if (part[0] == "BSplineTransformSplineOrder") {
      if (part.size() != 2) {
        Warning("Elastix input: Expected BSplineTransformSplineOrder to have exactly one value");
        return nullptr;
      }
      if (!FromString(part[1], order)) {
        Warning("Elastix input: Failed to parse BSplineTransformSplineOrder value");
        return nullptr;
      }
    } else if (part[0] == "UseCyclicTransform") {
      if (!FromString(part[1], cyclic)) {
        Warning("Elastix input: Failed to parse UseCyclicTransform value");
        return nullptr;
      }
    }
  }
  UniquePtr<FreeFormTransformation> dof;
  if (ndofs == 0) {
    Warning("Elastix input: NumberOfParameters not specified");
    return nullptr;
  }
  if (ndims < 2 || ndims > 4) {
    Warning("Elastix input: Can only read transformation with 2, 3, or 4 dimensions");
    return nullptr;
  }
  if (!grid) {
    Warning("Elastix input: Invalid BSplineTransform grid");
    return nullptr;
  }
  int ncps = grid.NumberOfPoints();
  if (ndims * ncps != ndofs) {
    Warning("Elastix input: NumberOfParameters does not match GridSize");
    return nullptr;
  }
  // Deal with itk::NiftiImageIO::SetNIfTIOrientationFromImageIO nonsense
  grid._xorigin *= -1.;
  grid._yorigin *= -1.;
  for (int i = 0; i < 2; ++i) {
    grid._xaxis[i] *= -1.;
    grid._yaxis[i] *= -1.;
    grid._zaxis[i] *= -1.;
  }
  // Move origin to grid center
  double tx = .5 * (grid._x - 1) * grid._dx;
  double ty = .5 * (grid._y - 1) * grid._dy;
  double tz = .5 * (grid._z - 1) * grid._dz;
  grid._xorigin += tx * grid._xaxis[0] + ty * grid._yaxis[0] + tz * grid._zaxis[0];
  grid._yorigin += tx * grid._xaxis[1] + ty * grid._yaxis[1] + tz * grid._zaxis[1];
  grid._zorigin += tx * grid._xaxis[2] + ty * grid._yaxis[2] + tz * grid._zaxis[2];
  // Create MIRTK transformation
  if (order == 1) {
    if (ndims == 4) {
      dof.reset(new LinearFreeFormTransformation4D(grid));
    } else {
      dof.reset(new LinearFreeFormTransformation3D(grid));
    }
  } else if (order == 3) {
    if (ndims == 4) {
      dof.reset(new BSplineFreeFormTransformation4D(grid));
    } else {
      dof.reset(new BSplineFreeFormTransformation3D(grid));
    }
  } else {
    Warning("Elastix input: Can only convert BSplineTransform with BSplineTransformSplineOrder 1 or 3");
    return nullptr;
  }
  double *x = dofs.get();
  double *y = x + ncps;
  double *z = y + ncps;
  for (int i = 0; i < ncps; ++i, ++x, ++y, ++z) {
    // x and y axes are reflected by ITK
    dof->Put(i, -(*x), -(*y), *z);
  }
  return dof.release();
}

// =============================================================================
// DRAMMS
// =============================================================================

// -----------------------------------------------------------------------------
/// Read transformation from displacement field image written by DRAMMS
Transformation *ReadDRAMMS(const char *fname)
{
  // DRAMMS data order is yxzyxzyxz... instead of xxx...yyy...zzz..., but the
  // header information still relates to the respective x, y, and z axes
  GenericImage<double> in(fname);
  GenericImage<double> disp(in.Attributes());
  double *data = in.Data();
  for (int k = 0; k < in.Z(); ++k)
  for (int j = 0; j < in.Y(); ++j)
  for (int i = 0; i < in.X(); ++i, data += 3) {
    disp(i, j, k, 0) = data[1];
    disp(i, j, k, 1) = data[0];
    disp(i, j, k, 2) = data[2];
  }
  in.Clear();
  // Convert to physical displacements
  ConvertVoxelToWorldDisplacement(disp);
  // Convert to linear FFD
  return ToLinearFFD(disp);
}

////////////////////////////////////////////////////////////////////////////////
// Write transformation
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Dense displacement field
// =============================================================================

// -----------------------------------------------------------------------------
/// Convert transformation to dense world space displacement field
bool ToWorldDisplacement(GenericImage<double> &disp, const Transformation  *dof,
                                                     const ImageAttributes &target,
                                                     double                 ts)
{
  if (target) {
    disp.Initialize(target, 3);
  } else {
    const MultiLevelTransformation *mffd = dynamic_cast<const MultiLevelTransformation *>(dof);
    const FreeFormTransformation   *ffd  = dynamic_cast<const FreeFormTransformation   *>(dof);
    if (mffd) ffd = mffd->GetLocalTransformation(-1);
    if (!ffd) {
      Warning("Cannot convert linear transformation to image without input -target image!");
      return false;
    }
    disp.Initialize(ffd->Attributes(), 3);
  }
  dof->Displacement(disp, ts, target._torigin);
  return true;
}

// -----------------------------------------------------------------------------
/// Convert displacements in world space to voxel space
void ConvertWorldToVoxelDisplacement(GenericImage<double> &disp)
{
  double x0, y0, z0, x1, y1, z1;
  for (int k = 0; k < disp.Z(); ++k)
  for (int j = 0; j < disp.Y(); ++j)
  for (int i = 0; i < disp.X(); ++i) {
    x0 = i, y0 = j, z0 = k;
    x1 = i + disp(i, j, k, 0);
    y1 = j + disp(i, j, k, 1);
    z1 = k + disp(i, j, k, 2);
    disp.ImageToWorld(x0, y0, z0);
    disp.ImageToWorld(x1, y1, z1);
    disp(i, j, k, 0) = x1 - x0;
    disp(i, j, k, 1) = y1 - y0;
    disp(i, j, k, 2) = z1 - z0;
  }
}

// -----------------------------------------------------------------------------
/// Write transformation to dense world space displacement field image
bool WriteWorldDisplacement(const char *fname, const Transformation *dof,
                            const ImageAttributes &target, double ts)
{
  GenericImage<double> disp;
  if (ToWorldDisplacement(disp, dof, target, ts)) {
    disp.Write(fname);
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Write transformation to dense world space displacement field component images
bool WriteWorldDisplacement(const char *dx_name, const char *dy_name, const char *dz_name,
                            const Transformation *dof, const ImageAttributes &target, double ts)
{
  GenericImage<double> disp;
  if (ToWorldDisplacement(disp, dof, target, ts)) {
    disp.GetFrame(0).Write(dx_name);
    disp.GetFrame(1).Write(dy_name);
    disp.GetFrame(2).Write(dz_name);
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Write transformation to dense voxel space displacement field image
bool WriteVoxelDisplacement(const char *fname, const Transformation *dof,
                            const ImageAttributes &target, double ts)
{
  GenericImage<double> disp;
  if (ToWorldDisplacement(disp, dof, target, ts)) {
    ConvertWorldToVoxelDisplacement(disp);
    disp.Write(fname);
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Write transformation to dense voxel space displacement field component images
bool WriteVoxelDisplacement(const char *dx_name, const char *dy_name, const char *dz_name,
                            const Transformation *dof, const ImageAttributes &target, double ts)
{
  GenericImage<double> disp;
  if (ToWorldDisplacement(disp, dof, target, ts)) {
    ConvertWorldToVoxelDisplacement(disp);
    disp.GetFrame(0).Write(dx_name);
    disp.GetFrame(1).Write(dy_name);
    disp.GetFrame(2).Write(dz_name);
    return true;
  }
  return false;
}

// =============================================================================
// MIRTK
// =============================================================================

// -----------------------------------------------------------------------------
/// Approximate a given MIRTK transformation by another
///
/// This function tries to approximate the input transformation by an output
/// transformation of a given type. When the output transformation is a
/// multi-level free-form deformation with global and local components,
/// the input transformation is approximated by a combination of the global
/// transformation and the active local transformations. Whenever possible,
/// this function attempts to copy the global and/or local transformation
/// parameters from the input to the output transformation in order to reduce
/// the approximation error and not require a more costly approximation step.
/// This generic function can be used to convert from one MIRTK transformation
/// type to any other transformation type. Another use case would be the
/// initialization of a registration output transformation given an initial
/// guess (i.e., manual user input transformation or result of a previous step).
///
/// \returns RMS error of approximation or NaN in case of a failure.
double ApproximateAsNew(const Transformation *dofin,
                        Transformation       *dofout,
                        ImageAttributes      *domain = nullptr,
                        FFDIMParams           ffdim  = FFDIMParams(),
                        BCHParams             bch    = BCHParams())
{
  UniquePtr<Transformation> itmp;

  // Just copy parameters whenever possible
  if (dofout->CopyFrom(dofin)) return .0;

  // Reset output transformation
  dofout->Reset();

  // Input...
  const HomogeneousTransformation *ilin  = nullptr; // ...linear transformation
  const FreeFormTransformation    *iffd  = nullptr; // or non-linear  FFD
  const MultiLevelTransformation  *imffd = nullptr; // or multi-level FFD

  ( ilin = dynamic_cast<const HomogeneousTransformation *>(dofin)) ||
  ( iffd = dynamic_cast<const FreeFormTransformation    *>(dofin)) ||
  (imffd = dynamic_cast<const MultiLevelTransformation  *>(dofin));

  if (imffd) {
    if (imffd->NumberOfLevels() == 0) {
      dofin = ilin = imffd->GetGlobalTransformation();
      imffd = nullptr;
    } else if (imffd->NumberOfLevels() == 1 && imffd->GetGlobalTransformation()->IsIdentity()) {
      dofin = iffd = imffd->GetLocalTransformation(0);
      imffd = nullptr;
    }
  }

  // Output...
  HomogeneousTransformation *olin  = nullptr; // ...linear transformation
  FreeFormTransformation    *offd  = nullptr; // or non-linear FFD
  MultiLevelTransformation  *omffd = nullptr; // or multi-level FFD

  ( olin = dynamic_cast<HomogeneousTransformation *>(dofout)) ||
  ( offd = dynamic_cast<FreeFormTransformation    *>(dofout)) ||
  (omffd = dynamic_cast<MultiLevelTransformation  *>(dofout));

  if (omffd && omffd->NumberOfActiveLevels() == 0) {
    dofout = olin = omffd->GetGlobalTransformation();
    omffd  = nullptr;
  }

  // When input is a single-level MFFD and the output is a free-form deformation
  // which has the same type as the FFD of the input transformation,
  // merge the global component of the input transformation into its local
  // component before proceeding.
  if (imffd && imffd->NumberOfLevels() == 1 && offd &&
      offd->TypeOfClass() == imffd->GetLocalTransformation(0)->TypeOfClass()) {
    MultiLevelTransformation *mffd;
    FreeFormTransformation   *affd;
    itmp.reset(Transformation::New(imffd));
    mffd = dynamic_cast<MultiLevelTransformation *>(itmp.get());
    mirtkAssert(mffd != nullptr, "copy must also be a multi-level transformation");
    mffd->MergeGlobalIntoLocalDisplacement();
    itmp.reset(affd = mffd->PopLocalTransformation());
    imffd = mffd = nullptr;
    dofin = iffd = affd;
  }

  // When both input and output transformation are a MFFD of the same type
  // and the output MFFD has only a single active FFD level, copy the global
  // component from the input MFFD to the output MFFD and approximate the local
  // component of the input MFFD by the output FFD.
  if (imffd && omffd && omffd->NumberOfActiveLevels() == 1 && imffd->TypeOfClass() == omffd->TypeOfClass()) {
    // Copy global component
    omffd->GetGlobalTransformation()->CopyFrom(imffd->GetGlobalTransformation());
    // Remaining local component of output MFFD
    offd = nullptr;
    for (int n = 0; n < omffd->NumberOfLevels(); ++n) {
      if (omffd->LocalTransformationIsActive(n)) {
        offd = omffd->GetLocalTransformation(n);
        break;
      }
    }
    mirtkAssert(offd != nullptr, "multi-level transformation claimed to have one active level");
    dofout = offd;
    omffd  = nullptr;
    // Remaining local component of input MFFD
    if (imffd->NumberOfLevels() == 1) {
      dofin = iffd = imffd->GetLocalTransformation(0);
      imffd = nullptr;
    } else {
      MultiLevelTransformation *mffd;
      FreeFormTransformation   *affd;
      itmp.reset(Transformation::New(imffd->TypeOfClass()));
      mffd = dynamic_cast<MultiLevelTransformation *>(itmp.get());
      mirtkAssert(mffd != nullptr, "copy must also be a multi-level transformation");
      for (int n = 0; n < imffd->NumberOfLevels(); ++n) {
        affd = const_cast<FreeFormTransformation *>(imffd->GetLocalTransformation(n));
        mffd->PushLocalTransformation(affd, /* transfer_ownership = */ false);
      }
      dofin = imffd = mffd;
    }
  }

  // Case 1: Input is linear homogeneous coordinate transformation
  if (ilin) {
    if (olin) {
      olin->CopyFrom(ilin);
      // No error when input transformation is subtype of output transformation
      if ( ilin->TypeOfClass() == TRANSFORMATION_RIGID ||
          (ilin->TypeOfClass() == TRANSFORMATION_SIMILARITY && olin->TypeOfClass() != TRANSFORMATION_RIGID) ||
          (ilin->TypeOfClass() == /*TRANSFORMATION_AFFINE == */ olin->TypeOfClass())) {
        return .0;
      }
      // Otherwise, compute RMS error after possible loss of linear component
      if (domain) {
        ImageAttributes lattice(*domain);
        Matrix   i2w = lattice.GetImageToWorldMatrix();
        lattice._i2w = &i2w;
        double x, y, z, x1, y1, z1, x2, y2, z2, error = .0;
        for (int k = 0; k < lattice._z; ++k)
        for (int j = 0; j < lattice._y; ++j)
        for (int i = 0; i < lattice._x; ++i) {
          x = i, y = j, z = k;
          lattice.LatticeToWorld(x, y, z);
          x1 = x2 = x;
          y1 = y2 = y;
          z1 = z2 = z;
          ilin->Transform(x1, y1, z1);
          olin->Transform(x2, y2, z2);
          error += pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2);
        }
        return sqrt(error);
      }
      return .0;
    }
    if (omffd) {
      omffd->GetGlobalTransformation()->CopyFrom(ilin);
      return .0;
    }
  // Case 2: Input is a single free-form deformation
  } else if (iffd) {
    if (offd && offd->CopyFrom(iffd)) return .0;
  }
  // Case 3: Input is multi-level free-form deformation with global and local components
  //         or parameters of input transformation in case 1 & 2 could not just
  //         be copied over to the output transformation...
  double error = numeric_limits<double>::quiet_NaN();

  // ...determine discrete lattice on which to perform approximation
  ImageAttributes lattice;
  if (domain) lattice = *domain;
  if (iffd && iffd->Attributes().NumberOfPoints() > lattice.NumberOfPoints()) {
    lattice = iffd->Attributes();
  }
  if (imffd) {
    const FreeFormTransformation *affd;
    for (int n = 0; n < imffd->NumberOfLevels(); ++n) {
      affd = imffd->GetLocalTransformation(n);
      if (affd && affd->Attributes().NumberOfPoints() > lattice.NumberOfPoints()) {
        lattice = affd->Attributes();
      }
    }
  }
  if (offd && offd->Attributes().NumberOfPoints() > lattice.NumberOfPoints()) {
    lattice = offd->Attributes();
  }
  if (omffd) {
    const FreeFormTransformation *affd;
    for (int n = 0; n < omffd->NumberOfLevels(); ++n) {
      affd = omffd->GetLocalTransformation(n);
      if (affd && affd->Attributes().NumberOfPoints() > lattice.NumberOfPoints()) {
        lattice = affd->Attributes();
      }
    }
  }
  mirtkAssert(lattice.NumberOfPoints() > 0, "approximation domain valid");
  if (lattice.NumberOfPoints() == 0) return error;

  // ...when possible, use specialized approximation methods
  if (offd) {
    LinearFreeFormTransformationTD  *ltdffd = nullptr;
    BSplineFreeFormTransformationSV *bsvffd = nullptr;
    BSplineFreeFormTransformationTD *btdffd = nullptr;

    (ltdffd = dynamic_cast<LinearFreeFormTransformationTD  *>(offd)) ||
    (bsvffd = dynamic_cast<BSplineFreeFormTransformationSV *>(offd)) ||
    (btdffd = dynamic_cast<BSplineFreeFormTransformationTD *>(offd));

    if (ltdffd || bsvffd || btdffd) {
      // Get input displacements at lattice points of approximation domain
      GenericImage<double> disp(lattice, 3);
      dofin->Displacement(disp);
      // Approximate displacements
      if (bsvffd) {
        error = bsvffd->ApproximateAsNew(disp, bch.smooth, bch.nterms, bch.nsteps);
      } else if (ltdffd) {
        GenericImage<double> *disps[1] = { &disp };
        error = ltdffd->ApproximateAsNew(disps, &ffdim.t1, &ffdim.t2, 1,
                                         bch.smooth, bch.nterms, bch.nsteps);
      } else if (btdffd) {
        GenericImage<double> *disps[1] = { &disp };
        error = btdffd->ApproximateAsNew(disps, &ffdim.t1, &ffdim.t2, 1,
                                         bch.smooth, bch.nterms, bch.nsteps);
      }
    }
  }

  // ...otherwise, use generic approximation interface
  if (IsNaN(error)) {
    error = dofout->ApproximateAsNew(lattice, dofin);
  }

  if (domain) *domain = lattice;
  return error;
}

// -----------------------------------------------------------------------------
/// Report approximation error
void PrintApproximationError(const Transformation  *dofin,
                             const Transformation  *dofout,
                             const ImageAttributes &domain,
                             double                 ts     = .0,
                             int                    margin = 0,
                             Indent                 indent = Indent())
{
  double x, y, z, error, mag;
  double dx1, dy1, dz1;
  double dx2, dy2, dz2;
  double avg_idisp = .0;
  double max_idisp = .0;
  double avg_odisp = .0;
  double max_odisp = .0;
  double avg_error = .0;
  double max_error = .0;

  // Lattice to world coordinate transformation matrix
  const Matrix i2w = domain.GetImageToWorldMatrix();

  // Whether to pre-compute displacements for entire spatial domain
  const bool use_idisp = dofin ->RequiresCachingOfDisplacements();
  const bool use_odisp = dofout->RequiresCachingOfDisplacements();

  // Initialize displacement caches
  GenericImage<double> idisp, odisp;
  if (use_idisp) idisp.Initialize(domain, 3);
  if (use_odisp) odisp.Initialize(domain, 3);

  // Time point of source image
  const double t0 = domain.LatticeToTime(0);

  // Pre-compute displacements for entire volume domain
  if (idisp) dofin ->Displacement(idisp, ts, t0);
  if (odisp) dofout->Displacement(odisp, ts, t0);

  // Evaluate approximation error for source image
  for (int k = margin; k < domain._z - margin; ++k)
  for (int j = margin; j < domain._y - margin; ++j)
  for (int i = margin; i < domain._x - margin; ++i) {
    x = i, y = j, z = k;
    Transform(i2w, x, y, z);

    if (idisp) {
      dx1 = idisp(i, j, k, 0);
      dy1 = idisp(i, j, k, 1);
      dz1 = idisp(i, j, k, 2);
    } else {
      dx1 = x, dy1 = y, dz1 = z;
      dofin->Displacement(dx1, dy1, dz1, ts, t0);
    }
    mag = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
    avg_idisp += mag;
    if (mag > max_idisp) max_idisp = mag;

    if (odisp) {
      dx2 = odisp(i, j, k, 0);
      dy2 = odisp(i, j, k, 1);
      dz2 = odisp(i, j, k, 2);
    } else {
      dx2 = x, dy2 = y, dz2 = z;
      dofout->Displacement(dx2, dy2, dz2, ts, t0);
    }
    mag = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
    avg_odisp += mag;
    if (mag > max_odisp) max_odisp = mag;

    error = sqrt((dx2 - dx1) * (dx2 - dx1) +
                 (dy2 - dy1) * (dy2 - dy1) +
                 (dz2 - dz1) * (dz2 - dz1));
    avg_error += error;
    if (error > max_error) max_error = error;
  }

  const int n = domain.NumberOfPoints();
  avg_idisp /= n;
  avg_odisp /= n;
  avg_error /= n;

  cout << indent << "Average input displacement:  " << avg_idisp << endl;
  cout << indent << "Maximum input displacement:  " << max_idisp << endl;
  cout << indent << "Average output displacement: " << avg_odisp << endl;
  cout << indent << "Maximum output displacement: " << max_odisp << endl;
  cout << indent << "Average RMS error:           " << avg_error << endl;
  cout << indent << "Maximum RMS error:           " << max_error << endl;
}

// -----------------------------------------------------------------------------
/// Write MIRTK transformation file
bool WriteMIRTK(const char *fname, Transformation *dof,
                ImageAttributes target_attr = ImageAttributes(), double ts = .0,
                double dx = .0, double dy = .0, double dz = .0,  double dt = .0,
                TransformationType type      = TRANSFORMATION_UNKNOWN,
                MFFDMode           mffd_type = MFFD_Default,
                FFDIMParams        ffdim     = FFDIMParams(),
                BCHParams          bch       = BCHParams())
{
  // Boundary margin to exclude from approximation error evaluation
  const int rms_excl_margin = 2;

  // Output transformation is a homogeneous coordinate transformation
  const bool type_is_linear = (type == TRANSFORMATION_RIGID      ||
                               type == TRANSFORMATION_SIMILARITY ||
                               type == TRANSFORMATION_AFFINE);

  // Determine type of input transformation
  TransformationType dof_type = dof->TypeOfClass();

  HomogeneousTransformation *aff  = dynamic_cast<HomogeneousTransformation *>(dof);
  FreeFormTransformation    *ffd  = dynamic_cast<FreeFormTransformation    *>(dof);
  MultiLevelTransformation  *mffd = dynamic_cast<MultiLevelTransformation  *>(dof);

  if (mffd) {
    aff = mffd->GetGlobalTransformation();
    // Case 1: MFFD only has global component
    if (mffd->NumberOfLevels() == 0) {
      dof      = aff;
      dof_type = TRANSFORMATION_AFFINE;
      mffd     = nullptr;
    } else {
      // Merge global transformation into local transformation when output type
      // is a local transformation without explicit global component
      if (!type_is_linear && mffd_type == MFFD_None) {
        mffd->MergeGlobalIntoLocalDisplacement();
      }
      // Get local transformation with highest resolution (i.e., most #DoFs)
      // and determine common type of all local transformations (if unique)
      ffd      = mffd->GetLocalTransformation(0);
      dof_type = ffd ->TypeOfClass();
      if (mffd->NumberOfLevels() > 1) {
        for (int n = 1; n < mffd->NumberOfLevels(); ++n) {
          FreeFormTransformation *affd = mffd->GetLocalTransformation(0);
          if (affd->TypeOfClass() != dof_type) {
            dof_type = TRANSFORMATION_UNKNOWN;
            break;
          }
          if (affd->NumberOfDOFs() > ffd->NumberOfDOFs()) ffd = affd;
        }
      // Case 2: MFFD has only a single local component
      } else if (mffd->GetGlobalTransformation()->IsIdentity()) {
        dof  = ffd;
        mffd = nullptr;
      }
      // Case 3: MFFD has both global and local components
    }
  }

  // Attributes of output FFD (if output is not only a homogeneous transformation)
  if (ffd) {
    if (!target_attr) target_attr = ffd->Attributes();
    if (dx <= .0) dx = ffd->GetXSpacing();
    if (dy <= .0) dy = ffd->GetYSpacing();
    if (dz <= .0) dz = ffd->GetZSpacing();
    if (dt <= .0) dt = ffd->GetTSpacing();
  } else if (target_attr) {
    if (dx <= .0) dx = target_attr._dx;
    if (dy <= .0) dy = target_attr._dy;
    if (dz <= .0) dz = target_attr._dz;
    if (dt <= .0) dt = target_attr._dt;
  }

  // Check if requested spacing of output FFD differs from input transformation
  // Always resample a MFFD with multiple levels such that output has single level
  bool resample_ffd;
  if (mffd && mffd->NumberOfLevels() > 1) {
    resample_ffd = true;
  } else if (!type_is_linear && ffd) {
    resample_ffd = ( !fequal(ffd->GetXSpacing(), dx) ||
                     !fequal(ffd->GetYSpacing(), dy) ||
                    (!fequal(ffd->GetZSpacing(), dz) && ffd->Z() > 1) ||
                    (!fequal(ffd->GetTSpacing(), dt) && ffd->T() > 1));
  } else {
    resample_ffd = false;
  }

  // Instantiate output MFFD (if any)
  UniquePtr<MultiLevelTransformation> omffd;
  switch (mffd_type) {
    case MFFD_None:
      break;
    case MFFD_Default:
      if (!type_is_linear) {
        omffd.reset(new MultiLevelFreeFormTransformation());
      }
      break;
    case MFFD_Sum:
      omffd.reset(new MultiLevelFreeFormTransformation());
      break;
    case MFFD_Fluid:
      omffd.reset(new FluidFreeFormTransformation());
      break;
    case MFFD_LogSum:
      if (type != TRANSFORMATION_BSPLINE_FFD_SV) {
        Warning("Multi-level mode " << ToString(MFFD_LogSum) << " only suitable for a B-spline SV FFD output transformation!");
        return false;
      }
      omffd.reset(new MultiLevelStationaryVelocityTransformation());
      break;
    default:
      Warning("The " << ToString(mffd_type) << " multi-level transformation mode is not supported!");
      return false;
  }

  // Write input transformation directly when possible
  if (type == TRANSFORMATION_UNKNOWN) type = dof_type;
  if (!resample_ffd && dof_type == type) {
    if (omffd) {
      if (mffd) {
        if (omffd->TypeOfClass() == mffd->TypeOfClass()) {
          dof->Write(fname);
          return true;
        }
      } else if (ffd) {
        const bool transfer_ownership = false;
        omffd->PushLocalTransformation(ffd, transfer_ownership);
        omffd->Write(fname);
        omffd->PopLocalTransformation();
        return true;
      } else if (aff) {
        omffd->GetGlobalTransformation()->CopyFrom(aff);
        omffd->Write(fname);
        return true;
      }
    } else {
      dof->Write(fname);
      return true;
    }
  }

  // Otherwise, approximate input transformation by new output transformation
  UniquePtr<Transformation> odof(Transformation::New(type));

  // When output type is a homogeneous coordinate transformation,
  // convert input transformation such that only global part remains
  // and write resulting transformation either as MFFD without local
  // component or a plain homogeneous coordinate transformation
  HomogeneousTransformation *oaff = dynamic_cast<HomogeneousTransformation *>(odof.get());
  if (oaff) {
    double rms = ApproximateAsNew(dof, oaff, &target_attr, ffdim, bch);
    if (IsNaN(rms)) return false;
    if (verbose > 1 && target_attr) {
      PrintApproximationError(dof, oaff, target_attr, ts, rms_excl_margin);
    } else if (verbose > 0) {
      cout << "RMS error of approximation = " << rms << endl;
    }
    if (omffd) {
      omffd->GetGlobalTransformation()->CopyFrom(oaff);
      omffd->Write(fname);
    } else {
      oaff->Write(fname);
    }
    return true;
  }

  // Otherwise, when output type is a FFD, set the requested lattice attributes and
  // insert FFD into desired output MFFD or write plain FFD without global component.
  // In case of a FFD parameterized by velocities, also set the parameters of
  // the integration method used to obtain the displacement vectors.
  FreeFormTransformation *offd = dynamic_cast<FreeFormTransformation *>(odof.get());
  if (offd) {
    if (!target_attr) {
      Warning("Cannot convert linear transformation to FFD without input -target image!");
      return false;
    }
    offd->Initialize(target_attr, dx, dy, dz, dt);
    if (omffd) {
      omffd->PushLocalTransformation(offd);
      odof.release();
      odof.reset(omffd.release());
    }

    BSplineFreeFormTransformationSV *bsvffd = nullptr;
    LinearFreeFormTransformationTD  *ltdffd = nullptr;
    BSplineFreeFormTransformationTD *btdffd = nullptr;

    (bsvffd = dynamic_cast<BSplineFreeFormTransformationSV *>(offd)) ||
    (ltdffd = dynamic_cast<LinearFreeFormTransformationTD  *>(offd)) ||
    (btdffd = dynamic_cast<BSplineFreeFormTransformationTD *>(offd));

    if (bsvffd) {
      if (ffdim.method != FFDIM_Unknown) {
        bsvffd->IntegrationMethod(ffdim.method);
      }
      if (ffdim.minsteps > 0) {
        bsvffd->NumberOfSteps(ffdim.minsteps);
      }
    } else if (ltdffd) {
      ltdffd->MinTimeStep(ffdim.MinTimeStep());
      ltdffd->MaxTimeStep(ffdim.MaxTimeStep());
    } else if (btdffd) {
      if (ffdim.method != FFDIM_Unknown) {
        btdffd->IntegrationMethod(ffdim.method);
      }
      btdffd->MinTimeStep(ffdim.MinTimeStep());
      btdffd->MaxTimeStep(ffdim.MaxTimeStep());
      btdffd->Tolerance  (ffdim.Tolerance());
    }
  }

  // Approximate input transformation by (M)FFD
  double rms = ApproximateAsNew(dof, odof.get(), &target_attr, ffdim, bch);
  if (IsNaN(rms)) return false;

  if (verbose > 1 && target_attr) {
    PrintApproximationError(dof, odof.get(), target_attr, ts, rms_excl_margin);
  } else if (verbose > 0) {
    cout << "RMS error of approximation = " << rms << endl;
  }

  // Write (M)FFD of requested FFD and MFFD types
  odof->Write(fname);
  return true;
}

// =============================================================================
// FSL
// =============================================================================

// -----------------------------------------------------------------------------
/// Convert linear transformation to FSL FLIRT matrix
///
/// @param[in] target Attributes of target image.
/// @param[in] source Attributes of source image.
/// @param[in] dof    Rigid/affine transformation.
///
/// @return FSL FLIRT transformation matrix.
///
/// @see fsl/src/newimage/newimagefns.cc : raw_affine_transform
Matrix ToFLIRTMatrix(ImageAttributes target,
                     ImageAttributes source,
                     const HomogeneousTransformation *dof)
{
  // x image axis must be mirrored if determinant of world to image orientation
  // matrix is positive, i.e., when the image is stored in neurological order
  if (target.GetWorldToLatticeOrientation().Det() > 0.) {
    if (verbose) cout << "qform determinant of target image is positive, reflecting x axis" << endl;
    target._xaxis[0] = -target._xaxis[0];
    target._xaxis[1] = -target._xaxis[1];
    target._xaxis[2] = -target._xaxis[2];
  } else {
    if (verbose) cout << "qform determinant of target image is negative, do not reflect x axis" << endl;
  }
  if (source.GetWorldToLatticeOrientation().Det() > 0.) {
    if (verbose) cout << "qform determinant of source image is positive, reflecting x axis" << endl;
    source._xaxis[0] = -source._xaxis[0];
    source._xaxis[1] = -source._xaxis[1];
    source._xaxis[2] = -source._xaxis[2];
  } else {
    if (verbose) cout << "qform determinant of source image is negative, do not reflect x axis" << endl;
  }

  Matrix t2w = target.GetImageToWorldMatrix();
  Matrix w2s = source.GetWorldToImageMatrix();

  Matrix S1(4, 4);
  S1(0, 0) = 1.0 / target._dx;
  S1(1, 1) = 1.0 / target._dy;
  S1(2, 2) = 1.0 / target._dz;
  S1(3, 3) = 1.0;

  Matrix S2(4, 4);
  S2(0, 0) = source._dx;
  S2(1, 1) = source._dy;
  S2(2, 2) = source._dz;
  S2(3, 3) = 1.0;

  Matrix A = S2 * w2s * dof->GetMatrix() * t2w * S1;
  return A.Inverse();
}

// -----------------------------------------------------------------------------
/// Write 4x4 transformation matrix to FLIRT ASCII file format
bool WriteFLIRTMatrix(const char *fname, const Matrix &m)
{
  ofstream os(fname);
  if (!os) {
    Warning("Failed to open file " << fname << " for writing");
    return false;
  }
  for (int i = 0; i < m.Rows(); ++i) {
    for (int j = 0; j < m.Cols(); ++j) {
      if (j > 0) os << " ";
      os << m(i, j);
    }
    os << "\n";
  }
  os.close();
  return bool(os);
}

// -----------------------------------------------------------------------------
/// Write FLIRT transformation file
bool WriteFLIRT(const char *fname, const ImageAttributes &target_attr,
                                   const ImageAttributes &source_attr,
                                   const Transformation  *dof)
{
  const HomogeneousTransformation *lin = dynamic_cast<const HomogeneousTransformation *>(dof);
  if (!lin) {
    const MultiLevelTransformation *mffd = dynamic_cast<const MultiLevelTransformation *>(dof);
    if (mffd) lin = mffd->GetGlobalTransformation();
    else {
      Warning("Cannot write non-linear transformation in Aladin format");
      return false;
    }
  }
  if (!target_attr) {
    Warning("Cannot convert linear transformation to Aladin format without input -target image!");
    return false;
  }
  if (!source_attr) {
    Warning("Cannot convert linear transformation to Aladin format without input -source image!");
    return false;
  }
  return WriteFLIRTMatrix(fname, ToFLIRTMatrix(target_attr, source_attr, lin));
}

// -----------------------------------------------------------------------------
/// Convert 3D+t displacement field to FSL warp
class ConvertToFSLWarp : public VoxelFunction
{
  const ImageAttributes &_Source;
  const BaseImage       &_Warp;
  bool                   _Relative;

  const Matrix _World2Source;

  static const int _x = 0;
  const int        _y, _z;

public:

  ConvertToFSLWarp(const ImageAttributes &source, const BaseImage &warp, bool relative)
  :
    _Source(source), _Warp(warp), _Relative(relative),
    _World2Source(source.GetWorldToImageMatrix()),
    _y(warp.NumberOfSpatialVoxels()), _z(_y + _y)
  {}

  template <class TIn, class TOut>
  void operator ()(int i, int j, int k, int, const TIn *din, TOut *dout) const
  {
    double x = i, y = j, z = k;
    _Warp.ImageToWorld(x, y, z);
    x += static_cast<double>(din[_x]);
    y += static_cast<double>(din[_y]);
    z += static_cast<double>(din[_z]);
    Transform(_World2Source, x, y, z);
    x *= _Source._dx, y *= _Source._dy, z *= _Source._dz;
    if (_Relative) {
      x -= i * _Warp.XSize();
      y -= j * _Warp.YSize();
      z -= k * _Warp.ZSize();
    }
    dout[_x] = static_cast<TOut>(x);
    dout[_y] = static_cast<TOut>(y);
    dout[_z] = static_cast<TOut>(z);
  }

  template <class TInOut>
  void operator ()(int i, int j, int k, int l, TInOut *d) const
  {
    this->operator()(i, j, k, l, d, d);
  }
};

// -----------------------------------------------------------------------------
/// Convert any transformation to FSL warp field
///
/// @param[in] target   Attributes of target image.
/// @param[in] source   Attributes of source image.
/// @param[in] dof      IRTK transformation.
/// @param[in] relative Return displacements instead of new voxel coordinates.
///
/// @return FSL warp field.
template <class TReal = double>
GenericImage<TReal> ToFSLWarpField(ImageAttributes target,
                                   ImageAttributes source,
                                   const Transformation *dof,
                                   bool relative = true)
{
  // x image axis must be mirrored if determinant of world to image orientation
  // matrix is positive, i.e., when the image is stored in neurological order
  if (target.GetWorldToLatticeOrientation().Det() > 0.) {
    if (verbose) cout << "qform determinant of target image is positive, reflecting x axis" << endl;
    target._xaxis[0] = -target._xaxis[0];
    target._xaxis[1] = -target._xaxis[1];
    target._xaxis[2] = -target._xaxis[2];
  }
  if (source.GetWorldToLatticeOrientation().Det() > 0.) {
    if (verbose) cout << "qform determinant of source image is positive, reflecting x axis" << endl;
    source._xaxis[0] = -source._xaxis[0];
    source._xaxis[1] = -source._xaxis[1];
    source._xaxis[2] = -source._xaxis[2];
  }

  // Get world displacement field
  GenericImage<TReal> warp(target, 3);
  dof->Displacement(warp);

  // Convert to FSL warp field
  ParallelForEachVoxel(ConvertToFSLWarp(source, warp, relative), warp.Attributes(), warp);

  // Set _dt != 0 such that warp field is saved as 3D+t time series instead of
  // 3D vector field, i.e., such that NIfTI dim[8] = {4, _x, _y, _z, 3, 1, 1, 1}
  warp.PutTSize(1.0);

  return warp;
}

// -----------------------------------------------------------------------------
/// Write FSL warp file
bool WriteFSLWarpField(const char *fname, ImageAttributes        target,
                                          const ImageAttributes &source,
                                          const Transformation  *dof,
                                          bool relative = true)
{
  if (!target) {
    // FIXME: Warp field not correct when voxel size/image orientation different
    //        from the target/reference image used for FSL's applywarp.
#if 1
    Warning("Cannot convert transformation to FSL warp field without input -target image!");
    return false;
#else
    const MultiLevelTransformation *mffd = dynamic_cast<const MultiLevelTransformation *>(dof);
    const FreeFormTransformation   *ffd  = dynamic_cast<const FreeFormTransformation   *>(dof);
    if (mffd) ffd = mffd->GetLocalTransformation(-1);
    if (!ffd) {
      Warning("Cannot convert linear transformation to FSL warp field without input -target image!");
      return false;
    }
    target = ffd->Attributes();
#endif
  }
  if (!source) {
    Warning("Cannot convert transformation to FSL warp field without input -source image!");
    return false;
  }
  ToFSLWarpField<double>(target, source, dof, relative).Write(fname);
  return true;
}

// =============================================================================
// Nifty Reg
// =============================================================================

// -----------------------------------------------------------------------------
/// Write Aladin transformation file
bool WriteAladin(const char *fname, const Transformation *dof)
{
  const HomogeneousTransformation *lin = dynamic_cast<const HomogeneousTransformation *>(dof);
  if (!lin) {
    const MultiLevelTransformation *mffd = dynamic_cast<const MultiLevelTransformation *>(dof);
    if (mffd) lin = mffd->GetGlobalTransformation();
    else {
      Warning("Cannot write non-linear transformation in Aladin format");
      return false;
    }
  }
  Matrix m = lin->GetMatrix();
  FILE *f = OpenFile(fname, "w");
  fprintf(f, "%lf %lf %lf %lf\n", m(0, 0), m(0, 1), m(0, 2), m(0, 3));
  fprintf(f, "%lf %lf %lf %lf\n", m(1, 0), m(1, 1), m(1, 2), m(1, 3));
  fprintf(f, "%lf %lf %lf %lf\n", m(2, 0), m(2, 1), m(2, 2), m(2, 3));
  fprintf(f, "0 0 0 1\n");
  fclose(f);
  return true;
}

// -----------------------------------------------------------------------------
/// Write F3D transformation file
bool WriteF3D(const char *fname, const Transformation *dof,
              F3DTransformationType type = F3D_TYPE_UNKNOWN)
{
  const HomogeneousTransformation        *lin  = NULL;
  const BSplineFreeFormTransformation3D  *ffd  = NULL;
  const MultiLevelFreeFormTransformation *mffd = NULL;

  (lin  = dynamic_cast<const HomogeneousTransformation        *>(dof)) ||
  (ffd  = dynamic_cast<const BSplineFreeFormTransformation3D  *>(dof)) ||
  (mffd = dynamic_cast<const MultiLevelFreeFormTransformation *>(dof));

  UniquePtr<MultiLevelFreeFormTransformation> new_mffd;
  if (lin) {
    new_mffd.reset(new MultiLevelFreeFormTransformation());
    new_mffd->GetGlobalTransformation()->PutMatrix(lin->GetMatrix());
    mffd = new_mffd.get();
  }

  AffineTransformation aff;
  if (mffd) {
    aff.PutMatrix(mffd->GetGlobalTransformation()->GetMatrix());
    if (mffd->NumberOfLevels() == 1) {
      ffd = dynamic_cast<const BSplineFreeFormTransformation3D *>(mffd->GetLocalTransformation(0));
    }
  }

  if (!ffd) {
    // TODO: - B-spline SV FFD --> SPLINE_VEL_GRID
    //       - Linear SV FFD   --> DISP_VEL_FIELD or DEF_VEL_FIELD
    //       - Linear FFD      --> DISP_FIELD or DEF_FIELD
    Warning("Input transformation must be 3D B-spline FFD or MFFD with more than one level.");
    Warning("Conversion of SV FFD and Linear FFD not implemented.");
    return false;
  }

  if (type != F3D_TYPE_UNKNOWN && type != F3D_SPLINE_GRID) {
    Warning("Cannot convert between F3D types! Can write B-spline (M)FFD as f3d_spline_grid only.");
    return false;
  }

  GenericImage<double> cpp(ffd->Attributes(), ffd->Z() == 1 ? 2 : 3);
  for (int k = 0; k < cpp.Z(); ++k)
  for (int j = 0; j < cpp.Y(); ++j)
  for (int i = 0; i < cpp.X(); ++i) {
    double x = i, y = j, z = k;
    cpp.ImageToWorld(x, y, z);
    aff.Transform(x, y, z);
    double dx, dy, dz;
    ffd->Get(i, j, k, dx, dy ,dz);
    cpp(i, j, k, 0) = x + dx;
    cpp(i, j, k, 1) = y + dy;
    if (cpp.T() > 2) cpp(i, j, k, 2) = z + dz;
  }
  // TODO: Set intent_code = NIFTI_INTENT_VECTOR, intent_name = NREG_TRANS, intent_p1 = F3D_SPLINE_GRID
  cpp.Write(fname);
  return true;
}

// =============================================================================
// STAR-CCM+
// =============================================================================

// -----------------------------------------------------------------------------
/// Write transformation for import into STAR-CCM+ GUI using XYZ Table format
bool WriteSTARCCMTable(const char *fname, const ImageAttributes &target, const Transformation *dof,
                       double      tmin          = -numeric_limits<double>::infinity(),
                       double      tmax          = +numeric_limits<double>::infinity(),
                       double      dt            = .0,
                       bool        displacements = false,
                       const char *points_name   = nullptr,
                       const char *delimiter     = " ",
                       int         precision     = 5)
{
  double t, t0 = target._torigin;

  // ---------------------------------------------------------------------------
  // Sampling grid
  ImageAttributes domain = target;
  if (!domain) {
    const FreeFormTransformation           *ffd  = nullptr;
    const MultiLevelFreeFormTransformation *mffd = nullptr;

    (ffd  = dynamic_cast<const FreeFormTransformation           *>(dof)) ||
    (mffd = dynamic_cast<const MultiLevelFreeFormTransformation *>(dof));
    if (mffd) ffd = mffd->GetLocalTransformation(-1);

    if (ffd) {
      domain = ffd->Attributes();
    } else if (!points_name) {
      Warning("Cannot convert linear transformation to STAR-CCM+ Table without input -target image or -points!");
      return false;
    }
  }

  if (IsInf(tmin)) tmin = domain.LatticeToTime(0);
  if (IsInf(tmax)) tmax = domain.LatticeToTime(domain._t - 1);
  if (dt == 0.)    dt   = domain._dt;

  domain._torigin = tmin;
  domain._t       = ifloor(abs((tmax - tmin) / dt));
  domain._dt      = dt;

  // Open text file
  ofstream table(fname);
  if (!table) {
    Warning("Failed to open file " << fname << " for writing!");
    return false;
  }

  table.precision(precision);
  table.flags(ios::fixed);

  // ---------------------------------------------------------------------------
  // Write header
  if (displacements) {
    table << "X" << delimiter << "Y" << delimiter << "Z";
    if (domain._t == 1) {
      table << delimiter << "DX" << delimiter << "DY" << delimiter << "DZ";
    } else {
      for (int l = 0; l < domain._t; ++l) {
        t = domain.LatticeToTime(l);
        table << delimiter << "DX[t=" << t << "ms]";
        table << delimiter << "DY[t=" << t << "ms]";
        table << delimiter << "DZ[t=" << t << "ms]";
      }
    }
  } else {
    if (domain._t == 1) {
      table << "X[t=0ms]" << delimiter << "Y[t=0ms]" << delimiter << "Z[t=0ms]";
      table << delimiter;
      table << "X[t=1ms]" << delimiter << "Y[t=1ms]" << delimiter << "Z[t=1ms]";
    } else {
      for (int l = 0, c = 0; l < domain._t; ++l) {
        t = domain.LatticeToTime(l);
        if (++c > 1) table << delimiter;
        table << "X[t=" << t << "ms]";
        table << delimiter << "Y[t=" << t << "ms]";
        table << delimiter << "Z[t=" << t << "ms]";
      }
    }
  }
  table << "\n";

  if (!table) return false;

  // ---------------------------------------------------------------------------
  // Either write displacement/new position for each input point
  if (points_name) {

    PointSet points;
    points.Read(points_name);

    double dx, dy, dz;
    for (int r = 0; r < points.Size(); ++r) {
      if (displacements) {
        table << points(r)._x << delimiter << points(r)._y << delimiter << points(r)._z;
        table << delimiter;
      }
      for (int l = 0, c = 0; l < domain._t; ++l) {
        t = domain.LatticeToTime(l);
        if (++c > 1) table << delimiter;
        dx = points(r)._x, dy = points(r)._y, dz = points(r)._z;
        if (displacements) {
          dof->Displacement(dx, dy, dz, t, t0);
        } else {
          dof->Transform(dx, dy, dz, t, t0);
        }
        table << dx << delimiter << dy << delimiter << dz;
      }
      table << "\n";

      if (!table) return false;
    }

  // ---------------------------------------------------------------------------
  // or write displacement/new position for each lattice point
  } else {

    // Evaluate displacements
    WorldCoordsImage wc;
    Array<RealImage> disp(domain._t);
    for (int l = 0; l < domain._t; ++l) {
      t = domain.LatticeToTime(l);
      disp[l].Initialize(domain, 3);
      disp[l].PutTOrigin(t);
      if (wc.IsEmpty()) disp[l].ImageToWorld(wc);
      dof->Displacement(disp[l], t0, &wc);
    }
    if (wc.IsEmpty()) {
      Warning("Invalid time interval [" << tmin << " " << tmax << "]");
      return false;
    }

    // Write rows
    double dx, dy, dz;
    const int nvox = domain.NumberOfSpatialPoints();
    const WorldCoordsImage::VoxelType *x = wc.Data();
    const WorldCoordsImage::VoxelType *y = x + nvox;
    const WorldCoordsImage::VoxelType *z = y + nvox;

    for (int r = 0; r < nvox; ++r, ++x, ++y, ++z) {
      if (displacements) {
        table << *x << delimiter << *y << delimiter << *z;
        table << delimiter;
      }
      for (int l = 0, c = 0; l < domain._t; ++l) {
        t = domain.LatticeToTime(l);
        if (++c > 1) table << delimiter;
        dx = *(disp[l].Data() + r           );
        dy = *(disp[l].Data() + r +     nvox);
        dz = *(disp[l].Data() + r + 2 * nvox);
        if (!displacements) {
          dx += *x, dy += *y, dz += *z;
        }
        table << dx << delimiter << dy << delimiter << dz;
      }
      table << "\n";

      if (!table) return false;
    }

  }

  // Close file
  table.close();
  return !table.fail();
}

////////////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // Parse arguments
  REQUIRES_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  if (NUM_POSARGS != 2 && NUM_POSARGS != 4) {
    FatalError("Invalid number of positional arguments!");
  }

  TransformationFileFormat format_in  = Format_Unknown;
  TransformationFileFormat format_out = Format_Unknown;
  MFFDMode                 mffd_type  = MFFD_Default;

  const char *dofin_name  = nullptr;
  const char *target_name = nullptr;
  const char *source_name = nullptr;
  const char *points_name = nullptr;
  int         xyz_units   = 0;
  double      dx          = .0;
  double      dy          = .0;
  double      dz          = .0;
  double      dt          = .0;
  double      t0          = numeric_limits<double>::quiet_NaN();
  double      ts          = numeric_limits<double>::quiet_NaN();
  double      tmin        = -numeric_limits<double>::infinity();
  double      tmax        = +numeric_limits<double>::infinity();
  const char *delimiter   = nullptr;
  int         precision   = -1;
  bool        velocities  = false;
  FFDIMParams ffdim;
  BCHParams   bchparam;

  #if MIRTK_IO_WITH_NIfTI
    xyz_units = NIFTI_UNITS_MM;
  #endif 

  for (ALL_OPTIONS) {
    if (OPTION("-input-format")) {
      PARSE_ARGUMENT(format_in);
    }
    else if (OPTION("-format") || OPTION("-output-format")) {
      PARSE_ARGUMENT(format_out);
    }
    else if (OPTION("-dofin")) {
      dofin_name = ARGUMENT;
    }
    else if (OPTION("-nomffd")) {
      mffd_type = MFFD_None;
    }
    else if (OPTION("-mffd")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(mffd_type);
      else              mffd_type = MFFD_Sum;
    }
    else if (OPTION("-fluid")) {
      mffd_type = MFFD_Fluid;
    }
    else if (OPTION("-msvffd")) {
      mffd_type  = MFFD_LogSum;
      format_out = Format_MIRTK_BSplineSVFFD;
    }
    else if (OPTION("-target")) target_name = ARGUMENT;
    else if (OPTION("-source")) source_name = ARGUMENT;
    else if (OPTION("-points")) points_name = ARGUMENT;
    else if (OPTION("-integration-method") || OPTION("-im") || OPTION("-ffdim")) {
      PARSE_ARGUMENT(ffdim.method);
    }
    else if (OPTION("-steps")) {
      PARSE_ARGUMENT(ffdim.minsteps);
      ffdim.maxsteps = ffdim.minsteps;
    }
    else if (OPTION("-min-steps") || OPTION("-minsteps")) {
      PARSE_ARGUMENT(ffdim.minsteps);
    }
    else if (OPTION("-max-steps") || OPTION("-maxsteps")) {
      PARSE_ARGUMENT(ffdim.maxsteps);
    }
    else if (OPTION("-Tt")   || OPTION("-t0")) PARSE_ARGUMENT(t0);
    else if (OPTION("-Ts")   || OPTION("-tT") || OPTION("-ts")) PARSE_ARGUMENT(ts);
    else if (OPTION("-tmin") || OPTION("-t1")) PARSE_ARGUMENT(tmin);
    else if (OPTION("-tmax") || OPTION("-t2")) PARSE_ARGUMENT(tmax);
    else if (OPTION("-ds")) {
      PARSE_ARGUMENT(dx);
      dy = dz = dx;
    }
    else if (OPTION("-dx")) PARSE_ARGUMENT(dx);
    else if (OPTION("-dy")) PARSE_ARGUMENT(dy);
    else if (OPTION("-dz")) PARSE_ARGUMENT(dz);
    else if (OPTION("-dt")) PARSE_ARGUMENT(dt);
    else if (OPTION("-bch-smooth")) bchparam.smooth = true;
    else if (OPTION("-bch-terms")) PARSE_ARGUMENT(bchparam.nterms);
    else if (OPTION("-bch-steps")) PARSE_ARGUMENT(bchparam.nsteps);
    else if (OPTION("-delimiter") || OPTION("-delim")) delimiter = ARGUMENT;
    else if (OPTION("-precision")) PARSE_ARGUMENT(precision);
    else if (OPTION("-xyz-units") || OPTION("-xyz_units")) {
      #if MIRTK_IO_WITH_NIfTI
        PARSE_ARGUMENT(xyz_units);
      #else
        ARGUMENT; // unused, but needs to be parsed
      #endif // MIRTK_IO_WITH_NIfTI
    }
    else HANDLE_BOOL_OPTION(velocities);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (precision < 0) precision = 5;

  const string ext_out       = Extension(output_name, EXT_LastWithoutGz);
  const bool   ext_out_nifti = (ext_out == ".nii" || ext_out == ".hdr" || ext_out == ".img");
  const bool   ext_out_mirtk = (ext_out == ".dof");

  // Choose default output format based on output file name extension
  if (format_out == Format_Unknown) {
    if      (ext_out_nifti) format_out = Format_WorldDisplacement;
    else if (ext_out_mirtk) format_out = Format_MIRTK;
    else if (ext_out == ".csv") format_out = Format_CSV;
    else if (ext_out == ".tsv" || ext_out == ".txt") format_out = Format_TSV;
    else {
      FatalError("No default output format available for extension " << ext_out << ", use [-output]-format option!");
    }
  }

  if (format_out == Format_IRTK_LinearFFD || format_out == Format_IRTK_BSplineFFD) {
    FatalError("Cannot yet write deformable IRTK transformation file!");
  }
  if (format_out == Format_MIRTK_LinearSVFFD) {
    FatalError("Cannot yet write linear SV FFD file!");
  }

  // Read target image attributes
  ImageAttributes target_attr;
  if (target_name) {
    InitializeIOLibrary();
    BinaryImage target(target_name);
    target_attr = target.Attributes();
  } else {
    if (IsNaN(t0)) t0 = .0;
  }
  if (IsNaN(t0)) t0 = target_attr._torigin;
  else           target_attr._torigin = t0;

  // Read source image attributes
  ImageAttributes source_attr;
  if (source_name) {
    InitializeIOLibrary();
    BinaryImage source(source_name);
    source_attr = source.Attributes();
  } else {
    if (IsNaN(ts)) ts = 1.0;
  }
  if (IsNaN(ts)) ts = source_attr._torigin;
  else           source_attr._torigin = ts;

  ffdim.t1 = t0;
  ffdim.t2 = ts;

  // Guess input file format
  const string ext_in = Extension(input_name, EXT_LastWithoutGz); 
  bool  ext_in_nifti  = (ext_in == ".nii" || ext_in == ".hdr" || ext_in == ".img");

  if (format_in == Format_Unknown) {
    if (Transformation::CheckHeader(input_name)) {
      format_in = Format_MIRTK;
    } else if (ext_in == ".xfm") {
      format_in = Format_MNI_XFM;
    } else if (ext_in == ".m3z") {
      format_in = Format_MNI_M3Z;
    } else if (ext_in == ".csv") {
      format_in = Format_CSV;
    } else if (ext_in == ".tsv" || ext_in == ".txt") {
      format_in = Format_TSV;
    }
  } else if (format_in == Format_FSL) {
    if (ext_in_nifti) {
      NiftiImageInfo info(input_name);
      switch (info.intent_code) {
        case FSL_FNIRT_DISPLACEMENT_FIELD:      format_in = Format_FSL_FNIRT_Displacement; break;
        case FSL_DCT_COEFFICIENTS:              format_in = Format_FSL_DctCoefficients; break;
        case FSL_CUBIC_SPLINE_COEFFICIENTS:     format_in = Format_FSL_CubicSpline; break;
        case FSL_QUADRATIC_SPLINE_COEFFICIENTS: format_in = Format_FSL_QuadraticSpline; break;
        default:
          FatalError("Unknown FSL NIfTI intent code: " << info.intent_code);
      }
    }
    else format_in = Format_FSL_FLIRT;
  } else if (format_in == Format_NREG) {
    if (ext_in_nifti) format_in = Format_F3D;
    else              format_in = Format_Aladin;
  } else if (format_in == Format_MNI) {
    if (ext_in == ".xfm") format_in = Format_MNI_XFM;
    else                  format_in = Format_MNI_M3Z;
  }

  const char *dx_input = nullptr;
  const char *dy_input = nullptr;
  const char *dz_input = nullptr;

  if (format_in == Format_WorldDisplacement || format_in == Format_VoxelDisplacement ||
      format_in == Format_WorldSVF || format_in == Format_VoxelSVF ||
      format_in == Format_ITKDisplacement || format_in == Format_ITKSVF) {
    if (NUM_POSARGS == 4 || NUM_POSARGS == 6) {
      dx_input    = POSARG(1);
      dy_input    = POSARG(2);
      dz_input    = POSARG(3);
      input_name  = nullptr;
      output_name = POSARG(4);
    }
  } else {
    if (NUM_POSARGS == 6) {
      FatalError("Invalid number of positional arguments!");
    }
  }

  // Read input transformation
  UniquePtr<Transformation> dof;
  switch (format_in) {

    // Dense displacement field
    case Format_WorldDisplacement: {
      if (input_name) {
        dof.reset(ReadWorldDisplacement(input_name));
      } else {
        dof.reset(ReadWorldDisplacement(dx_input, dy_input, dz_input));
      }
    } break;
    case Format_VoxelDisplacement: {
      if (input_name) {
        dof.reset(ReadVoxelDisplacement(input_name));
      } else {
        dof.reset(ReadVoxelDisplacement(dx_input, dy_input, dz_input));
      }
    } break;
    case Format_ITKDisplacement: {
      if (input_name) {
        dof.reset(ReadITKDisplacement(input_name));
      } else {
        dof.reset(ReadITKDisplacement(dx_input, dy_input, dz_input));
      }
    } break;

    // Dense stationary velocity field
    case Format_WorldSVF: {
      if (input_name) {
        dof.reset(ReadWorldSVF(input_name));
      } else {
        dof.reset(ReadWorldSVF(dx_input, dy_input, dz_input));
      }
    } break;
    case Format_VoxelSVF: {
      if (input_name) {
        dof.reset(ReadVoxelSVF(input_name));
      } else {
        dof.reset(ReadVoxelSVF(dx_input, dy_input, dz_input));
      }
    } break;
    case Format_ITKSVF: {
      if (input_name) {
        dof.reset(ReadITKSVF(input_name));
      } else {
        dof.reset(ReadITKSVF(dx_input, dy_input, dz_input));
      }
    } break;

    // MIRTK
    case Format_MIRTK:
    case Format_MIRTK_Rigid:
    case Format_MIRTK_Similarity:
    case Format_MIRTK_Affine:
    case Format_MIRTK_LinearFFD:
    case Format_MIRTK_LinearSVFFD:
    case Format_MIRTK_LinearTDFFD:
    case Format_MIRTK_BSplineFFD:
    case Format_MIRTK_BSplineSVFFD:
    case Format_MIRTK_BSplineTDFFD: {
      dof.reset(ReadMIRTK(input_name));
      if (format_in != Format_MIRTK && dof && dof->TypeOfClass() != ToMIRTKTransformationType(format_in)) {
        Warning("Type of input transformation differs from specified -input-format! Ignoring option.");
      }
    } break;

    // IRTK
    case Format_IRTK: {
      dof.reset(ReadIRTK(input_name));
    } break;

    // FSL
    case Format_FSL_FLIRT: {
      dof.reset(ReadFLIRT(input_name, target_attr, source_attr));
    } break;

    case Format_FSL_FNIRT_Displacement: {
      dof.reset(ReadFNIRTDisplacement(input_name));
    } break;

    // Nifty Reg
    case Format_Aladin: {
      dof.reset(ReadAladin(input_name));
    } break;

    case Format_F3D:
    case Format_F3D_DEF_FIELD:
    case Format_F3D_DISP_FIELD:
    case Format_F3D_SPLINE_GRID:
    case Format_F3D_DEF_VEL_FIELD:
    case Format_F3D_DISP_VEL_FIELD:
    case Format_F3D_SPLINE_VEL_GRID: {
      dof.reset(ReadF3D(input_name, dofin_name, xyz_units, ffdim.minsteps,
                        ToF3DTransformationType(format_in)));
    } break;

    // Elastix
    case Format_Elastix:
    case Format_Elastix_FFD: {
      dof.reset(ReadElastix(input_name, target_attr));
    } break;

    // FreeSurfer
    case Format_MNI_XFM: {
      dof.reset(ReadXFM(input_name));
    } break;

    // DRAMMS
    case Format_DRAMMS: {
      dof.reset(ReadDRAMMS(input_name));
    } break;

    // STAR-CCM+
    case Format_STAR_CCM:
    case Format_STAR_CCM_Table:
    case Format_STAR_CCM_Table_XYZ: {
      FatalError("Cannot read transformation from STAR-CCM+ file!");
    } break;

    // CSV/TSV table
    case Format_CSV:
    case Format_CSV_XYZ:
    case Format_TSV:
    case Format_TSV_XYZ: {
      FatalError("Cannot read transformation from CSV/TSV file!");
    } break;

    // Still unknown...
    default:
      FatalError("Unknown input transformation file format! Use -input-format option.");
  }

  if (!dof) {
    FatalError("Failed to convert transformation! Either required additional input"
               " format options missing or memory allocation failed.");
  }

  // Guess output file format
  bool is_linear = (dynamic_cast<HomogeneousTransformation *>(dof.get()) != nullptr);

  if (format_out == Format_FSL) {
    if (ext_out_nifti || !is_linear) format_out = Format_FSL_WarpRelative;
    else                             format_out = Format_FSL_FLIRT;
  } else if (format_out == Format_NREG) {
    if (ext_out_nifti || !is_linear) format_out = Format_F3D;
    else                             format_out = Format_Aladin;
  } else if (format_out == Format_STAR_CCM) {
    format_out = Format_STAR_CCM_Table;
  }

  const char *dx_output = nullptr;
  const char *dy_output = nullptr;
  const char *dz_output = nullptr;

  if ((format_in  == Format_WorldDisplacement || format_in  == Format_VoxelDisplacement || format_out == Format_WorldSVF || format_out == Format_VoxelSVF) &&
      (format_out == Format_WorldDisplacement || format_out == Format_VoxelDisplacement || format_out == Format_WorldSVF || format_out == Format_VoxelSVF)) {
    if (NUM_POSARGS == 4) {
      FatalError("Ambiguous usage when both input and output formats are dense vector fields\n"
                 " and 4 positional arguments given. Either read *and* write components from/to\n"
                 " separate image files (6 arguments) or from/to a single 3D+3 image.\n");
    }
  }
  if (format_out == Format_WorldDisplacement || format_out == Format_VoxelDisplacement || format_out == Format_WorldSVF || format_out == Format_VoxelSVF) {
    if (NUM_POSARGS == 4) {
      if (dx_input || dy_input || dz_input) {
        FatalError("Ambiguous usage when both input and output formats are dense vector fields\n"
                   " and 4 positional arguments given. Either read *and* write components from/to\n"
                   " separate image files (6 arguments) or from/to a single 3D+3 image.\n");
      }
      dx_output = POSARG(2);
      dy_output = POSARG(3);
      dz_output = POSARG(4);
      output_name = nullptr;
    } else if (NUM_POSARGS == 6) {
      dx_output = POSARG(4);
      dy_output = POSARG(5);
      dz_output = POSARG(6);
      output_name = nullptr;
    }
  }

  // Cast velocity based MIRTK transformation to displacement based MIRTK transformation
  if (format_out == Format_WorldSVF || format_out == Format_VoxelSVF) {
    auto svffd = dynamic_cast<BSplineFreeFormTransformationSV *>(dof.get());
    if (svffd == nullptr) {
      auto mffd = dynamic_cast<MultiLevelFreeFormTransformation *>(dof.get());
      if (mffd && mffd->NumberOfLevels() == 1) {
        svffd = dynamic_cast<BSplineFreeFormTransformationSV *>(mffd->GetLocalTransformation(0));
      }
      if (svffd == nullptr) {
        FatalError("Output format 'svf_world' and 'svf_voxel' requires input SVFFD!");
      }
    }
    dof.reset(new BSplineFreeFormTransformation3D(*svffd));
    if (format_out == Format_WorldSVF) format_out = Format_WorldDisplacement;
    else                               format_out = Format_VoxelDisplacement;
  }

  // Write transformation in requested output format
  bool success;
  switch (format_out) {

    // Dense displacement field or stationary velocity field, respectively
    case Format_WorldDisplacement: {
      if (output_name) {
        success = WriteWorldDisplacement(output_name, dof.get(), target_attr, ts);
      } else {
        success = WriteWorldDisplacement(dx_output, dy_output, dz_output, dof.get(), target_attr, ts);
      }
    } break;
    case Format_VoxelDisplacement: {
      if (output_name) {
        success = WriteVoxelDisplacement(output_name, dof.get(), target_attr, ts);
      } else {
        success = WriteVoxelDisplacement(dx_output, dy_output, dz_output, dof.get(), target_attr, ts);
      }
    } break;

    // MIRTK
    case Format_MIRTK:
    case Format_MIRTK_Rigid:
    case Format_MIRTK_Similarity:
    case Format_MIRTK_Affine:
    case Format_MIRTK_LinearFFD:
    case Format_MIRTK_LinearTDFFD:
    case Format_MIRTK_BSplineFFD:
    case Format_MIRTK_BSplineSVFFD:
    case Format_MIRTK_BSplineTDFFD: {
      success = WriteMIRTK(output_name, dof.get(), target_attr, ts, dx, dy, dz, dt,
                           ToMIRTKTransformationType(format_out), mffd_type,
                           ffdim, bchparam);
    } break;

    case Format_MIRTK_LinearSVFFD: {
      // FIXME: Once LinearFreeFormTransformationSV is implemented!
      FatalError("Cannot create linear MIRTK SV FFD transformation file!");
    } break;

    // IRTK
    case Format_IRTK_Rigid: {
      success = WriteMIRTK(output_name, dof.get(), target_attr, ts, dx, dy, dz, dt,
                           TRANSFORMATION_RIGID, MFFD_None);
    } break;
    case Format_IRTK_Affine: {
      success = WriteMIRTK(output_name, dof.get(), target_attr, ts, dx, dy, dz, dt,
                           TRANSFORMATION_AFFINE, MFFD_None);
    } break;
    case Format_IRTK: {
      const RigidTransformation      *rig;
      const SimilarityTransformation *sim;
      const AffineTransformation     *aff;
      rig = dynamic_cast<const RigidTransformation      *>(dof.get());
      sim = dynamic_cast<const SimilarityTransformation *>(dof.get());
      aff = dynamic_cast<const AffineTransformation     *>(dof.get());
      if (aff || sim) {
        success = WriteMIRTK(output_name, dof.get(), target_attr, ts, dx, dy, dz, dt,
                             TRANSFORMATION_AFFINE, MFFD_None);
      } else if (rig) {
        success = WriteMIRTK(output_name, dof.get(), target_attr, ts, dx, dy, dz, dt,
                             TRANSFORMATION_RIGID, MFFD_None);
      } else {
        // TODO: Write (M)FFD in old IRTK format
        FatalError("Cannot write deformable IRTK transformation!");
      }
    } break;

    case Format_IRTK_LinearFFD:
    case Format_IRTK_BSplineFFD: {
      // TODO: Write (M)FFD in old IRTK format
      FatalError("Cannot write deformable IRTK transformation!");
    } break;

    // FSL
    case Format_FSL_FLIRT: {
      success = WriteFLIRT(output_name, target_attr, source_attr, dof.get());
    } break;

    case Format_FSL_WarpAbsolute:
    case Format_FSL_WarpRelative: {
      success = WriteFSLWarpField(output_name, target_attr, source_attr, dof.get(),
                                  format_out == Format_FSL_WarpRelative);
    } break;

    case Format_FSL_FNIRT_Displacement: {
      success = false;
    } break;

    // NiftyReg
    case Format_Aladin: {
      success = WriteAladin(output_name, dof.get());
    } break;

    case Format_F3D:
    case Format_F3D_DEF_FIELD:
    case Format_F3D_DISP_FIELD:
    case Format_F3D_SPLINE_GRID:
    case Format_F3D_DEF_VEL_FIELD:
    case Format_F3D_DISP_VEL_FIELD:
    case Format_F3D_SPLINE_VEL_GRID: {
      success = WriteF3D(output_name, dof.get(), ToF3DTransformationType(format_out));
    } break;

    // STAR-CCM+
    case Format_STAR_CCM_Table: {
      const bool disps = true;
      if (delimiter == nullptr) delimiter = " ";
      success = WriteSTARCCMTable(output_name, target_attr, dof.get(),
                                  tmin, tmax, dt, disps, points_name,
                                  delimiter, precision);
    } break;

    case Format_STAR_CCM_Table_XYZ: {
      const bool disps = false;
      if (delimiter == nullptr) delimiter = " ";
      success = WriteSTARCCMTable(output_name, target_attr, dof.get(),
                                  tmin, tmax, dt, disps, points_name,
                                  delimiter, precision);
    } break;

    // CSV/TSV table
    case Format_CSV: {
      const bool disps = true;
      if (delimiter == nullptr) delimiter = ",";
      success = WriteSTARCCMTable(output_name, target_attr, dof.get(),
                                  tmin, tmax, dt, disps, points_name,
                                  delimiter, precision);
    } break;

    case Format_CSV_XYZ: {
      const bool disps = false;
      if (delimiter == nullptr) delimiter = ",";
      success = WriteSTARCCMTable(output_name, target_attr, dof.get(),
                                  tmin, tmax, dt, disps, points_name,
                                  delimiter, precision);
    } break;

    case Format_TSV: {
      const bool disps = true;
      if (delimiter == nullptr) delimiter = "\t";
      success = WriteSTARCCMTable(output_name, target_attr, dof.get(),
                                  tmin, tmax, dt, disps, points_name,
                                  delimiter, precision);
    } break;

    case Format_TSV_XYZ: {
      const bool disps = false;
      if (delimiter == nullptr) delimiter = "\t";
      success = WriteSTARCCMTable(output_name, target_attr, dof.get(),
                                  tmin, tmax, dt, disps, points_name,
                                  delimiter, precision);
    } break;

    // Unknown?
    default: success = false;
  }
  if (!success) {
    FatalError("Cannot write transformation in file format: " << ToString(format_out));
  }

  return 0;
}

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

#include <mirtkCommon.h>
#include <mirtkOptions.h>

#include <mirtkImageIOConfig.h>
#include <mirtkGenericImage.h>
#include <mirtkVoxelFunction.h>
#include <mirtkTransformations.h>

#if MIRTK_ImageIO_WITH_NIfTI
#  include <mirtkNiftiImageInfo.h>
#  include <mirtkNiftiImageReader.h>
#endif // MIRTK_ImageIO_WITH_NIfTI


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
  cout << "  =====================  =================================================================================\n";
  cout << "  unknown                Unknown, try to guess it from file header/type.\n";
  cout << "  disp_world|disp|image  Dense displacement field image with world space displacement vectors [mm].\n";
  cout << "  disp_voxel             Dense displacement field image with target space displacement vectors [voxel].\n";
  cout << "  mirtk|irtk2            MIRTK/IRTK 2 transformation file format.\n";
  cout << "  mirtk_rigid|rigid      Rigid MIRTK/IRTK 2 transformation file format (6 DoFs).\n";
  cout << "  mirtk_similarity       Similarity MIRTK/IRTK 2 transformation file format (7 DoFs).\n";
  cout << "  mirtk_affine|affine    Affine MIRTK/IRTK 2 transformation file format (12 DoFs).\n";
  cout << "  mirtk_linear_ffd       Linear free-form deformation.\n";
  cout << "  mirtk_bspline_ffd      Cubic B-spline free-form deformation.\n";
  cout << "  mni_xfm|xfm            Linear FreeSurfer transformation (.xfm file).\n";
  cout << "  fsl                    Guess/choose FSL output file format.\n";
  cout << "  flirt                  FSL FLIRT output file format.\n";
  cout << "  fnirt                  FSL FNIRT output file format.\n";
  cout << "  nreg                   Guess/choose Nifty Reg transformation output file format.\n";
  cout << "  aladin                 Nifty Reg Aladin output file format.\n";
  cout << "  f3d                    Nifty Reg reg_f3d output file format with nifti1.intent_p1 code.\n";
  cout << "  f3d_def_field          Nifty Reg reg_f3d output image deformation  field.\n";
  cout << "  f3d_disp_field         Nifty Reg reg_f3d output image displacement field.\n";
  cout << "  f3d_spline_grid        Nifty Reg reg_f3d output control point displacement field.\n";
  cout << "  f3d_def_vel_field      Nifty Reg reg_f3d output image deformation  field as stationary velocity field.\n";
  cout << "  f3d_disp_vel_field     Nifty Reg reg_f3d output image displacement field as stationary velocity field.\n";
  cout << "  f3d_spline_vel_grid    Nifty Reg reg_f3d output control point velocity field.\n";
  cout << "  =====================  =================================================================================\n";
#if !MIRTK_ImageIO_WITH_NIfTI
  cout << "\n";
  cout << "  Cannot convert from/to the following formats because the ImageIO module is missing NIfTI support:\n";
  cout << "\n";
  cout << "    f3d*, fnirt\n";
  cout << "\n";
#endif // !MIRTK_ImageIO_WITH_NIfTI
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Input transformation file.\n";
  cout << "  output   Output transformation file.\n";
  cout << "\n";
  cout << "Options:\n";
  cout << "  -input-format  <format>    Format of input file. (default: unknown)\n";
  cout << "  -output-format <format>    Format of output file. (default: mirtk)\n";
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
  cout << "  -Ts <time>                 Time point of source image. Used by 3D+t, TD, and SV FFDs.\n";
  cout << "  -ds <value>                Output control point spacing. (default: input spacing)\n";
  cout << "  -dx <value>                Output control point spacing in x dimension. (default: input spacing)\n";
  cout << "  -dy <value>                Output control point spacing in y dimension. (default: input spacing)\n";
  cout << "  -dz <value>                Output control point spacing in z dimension. (default: input spacing)\n";
  cout << "  -scaling-steps <int>       Number of scaling and squaring steps.\n";
  cout << "                             In case of f3d_*_vel_*, use nifti1.intent_p2 by default.\n";
  cout << "  -xyz_units (m|mm|mu)       Spatial units of original target NIfTI header\n";
  cout << "                             if ignored by Nifty Reg's reg_f3d. (default: mm)\n";
  PrintStandardOptions(cout);
  cout << endl;
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
  // (M)IRTK
  Format_MIRTK,
  Format_MIRTK_Rigid,
  Format_MIRTK_Similarity,
  Format_MIRTK_Affine,
  Format_MIRTK_LinearFFD,
  Format_MIRTK_BSplineFFD,
  Format_LegacyIRTK,
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
  // FreeSurfer
  Format_MNI,
  Format_MNI_XFM,
  Format_MNI_M3Z,
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
    case Format_MIRTK:                    str = "mirtk"; break;
    case Format_MIRTK_Rigid:              str = "mirtk_rigid"; break;
    case Format_MIRTK_Similarity:         str = "mirtk_similarity"; break;
    case Format_MIRTK_Affine:             str = "mirtk_affine"; break;
    case Format_MIRTK_LinearFFD:          str = "mirtk_linear_ffd"; break;
    case Format_MIRTK_BSplineFFD:         str = "mirtk_bspline_ffd"; break;
    case Format_LegacyIRTK:               str = "legacy_irtk"; break;
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
    default:                              str = "unknown"; break;
  }
  return ToString(str, w, c, left);
}

// -----------------------------------------------------------------------------
/// Convert transformation file format string to enumeration value
template <>
inline bool FromString(const char *str, TransformationFileFormat &format)
{
  string format_name = ToLower(str);
  format = Format_Unknown;

  // Alternative format names
  if      (format_name == "irtk1") format = Format_LegacyIRTK;
  else if (format_name.compare(0, 5, "irtk2") == 0) format_name = "mirtk" + format_name.substr(5);
  else if (format_name == "disp" || format_name == "image") format = Format_WorldDisplacement;
  else if (format_name == "rigid") format = Format_MIRTK_Rigid;
  else if (format_name == "similarity") format = Format_MIRTK_Similarity;
  else if (format_name == "affine") format = Format_MIRTK_Affine;
  else if (format_name == "linear_ffd") format = Format_MIRTK_LinearFFD;
  else if (format_name == "bspline_ffd") format = Format_MIRTK_BSplineFFD;
  else if (format_name == "mirtk_ffd" || format_name == "ffd") format = Format_MIRTK_BSplineFFD;
  else if (format_name == "fsl_warp" || format_name == "warp") format = Format_FSL_WarpRelative;
  else if (format_name == "flirt") format = Format_FSL_FLIRT;
  else if (format_name == "fnirt") format = Format_FSL_FNIRT_Displacement;
  else if (format_name == "niftk" || format_name == "niftyreg") format = Format_NREG;
  else if (format_name == "freesurfer") format = Format_MNI;
  else if (format_name == "freesurfer_xfm" || format_name == "xfm") format = Format_MNI_XFM;
  else if (format_name == "freesurfer_m3z" || format_name == "m3z") format = Format_MNI_M3Z;
  else if (format_name.compare(0, 4, "reg_") == 0) format_name = format_name.substr(4);

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
    case Format_MIRTK_BSplineFFD:    return TRANSFORMATION_BSPLINE_FFD_3D;
    case Format_MIRTK_LinearFFD:     return TRANSFORMATION_LINEAR_FFD_3D;
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


} // namespace mirtk

////////////////////////////////////////////////////////////////////////////////
// Read transformation
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// (M)IRTK
// =============================================================================

// -----------------------------------------------------------------------------
/// Read transformation from MIRTK transformation file
Transformation *ReadMIRTK(const char *fname)
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
  if (dx <= .0) dx = disp.GetXSize();
  if (dy <= .0) dy = disp.GetYSize();
  if (dz <= .0) dz = disp.GetZSize();

  double xaxis[3], yaxis[3], zaxis[3];
  disp.GetOrientation(xaxis, yaxis, zaxis);
  unique_ptr<LinearFreeFormTransformation> ffd;
  ffd.reset(new LinearFreeFormTransformation3D(0, 0, 0, disp.X() - 1, disp.Y() - 1, disp.Z() - 1,
                                               dx, dy, dz, xaxis, yaxis, zaxis));

  if (fequal(dx, disp.GetXSize()) && fequal(dy, disp.GetYSize()) && fequal(dz, disp.GetZSize())) {
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
  } else {
    for (int k = 0; k < ffd->Z(); ++k)
    for (int j = 0; j < ffd->Y(); ++j)
    for (int i = 0; i < ffd->X(); ++i) {
      ffd->Put(i, j, k, disp(i, j, k, 0), disp(i, j, k, 1), disp(i, j, k, 2));
    }
  }

  unique_ptr<MultiLevelFreeFormTransformation> dof;
  dof.reset(new MultiLevelFreeFormTransformation());
  dof->PushLocalTransformation(ffd.release());
  return dof.release();
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
Transformation *ReadFLIRT(const char *fname, const ImageAttributes &target,
                                             const ImageAttributes &source)
{
  Matrix A   = ReadFLIRTMatrix(fname);
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

  unique_ptr<AffineTransformation> dof(new AffineTransformation());
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
  FILE *f = fopen(fname, "r");
  if (fscanf(f, "%lf %lf %lf %lf\n", &m(0, 0), &m(0, 1), &m(0, 2), &m(0, 3)) != 4 ||
      fscanf(f, "%lf %lf %lf %lf\n", &m(1, 0), &m(1, 1), &m(1, 2), &m(1, 3)) != 4 ||
      fscanf(f, "%lf %lf %lf %lf\n", &m(2, 0), &m(2, 1), &m(2, 2), &m(2, 3)) != 4) {
    fclose(f);
    FatalError("File does not appear to be a valid Aladin output file: " << fname);
    exit(1);
  }
  m(3, 0) = m(3, 1) = m(3, 2) = .0; m(3, 3) = 1.0;
  fclose(f);
  unique_ptr<AffineTransformation> dof(new AffineTransformation);
  dof->PutMatrix(m);
  return dof.release();
}

// -----------------------------------------------------------------------------
/// Read transformation from reg_f3d output file
Transformation *ReadF3D(const char *fname, const char *dofin_name = NULL,
                        int xyz_units = 0, int steps = 0,
                        F3DTransformationType type = F3D_TYPE_UNKNOWN)
{
  #if MIRTK_ImageIO_WITH_NIfTI
    if (NiftiImageReader::CheckHeader(fname)) {
      FatalError("Input file is not a F3D output NIfTI image file!");
    }

    // Read NIfTI header
    NiftiImageInfo hdr(fname);
    if (type == F3D_TYPE_UNKNOWN) {
      if (hdr.intent_code == NIFTI_INTENT_VECTOR &&
          hdr.intent_name == "NREG_TRANS") {
        type = static_cast<F3DTransformationType>(hdr.intent_p1);
      } else {
        FatalError("Cannot determine format of input F3D vector field from NIfTI intent code!");
      }
    }
    if (steps == 0) {
      steps = static_cast<int>(hdr.intent_p2);
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
    unique_ptr<MultiLevelFreeFormTransformation> mffd(new MultiLevelFreeFormTransformation());

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
        ffd->NumberOfSteps(pow(2.0, steps));
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
  #else // MIRTK_ImageIO_WITH_NIfTI
    Warning("Cannot read F3D output file without NIfTI module");
    return NULL;
  #endif // MIRTK_ImageIO_WITH_NIfTI
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
  unique_ptr<AffineTransformation> dof(new AffineTransformation);
  dof->PutMatrix(m);
  return dof.release();
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
// (M)IRTK
// =============================================================================

// -----------------------------------------------------------------------------
/// Convert between MIRTK transformation types
bool ConvertMIRTKTransformation(const Transformation  *dofin,
                                Transformation        *dofout,
                                const ImageAttributes *target = nullptr)
{
  // Just copy parameters whenever possible
  if (dofout->CopyFrom(dofin)) return true;

  // Input...
  const HomogeneousTransformation        *ilin  = nullptr; // ...linear transformation
  const FreeFormTransformation           *iffd  = nullptr; // or non-linear  FFD
  const MultiLevelFreeFormTransformation *imffd = nullptr; // or multi-level FFD

  ( ilin = dynamic_cast<const HomogeneousTransformation        *>(dofin)) ||
  ( iffd = dynamic_cast<const FreeFormTransformation           *>(dofin)) ||
  (imffd = dynamic_cast<const MultiLevelFreeFormTransformation *>(dofin));

  if (imffd && imffd->NumberOfLevels() == 0) {
    ilin  = imffd->GetGlobalTransformation();
    imffd = nullptr;
  }

  // Output...
  HomogeneousTransformation        *olin  = nullptr; // ...linear transformation
  FreeFormTransformation           *offd  = nullptr; // or non-linear FFD
  MultiLevelTransformation         *omffd = nullptr; // or multi-level FFD
  MultiLevelFreeFormTransformation *osum  = nullptr; // (i.e., additive MFFD)

  ( olin = dynamic_cast<HomogeneousTransformation *>(dofout)) ||
  ( offd = dynamic_cast<FreeFormTransformation    *>(dofout)) ||
  (omffd = dynamic_cast<MultiLevelTransformation  *>(dofout));

  if (omffd) {
    const int nactive = omffd->NumberOfActiveLevels();
    if (nactive == 0) {
      FatalError("ConvertMIRTKTransformation: Expected output MFFD to have at least one active level!");
    } else if (nactive == 1) {
      for (int l = omffd->NumberOfLevels(); l >= 0; --l) {
        if (!omffd->LocalTransformationIsActive(l)) continue;
        offd = omffd->GetLocalTransformation(l);
      }
    }
    osum = dynamic_cast<MultiLevelFreeFormTransformation *>(omffd);
  }

  // Copy global transformation
  if (ilin) {
    if (olin) {
      olin->CopyFrom(ilin);
      return true;
    } else if (omffd) {
      omffd->GetGlobalTransformation()->CopyFrom(ilin);
      return true;
    }
  // Copy local transformation
  } else if (iffd && offd) {
    if (offd->CopyFrom(iffd)) return true;
  // Copy global and local transformation (additive MFFD only!)
  } else if (imffd && imffd->NumberOfLevels() == 1 && osum && offd) {
    osum->GetGlobalTransformation()->CopyFrom(imffd->GetGlobalTransformation());
    if (offd->CopyFrom(imffd->GetLocalTransformation(0))) return true;
  }

  // Domain for approximation
  ImageAttributes domain;
  if (target) {
    domain = *target;
  }
  if (!domain) {
    if      (iffd) domain = iffd->Attributes();
    else if (offd) domain = offd->Attributes();
    else {
      Warning("Cannot convert transformation without input -target image.");
      return false;
    }
  }

  // Otherwise, approximate the input transformation
  double error = dofout->ApproximateAsNew(domain, dofin);
  if (verbose) cout << "RMS error of approximation = " << error << endl;
  return true;
}

// -----------------------------------------------------------------------------
/// Write (M)IRTK transformation file
bool WriteMIRTK(const char *fname, const Transformation *dof,
                                   ImageAttributes target = ImageAttributes(),
                                   double dx = .0, double dy = .0, double dz = .0,
                                   TransformationType type = TRANSFORMATION_UNKNOWN)
{
  const FreeFormTransformation   *ffd  = dynamic_cast<const FreeFormTransformation    *>(dof);
  const MultiLevelTransformation *mffd = dynamic_cast<const MultiLevelTransformation  *>(dof);
  if (mffd) ffd = mffd->GetLocalTransformation(-1);

  if (ffd) {
    if (dx <= .0) dx = ffd->GetXSpacing();
    if (dy <= .0) dy = ffd->GetYSpacing();
    if (dz <= .0) dz = ffd->GetZSpacing();
    if (!target) target = ffd->Attributes();
  }

  bool type_is_linear = (type == TRANSFORMATION_RIGID ||
                         type == TRANSFORMATION_SIMILARITY ||
                         type == TRANSFORMATION_AFFINE);

  bool resample = ffd && !type_is_linear && (!fequal(ffd->GetXSpacing(), dx) ||
                                             !fequal(ffd->GetYSpacing(), dy) ||
                                             !fequal(ffd->GetZSpacing(), dz));

  if (resample || (type != TRANSFORMATION_UNKNOWN && dof->TypeOfClass() != type)) {
    unique_ptr<Transformation> dofout(Transformation::New(type));
    FreeFormTransformation *ffdout = dynamic_cast<FreeFormTransformation *>(dofout.get());
    if (ffdout) {
      if (!target) {
        Warning("Cannot convert linear transformation to FFD without input -target image!");
        return false;
      }
      ffdout->Initialize(target, dx, dy, dz);
    }
    if (ConvertMIRTKTransformation(dof, dofout.get(), bool(target) ? &target : nullptr)) {
      dofout->Write(fname);
      return true;
    }
  } else {
    dof->Write(fname);
    return true;
  }

  return false;
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
Matrix ToFLIRTMatrix(const ImageAttributes &target,
                     const ImageAttributes &source,
                     const HomogeneousTransformation *dof)
{
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
  return WriteFLIRTMatrix(fname, ToFLIRTMatrix(target_attr, source_attr, lin));
}

// -----------------------------------------------------------------------------
/// Convert 3D+t displacement field to FSL warp
class ConvertToFSLWarp : public VoxelFunction
{
  const Matrix    *_World2Source;
  const BaseImage *_Warp;
  bool             _Relative;

  static const int _x = 0;
  const int        _y, _z;

public:

  ConvertToFSLWarp(const Matrix *w2s, const BaseImage &warp, bool relative)
  :
    _World2Source(w2s), _Warp(&warp), _Relative(relative),
    _y(warp.NumberOfSpatialVoxels()), _z(_y + _y)
  {}

  template <class TIn, class TOut>
  void operator ()(int i, int j, int k, int, const TIn *din, TOut *dout) const
  {
    double x, y, z;
    x = i, y = j, z = k;
    _Warp->ImageToWorld(x, y, z);
    x += static_cast<double>(din[_x]);
    y += static_cast<double>(din[_y]);
    z += static_cast<double>(din[_z]);
    Transform(*_World2Source, x, y, z);
    if (_Relative) x -= i, y -= j, z -= k;
    dout[_x] = static_cast<TOut>(x * _Warp->GetXSize());
    dout[_y] = static_cast<TOut>(y * _Warp->GetYSize());
    dout[_z] = static_cast<TOut>(z * _Warp->GetZSize());
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
GenericImage<TReal> ToFSLWarpField(const ImageAttributes &target,
                                   const ImageAttributes &source,
                                   const Transformation  *dof,
                                   bool relative = true)
{
  GenericImage<TReal> warp(target, 3);

  // Get world displacement field
  dof->Displacement(warp);

  // Convert to FSL warp field
  const Matrix i2s = source.GetWorldToImageMatrix();
  ParallelForEachVoxel(ConvertToFSLWarp(&i2s, warp, relative), warp.Attributes(), warp);

  // Set _dt != 0 such that warp field is saved as 3D+t time series instead of
  // 3D vector field, i.e., such that NIfTI dim[8] = {4, _x, _y, _z, 3, 1, 1, 1}
  warp.PutTSize(1.0);

  // x image axis must be mirrored if determinant of world to image orientation
  // matrix is positive, i.e., when the image is stored in neurological order
  if (warp.Attributes().GetWorldToLatticeOrientation().Det() > .0) {
    if (verbose) cout << "qform determinant is positive, reflecting x axis of FSL warp field" << endl;
    warp.ReflectX();
  }

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
    const MultiLevelTransformation *mffd = dynamic_cast<const MultiLevelTransformation *>(dof);
    const FreeFormTransformation   *ffd  = dynamic_cast<const FreeFormTransformation   *>(dof);
    if (mffd) ffd = mffd->GetLocalTransformation(-1);
    if (!ffd) {
      Warning("Cannot convert linear transformation to FSL warp field without input -target image!");
      return false;
    }
    target = ffd->Attributes();
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
  FILE *f = fopen(fname, "w");
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

  unique_ptr<MultiLevelFreeFormTransformation> new_mffd;
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
  TransformationFileFormat format_out = Format_MIRTK;

  const char *dofin_name  = nullptr;
  const char *target_name = nullptr;
  const char *source_name = nullptr;
  int         steps       = 0;
  int         xyz_units   = 0;
  double      dx = .0, dy = .0, dz = .0;
  double      ts = .0;

  #if MIRTK_ImageIO_WITH_NIfTI
    xyz_units = NIFTI_UNITS_MM;
  #endif 

  for (ALL_OPTIONS) {
    if (OPTION("-input-format")) {
      const char *arg = ARGUMENT;
      if (!FromString(arg, format_in)) {
        FatalError("Invalid -input-format file format argument: " << arg);
        exit(1);
      }
    }
    else if (OPTION("-format") || OPTION("-output-format")) {
      const char *arg = ARGUMENT;
      if (!FromString(arg, format_out)) {
        FatalError("Invalid [-output]-format file format argument: " << arg);
        exit(1);
      }
    }
    else if (OPTION("-dofin")) dofin_name = ARGUMENT;
    else if (OPTION("-target")) target_name = ARGUMENT;
    else if (OPTION("-source")) source_name = ARGUMENT;
    else if (OPTION("-steps")) steps = atoi(ARGUMENT);
    else if (OPTION("-Ts")) PARSE_ARGUMENT(ts);
    else if (OPTION("-ds")) {
      PARSE_ARGUMENT(dx);
      dy = dz = dx;
    }
    else if (OPTION("-dx")) PARSE_ARGUMENT(dx);
    else if (OPTION("-dy")) PARSE_ARGUMENT(dy);
    else if (OPTION("-dz")) PARSE_ARGUMENT(dz);
    else if (OPTION("-xyz_units")) {
      const char *arg = ARGUMENT;
      #if MIRTK_ImageIO_WITH_NIfTI
        if (!FromString(arg, xyz_units)) {
          FatalError("Invalid argument for option -xyz_units: " << arg);
          exit(1);
        }
      #endif // MIRTK_ImageIO_WITH_NIfTI
    }
    else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
  }

  if (format_out == Format_Unknown) {
    FatalError("Output file format cannot be unknown");
  }

  // Read target/source attributes
  ImageAttributes target_attr, source_attr;
  if (target_name) {
    InitializeImageIOLibrary();
    BinaryImage target(target_name);
    target_attr = target.Attributes();
  }
  if (source_name) {
    InitializeImageIOLibrary();
    BinaryImage source(source_name);
    source_attr = source.Attributes();
  }

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

  if (format_in == Format_WorldDisplacement || format_in == Format_VoxelDisplacement) {
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
  unique_ptr<Transformation> dof;
  switch (format_in) {

    // Image
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

    // (M)IRTK
    case Format_MIRTK: {
      dof.reset(ReadMIRTK(input_name));
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
      dof.reset(ReadF3D(input_name, dofin_name, xyz_units, steps, ToF3DTransformationType(format_in)));
    } break;

    // FreeSurfer
    case Format_MNI_XFM: {
      dof.reset(ReadXFM(input_name));
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
  const string ext_out = Extension(input_name, EXT_LastWithoutGz); 
  bool  ext_out_nifti  = (ext_out == ".nii" || ext_out == ".hdr" || ext_out == ".img");
  bool  is_linear      = (dynamic_cast<HomogeneousTransformation *>(dof.get()) != nullptr);

  if (format_out == Format_FSL) {
    if (ext_out_nifti || !is_linear) format_out = Format_FSL_WarpRelative;
    else                             format_out = Format_FSL_FLIRT;
  } else if (format_out == Format_NREG) {
    if (ext_out_nifti || !is_linear) format_out = Format_F3D;
    else                             format_out = Format_Aladin;
  }

  const char *dx_output = nullptr;
  const char *dy_output = nullptr;
  const char *dz_output = nullptr;

  if ((format_in  == Format_WorldDisplacement || format_in  == Format_VoxelDisplacement) &&
      (format_out == Format_WorldDisplacement || format_out == Format_VoxelDisplacement)) {
    if (NUM_POSARGS == 4) {
      FatalError("Ambiguous usage when both input and output formats are dense displacement\n"
                 " fields and 4 positional arguments given. Either read *and* write components\n"
                 " from/to separate image files (6 arguments) or from/to a single 3D+3 image.\n");
    }
  }
  if (format_out == Format_WorldDisplacement || format_out == Format_VoxelDisplacement) {
    if (NUM_POSARGS == 4) {
      if (dx_input || dy_input || dz_input) {
        FatalError("Ambiguous usage when both input and output formats are dense displacement\n"
                   " fields and 4 positional arguments given. Either read *and* write components\n"
                   " from/to separate image files (6 arguments) or from/to a single 3D+3 image.\n");
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

  // Write transformation in requested output format
  bool success;
  switch (format_out) {

    // Image
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

    // (M)IRTK
    case Format_MIRTK:
    case Format_MIRTK_Rigid:
    case Format_MIRTK_Similarity:
    case Format_MIRTK_Affine:
    case Format_MIRTK_LinearFFD:
    case Format_MIRTK_BSplineFFD: {
      success = WriteMIRTK(output_name, dof.get(), target_attr, dx, dy, dz,
                           ToMIRTKTransformationType(format_out));
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

    // Unknown?
    default: success = false;
  }
  if (!success) {
    FatalError("Cannot write transformation in file format: " << ToString(format_out));
  }

  return 0;
}

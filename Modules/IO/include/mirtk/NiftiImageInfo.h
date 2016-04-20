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

#ifndef MIRTK_NiftiImageInfo_H
#define MIRTK_NiftiImageInfo_H

#include "mirtk/String.h"
#include "mirtk/Matrix.h"


namespace mirtk {


/// NIfTI datatype codes
enum NiftiDataType
{
  NIFTI_TYPE_UINT8      =    2,
  NIFTI_TYPE_INT16      =    4,
  NIFTI_TYPE_INT32      =    8,
  NIFTI_TYPE_FLOAT32    =   16,
  NIFTI_TYPE_COMPLEX64  =   32,
  NIFTI_TYPE_FLOAT64    =   64,
  NIFTI_TYPE_RGB24      =  128,
  NIFTI_TYPE_INT8       =  256,
  NIFTI_TYPE_UINT16     =  512,
  NIFTI_TYPE_UINT32     =  768,
  NIFTI_TYPE_INT64      = 1024,
  NIFTI_TYPE_UINT64     = 1280,
  NIFTI_TYPE_FLOAT128   = 1536,
  NIFTI_TYPE_COMPLEX128 = 1792,
  NIFTI_TYPE_COMPLEX256 = 2048,
  NIFTI_TYPE_RGBA32     = 2304
};

/// NIfTI intent codes, to describe intended meaning of dataset contents
enum NiftiIntent
{
  NIFTI_INTENT_NONE         =    0,
  NIFTI_FIRST_STATCODE      =    2,
  NIFTI_INTENT_CORREL       =    2,
  NIFTI_INTENT_TTEST        =    3,
  NIFTI_INTENT_FTEST        =    4,
  NIFTI_INTENT_ZSCORE       =    5,
  NIFTI_INTENT_CHISQ        =    6,
  NIFTI_INTENT_BETA         =    7,
  NIFTI_INTENT_BINOM        =    8,
  NIFTI_INTENT_GAMMA        =    9,
  NIFTI_INTENT_POISSON      =   10,
  NIFTI_INTENT_NORMAL       =   11,
  NIFTI_INTENT_FTEST_NONC   =   12,
  NIFTI_INTENT_CHISQ_NONC   =   13,
  NIFTI_INTENT_LOGISTIC     =   14,
  NIFTI_INTENT_LAPLACE      =   15,
  NIFTI_INTENT_UNIFORM      =   16,
  NIFTI_INTENT_TTEST_NONC   =   17,
  NIFTI_INTENT_WEIBULL      =   18,
  NIFTI_INTENT_CHI          =   19,
  NIFTI_INTENT_INVGAUSS     =   20,
  NIFTI_INTENT_EXTVAL       =   21,
  NIFTI_INTENT_PVAL         =   22,
  NIFTI_INTENT_LOGPVAL      =   23,
  NIFTI_INTENT_LOG10PVAL    =   24,
  NIFTI_LAST_STATCODE       =   24,
  NIFTI_INTENT_ESTIMATE     = 1001,
  NIFTI_INTENT_LABEL        = 1002,
  NIFTI_INTENT_NEURONAME    = 1003,
  NIFTI_INTENT_GENMATRIX    = 1004,
  NIFTI_INTENT_SYMMATRIX    = 1005,
  NIFTI_INTENT_DISPVECT     = 1006,
  NIFTI_INTENT_VECTOR       = 1007,
  NIFTI_INTENT_POINTSET     = 1008,
  NIFTI_INTENT_TRIANGLE     = 1009,
  NIFTI_INTENT_QUATERNION   = 1010,
  NIFTI_INTENT_DIMLESS      = 1011,
  NIFTI_INTENT_TIME_SERIES  = 2001,
  NIFTI_INTENT_NODE_INDEX   = 2002,
  NIFTI_INTENT_RGB_VECTOR   = 2003,
  NIFTI_INTENT_RGBA_VECTOR  = 2004,
  NIFTI_INTENT_SHAPE        = 2005
};

/// Convert string to NIfTI intent code
template <> bool FromString(const char *str, NiftiIntent &value);

/// Convert NIfTI intent code to string
template <> string ToString(const NiftiIntent &value, int w, char c, bool left);

/// NIfTI xform codes to describe the "standard" coordinate system
enum NiftiXForm
{
  NIFTI_XFORM_UNKNOWN      = 0, ///< Arbitrary coordinates (Method 1)
  NIFTI_XFORM_SCANNER_ANAT = 1, ///< Scanner-based anatomical coordinates
  NIFTI_XFORM_ALIGNED_ANAT = 2, ///< Coordinates aligned to another file's,
                                ///< or to anatomical "truth"
  NIFTI_XFORM_TALAIRACH    = 3, ///< Coordinates aligned to Talairach-
                                ///< Tournoux Atlas; (0,0,0)=AC, etc.
  NIFTI_XFORM_MNI_152      = 4  ///< MNI 152 normalized coordinates
};

/// NIfTI units codes to describe the unit of measurement for
/// each dimension of the dataset
enum NiftiUnits
{
  NIFTI_UNITS_UNKNOWN =  0, ///< NIFTI code for unspecified units
  NIFTI_UNITS_METER   =  1, ///< NIFTI code for meters
  NIFTI_UNITS_MM      =  2, ///< NIFTI code for millimeters
  NIFTI_UNITS_MICRON  =  3, ///< NIFTI code for micrometers
  NIFTI_UNITS_SEC     =  8, ///< NIFTI code for seconds
  NIFTI_UNITS_MSEC    = 16, ///< NIFTI code for milliseconds
  NIFTI_UNITS_USEC    = 24, ///< NIFTI code for microseconds
  NIFTI_UNITS_HZ      = 32, ///< NIFTI code for Hertz
  NIFTI_UNITS_PPM     = 40, ///< NIFTI code for ppm
  NIFTI_UNITS_RADS    = 48  ///< NIFTI code for radians per second
};

/// Convert string to NIfTI units code
bool FromString(const char *, NiftiUnits &);

/**
 * NIfTI image header information.
 *
 * \note Only contains most commonly used entries of nifti_image struct
 *       and uses MIRTK types for non-builtin types.
 */
struct NiftiImageInfo
{
  int    ndim;           ///< last dimension greater than 1 (1..7)
  int    nx;             ///< dimensions of grid array
  int    ny;             ///< dimensions of grid array
  int    nz;             ///< dimensions of grid array
  int    nt;             ///< dimensions of grid array
  int    nu;             ///< dimensions of grid array
  int    nv;             ///< dimensions of grid array
  int    nw;             ///< dimensions of grid array
  size_t nvox;           ///< number of voxels = nx*ny*nz*...*nw
  int    nbyper;         ///< bytes per voxel, matches datatype
  int    datatype;       ///< type of data in voxels: DT_* code
  double dx;             ///< grid spacings
  double dy;             ///< grid spacings
  double dz;             ///< grid spacings
  double dt;             ///< grid spacings
  double du;             ///< grid spacings
  double dv;             ///< grid spacings
  double dw;             ///< grid spacings
  double scl_slope;      ///< scaling parameter - slope
  double scl_inter;      ///< scaling parameter - intercept
  double cal_min;        ///< calibration parameter, minimum
  double cal_max;        ///< calibration parameter, maximum
  int    slice_code;     ///< code for slice timing pattern
  int    slice_start;    ///< index for start of slices
  int    slice_end;      ///< index for end of slices
  double slice_duration; ///< time between individual slices
  int    qform_code;     ///< codes for (x,y,z) space meaning
  Matrix qto_xyz;        ///< qform: transform (i,j,k) to (x,y,z)
  Matrix qto_ijk;        ///< qform: transform (x,y,z) to (i,j,k)
  int    sform_code;     ///< codes for (x,y,z) space meaning
  Matrix sto_xyz;        ///< sform: transform (i,j,k) to (x,y,z)
  Matrix sto_ijk;        ///< sform: transform (x,y,z) to (i,j,k)
  double toffset;        ///< time coordinate offset
  int    xyz_units;      ///< dx,dy,dz units: NIFTI_UNITS_* code
  int    time_units;     ///< dt       units: NIFTI_UNITS_* code
  int    nifti_type;     ///< 0==ANALYZE, 1==NIFTI-1 (1 file),
                         ///              2==NIFTI-1 (2 files),
                         ///              3==NIFTI-ASCII (1 file)
  int    intent_code;    ///< statistic type (or something)
  double intent_p1;      ///< intent parameters
  double intent_p2;      ///< intent parameters
  double intent_p3;      ///< intent parameters
  string intent_name;    ///< optional description of intent data
  string descrip;        ///< optional text to describe dataset
  string aux_file;       ///< auxiliary filename
  string fname;          ///< header filename (.hdr or .nii)
  string iname;          ///< image filename  (.img or .nii)
  int    iname_offset;   ///< offset into iname where data starts
  int    swapsize;       ///< swap unit in image data (might be 0)
  int    byteorder;      ///< byte order on disk (MSB_ or LSB_FIRST)

  /// Constructor, reads header from NIfTI image file
  NiftiImageInfo(const char * = nullptr);
};


} // namespace mirtk

#endif // MIRTK_NiftiImageInfo_H

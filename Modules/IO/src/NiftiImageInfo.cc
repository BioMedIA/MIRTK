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

#include "mirtk/NiftiImageInfo.h"


// =============================================================================
// Enumerations -- **BEFORE** including nifti2_io.h which defines macros!!!
// =============================================================================

namespace mirtk {


// -----------------------------------------------------------------------------
template <> bool FromString(const char *str, NiftiIntent &value)
{
  string lstr = Trim(ToLower(str));
  if (lstr.rfind(" distribution") == lstr.length() - 13) {
    lstr = lstr.substr(0, lstr.length() - 13);
  } else if (lstr.rfind(" statistic") == lstr.length() - 10) {
    lstr = lstr.substr(0, lstr.length() - 10);
  }

  if      (lstr == "unknown" || lstr == "none") value = NIFTI_INTENT_NONE;
  else if (lstr == "correlation") value = NIFTI_INTENT_CORREL;
  else if (lstr == "t-statistic" || lstr == "t-test") value = NIFTI_INTENT_TTEST;
  else if (lstr == "f-statistic" || lstr == "f-test") value = NIFTI_INTENT_FTEST;
  else if (lstr == "z-score") value = NIFTI_INTENT_ZSCORE;
  else if (lstr == "f-statistic noncentral") value = NIFTI_INTENT_FTEST_NONC;
  else if (lstr == "chi-squared noncentral") value = NIFTI_INTENT_CHISQ_NONC;
  else if (lstr == "t-statistic noncentral") value = NIFTI_INTENT_TTEST_NONC;
  else if (lstr == "chi-squared") value = NIFTI_INTENT_CHISQ;
  else if (lstr == "beta") value = NIFTI_INTENT_BETA;
  else if (lstr == "gamma") value = NIFTI_INTENT_GAMMA;
  else if (lstr == "poisson") value = NIFTI_INTENT_POISSON;
  else if (lstr == "normal") value = NIFTI_INTENT_NORMAL;
  else if (lstr == "logistic") value = NIFTI_INTENT_LOGISTIC;
  else if (lstr == "laplace") value = NIFTI_INTENT_LAPLACE;
  else if (lstr == "uniform") value = NIFTI_INTENT_UNIFORM;
  else if (lstr == "weibull") value = NIFTI_INTENT_WEIBULL;
  else if (lstr == "chi") value = NIFTI_INTENT_CHI;
  else if (lstr == "inverse gaussian") value = NIFTI_INTENT_INVGAUSS;
  else if (lstr == "extreme value") value = NIFTI_INTENT_EXTVAL;
  else if (lstr == "p-value") value = NIFTI_INTENT_PVAL;
  else if (lstr == "log p-value") value = NIFTI_INTENT_LOGPVAL;
  else if (lstr == "log10 p-value") value = NIFTI_INTENT_LOG10PVAL;
  else if (lstr == "estimate") value = NIFTI_INTENT_ESTIMATE;
  else if (lstr == "label index" || lstr == "label") value = NIFTI_INTENT_LABEL;
  else if (lstr == "neuronames index" || lstr == "neuroname") value = NIFTI_INTENT_NEURONAME;
  else if (lstr == "general matrix") value = NIFTI_INTENT_GENMATRIX;
  else if (lstr == "symmetric matrix") value = NIFTI_INTENT_SYMMATRIX;
  else if (lstr == "displacement vector") value = NIFTI_INTENT_DISPVECT;
  else if (lstr == "vector") value = NIFTI_INTENT_VECTOR;
  else if (lstr == "pointset") value = NIFTI_INTENT_POINTSET;
  else if (lstr == "triangle") value = NIFTI_INTENT_TRIANGLE;
  else if (lstr == "quaternion") value = NIFTI_INTENT_QUATERNION;
  else if (lstr == "time series") value = NIFTI_INTENT_TIME_SERIES;
  else if (lstr == "node index") value = NIFTI_INTENT_NODE_INDEX;
  else if (lstr == "shape") value = NIFTI_INTENT_SHAPE;
  else if (lstr == "rgb") value = NIFTI_INTENT_RGB_VECTOR;
  else if (lstr == "rgba") value = NIFTI_INTENT_RGBA_VECTOR;
  else if (lstr == "dimensionless number" || lstr == "dimensionless") {
    value = NIFTI_INTENT_DIMLESS;
  } else {
    return false;
  }

  return true;
}

// -----------------------------------------------------------------------------
template <> string ToString(const NiftiIntent &value, int w, char c, bool left)
{
  const char *str = "Unknown";
  switch (value) {
    case NIFTI_INTENT_NONE:       str = "Unknown"; break;
    case NIFTI_INTENT_CORREL:     str = "Correlation statistic"; break;
    case NIFTI_INTENT_TTEST:      str = "T-statistic"; break;
    case NIFTI_INTENT_FTEST:      str = "F-statistic"; break;
    case NIFTI_INTENT_ZSCORE:     str = "Z-score"; break;
    case NIFTI_INTENT_CHISQ:      str = "Chi-squared distribution"; break;
    case NIFTI_INTENT_BETA:       str = "Beta distribution"; break;
    case NIFTI_INTENT_BINOM:      str = "Binomial distribution"; break;
    case NIFTI_INTENT_GAMMA:      str = "Gamma distribution"; break;
    case NIFTI_INTENT_POISSON:    str = "Poisson distribution"; break;
    case NIFTI_INTENT_NORMAL:     str = "Normal distribution"; break;
    case NIFTI_INTENT_FTEST_NONC: str = "F-statistic noncentral"; break;
    case NIFTI_INTENT_CHISQ_NONC: str = "Chi-squared noncentral"; break;
    case NIFTI_INTENT_LOGISTIC:   str = "Logistic distribution"; break;
    case NIFTI_INTENT_LAPLACE:    str = "Laplace distribution"; break;
    case NIFTI_INTENT_UNIFORM:    str = "Uniform distribution"; break;
    case NIFTI_INTENT_TTEST_NONC: str = "T-statistic noncentral"; break;
    case NIFTI_INTENT_WEIBULL:    str = "Weibull distribution"; break;
    case NIFTI_INTENT_CHI:        str = "Chi distribution"; break;
    case NIFTI_INTENT_INVGAUSS:   str = "Inverse Gaussian distribution"; break;
    case NIFTI_INTENT_EXTVAL:     str = "Extreme Value distribution"; break;
    case NIFTI_INTENT_PVAL:       str = "P-value"; break;

    case NIFTI_INTENT_LOGPVAL:    str = "Log P-value"; break;
    case NIFTI_INTENT_LOG10PVAL:  str = "Log10 P-value"; break;

    case NIFTI_INTENT_ESTIMATE:   str = "Estimate"; break;
    case NIFTI_INTENT_LABEL:      str = "Label index"; break;
    case NIFTI_INTENT_NEURONAME:  str = "NeuroNames index"; break;
    case NIFTI_INTENT_GENMATRIX:  str = "General matrix"; break;
    case NIFTI_INTENT_SYMMATRIX:  str = "Symmetric matrix"; break;
    case NIFTI_INTENT_DISPVECT:   str = "Displacement vector"; break;
    case NIFTI_INTENT_VECTOR:     str = "Vector"; break;
    case NIFTI_INTENT_POINTSET:   str = "Pointset"; break;
    case NIFTI_INTENT_TRIANGLE:   str = "Triangle"; break;
    case NIFTI_INTENT_QUATERNION: str = "Quaternion"; break;

    case NIFTI_INTENT_DIMLESS:    str = "Dimensionless number"; break;

    case NIFTI_INTENT_TIME_SERIES: str = "Time series"; break;
    case NIFTI_INTENT_NODE_INDEX:  str = "Node index"; break;
    case NIFTI_INTENT_SHAPE:       str = "Shape"; break;
    case NIFTI_INTENT_RGB_VECTOR:  str = "RGB"; break;
    case NIFTI_INTENT_RGBA_VECTOR: str = "RGBA"; break;
  }
  return ToString(str, w, c, left);
}

// -----------------------------------------------------------------------------
bool FromString(const char *str, NiftiUnits &units)
{
  const string lstr = Trim(ToLower(str));
  if (lstr == "m" || lstr == "meter" || lstr == ToString(int(NIFTI_UNITS_METER))) {
    units = NIFTI_UNITS_METER;
  } else if (lstr == "mm" || lstr == "millimeter" || lstr == ToString(int(NIFTI_UNITS_MM))) {
    units = NIFTI_UNITS_MM;
  } else if (lstr == "mu" || lstr == "micrometer" || lstr == "micron" || lstr == ToString(int(NIFTI_UNITS_MICRON))) {
    units = NIFTI_UNITS_MICRON;
  } else {
    return false;
  }
  return true;
}


} // namespace mirtk

// =============================================================================
// Include NIfTI I/O library header
// =============================================================================

#include "nifti/nifti2_io.h"


namespace mirtk {


// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
static inline Matrix ToMatrix(const nifti_dmat44 &mat)
{
  Matrix m;
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j) {
    m(i, j) = mat.m[i][j];
  }
  return m;
}

// -----------------------------------------------------------------------------
inline NiftiImageInfo ToNiftiImageInfo(const nifti_image *nim)
{
  NiftiImageInfo info;
  info.ndim           = static_cast<int>(nim->ndim);
  info.nx             = static_cast<int>(nim->nx);
  info.ny             = static_cast<int>(nim->ny);
  info.nz             = static_cast<int>(nim->nz);
  info.nt             = static_cast<int>(nim->nt);
  info.nu             = static_cast<int>(nim->nu);
  info.nv             = static_cast<int>(nim->nv);
  info.nw             = static_cast<int>(nim->nw);
  info.nvox           = static_cast<int>(nim->nvox);
  info.nbyper         = nim->nbyper;
  info.datatype       = nim->datatype;
  info.dx             = nim->dx;
  info.dy             = nim->dy;
  info.dz             = nim->dz;
  info.dt             = nim->dt;
  info.du             = nim->du;
  info.dv             = nim->dv;
  info.dw             = nim->dw;
  info.scl_slope      = nim->scl_slope;
  info.scl_inter      = nim->scl_inter;
  info.cal_min        = nim->cal_min;
  info.cal_max        = nim->cal_max;
  info.slice_code     = nim->slice_code;
  info.slice_start    = static_cast<int>(nim->slice_start);
  info.slice_end      = static_cast<int>(nim->slice_end);
  info.slice_duration = nim->slice_duration;
  info.qform_code     = nim->qform_code;
  info.qto_xyz        = ToMatrix(nim->qto_xyz);
  info.qto_ijk        = ToMatrix(nim->qto_ijk);
  info.sform_code     = nim->sform_code;
  info.sto_xyz        = ToMatrix(nim->sto_xyz);
  info.sto_ijk        = ToMatrix(nim->sto_ijk);
  info.toffset        = nim->toffset;
  info.xyz_units      = nim->xyz_units;
  info.time_units     = nim->time_units;
  info.nifti_type     = nim->nifti_type;
  info.intent_code    = nim->intent_code;
  info.intent_p1      = nim->intent_p1;
  info.intent_p2      = nim->intent_p2;
  info.intent_p3      = nim->intent_p3;
  info.intent_name    = nim->intent_name;
  info.descrip        = nim->descrip;
  info.aux_file       = nim->aux_file;
  info.fname          = nim->fname;
  info.iname          = nim->iname;
  info.iname_offset   = static_cast<int>(nim->iname_offset);
  info.swapsize       = nim->swapsize;
  info.byteorder      = nim->byteorder;
  return info;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
NiftiImageInfo::NiftiImageInfo(const char *fname)
{
  nifti_image *nim;
  if (fname) {
    const int read_data = 0;
    nim = nifti_image_read(fname, read_data);
  } else {
    nim = nifti_simple_init_nim();
  }
  *this = ToNiftiImageInfo(nim);
  nifti_image_free(nim);
}


} // namespace mirtk

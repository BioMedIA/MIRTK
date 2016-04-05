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

#include "nifti/nifti2_io.h"

#undef NIFTI_UNITS_METER
#undef NIFTI_UNITS_MM
#undef NIFTI_UNITS_MICRON


// =============================================================================
// Enumerations -- **BEFORE** including nifti1_io.h which defines macros!!!
// =============================================================================

namespace mirtk {


// -----------------------------------------------------------------------------
template <> string ToString(const NiftiIntent &value, int w, char c, bool left)
{
  return ToString(nifti_intent_string(value), w, c, left);
}

// -----------------------------------------------------------------------------
bool FromString(const char *str, NiftiUnits &units)
{
  const string lstr = ToLower(str);
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
  info.ndim           = nim->ndim;
  info.nx             = nim->nx;
  info.ny             = nim->ny;
  info.nz             = nim->nz;
  info.nt             = nim->nt;
  info.nu             = nim->nu;
  info.nv             = nim->nv;
  info.nw             = nim->nw;
  info.nvox           = nim->nvox;
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
  info.slice_start    = nim->slice_start;
  info.slice_end      = nim->slice_end;
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
  info.iname_offset   = nim->iname_offset;
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

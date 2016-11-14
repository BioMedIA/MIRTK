/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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

#include "NiftiImage.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Matrix.h"


namespace mirtk {


// -----------------------------------------------------------------------------
NiftiImage::NiftiImage(const char *fname)
:
  nim(nullptr)
{
  if (fname) Read(fname);
}

// -----------------------------------------------------------------------------
NiftiImage::~NiftiImage()
{
  if (nim != nullptr) {
    nim->data = nullptr; // (de)allocated elsewhere
    nifti_image_free(nim);
  }
  nim = nullptr;
}

// -----------------------------------------------------------------------------
void NiftiImage::Read(const char *filename)
{
  if (nim != nullptr) {
    nim->data = nullptr; // (de)allocated elsewhere
    nifti_image_free(nim);
  }
  const int read_data = 0;
  nim = nifti_image_read(filename, read_data);
}

// -----------------------------------------------------------------------------
void NiftiImage::Initialize(int x, int y, int z, int t, int u,
                            double xsize, double ysize, double zsize, double tsize,
                            int datatype,
                            Matrix &qmat, Matrix *smat, double torigin,
                            const void *data)
{
  nifti_dmat44 mat_44;
  int nbytepix = 0, ss = 0;
  int i, j;

  // Initialize NIfTI I/O data structure
  if (nim != nullptr) {
    nim->data = nullptr; // (de)allocated elsewhere
    nifti_image_free(nim);
  }
  nim = nifti_simple_init_nim();

  // Spatial dimensions
  nim->nifti_type = 1;               // 1==nii (1 file) - should set the magic in nhdr in Write()
  nim->datatype   = datatype;        // Will be NIFTI_TYPE_UINT8 | NIFTI_TYPE_INT16 | NIFTI_TYPE_FLOAT32
  nim->ndim       = ((u > 1) ? 5 : (t > 1 ? 4 : 3));
  nim->nx         = (x > 1 ? x : 1); // So that nvox can be computed correctly, see below
  nim->ny         = (y > 1 ? y : 1); // dito
  nim->nz         = (z > 1 ? z : 1); // ...
  nim->nt         = (t > 1 ? t : 1); // ...
  nim->nu         = (u > 1 ? u : 1); // ...
  nim->nv         = 1;               // ...
  nim->nw         = 1;               // ...
  nim->dx         = abs(xsize);      // Store only absolute pixel size values
  nim->dy         = abs(ysize);      // dito
  nim->dz         = abs(zsize);      // ...
  nim->dt         = abs(tsize);
  nim->toffset    = torigin;
  // Redundant fields ndim->dim/pixdim will be filled by nim2nhdr/nhdr2nim conversions in Write()

  // Derived values
  nifti_datatype_sizes(datatype, &nbytepix, &ss);
  nim->nbyper = nbytepix;
  nim->nvox   = nim->nx * nim->ny * nim->nz * nim->nt * nim->nu * nim->nv * nim->nw;

  // Compute qform
  for (i = 0; i < 4; ++i)
  for (j = 0; j < 4; ++j) {
    mat_44.m[i][j] = qmat(i, j);
  }

  nim->qform_code = 1;
  nim->qto_xyz = mat_44;
  nifti_dmat44_to_quatern(
    mat_44,
    &(nim->quatern_b), &(nim->quatern_c), &(nim->quatern_d),
    &(nim->qoffset_x), &(nim->qoffset_y), &(nim->qoffset_z),
    &(nim->dx)       , &(nim->dy)       , &(nim->dz)       ,
    &(nim->qfac));

  // Set sform
  if (smat) {
    for (i = 0; i < 4; ++i)
    for (j = 0; j < 4; ++j) {
      mat_44.m[i][j] = (*smat)(i, j);
    }
    nim->sform_code = NIFTI_XFORM_ALIGNED_ANAT;
    nim->sto_xyz    = mat_44;
    nim->sto_ijk    = nifti_dmat44_inverse(mat_44);
  } else {
    nim->sform_code = NIFTI_XFORM_UNKNOWN;
  }

  // Units
  nim->xyz_units  = NIFTI_UNITS_MM;
  nim->time_units = NIFTI_UNITS_MSEC;

  // Intensity rescaling
  nim->scl_slope = 1;
  nim->scl_inter = 0;

  // Intent
  nim->intent_code = ((u > 1) ? NIFTI_INTENT_VECTOR : NIFTI_INTENT_NONE);
  nim->intent_p1 = nim->intent_p2 = nim->intent_p3 = 0;

  // Data will be written out by ImageToFile class
  nim->data = const_cast<void *>(data);
}

// -----------------------------------------------------------------------------
void NiftiImage::Initialize(int x, int y, int z, int t,
                            double xsize, double ysize, double zsize, double tsize,
                            int datatype,
                            Matrix &qmat, Matrix *smat, double torigin, const void *data)
{
  int u = 1;
  if (t > 1 && !tsize) swap(t, u);
  Initialize(x, y, z, t, u, xsize, ysize, zsize, tsize, datatype, qmat, smat, torigin, data);
}

// -----------------------------------------------------------------------------
void NiftiImage::Print()
{
  cerr << "NiftiImage::Print() : \n";
  if (nim != nullptr) nifti_image_infodump(nim);
}


} // namespace mirtk

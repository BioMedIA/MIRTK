/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
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

#define _NIFTI2_IO_C_
#include "NiftiImage.h"

#include "mirtk/NiftiImageReader.h"
#include "mirtk/ImageReaderFactory.h"

#include "mirtk/Array.h"
#include "mirtk/Voxel.h"
#include "mirtk/Path.h"
#include "mirtk/Matrix.h"

#include <sys/types.h>
#include <sys/stat.h>


namespace mirtk {


// Register image reader with object factory during static initialization
mirtkAutoRegisterImageReaderMacro(NiftiImageReader);


// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
static inline nifti_dmat33 nifti_dmat44_to_dmat33(const nifti_dmat44 &x)
{
  nifti_dmat33 y;
  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j) {
    y.m[i][j] = x.m[i][j];
  }
  return y;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
NiftiImageReader::NiftiImageReader()
:
  _Nifti(new NiftiImage())
{
}

// -----------------------------------------------------------------------------
NiftiImageReader::~NiftiImageReader()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
bool NiftiImageReader::CheckHeader(const char *fname)
{
  // Get name of corresponding image header file
  char *hname = nifti_findhdrname(fname);
  if (hname == nullptr) return false;
  // Open header file stream
  Cifstream from(hname);
  // Check if it as a NIfTI image in ASCII format
  // (cf. has_ascii_header function defined in nifti2_io.cc)
  char nia_magic[16];
  if (from.ReadAsChar(nia_magic, 12)) {
    if (strncmp(nia_magic, "<nifti_image", 12) == 0) {
      free(hname);
      return true;
    }
    from.Seek(0);
  }
  // Try to read NIfTI-1 header
  nifti_1_header n1hdr;
  const size_t h1size = sizeof(nifti_1_header);
  from.ReadAsChar((char *)&n1hdr, h1size);
  // Check if it is a NIfTI header and which version
  const int ni_ver = nifti_header_version((char *)&n1hdr, h1size);
  if (ni_ver == -1) return false;
  // In case of NIfTI-2, the file may be a CIFTI file instead
  if (ni_ver == 2) {
    bool maybe_cifti = false, is_cifti = false;
    // Fill remaining bytes of NIfTI-2 header
    nifti_2_header n2hdr;
    const size_t h2size = sizeof(nifti_2_header);
    memcpy(&n2hdr, &n1hdr, h1size);
    if (!from.ReadAsChar((char *)&n2hdr + h1size, static_cast<int>(h2size - h1size))) {
      free(hname);
      return false;
    }
    // Check if header has extensions
    bool has_extensions = false;
    nifti1_extender extdr;
    if (from.ReadAsChar((char *)&extdr, 4)) {
      has_extensions = (extdr.extension[0] == 1 &&
                        extdr.extension[1] == 0 &&
                        extdr.extension[2] == 0 &&
                        extdr.extension[3] == 0);
    }
    // Check dimension sizes if it could be a CIFTI before checking the extensions
    if (has_extensions) {
      nifti_image *nim = nifti_convert_n2hdr2nim(n2hdr, hname);
      maybe_cifti = (nim->nx < 2 && nim->ny < 2 && nim->nz < 2 && nim->nt < 2);
      if (nim->nu < 2 && nim->nt < 2) maybe_cifti = false;
      nifti_image_free(nim);
    }
    // Check presence of CIFTI extension
    if (maybe_cifti) {
      int size, code;
      const bool needs_swap = NIFTI_NEEDS_SWAP(n2hdr);
      while (from.ReadAsChar((char *)&size, 4) &&
             from.ReadAsChar((char *)&code, 4)) {
        if (needs_swap) {
          nifti_swap_4bytes(1, &size);
          nifti_swap_4bytes(1, &code);
        }
        if (size < 16 || size & 0xf) break;     // size must be multiple of 16
        if (!nifti_is_valid_ecode(code)) break; // extension code must be valid
        if (code == NIFTI_ECODE_CIFTI) {        // found CIFTI extension
          is_cifti = true;
          break;
        }
        from.Seek(from.Tell() + size - 8);
      }
    }
    // CIFTI file is no image file
    if (is_cifti) {
      free(hname);
      return false;
    }
  }
  free(hname);
  return true;
}

// -----------------------------------------------------------------------------
bool NiftiImageReader::CanRead(const char *fname) const
{
  return CheckHeader(fname);
}

// -----------------------------------------------------------------------------
void NiftiImageReader::Initialize()
{
  // Read header
  this->ReadHeader();

  // Open image data file
  char *iname = nifti_findimgname(_FileName.c_str(), _Nifti->nim->nifti_type);
  if (iname == nullptr) {
    cerr << this->NameOfClass() << "::ReadHeader: Could not find image data file for " << _FileName << endl;
    exit(1);
  }
  _ImageName = iname;
  free(iname);
  this->Open(_ImageName.c_str());
}

// -----------------------------------------------------------------------------
void NiftiImageReader::ReadHeader()
{
  int          order;
  double       det;
  nifti_dmat44 mat_44, mat_inv_44;
  nifti_dmat33 mat_33;

  // Read header
  char *hname = nifti_findhdrname(_FileName.c_str());
  if (hname == nullptr) {
    cerr << this->NameOfClass() << "::ReadHeader: Could not find header file for " << _FileName << endl;
    exit(1);
  }
  _Nifti->Read(hname);
  free(hname);

  // Check dimension
  if (_Nifti->nim->dim[0] > 5) {
    cerr << this->NameOfClass() << "::ReadHeader: Number of dimensions > 5 (Number of dimensions = " << _Nifti->nim->dim[0] << ") \n";
    exit(1);
  }
  if ((_Nifti->nim->dim[0] == 5) && (_Nifti->nim->dim[4] != 1)) {
    cerr << this->NameOfClass() << "::ReadHeader: Number of dimensions == 5 (Number of dimensions = " << _Nifti->nim->dim[0] << ") \n";
    exit(1);
  }

  // Check for swapping
  _Swapped = (_Nifti->nim->byteorder != nifti_short_order());

  // Check data scaling
  if (_Nifti->nim->scl_slope != 0) {
    _Slope = _Nifti->nim->scl_slope;
  } else {
    _Slope = 1.0;
  }
  _Intercept = _Nifti->nim->scl_inter;

  // Copy header information
  // (from now on force all vox dims to be positive - LR info is in sform)
  _Attributes._x  = static_cast<int>(_Nifti->nim->dim[1]);
  _Attributes._y  = static_cast<int>(_Nifti->nim->dim[2]);
  _Attributes._z  = static_cast<int>(_Nifti->nim->dim[3]);
  _Attributes._dx = abs(_Nifti->nim->pixdim[1]);
  _Attributes._dy = abs(_Nifti->nim->pixdim[2]);
  _Attributes._dz = abs(_Nifti->nim->pixdim[3]);
  if (_Nifti->nim->dim[0] == 4) {
    _Attributes._t  = static_cast<int>(_Nifti->nim->dim[4]);
    _Attributes._dt = abs(_Nifti->nim->pixdim[4]);
  } else if (_Nifti->nim->dim[0] == 5) {
    _Attributes._t  = static_cast<int>(_Nifti->nim->dim[5]);
    _Attributes._dt = (_Nifti->nim->intent_code == NIFTI_INTENT_VECTOR ? 0 : 1);
  } else {
    _Attributes._t  = 1;
    _Attributes._dt = 1;
  }
  if (_Attributes._x < 1) _Attributes._x = 1;
  if (_Attributes._y < 1) _Attributes._y = 1;
  if (_Attributes._z < 1) _Attributes._z = 1;
  if (_Attributes._t < 1) _Attributes._t = 1;

  // Check which coordinate system to use
  if (_Nifti->nim->qform_code > 0) {
    // Access qform
    mat_44 = _Nifti->nim->qto_xyz;
    mat_inv_44 = _Nifti->nim->qto_ijk;
  } else if (_Nifti->nim->sform_code > 0) {
    mat_44     = _Nifti->nim->sto_xyz;
    mat_inv_44 = _Nifti->nim->sto_ijk;
  } else {

    // What's below should be in _Nifti->nim->qto_xyz

    // grid spacings along diagonal (default radiological convention)
    mat_44.m[0][0] = -abs(_Nifti->nim->pixdim[1]);
    mat_44.m[1][1] =  abs(_Nifti->nim->pixdim[2]);
    mat_44.m[2][2] =  abs(_Nifti->nim->pixdim[3]);

    // off diagonal is zero
    mat_44.m[0][1] = mat_44.m[0][2] = .0;
    mat_44.m[1][0] = mat_44.m[1][2] = .0;
    mat_44.m[2][0] = mat_44.m[2][1] = .0;

    // Apart from (inverted) offset in last column for origin to become (0,0,0) below
    mat_44.m[0][3] =   (_Attributes._dx * (_Nifti->nim->dim[1] - 1)) / 2.0;
    mat_44.m[1][3] = - (_Attributes._dy * (_Nifti->nim->dim[2] - 1)) / 2.0;
    mat_44.m[2][3] = - (_Attributes._dz * (_Nifti->nim->dim[3] - 1)) / 2.0;

    // last row is always [ 0 0 0 1 ]
    mat_44.m[3][0] = mat_44.m[3][1] = mat_44.m[3][2] = .0, mat_44.m[3][3 ]= 1.0;

    // Invert
    mat_inv_44 = nifti_dmat44_inverse(mat_44);
  }

  // Get left/right order (no check for inconsistency between qform/sform since only one is used)
  mat_33 = nifti_dmat44_to_dmat33(mat_44);
  det    = nifti_dmat33_determ(mat_33);
  if (det < .0) {
    order = NiftiImage::RADIOLOGICAL;
  } else {
    order = NiftiImage::NEUROLOGICAL;
  }

  // Set axis orientation, including zaxis.
  // Need to preserve sign of axis, hence use absolute pixel sizes for descaling.
  for (int i = 0; i < 3; i++) {
    _Attributes._xaxis[i] = mat_44.m[i][0] / _Attributes._dx;
    _Attributes._yaxis[i] = mat_44.m[i][1] / _Attributes._dy;
    _Attributes._zaxis[i] = mat_44.m[i][2] / _Attributes._dz;
  }

  // Convert between nifti and  coordinate systems
  // See https://www.fmrib.ox.ac.uk/ibim/uploads/coordtransforms.pdf
  Matrix D(4, 4), D_inv(4, 4), M(4, 4), R;
  for (int j = 0; j < 4; j++) {
    M(j, j) = 1.0;
    for (int i = 0; i < 4; i++) {
      D(i, j)     = mat_44    .m[i][j];
      D_inv(i, j) = mat_inv_44.m[i][j];
    }
  }
  for (int i = 0; i < 3; i++) {
    M(i, 3) = static_cast<double>(_Nifti->nim->dim[i+1] - 1) / 2.0;
  }
  R = D * M * D_inv;

  // Set image origin by adding q/sform offset to third column of R:
  _Attributes._xorigin = R(0, 3) + mat_44.m[0][3];
  _Attributes._yorigin = R(1, 3) + mat_44.m[1][3];
  _Attributes._zorigin = R(2, 3) + mat_44.m[2][3];

  // Set temporal origin
  _Attributes._torigin = _Nifti->nim->toffset;

  // Convert spatial units to mm
  double xyz_scale = 1.0;
  if      (_Nifti->nim->xyz_units == NIFTI_UNITS_METER)  xyz_scale = 1.0e3;
  else if (_Nifti->nim->xyz_units == NIFTI_UNITS_MICRON) xyz_scale = 1.0e-3;
  if (xyz_scale != 1.0) {
    _Attributes._dx      *= xyz_scale;
    _Attributes._dy      *= xyz_scale;
    _Attributes._dz      *= xyz_scale;
    _Attributes._xorigin *= xyz_scale;
    _Attributes._yorigin *= xyz_scale;
    _Attributes._zorigin *= xyz_scale;
  }

  // Convert temporal units to ms
  double t_scale = 1.0;
  if      (_Nifti->nim->time_units == NIFTI_UNITS_SEC)  t_scale = 1.0e3;
  else if (_Nifti->nim->time_units == NIFTI_UNITS_USEC) t_scale = 1.0e-3;
  if (t_scale != 1.0) {
    _Attributes._dt      *= t_scale;
    _Attributes._torigin *= t_scale;
  }

  // If both qform and sform are given and sform_code is NIFTI_XFORM_ALIGNED_ANAT
  // use the qform for the imaging geometry (i.e., attributes) and the sform as
  // actual image to world coordinate transformation with additional affine transformation
  if (_Nifti->nim->qform_code > 0 && _Nifti->nim->sform_code == NIFTI_XFORM_ALIGNED_ANAT) {
    for (int j = 0; j < 4; ++j) {
      for (int i = 0; i < 4; ++i) {
        _Attributes._smat(i, j) = _Nifti->nim->sto_xyz.m[i][j]; // sform = T * qform
        D_inv            (i, j) = _Nifti->nim->qto_ijk.m[i][j]; // inverse of qform
      }
    }
    _Attributes._smat *= D_inv; // isolate affine transformation T
  } else {
    _Attributes._smat.Ident();
  }

  // Set data type and number of bytes per voxels
  switch (_Nifti->nim->datatype) {
    case NIFTI_TYPE_INT8:
       _DataType = MIRTK_VOXEL_CHAR;
       _Bytes    = 1;
       break;
    case NIFTI_TYPE_UINT8:
       _DataType = MIRTK_VOXEL_UNSIGNED_CHAR;
       _Bytes    = 1;
       break;
    case NIFTI_TYPE_INT16:
      _DataType = MIRTK_VOXEL_SHORT;
      _Bytes    = 2;
      break;
    case NIFTI_TYPE_UINT16:
      _DataType = MIRTK_VOXEL_UNSIGNED_SHORT;
      _Bytes    = 2;
      break;
    case NIFTI_TYPE_INT32:
      _DataType = MIRTK_VOXEL_INT;
      _Bytes    = 4;
      break;
    case NIFTI_TYPE_UINT32:
      _DataType = MIRTK_VOXEL_UNSIGNED_INT;
      _Bytes    = 4;
      break;
    case NIFTI_TYPE_FLOAT32:
      _DataType = MIRTK_VOXEL_FLOAT;
      _Bytes    = 4;
      break;
    case NIFTI_TYPE_FLOAT64:
      _DataType = MIRTK_VOXEL_DOUBLE;
      _Bytes    = 8;
      break;
    default:
      cerr << this->NameOfClass() << ": Data type " << _Nifti->nim->datatype << " not supported" << endl;
      exit(1);
  }

  // Data starts here 
  _Start = static_cast<int>(_Nifti->nim->iname_offset);
}

// -----------------------------------------------------------------------------
void NiftiImageReader::Print() const
{
  Matrix mat(3, 3);
  for (int i = 0; i < 3; i++) {
    mat(0, i) = _Attributes._xaxis[i] * _Attributes._dx;
    mat(1, i) = _Attributes._yaxis[i] * _Attributes._dy;
    mat(2, i) = _Attributes._zaxis[i] * _Attributes._dz;
  }

  cout << "Name of class is " << this->NameOfClass() << "\n";

  // Check which coordinate system is used.
  cout << "Transformation specified in NIFTI header using ";
  if (_Nifti->nim->qform_code > 0) {
    cout << "quaternion." << endl;
  } else if (_Nifti->nim->sform_code > 0) {
    cout << "affine transformation.\n";
  } else {
    cout << "default transformation.\n";
  }

  if (mat.Det() < .0) {
    cout << "Order is radiological\n";
  } else {
    cout << "Order is neurological\n";
  }

  cout << "Scale slope = " << _Slope << "\n";
  cout << "Intercept   = " << _Intercept << "\n";

  ImageReader::Print();
}


} // namespace mirtk

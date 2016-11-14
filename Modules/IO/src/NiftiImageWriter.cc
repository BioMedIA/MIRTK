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

#include "mirtk/NiftiImageWriter.h"
#include "mirtk/ImageWriterFactory.h"

#include "mirtk/Voxel.h"


namespace mirtk {


// Register image reader with object factory during static initialization
mirtkAutoRegisterImageWriterMacro(NiftiImageWriter);


// -----------------------------------------------------------------------------
Array<string> NiftiImageWriter::Extensions()
{
  Array<string> exts(4);
  exts[0] = ".nii";
  exts[1] = ".nii.gz";
  exts[2] = ".hdr";
  exts[3] = ".img";
  return exts;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
// Change memory layout from xyz...xyz...xyz... to xxx...yyy...zzz......
template <class VectorType>
inline void *ReformatImageData(const void *input, int nx, int ny, int nz, int nt, int nu)
{
  typedef typename voxel_info<VectorType>::ScalarType ScalarType;
  void *output = malloc(nx * ny * nz * nt * nu * sizeof(ScalarType));
  const ScalarType *in  = reinterpret_cast<const ScalarType *>(input);
  ScalarType       *out = reinterpret_cast<      ScalarType *>(output);
  const int offset = nx * ny * nz * nt;
  for (int l = 0; l < nt; ++l)
  for (int k = 0; k < nz; ++k)
  for (int j = 0; j < ny; ++j)
  for (int i = 0; i < nx; ++i) {
    ScalarType *tmp = out;
    for (int m = 0; m < nu; ++m) {
      *tmp = (*in++);
      tmp += offset;
    }
    ++out;
  }
  return output;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
NiftiImageWriter::NiftiImageWriter()
:
  _Nifti(new NiftiImage())
{
}

// -----------------------------------------------------------------------------
NiftiImageWriter::~NiftiImageWriter()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void NiftiImageWriter::Initialize()
{
  // Get image info
  const int nx = _Input->X();
  const int ny = _Input->Y();
  const int nz = _Input->Z();
  const int nt = _Input->T();

  double xsize, ysize, zsize, tsize;
  _Input->GetPixelSize(&xsize, &ysize, &zsize, &tsize);

  Matrix qmat = _Input->GetImageToWorldMatrix();
  UniquePtr<Matrix> smat;
  if (!_Input->GetAffineMatrix().IsIdentity()) {
    smat.reset(new Matrix(qmat));
    qmat = _Input->GetAffineMatrix().Inverse() * qmat;
  }
  const double torigin = _Input->ImageToTime(0);

  // Init header
  const void * const data = _Input->GetDataPointer();
  switch (_Input->GetDataType()) {
    case MIRTK_VOXEL_CHAR: {
      _Nifti->Initialize(nx, ny, nz, nt, xsize, ysize, zsize, tsize,
                         NIFTI_TYPE_INT8, qmat, smat.get(), torigin, data);
      break;
    }
    case MIRTK_VOXEL_UNSIGNED_CHAR: {
      _Nifti->Initialize(nx, ny, nz, nt, xsize, ysize, zsize, tsize,
                         NIFTI_TYPE_UINT8, qmat, smat.get(), torigin, data);
      break;
    }
    case MIRTK_VOXEL_SHORT: {
      _Nifti->Initialize(nx, ny, nz, nt, xsize, ysize, zsize, tsize,
                         NIFTI_TYPE_INT16, qmat, smat.get(), torigin, data);
      break;
    }
    case MIRTK_VOXEL_UNSIGNED_SHORT: {
      _Nifti->Initialize(nx, ny, nz, nt, xsize, ysize, zsize, tsize,
                         NIFTI_TYPE_UINT16, qmat, smat.get(), torigin, data);
      break;
    }
    case MIRTK_VOXEL_INT: {
      _Nifti->Initialize(nx, ny, nz, nt, xsize, ysize, zsize, tsize,
                         NIFTI_TYPE_INT32, qmat, smat.get(), torigin, data);
      break;
    }
    case MIRTK_VOXEL_UNSIGNED_INT: {
      _Nifti->Initialize(nx, ny, nz, nt, xsize, ysize, zsize, tsize,
                         NIFTI_TYPE_UINT32, qmat, smat.get(), torigin, data);
      break;
    }
    case MIRTK_VOXEL_FLOAT: {
      _Nifti->Initialize(nx, ny, nz, nt, xsize, ysize, zsize, tsize,
                         NIFTI_TYPE_FLOAT32, qmat, smat.get(), torigin, data);
      break;
    }
//    case MIRTK_VOXEL_FLOAT2: {
//      _Nifti->Initialize(nx, ny, nz, nt, 2, xsize, ysize, zsize, tsize,
//                         NIFTI_TYPE_FLOAT32, qmat, smat.get(), torigin,
//                         ReformatImageData<Float2>(data, nx, ny, nz, nt, 2));
//      break;
//    }
    case MIRTK_VOXEL_FLOAT3: {
      _Nifti->Initialize(nx, ny, nz, nt, 3, xsize, ysize, zsize, tsize,
                         NIFTI_TYPE_FLOAT32, qmat, smat.get(), torigin,
                         ReformatImageData<Float3>(data, nx, ny, nz, nt, 3));
      break;
    }
    case MIRTK_VOXEL_FLOAT4: {
      _Nifti->Initialize(nx, ny, nz, nt, 4, xsize, ysize, zsize, tsize,
                         NIFTI_TYPE_FLOAT32, qmat, smat.get(), torigin,
                         ReformatImageData<Float4>(data, nx, ny, nz, nt, 4));
      break;
    }
    case MIRTK_VOXEL_DOUBLE: {
      _Nifti->Initialize(nx, ny, nz, nt, xsize, ysize, zsize, tsize,
                         NIFTI_TYPE_FLOAT64, qmat, smat.get(), torigin, data);
      break;
    }
//    case MIRTK_VOXEL_DOUBLE2: {
//      _Nifti->Initialize(nx, ny, nz, nt, 2, xsize, ysize, zsize, tsize,
//                         NIFTI_TYPE_FLOAT64, qmat, smat.get(), torigin,
//                         ReformatImageData<Double2>(data, nx, ny, nz, nt, 2));
//      break;
//    }
    case MIRTK_VOXEL_DOUBLE3: {
      _Nifti->Initialize(nx, ny, nz, nt, 3, xsize, ysize, zsize, tsize,
                         NIFTI_TYPE_FLOAT64, qmat, smat.get(), torigin,
                         ReformatImageData<Double3>(data, nx, ny, nz, nt, 3));
      break;
    }
    case MIRTK_VOXEL_DOUBLE4: {
      _Nifti->Initialize(nx, ny, nz, nt, 4, xsize, ysize, zsize, tsize,
                         NIFTI_TYPE_FLOAT64, qmat, smat.get(), torigin,
                         ReformatImageData<Double4>(data, nx, ny, nz, nt, 4));
      break;
    }
    default:
      cerr << this->NameOfClass() << "::Initialize(): Unsupported voxel type: " << _Input->GetDataType() << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
void NiftiImageWriter::Run()
{
  // Initialize filter
  this->Initialize();

  // Set filename in _hdr
  struct nifti_2_header nhdr;
  nifti_convert_nim2n2hdr(_Nifti->nim, &nhdr);

  void * const data = _Nifti->nim->data; // Keep copy of data pointer
  _Nifti->nim->data = nullptr;           // Free nifti_image, but not data
  nifti_image_free(_Nifti->nim);

  _Nifti->nim               = nifti_convert_n2hdr2nim(nhdr, _FileName.c_str()); // This sets fname and iname
  _Nifti->nim->iname_offset = 352;       // Some nifti versions lose this on the way!
  _Nifti->nim->data         = data;      // Restore data pointer

  // Write hdr and data
  nifti_image_write(_Nifti->nim);

  // Finalize filter
  this->Finalize();
}

// -----------------------------------------------------------------------------
void NiftiImageWriter::Finalize()
{
  if (_Nifti->nim->data != _Input->GetDataPointer()) {
    free(_Nifti->nim->data);
  }
  _Nifti->nim->data = nullptr;
}


} // namespace mirtk

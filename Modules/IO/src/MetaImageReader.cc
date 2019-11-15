/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2019 Imperial College London
 * Copyright 2019 Andreas Schuh
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

#include "metaImage.h"

#include "mirtk/MetaImageReader.h"
#include "mirtk/ImageReaderFactory.h"


namespace mirtk {


// Register image reader with object factory during static initialization
mirtkAutoRegisterImageReaderMacro(MetaImageReader);


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
MetaImageReader::MetaImageReader()
:
  _MetaImage(new MetaImage())
{
}

// -----------------------------------------------------------------------------
MetaImageReader::~MetaImageReader()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
bool MetaImageReader::CheckHeader(const char *fname)
{
  return MetaImage().CanRead(fname);
}

// -----------------------------------------------------------------------------
bool MetaImageReader::CanRead(const char *fname) const
{
  return CheckHeader(fname);
}

// -----------------------------------------------------------------------------
void MetaImageReader::Initialize()
{
  // Close (unused) file stream
  this->Close();

  // Read image header
  this->ReadHeader();
}

// -----------------------------------------------------------------------------
void MetaImageReader::CopyHeader(const MetaImage &meta_image)
{
  const int ndims = meta_image.NDims();
  const int channels = meta_image.ElementNumberOfChannels();

  if (ndims < 2 || ndims > 4) {
    cerr << this->NameOfClass() << ": Only 2D, 3D, and 4D images supported, got ndims=" << ndims << endl;
    exit(1);
  }
  if (channels > 1) {
    cerr << this->NameOfClass() << ": Multi-channel image are currently not supported" << endl;
    exit(1);
  }

  const int * const size = meta_image.DimSize();
  _Attributes._x = size[0];
  _Attributes._y = size[1];
  _Attributes._z = (ndims > 2 ? size[2] : 1);
  _Attributes._t = (ndims > 3 ? size[3] : channels);

  const double * const spacing = meta_image.ElementSpacing();
  _Attributes._dx = spacing[0];
  _Attributes._dy = spacing[1];
  _Attributes._dz = (ndims > 2 ? spacing[2] : 1.);
  _Attributes._dt = (ndims > 3 ? spacing[3] : (channels > 1 ? 0. : 1.));

  const double * const offset = meta_image.Offset();
  _Attributes._xorigin = offset[0];
  _Attributes._yorigin = offset[1];
  _Attributes._zorigin = (ndims > 2 ? offset[2] : 0.);
  _Attributes._torigin = (ndims > 3 ? offset[3] : 0.);

  _Attributes._xaxis[0] = 1.;
  _Attributes._xaxis[1] = 0.;
  _Attributes._xaxis[2] = 0.;

  _Attributes._yaxis[0] = 0.;
  _Attributes._yaxis[1] = 1.;
  _Attributes._yaxis[2] = 0.;

  _Attributes._zaxis[0] = 0.;
  _Attributes._zaxis[1] = 0.;
  _Attributes._zaxis[2] = 1.;

  const double * const matrix = meta_image.TransformMatrix();
  memcpy(_Attributes._xaxis, matrix + 0 * ndims, ndims);
  memcpy(_Attributes._yaxis, matrix + 1 * ndims, ndims);
  if (ndims > 2) {
    memcpy(_Attributes._zaxis, matrix + 2 * ndims, ndims);
  }

  _Slope     = meta_image.ElementToIntensityFunctionSlope();
  _Intercept = meta_image.ElementToIntensityFunctionOffset();

  // Set data type and number of bytes per voxels
  switch (meta_image.ElementType()) {
    case MET_CHAR:
       _DataType = MIRTK_VOXEL_CHAR;
       _Bytes    = 1;
       break;
    case MET_UCHAR:
       _DataType = MIRTK_VOXEL_UNSIGNED_CHAR;
       _Bytes    = 1;
       break;
    case MET_SHORT:
      _DataType = MIRTK_VOXEL_SHORT;
      _Bytes    = 2;
      break;
    case MET_USHORT:
      _DataType = MIRTK_VOXEL_UNSIGNED_SHORT;
      _Bytes    = 2;
      break;
    case MET_INT:
      _DataType = MIRTK_VOXEL_INT;
      _Bytes    = 4;
      break;
    case MET_UINT:
      _DataType = MIRTK_VOXEL_UNSIGNED_INT;
      _Bytes    = 4;
      break;
    case MET_FLOAT:
      _DataType = MIRTK_VOXEL_FLOAT;
      _Bytes    = 4;
      break;
    case MET_DOUBLE:
      _DataType = MIRTK_VOXEL_DOUBLE;
      _Bytes    = 8;
      break;
    default:
      char type_string[sizeof(MET_ValueTypeName[0])];
      if (MET_TypeToString(meta_image.ElementType(), type_string)) {
        cerr << this->NameOfClass() << ": Data type " << type_string << " not supported" << endl;
      } else {
        cerr << this->NameOfClass() << ": Data type " << meta_image.ElementType() << " not supported" << endl;
      }
      exit(1);
  }

  // Data starts here 
  _Start = static_cast<int>(meta_image.HeaderSize());
}

// -----------------------------------------------------------------------------
void MetaImageReader::ReadHeader()
{
  if (!_MetaImage->Read(_FileName.c_str(), false)) {
    cerr << this->NameOfClass() << ": Failed to read MetaImage file " << _FileName << endl;
    exit(1);
  }

  CopyHeader(*_MetaImage);
}

// -----------------------------------------------------------------------------
BaseImage *MetaImageReader::Run()
{
  if (!_MetaImage->Read(_FileName.c_str(), true)) {
    cerr << this->NameOfClass() << ": Failed to read MetaImage file " << _FileName << endl;
    exit(1);
  }

  CopyHeader(*_MetaImage);

  const int n = _Attributes.NumberOfLatticePoints();
  UniquePtr<BaseImage> output(BaseImage::New(_DataType));
  output->Initialize(_Attributes);
  memcpy(output->GetScalarPointer(), _MetaImage->ElementData(), n * _Bytes);
  this->Finalize(output.get());

  return output.release();
}


} // namespace mirtk

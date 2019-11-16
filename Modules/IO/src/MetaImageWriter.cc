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

#include "mirtk/MetaImageWriter.h"
#include "mirtk/ImageWriterFactory.h"

#include "mirtk/Memory.h"
#include "mirtk/Matrix.h"


namespace mirtk {


// Register image reader with object factory during static initialization
mirtkAutoRegisterImageWriterMacro(MetaImageWriter);


// -----------------------------------------------------------------------------
Array<string> MetaImageWriter::Extensions()
{
  Array<string> exts(2);
  exts[0] = ".mha"; // combined header and data file
  exts[1] = ".mhd"; // separate header file
  return exts;
}


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
MetaImageWriter::MetaImageWriter()
{
}

// -----------------------------------------------------------------------------
MetaImageWriter::~MetaImageWriter()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void MetaImageWriter::Run()
{
  this->Initialize();
  this->Close();

  const int channels = IsZero(_Input->TSize()) ? _Input->T() : 1;

  int ndims;
  if (_Input->T() > 1 && channels == 1) {
    ndims = 4;
  } else if (_Input->Z() > 1) {
    ndims = 3;
  } else {
    ndims = 2;
  }
  const int spatial_dims = max(ndims, 3);

  const int size[] = {
    _Input->X(),
    _Input->Y(),
    _Input->Z(),
    _Input->T(),
  };
  const double spacing[] = {
    _Input->XSize(),
    _Input->YSize(),
    _Input->ZSize(),
    _Input->TSize(),
  };

  MET_ValueEnumType element_type;
  const auto data_type = ImageDataType(_Input->GetDataType());
  switch (data_type) {
    case MIRTK_VOXEL_CHAR:
      element_type = MET_CHAR;
      break;
    case MIRTK_VOXEL_UNSIGNED_CHAR:
      element_type = MET_UCHAR;
      break;
    case MIRTK_VOXEL_SHORT:
      element_type = MET_SHORT;
      break;
    case MIRTK_VOXEL_UNSIGNED_SHORT:
      element_type = MET_USHORT;
      break;
    case MIRTK_VOXEL_INT:
      element_type = MET_INT;
      break;
    case MIRTK_VOXEL_UNSIGNED_INT:
      element_type = MET_UINT;
      break;
    case MIRTK_VOXEL_FLOAT:
      element_type = MET_FLOAT;
      break;
    case MIRTK_VOXEL_DOUBLE:
      element_type = MET_DOUBLE;
      break;
    default:
      cerr << this->NameOfClass() << ": Unsupported data type " << ToString(data_type) << endl;
      exit(1);
  }

  UniquePtr<MetaImage> meta_image;
  if (channels == 1) {
    void * const data = const_cast<void *>(_Input->GetDataPointer());
    meta_image.reset(new MetaImage(ndims, size, spacing, element_type, channels, data));
  } else {
    meta_image.reset(new MetaImage(ndims, size, spacing, element_type, channels));
    const size_t element_bytes = MET_ValueTypeSize[element_type];
    char *ptr = reinterpret_cast<char *>(meta_image->ElementData());
    for (int k = 0; k < _Input->Z(); ++k)
    for (int j = 0; j < _Input->Y(); ++j)
    for (int i = 0; i < _Input->X(); ++i) {
      for (int l = 0; l < _Input->T(); ++l) {
        memcpy(ptr, _Input->GetDataPointer(i, j, k, l), element_bytes);
        ptr += element_bytes;
      }
    }
  }

  meta_image->BinaryData(true);
  meta_image->CompressedData(true);

  double offset[3] = {0., 0., 0.};
  _Input->ImageToWorld(offset[0], offset[1], offset[2]);
  for (int i = 0; i < spatial_dims; ++i) {
    meta_image->Offset(i, offset[i]);
  }
  if (ndims > 3) {
    meta_image->Offset(3, _Input->GetTOrigin());
  }

  const auto matrix = _Input->Attributes().GetLatticeToWorldOrientation();
  for (int c = 0; c < spatial_dims; ++c)
  for (int r = 0; r < spatial_dims; ++r) {
    meta_image->TransformMatrix(r, c, matrix(r, c));
  }

  BaseImage::OrientationCode code[3];
  _Input->Orientation(code[0], code[1], code[2]);
  for (int i = 0; i < spatial_dims; ++i) {
    MET_OrientationEnumType met_orientation;
    switch (code[i]) {
      case BaseImage::L2R:
        met_orientation = MET_ORIENTATION_RL;
        break;
      case BaseImage::R2L:
        met_orientation = MET_ORIENTATION_LR;
        break;
      case BaseImage::A2P:
        met_orientation = MET_ORIENTATION_AP;
        break;
      case BaseImage::P2A:
        met_orientation = MET_ORIENTATION_PA;
        break;
      case BaseImage::S2I:
        met_orientation = MET_ORIENTATION_SI;
        break;
      case BaseImage::I2S:
        met_orientation = MET_ORIENTATION_IS;
        break;
      case BaseImage::Orientation_Unknown:
        met_orientation = MET_ORIENTATION_UNKNOWN;
        break;
    };
    meta_image->AnatomicalOrientation(i, met_orientation);
  }

  if (!meta_image->Write(_FileName.c_str())) {
    cerr << this->NameOfClass() << ": Failed to write image to file '" << _FileName << "'" << endl;
    exit(1);
  }

  this->Finalize();
}


} // namespace mirtk

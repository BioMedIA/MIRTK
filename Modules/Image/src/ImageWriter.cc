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

#include "mirtk/ImageWriter.h"

#include "mirtk/ImageWriterFactory.h"


namespace mirtk {


// =============================================================================
// Factory method
// =============================================================================

// -----------------------------------------------------------------------------
ImageWriter *ImageWriter::New(const char *fname)
{
  ImageWriter *writer = ImageWriterFactory::Instance().New(fname);
  if (!writer) {
    cerr << NameOfType() << "::New: Unsupported file name extension/format: " << fname << endl;
    exit(1);
  }
  writer->FileName(fname);
  return writer;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
ImageWriter::ImageWriter()
:
  _Input(NULL),
  _Start(0),
  _ReflectX(false),
  _ReflectY(false),
  _ReflectZ(false)
{
}

// -----------------------------------------------------------------------------
ImageWriter::~ImageWriter()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void ImageWriter::Initialize()
{
  // Check inputs and outputs
  if (!_Input) {
    cerr << this->NameOfClass() << "::Initialize: Filter has no input image" << endl;
    exit(1);
  }
  if (_FileName.empty()) {
    cerr << this->NameOfClass() << "::Initialize: Output file name not set" << endl;
    exit(1);
  }

  // Open file for writing
  this->Open(_FileName.c_str());

  // Reflect if necessary
  BaseImage *input = const_cast<BaseImage *>(_Input);
  if (_ReflectX) input->ReflectX();
  if (_ReflectY) input->ReflectY();
  if (_ReflectZ) input->ReflectZ();
}

// -----------------------------------------------------------------------------
void ImageWriter::Run()
{
  this->Initialize();

  const void *data = _Input->GetDataPointer();
  const int   n    = _Input->NumberOfVoxels();

  switch (_Input->GetDataType()) {
    case MIRTK_VOXEL_CHAR:           WriteAsChar  ((char           *)data, n, _Start); break;
    case MIRTK_VOXEL_UNSIGNED_CHAR:  WriteAsUChar ((unsigned char  *)data, n, _Start); break;
    case MIRTK_VOXEL_SHORT:          WriteAsShort ((short          *)data, n, _Start); break;
    case MIRTK_VOXEL_UNSIGNED_SHORT: WriteAsUShort((unsigned short *)data, n, _Start); break;
    case MIRTK_VOXEL_FLOAT:          WriteAsFloat ((float          *)data, n, _Start); break;
    case MIRTK_VOXEL_DOUBLE:         WriteAsDouble((double         *)data, n, _Start); break;
    default:
      cerr << this->NameOfClass() << "::Run(): Unsupported voxel type" << endl;
      exit(1);
  }

  this->Finalize();
}

// -----------------------------------------------------------------------------
void ImageWriter::Finalize()
{
  // Close file
  this->Close();

  // Reflect back if necessary
  BaseImage *input = const_cast<BaseImage *>(_Input);
  if (_ReflectX) input->ReflectX();
  if (_ReflectY) input->ReflectY();
  if (_ReflectZ) input->ReflectZ();
}


} // namespace mirtk

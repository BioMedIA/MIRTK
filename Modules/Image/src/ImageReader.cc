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

#include "mirtk/ImageReader.h"
#include "mirtk/ImageReaderFactory.h"

#include "mirtk/Voxel.h"
#include "mirtk/GenericImage.h"


namespace mirtk {


// =============================================================================
// Factory method
// =============================================================================

// -----------------------------------------------------------------------------
ImageReader *ImageReader::TryNew(const char *fname)
{
  ImageReader *reader = ImageReaderFactory::Instance().New(fname);
  if (reader) {
    reader->FileName(fname);
    reader->Initialize();
  }
  return reader;
}

// -----------------------------------------------------------------------------
ImageReader *ImageReader::New(const char *fname)
{
  ImageReader *reader = TryNew(fname);
  if (!reader) {
    cerr << NameOfType() << "::New: Unsupported file format or image reader not available for: " << fname << endl;
    exit(1);
  }
  return reader;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
ImageReader::ImageReader()
{
  _DataType  = MIRTK_VOXEL_UNKNOWN;
  _Slope     = 1.;
  _Intercept = 0.;
  _ReflectX  = false;
  _ReflectY  = false;
  _ReflectZ  = false;
  _Start     = 0;
  _Bytes     = 0;
}

// -----------------------------------------------------------------------------
ImageReader::~ImageReader()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void ImageReader::Initialize()
{
  // Close old file
  this->Close();

  // Open new file for reading
  this->Open(_FileName.c_str());

  // Read header
  this->ReadHeader();
}

// -----------------------------------------------------------------------------
BaseImage *ImageReader::Run()
{
  const int n = _Attributes.NumberOfLatticePoints();
  BaseImage *output = NULL;

  switch (_DataType) {
    case MIRTK_VOXEL_CHAR: {
        output = new GenericImage<char>(_Attributes);
        this->ReadAsChar((char *)output->GetDataPointer(), n, _Start);
      } break;

    case MIRTK_VOXEL_UNSIGNED_CHAR: {
        output = new GenericImage<unsigned char>(_Attributes);
        this->ReadAsUChar((unsigned char *)output->GetDataPointer(), n, _Start);
      } break;

    case MIRTK_VOXEL_SHORT: {
        output = new GenericImage<short>(_Attributes);
        this->ReadAsShort((short *)output->GetDataPointer(), n, _Start);
      } break;

    case MIRTK_VOXEL_UNSIGNED_SHORT: {
        output = new GenericImage<unsigned short>(_Attributes);
        this->ReadAsUShort((unsigned short *)output->GetDataPointer(), n, _Start);
      } break;

    case MIRTK_VOXEL_INT: {
        output = new GenericImage<int>(_Attributes);
        this->ReadAsInt((int *)output->GetDataPointer(), n, _Start);
      } break;

    case MIRTK_VOXEL_FLOAT: {
        output = new GenericImage<float>(_Attributes);
        this->ReadAsFloat((float *)output->GetDataPointer(), n, _Start);

        // Set background value to NaN if image contains NaNs
        float *p = reinterpret_cast<float *>(output->GetDataPointer());
        for (int i = 0; i < n; ++i, ++p) {
          if (IsNaN(*p)) {
            output->PutBackgroundValueAsDouble(numeric_limits<float>::quiet_NaN());
            break;
          }
        }
      } break;

    case MIRTK_VOXEL_DOUBLE: {
        output = new GenericImage<double>(_Attributes);
        this->ReadAsDouble((double *)output->GetDataPointer(), n, _Start);

        // Set background value to NaN if image contains NaNs
        double *p = reinterpret_cast<double *>(output->GetDataPointer());
        for (int i = 0; i < n; ++i, ++p) {
          if (IsNaN(*p)) {
            output->PutBackgroundValueAsDouble(numeric_limits<double>::quiet_NaN());
            break;
          }
        }
      } break;

    default:
      cerr << "FileToImage::Run: Unknown voxel type" << endl;
      exit(1);
  }

  if (_ReflectX) output->ReflectX();
  if (_ReflectY) output->ReflectY();
  if (_ReflectZ) output->ReflectZ();

  return output;
}

// -----------------------------------------------------------------------------
void ImageReader::Print() const
{
  cout << "Name of class is " << this->NameOfClass() << "\n";
  cout << "File name is " << _FileName << "\n";
  cout << "Image dimensions are " << _Attributes._x << " " << _Attributes._y << " " << _Attributes._z << " " << _Attributes._t << "\n";
  cout << "Image has " << _Bytes << " bytes per voxel" << "\n";
  cout << "Voxel dimensions are " << _Attributes._dx << " " << _Attributes._dy << " "
                                  << _Attributes._dz << " " << _Attributes._dt << "\n";
  cout << "Voxel type is ";
  switch (_DataType) {
    case MIRTK_VOXEL_CHAR:
      cout << "char";
      break;
    case MIRTK_VOXEL_UNSIGNED_CHAR:
      cout << "unsigned char";
      break;
    case MIRTK_VOXEL_SHORT:
      cout << "short";
      break;
    case MIRTK_VOXEL_UNSIGNED_SHORT:
      cout << "unsigned short";
      break;
    case MIRTK_VOXEL_INT:
      cout << "int";
      break;
    case MIRTK_VOXEL_UNSIGNED_INT:
      cout << "unsigned int";
      break;
    case MIRTK_VOXEL_FLOAT:
      cout << "float";
      break;
    case MIRTK_VOXEL_DOUBLE:
      cout << "double";
      break;
    default:
      cout << "unknown";
  }
  cout << "\n";
}


} // namespace mirtk

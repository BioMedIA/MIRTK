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

#include "PGM.h"

#include "mirtk/PGMImageWriter.h"
#include "mirtk/ImageWriterFactory.h"

#include "mirtk/Math.h"

#include <cstdio>


namespace mirtk {


// Register image reader with object factory during static initialization
mirtkAutoRegisterImageWriterMacro(PGMImageWriter);


// -----------------------------------------------------------------------------
Array<string> PGMImageWriter::Extensions()
{
  Array<string> exts(1);
  exts[0] = ".pgm";
  return exts;
}

// -----------------------------------------------------------------------------
void PGMImageWriter::Initialize()
{
  // Initialize base class
  ImageWriter::Initialize();

  if (_Input->Z() != 1) {
    cerr << this->NameOfClass() << " supports only 2D images" << endl;
    exit(1);
  }

  // Construct header
  const int sz = 256;
  char      header[sz];

  snprintf(header, sz, "%s\n# Created by %s\n%d %d\n255\n", PGM_MAGIC,
          this->NameOfClass(), _Input->X(), _Input->Y());
 
  // Data address/length of header
  _Start = static_cast<int>(strlen(header));

  // Write header
  this->Write(header, 0, _Start * sizeof(char));
}

// -----------------------------------------------------------------------------
void PGMImageWriter::Run()
{
  // Initialize filter
  this->Initialize();

  // Find dynamic range in image
  double min, max;
  _Input->GetMinMaxAsDouble(min, max);

  // Write image data
  const void * const data = _Input->GetDataPointer();
  const int          n    = _Input->NumberOfVoxels();

  switch (_Input->GetDataType()) {
    case MIRTK_VOXEL_CHAR: {
      unsigned char *buffer = new unsigned char[n];
      const char    *ptr    = reinterpret_cast<const char *>(data);
      for (int i = 0; i < n; ++i, ++ptr) {
        buffer[i] = static_cast<char>(round(255 * (static_cast<double>(*ptr) - min) / (max - min)));
      }
      this->WriteAsUChar(buffer, n, _Start);
      delete[] buffer;
    } break;

    case MIRTK_VOXEL_UNSIGNED_CHAR: {
      this->WriteAsUChar(reinterpret_cast<const unsigned char *>(data), n, _Start);
    } break;

    case MIRTK_VOXEL_SHORT: {
      unsigned char *buffer = new unsigned char[n];
      const short   *ptr    = reinterpret_cast<const short *>(data);
      for (int i = 0; i < n; ++i, ++ptr) {
        buffer[i] = static_cast<char>(round(255 * (static_cast<double>(*ptr) - min) / (max - min)));
      }
      this->WriteAsUChar(buffer, n, _Start);
      delete[] buffer;
    } break;

    case MIRTK_VOXEL_FLOAT: {
      unsigned char *buffer = new unsigned char[n];
      const float   *ptr    = reinterpret_cast<const float *>(data);
      for (int i = 0; i < n; ++i, ++ptr) {
        buffer[i] = static_cast<char>(round(255 * (static_cast<double>(*ptr) - min) / (max - min)));
      }
      this->WriteAsUChar(buffer, n, _Start);
      delete[] buffer;
    } break;

    case MIRTK_VOXEL_DOUBLE: {
      unsigned char *buffer = new unsigned char[n];
      const double  *ptr    = reinterpret_cast<const double *>(data);
      for (int i = 0; i < n; ++i, ++ptr) {
        buffer[i] = static_cast<char>(round(255 * (*ptr - min) / (max - min)));
      }
      this->WriteAsUChar(buffer, n, _Start);
      delete[] buffer;
    }

    default:
      cerr << this->NameOfClass() << "::Run(): Unsupported voxel type" << endl;
      exit(1);
  }

  // Finalize filter
  this->Finalize();
}


} // namespace mirtk


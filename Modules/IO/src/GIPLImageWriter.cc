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

#include "GIPL.h"

#include "mirtk/GIPLImageWriter.h"
#include "mirtk/ImageWriterFactory.h"


namespace mirtk {


// Register image reader with object factory during static initialization
mirtkAutoRegisterImageWriterMacro(GIPLImageWriter);


// -----------------------------------------------------------------------------
Array<string> GIPLImageWriter::Extensions()
{
  Array<string> exts(1);
  exts[0] = ".gipl";
  return exts;
}

// -----------------------------------------------------------------------------
void GIPLImageWriter::Initialize()
{
  // Initialize base class
  ImageWriter::Initialize();

  // Write image dimensions
  this->WriteAsShort(_Input->X(), 0);
  this->WriteAsShort(_Input->Y(), 2);
  this->WriteAsShort(_Input->Z(), 4);
  this->WriteAsShort(_Input->T(), 6);

  // Write type
  switch (_Input->GetScalarType()) {
    case MIRTK_VOXEL_CHAR: {
        this->WriteAsShort(GIPL_CHAR, 8);
        break;
      }
    case MIRTK_VOXEL_UNSIGNED_CHAR: {
        this->WriteAsShort(GIPL_U_CHAR, 8);
        break;
     }
    case MIRTK_VOXEL_SHORT: {
        this->WriteAsShort(GIPL_SHORT, 8);
        break;
      }
    case MIRTK_VOXEL_UNSIGNED_SHORT: {
        this->WriteAsShort(GIPL_U_SHORT, 8);
        break;
      }
    case MIRTK_VOXEL_FLOAT: {
        this->WriteAsShort(GIPL_FLOAT, 8);
        break;
      }
    case MIRTK_VOXEL_DOUBLE: {
        this->WriteAsShort(GIPL_DOUBLE, 8);
        break;
      }
    default:
      cerr << this->NameOfClass() << "::Initialize: Not supported for this image type" << endl;
      exit(1);
  }

  // Write voxel dimensions
  double xsize, ysize, zsize, tsize;
  _Input->GetPixelSize(&xsize, &ysize, &zsize, &tsize);
  this->WriteAsFloat(static_cast<float>(xsize), 10);
  this->WriteAsFloat(static_cast<float>(ysize), 14);
  this->WriteAsFloat(static_cast<float>(zsize), 18);
  this->WriteAsFloat(static_cast<float>(tsize), 22);

  // Write patient description
  for (int i = 0; i < 80; i++) {
    this->WriteAsChar(char(0), 26 + i*sizeof(char));
  }

  // Write image orientation
  double xaxis[3], yaxis[3], zaxis[3];
  _Input->GetOrientation(xaxis, yaxis, zaxis);
  this->WriteAsFloat(static_cast<float>(xaxis[0]), 106);
  this->WriteAsFloat(static_cast<float>(xaxis[1]), 110);
  this->WriteAsFloat(static_cast<float>(xaxis[2]), 114);
  this->WriteAsFloat(static_cast<float>(yaxis[0]), 122);
  this->WriteAsFloat(static_cast<float>(yaxis[1]), 126);
  this->WriteAsFloat(static_cast<float>(yaxis[2]), 130);

  // Write identifier
  this->WriteAsInt(0, 154);

  // Write spare bytes
  for (int i = 0; i < 28; i++) {
    this->WriteAsChar(char(0), 158 + i*sizeof(char));
  }

  // Write flag and orientation
  this->WriteAsChar(char(0), 186);
  this->WriteAsChar(char(0), 187);

  // Write min and max values
  double min, max;
  _Input->GetMinMaxAsDouble(&min, &max);
  this->WriteAsDouble(min, 188);
  this->WriteAsDouble(max, 196);

  // Write origin
  double pos[3];
  _Input->GetOrigin(pos[0], pos[1], pos[2]);
  this->WriteAsDouble(pos[0], 204);
  this->WriteAsDouble(pos[1], 212);
  this->WriteAsDouble(pos[2], 220);
  this->WriteAsDouble(.0, 228);

  // Write magic number
  this->WriteAsFloat(.0f, 236);
  this->WriteAsFloat(.0f, 240);
  this->WriteAsUInt(GIPL_MAGIC_EXT, 244);
  this->WriteAsFloat(.0f, 248);
  this->WriteAsUInt(GIPL_MAGIC1, 252);

  // Calculate data address
  _Start = GIPL_HEADERSIZE;
}


} // namespace mirtk

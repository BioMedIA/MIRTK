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

#include "mirtk/GIPLImageReader.h"
#include "mirtk/ImageReaderFactory.h"


namespace mirtk {


// Register image reader with object factory during static initialization
mirtkAutoRegisterImageReaderMacro(GIPLImageReader);


// -----------------------------------------------------------------------------
bool GIPLImageReader::CheckHeader(const char *fname)
{
  unsigned int magic_number;
  Cifstream from(fname);
  from.ReadAsUInt(&magic_number, 1, 252);
  from.Close();
  return magic_number == GIPL_MAGIC1 || magic_number == GIPL_MAGIC2;
}

// -----------------------------------------------------------------------------
bool GIPLImageReader::CanRead(const char *fname) const
{
  return CheckHeader(fname);
}

// -----------------------------------------------------------------------------
void GIPLImageReader::ReadHeader()
{
  float          size;
  double         origin;
  unsigned short dim, type;
  unsigned int   magic_number, magic_ext;

  // Read header
  this->ReadAsUInt(&magic_number, 1, 252);

  // Check header
  if (magic_number != GIPL_MAGIC1 && magic_number != GIPL_MAGIC2) {
    cerr << this->NameOfClass() << "::ReadHeader: Can't read magic number" << endl;
    exit(1);
  }

  // Read image dimensions
  this->ReadAsUShort(&dim, 1, 0);
  _Attributes._x = dim;
  this->ReadAsUShort(&dim, 1, 2);
  _Attributes._y = dim;
  this->ReadAsUShort(&dim, 1, 4);
  _Attributes._z = dim;
  this->ReadAsUShort(&dim, 1, 6);
  _Attributes._t = dim;
  if(_Attributes._t <= 0) _Attributes._t = 1;

  // Read voxel dimensions
  this->ReadAsFloat(&size, 1, 10);
  _Attributes._dx = size;
  this->ReadAsFloat(&size, 1, 14);
  _Attributes._dy = size;
  this->ReadAsFloat(&size, 1, 18);
  _Attributes._dz = size;
  this->ReadAsFloat(&size, 1, 22);
  _Attributes._dt = size;
  if(_Attributes._dt <= 0) _Attributes._dt = 1;

  // Read extension flag
  this->ReadAsUInt(&magic_ext, 1, 244);

  // Check if additional information about orientation and origin is available
  if (magic_ext == GIPL_MAGIC_EXT) {

    // Read image orientation if available
    this->ReadAsFloat(&size, 1, 106);
    _Attributes._xaxis[0] = size;
    this->ReadAsFloat(&size, 1, 110);
    _Attributes._xaxis[1] = size;
    this->ReadAsFloat(&size, 1, 114);
    _Attributes._xaxis[2] = size;
    this->ReadAsFloat(&size, 1, 122);
    _Attributes._yaxis[0] = size;
    this->ReadAsFloat(&size, 1, 126);
    _Attributes._yaxis[1] = size;
    this->ReadAsFloat(&size, 1, 130);
    _Attributes._yaxis[2] = size;

    // Construct the z-axis.
    // Assume a right handed (`neurological') set of axes!
    _Attributes._zaxis[0] = _Attributes._xaxis[1]*_Attributes._yaxis[2] - _Attributes._xaxis[2]*_Attributes._yaxis[1];
    _Attributes._zaxis[1] = _Attributes._xaxis[2]*_Attributes._yaxis[0] - _Attributes._xaxis[0]*_Attributes._yaxis[2];
    _Attributes._zaxis[2] = _Attributes._xaxis[0]*_Attributes._yaxis[1] - _Attributes._xaxis[1]*_Attributes._yaxis[0];

    // Read image origin if available
    this->ReadAsDouble(&origin, 1, 204);
    _Attributes._xorigin = origin;
    this->ReadAsDouble(&origin, 1, 212);
    _Attributes._yorigin = origin;
    this->ReadAsDouble(&origin, 1, 220);
    _Attributes._zorigin = origin;
  }

  // Read type of voxels
  this->ReadAsUShort(&type, 1, 8);

  // Calculate type of voxels and number of bytes per voxel
  switch (type) {
    case GIPL_CHAR:
      _DataType = MIRTK_VOXEL_CHAR;
      _Bytes    = 1;
      break;
    case GIPL_U_CHAR:
      _DataType = MIRTK_VOXEL_UNSIGNED_CHAR;
      _Bytes    = 1;
      break;
    case GIPL_SHORT:
      _DataType = MIRTK_VOXEL_SHORT;
      _Bytes    = 2;
      break;
    case GIPL_U_SHORT:
      _DataType = MIRTK_VOXEL_UNSIGNED_SHORT;
      _Bytes    = 2;
      break;
    case GIPL_INT:
      _DataType = MIRTK_VOXEL_INT;
      _Bytes    = 4;
      break;
    case GIPL_U_INT:
      _DataType = MIRTK_VOXEL_UNSIGNED_INT;
      _Bytes    = 4;
      break;
    case GIPL_FLOAT:
      _DataType = MIRTK_VOXEL_FLOAT;
      _Bytes    = 4;
      break;
    case GIPL_DOUBLE:
      _DataType = MIRTK_VOXEL_DOUBLE;
      _Bytes    = 8;
      break;
    default:
      cerr << this->NameOfClass() << "::ReadHeader(): Unknown voxel type" << endl;
      exit(1);
  }

  // Data starts here
  _Start = GIPL_HEADERSIZE;
}


} // namespace mirtk

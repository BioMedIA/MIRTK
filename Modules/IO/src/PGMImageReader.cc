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

#include "mirtk/PGMImageReader.h"
#include "mirtk/ImageReaderFactory.h"

#include "mirtk/Config.h" // WINDOWS

#include <cstdio>


namespace mirtk {


// Register image reader with object factory during static initialization
mirtkAutoRegisterImageReaderMacro(PGMImageReader);


// -----------------------------------------------------------------------------
bool PGMImageReader::CheckHeader(const char *fname)
{
  char buffer[255];
  ifstream from(fname);
  do {
    from.get(buffer, 255);
    from.seekg(1, ios::cur);
  } while (buffer[0] == '#');
  from.close();
  return (strcmp(buffer, PGM_MAGIC) == 0);
}

// -----------------------------------------------------------------------------
bool PGMImageReader::CanRead(const char *fname) const
{
  return CheckHeader(fname);
}

// -----------------------------------------------------------------------------
void PGMImageReader::ReadHeader()
{
  const long bufsz = 256l;
  char buffer[bufsz];

  // Read header, skip comments
  do {
    this->ReadAsString(buffer, bufsz-1l);
  } while (buffer[0] == '#');

  // Check header
  if (strcmp(buffer, PGM_MAGIC) != 0) {
    cerr << this->NameOfClass() << "::ReadHeader: Can't read magic number: " << buffer << endl;
    exit(1);
  }

  // Read voxel dimensions, skip comments
  do {
    this->ReadAsString(buffer, bufsz-1l);
  } while (buffer[0] == '#');

  // Parse voxel dimensions
#ifdef WINDOWS
  sscanf_s(buffer, "%d %d", &_Attributes._x, &_Attributes._y);
#else
  sscanf(buffer, "%d %d", &_Attributes._x, &_Attributes._y);
#endif

  // Ignore maximum greyvalue, skip comments
  do {
    this->ReadAsString(buffer, bufsz - 1l);
  } while (buffer[0] == '#');

  // PGM files support only 2D images, so set z and t to 1
  _Attributes._z = 1;
  _Attributes._t = 1;

  // PGM files do not have voxel dimensions, so set them to default values
  _Attributes._dx = 1;
  _Attributes._dy = 1;
  _Attributes._dz = 1;
  _Attributes._dt = 1;

  // PGM files have voxels which are unsigned char
  _DataType = MIRTK_VOXEL_UNSIGNED_CHAR;
  _Bytes    = 1;

  // Data starts here
  _Start = this->Tell();
}


} // namespace mirtk

/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2016 Imperial College London
 * Copyright 2008-2016 Daniel Rueckert, Julia Schnabel
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

#include "mirtk/Cofstream.h"

#include "mirtk/Config.h"       // WINDOWS
#include "mirtk/CommonConfig.h" // MIRTK_Common_WITH_ZLIB
#include "mirtk/Memory.h"       // swap16, swap32

#if MIRTK_Common_WITH_ZLIB
#  include <zlib.h>
#endif


namespace mirtk {


// -----------------------------------------------------------------------------
Cofstream::Cofstream(const char *fname)
:
  _File(nullptr),
  #if MIRTK_Common_WITH_ZLIB
    _ZFile(nullptr),
  #endif
  _Compressed(false),
  _Swapped(GetByteOrder() == LittleEndian)
{
  if (fname) Open(fname);
}

// -----------------------------------------------------------------------------
Cofstream::~Cofstream()
{
  Close();
}

// -----------------------------------------------------------------------------
void Cofstream::Open(const char *fname)
{
  Close();
  size_t len = strlen(fname);
  if (len == 0) {
    Throw(ERR_InvalidArgument, __func__, "Filename is empty");
  }
  if (len > 3 && (strncmp(fname + len-3, ".gz", 3) == 0 || strncmp(fname + len-3, ".GZ", 3) == 0)) {
    #if MIRTK_Common_WITH_ZLIB
      gzFile fp = gzopen(fname, "wb");
      if (fp == nullptr) {
        Throw(ERR_IOError, __func__, "Failed to open file ", fname);
      }
      _ZFile = fp;
      _Compressed = true;
    #else // MIRTK_Common_WITH_ZLIB
      Throw(ERR_IOError, __func__, "Cannot write compressed file when Common module not built WITH_ZLIB");
    #endif // MIRTK_Common_WITH_ZLIB
  } else {
    #ifdef WINDOWS
      errno_t err = fopen_s(&_File, fname, "wb");
      if (err != 0) _File = nullptr;
    #else
      _File = fopen(fname, "wb");
    #endif
    if (_File == nullptr) {
      Throw(ERR_IOError, __func__, "Failed to open file ", fname);
    }
    _Compressed = false;
  }
}

// -----------------------------------------------------------------------------
void Cofstream::Close()
{
  #if MIRTK_Common_WITH_ZLIB
    if (_ZFile != nullptr) {
      gzclose(reinterpret_cast<gzFile>(_ZFile));
      _ZFile = nullptr;
    }
  #endif // MIRTK_Common_WITH_ZLIB
  if (_File != nullptr) {
    fclose(_File);
    _File = nullptr;
  }
}

// -----------------------------------------------------------------------------
int Cofstream::IsCompressed() const
{
  return _Compressed ? 1 : 0;
}

// -----------------------------------------------------------------------------
void Cofstream::IsCompressed(int compressed)
{
  _Compressed = (compressed != 0);
}

// -----------------------------------------------------------------------------
int Cofstream::IsSwapped() const
{
  return _Swapped ? 1 : 0;
}

// -----------------------------------------------------------------------------
void Cofstream::IsSwapped(int swapped)
{
  _Swapped = (swapped != 0);
}

// -----------------------------------------------------------------------------
bool Cofstream::Write(const char *data, long offset, long length)
{
  #if MIRTK_Common_WITH_ZLIB
    if (_Compressed) {
      gzFile fp = reinterpret_cast<gzFile>(_ZFile);
      if (offset != -1) {
        if (gztell(fp) > offset) {
          Throw(ERR_IOError, __func__, "Writing compressed files only supports"
                " forward seek (pos=", gztell(fp), ", offset=", offset, ")");
        }
        gzseek(fp, offset, SEEK_SET);
      }
      return (gzwrite(fp, data, length) == length);
    }
  #endif // MIRTK_Common_WITH_ZLIB
  if (offset != -1) fseek(_File, offset, SEEK_SET);
  return (fwrite(data, length, 1, _File) == 1);
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsChar(const char data, long offset)
{
  return Write(&data, offset, sizeof(const char));
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsChar(const char *data, long length, long offset)
{
  return Write(data, offset, length*sizeof(const char));
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsUChar(unsigned char data, long offset)
{
  return Write((const char *)&data, offset, sizeof(unsigned char));
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsUChar(const unsigned char *data, long length, long offset)
{
  return Write((const char *)data, offset, length*sizeof(unsigned char));
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsShort(short data, long offset)
{
  if (_Swapped) swap16((char *)&data, (char *)&data, 1);
  return Write((const char *)&data, offset, sizeof(short));
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsShort(const short *data, long length, long offset)
{
  if (_Swapped) swap16((char *)data, (char *)data, length);
  if (!Write((const char *)data, offset, length*sizeof(short))) return false;
  if (_Swapped) swap16((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsUShort(unsigned short data, long offset)
{
  if (_Swapped) swap16((char *)&data, (char *)&data, 1);
  return Write((const char *)&data, offset, sizeof(unsigned short));
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsUShort(const unsigned short *data, long length, long offset)
{
  if (_Swapped) swap16((char *)data, (char *)data, length);
  if (!Write((const char *)data, offset, length*sizeof(unsigned short))) return false;
  if (_Swapped) swap16((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsInt(int data, long offset)
{
  if (_Swapped) swap32((char *)&data, (char *)&data, 1);
  return Write((const char *)&data, offset, sizeof(int));
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsInt(const int *data, long length, long offset)
{
  if (_Swapped) swap32((char *)data, (char *)data, length);
  if (!Write((const char *)data, offset, length*sizeof(int))) return false;
  if (_Swapped) swap32((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsUInt(unsigned int data, long offset)
{
  if (_Swapped) swap32((char *)&data, (char *)&data, 1);
  return Write((const char *)&data, offset, sizeof(unsigned int));
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsUInt(const unsigned int *data, long length, long offset)
{
  if (_Swapped) swap32((char *)data, (char *)data, length);
  if (!Write((const char *)data, offset, length*sizeof(unsigned int))) return false;
  if (_Swapped) swap32((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsFloat(float data, long offset)
{
  if (_Swapped) swap32((char *)&data, (char *)&data, 1);
  return Write((const char *)&data, offset, sizeof(float));
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsFloat(const float *data, long length, long offset)
{
  if (_Swapped) swap32((char *)data, (char *)data, length);
  if (!Write((const char *)data, offset, length*sizeof(float))) return false;
  if (_Swapped) swap32((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsDouble(const double *data, long length, long offset)
{
  if (_Swapped) swap64((char *)data, (char *)data, length);
  if (!Write((const char *)data, offset, length*sizeof(double))) return false;
  if (_Swapped) swap64((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsDouble(double data, long offset)
{
  if (_Swapped) swap64((char *)&data, (char *)&data, 1);
  return Write((const char *)&data, offset, sizeof(double));
}

// -----------------------------------------------------------------------------
bool Cofstream::WriteAsString(const char *data, long offset)
{
  #if MIRTK_Common_WITH_ZLIB
    if (_Compressed) {
      gzFile fp = reinterpret_cast<gzFile>(_ZFile);
      if (offset != -1) gzseek(fp, offset, SEEK_SET);
      return (gzputs(fp, data) > 0);
    }
  #endif // MIRTK_Common_WITH_ZLIB
  if (offset != -1) fseek(_File, offset, SEEK_SET);
  return (fputs(data, _File) != EOF);
}


} // namespace mirtk

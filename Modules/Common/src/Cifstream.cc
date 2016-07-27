/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#include "mirtk/Cifstream.h"

#include "mirtk/Config.h"       // WINDOWS
#include "mirtk/CommonConfig.h" // MIRTK_Common_WITH_ZLIB
#include "mirtk/Memory.h"       // swap16, swap32

#if MIRTK_Common_WITH_ZLIB
#  include <zlib.h>
#else
#  include <cstdio>
#endif


namespace mirtk {


// -----------------------------------------------------------------------------
Cifstream::Cifstream(const char *fname)
:
  _File(nullptr),
  _Swapped(GetByteOrder() == LittleEndian)
{
  if (fname) Open(fname);
}

// -----------------------------------------------------------------------------
Cifstream::~Cifstream()
{
  Close();
}

// -----------------------------------------------------------------------------
void Cifstream::Open(const char *fname)
{
  #if MIRTK_Common_WITH_ZLIB
    _File = gzopen(fname, "rb");
  #elif defined(WINDOWS)
    FILE *fp;
    errno_t err = fopen_s(&fp, fname, "rb");
    if (err != 0) _File = nullptr;
    _File = fp;
  #else
    _File = fopen(fname, "rb");
  #endif
  if (_File == nullptr) {
    Throw(ERR_IOError, __func__, "Failed to open file ", fname);
  }
}

// -----------------------------------------------------------------------------
void Cifstream::Close()
{
  if (_File != nullptr) {
    #if MIRTK_Common_WITH_ZLIB
      gzclose(reinterpret_cast<gzFile>(_File));
    #else
      fclose(reinterpret_cast<FILE *>(_File));
    #endif
    _File = nullptr;
  }
}

// -----------------------------------------------------------------------------
int Cifstream::IsSwapped() const
{
  return _Swapped ? 1 : 0;
}

// -----------------------------------------------------------------------------
void Cifstream::IsSwapped(int swapped)
{
  _Swapped = (swapped != 0);
}

// -----------------------------------------------------------------------------
long Cifstream::Tell() const
{
#if MIRTK_Common_WITH_ZLIB
  return gztell(reinterpret_cast<gzFile>(_File));
#else
  return ftell(reinterpret_cast<FILE *>(_File));
#endif
}

// -----------------------------------------------------------------------------
void Cifstream::Seek(long offset)
{
#if MIRTK_Common_WITH_ZLIB
  gzseek(reinterpret_cast<gzFile>(_File), offset, SEEK_SET);
#else
  fseek(reinterpret_cast<FILE *>(_File), offset, SEEK_SET);
#endif
}

// -----------------------------------------------------------------------------
bool Cifstream::Read(char *mem, long start, long num)
{
#if MIRTK_Common_WITH_ZLIB
  gzFile fp = reinterpret_cast<gzFile>(_File);
  if (start != -1) gzseek(fp, start, SEEK_SET);
  return (gzread(fp, mem, num) == num);
#else
  FILE * fp = reinterpret_cast<FILE *>(_File);
  if (start != -1) fseek(fp, start, SEEK_SET);
  return (fread(mem, num, 1, fp) == 1);
#endif
}

// -----------------------------------------------------------------------------
bool Cifstream::ReadAsChar(char *data, long length, long offset)
{
  return Read((char *)data, offset, length * sizeof(char));
}

// -----------------------------------------------------------------------------
bool Cifstream::ReadAsUChar(unsigned char *data, long length, long offset)
{
  return Read((char *)data, offset, length * sizeof(unsigned char));
}

// -----------------------------------------------------------------------------
bool Cifstream::ReadAsShort(short *data, long length, long offset)
{
  if (!Read((char *)data, offset, length * sizeof(short))) return false;
  if (_Swapped) swap16((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool Cifstream::ReadAsUShort(unsigned short *data, long length, long offset)
{
  if (!Read((char *)data, offset, length * sizeof(unsigned short))) return false;
  if (_Swapped) swap16((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool Cifstream::ReadAsInt(int *data, long length, long offset)
{
  if (!Read((char *)data, offset, length * sizeof(int))) return false;
  if (_Swapped) swap32((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool Cifstream::ReadAsUInt(unsigned int *data, long length, long offset)
{
  if (!Read((char *)data, offset, length * sizeof(unsigned int))) return false;
  if (_Swapped) swap32((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool Cifstream::ReadAsFloat(float *data, long length, long offset)
{
  if (!Read((char *)data, offset, length * sizeof(float))) return false;
  if (_Swapped) swap32((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool Cifstream::ReadAsDouble(double *data, long length, long offset)
{
  if (!Read((char *)data, offset, length * sizeof(double))) return false;
  if (_Swapped) swap64((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool Cifstream::ReadAsString(char *data, long length, long offset)
{
  // Read string
#if MIRTK_Common_WITH_ZLIB
  gzFile fp = reinterpret_cast<gzFile>(_File);
  if (offset != -1) gzseek(fp, offset, SEEK_SET);
  if (gzgets(fp, data, length) == Z_NULL) return false;
#else
  FILE *fp = reinterpret_cast<FILE *>(_File);
  if (offset != -1) fseek(fp, offset, SEEK_SET);
  if (fgets(data, length, fp) == nullptr) return false;
#endif // MIRTK_Common_WITH_ZLIB

  // Discard end-of-line character(s)
  const size_t len = strlen(data);
  if (len > 0 && (data[len-1] == '\n' || data[len-1] == '\r')) {
    data[len-1] = '\0';
    if (len > 1 && data[len-2] == '\r') {
      data[len-2] = '\0';
    }
  }
  return true;
}


} // namespace mirtk

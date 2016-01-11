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

#include <mirtkCofstream.h>

#include <mirtkCommonConfig.h>
#include <mirtkMemory.h> // swap16, swap32

#if MIRTK_Common_WITH_ZLIB
#  include <zlib.h>
#endif


namespace mirtk {


// -----------------------------------------------------------------------------
Cofstream::Cofstream(const char *fname)
{
  _File = NULL;
#if MIRTK_Common_WITH_ZLIB
  _ZFile = NULL;
#endif
#if MIRTK_Common_BIG_ENDIAN
  _Swapped = false;
#else
  _Swapped = true;
#endif // MIRTK_Common_BIG_ENDIAN
  _Compressed = false;
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
  size_t len = strlen(fname);
  if (len == 0) {
    cerr << "Cofstream::Open: Filename is empty!" << endl;
    exit(1);
  }
  if (len > 3 && (strncmp(fname + len-3, ".gz", 3) == 0 || strncmp(fname + len-3, ".GZ", 3) == 0)) {
#if MIRTK_Common_WITH_ZLIB
    _ZFile = gzopen(fname, "wb");
    if (_ZFile == NULL) {
      cerr << "Cofstream::Open: Cannot open file " << fname << endl;
      exit(1);
    }
    _Compressed = true;
#else // MIRTK_Common_WITH_ZLIB
    cerr << "Cofstream::Open: Cannot write compressed file when Common module not built WITH_ZLIB" << endl;
    exit(1);
#endif // MIRTK_Common_WITH_ZLIB
  } else {
    _File = fopen(fname, "wb");
    if (_File == NULL) {
      cerr << "Cofstream::Open: Cannot open file " << fname << endl;
      exit(1);
    }
    _Compressed = false;
  }
}

// -----------------------------------------------------------------------------
void Cofstream::Close()
{
#if MIRTK_Common_WITH_ZLIB
  if (_ZFile) {
    gzclose(_ZFile);
    _ZFile = NULL;
  }
#endif // MIRTK_Common_WITH_ZLIB
  if (_File) {
    fclose(_File);
    _File = NULL;
  }
}

// -----------------------------------------------------------------------------
int Cofstream::IsCompressed() const
{
  return static_cast<int>(_Compressed);
}

// -----------------------------------------------------------------------------
void Cofstream::IsCompressed(int compressed)
{
  _Compressed = static_cast<bool>(compressed);
}

// -----------------------------------------------------------------------------
int Cofstream::IsSwapped() const
{
  return static_cast<int>(_Swapped);
}

// -----------------------------------------------------------------------------
void Cofstream::IsSwapped(int swapped)
{
  _Swapped = static_cast<bool>(swapped);
}

// -----------------------------------------------------------------------------
bool Cofstream::Write(const char *data, long offset, long length)
{
#if MIRTK_Common_WITH_ZLIB
  if (_Compressed) {
    if (offset != -1) {
      if (gztell(_ZFile) > offset) {
        cerr << "Error: Writing compressed files only supports forward seek (pos="
                  << gztell(_ZFile) << ", offset=" << offset << ")" << endl;
        exit(1);
      }
      gzseek(_ZFile, offset, SEEK_SET);
    }
    return (gzwrite(_ZFile, data, length) == length);
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
    if (offset != -1) gzseek(_ZFile, offset, SEEK_SET);
    return (gzputs(_ZFile, data) > 0);
  }
#endif // MIRTK_Common_WITH_ZLIB
  if (offset != -1) fseek(_File, offset, SEEK_SET);
  return (fputs(data, _File) != EOF);
}


} // namespace mirtk

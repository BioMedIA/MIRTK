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

#ifndef MIRTK_Cofstream_H
#define MIRTK_Cofstream_H

#include "mirtk/CommonConfig.h"
#include "mirtk/CommonExport.h"

#include "mirtk/Object.h"


namespace mirtk {


/**
 * Class for writing (compressed) file streams.
 *
 * This class defines and implements functions for writing compressed file
 * streams. At the moment only the writing of uncompressed file streams is
 * supported.
 */

class Cofstream : public Object
{
  mirtkObjectMacro(Cofstream);

  /// File pointer to uncompressed file
  FILE *_File;

  /// File pointer to compressed file
#if MIRTK_Common_WITH_ZLIB
  void *_ZFile;
#endif

  /// Flag whether file is compressed
  mirtkPublicAttributeMacro(bool, Compressed);

  /// Flag whether file is swapped
  mirtkPublicAttributeMacro(bool, Swapped);

public:

  /// Constructor
  Cofstream(const char * = NULL);

  /// Destructor
  ~Cofstream();

  /// Write n data as array from offset
  bool Write(const char *, long, long);

  /// Write data as char (possibly compressed) from offset
  bool WriteAsChar(char, long = -1);
  /// Write data as char (possibly compressed) from offset
  bool WriteAsChar(const char *, long, long = -1);

  /// Write data as unsigned char (possibly compressed) from offset
  bool WriteAsUChar(unsigned char, long = -1);
  /// Write data as unsigned char (possibly compressed) from offset
  bool WriteAsUChar(const unsigned char *, long, long = -1);

  /// Write data as short (possibly compressed) from offset
  bool WriteAsShort(short, long = -1);
  /// Write data as short (possibly compressed) from offset
  bool WriteAsShort(const short *, long, long = -1);

  /// Write data as unsigned short (possibly compressed) from offset
  bool WriteAsUShort(unsigned short, long = -1);
  /// Write data as unsigned short (possibly compressed) from offset
  bool WriteAsUShort(const unsigned short *, long, long = -1);

  /// Write data as int (possibly compressed) from offset
  bool WriteAsInt(int, long = -1);
  /// Write data as int (possibly compressed) from offset
  bool WriteAsInt(const int *, long, long = -1);

  /// Write data as unsigned int (possibly compressed) from offset
  bool WriteAsUInt(unsigned int, long = -1);
  /// Write data as unsigned int (possibly compressed) from offset
  bool WriteAsUInt(const unsigned int *, long, long = -1);

  /// Write data as float (possibly compressed) from offset
  bool WriteAsFloat(float, long = -1);
  /// Write data as float (possibly compressed) from offset
  bool WriteAsFloat(const float *, long, long = -1);

  /// Write data as double (possibly compressed) from offset
  bool WriteAsDouble(double, long = -1);
  /// Write data as double (possibly compressed) from offset
  bool WriteAsDouble(const double *, long, long = -1);

  /// Write data as string (possibly compressed)
  bool WriteAsString(const char *, long = -1);

  /// Open file
  void Open(const char *);

  /// Close file
  void Close();

  /// Returns whether file is compressed
  /// \deprecated Used Compressed() instead.
  MIRTK_Common_DEPRECATED int IsCompressed() const;

  /// Sets whether file is compressed
  /// \deprecated Used Compressed(bool) instead.
  MIRTK_Common_DEPRECATED void IsCompressed(int);

  /// Returns whether file is swapped
  /// \deprecated Used Swapped() instead.
  MIRTK_Common_DEPRECATED int IsSwapped() const;

  /// Sets whether file is swapped
  /// \deprecated Used Swapped(bool) instead.
  MIRTK_Common_DEPRECATED void IsSwapped(int);

};


} // namespace mirtk

#endif // MIRTK_Cofstream_H

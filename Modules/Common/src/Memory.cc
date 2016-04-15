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

#include "mirtk/Memory.h"


namespace mirtk {


// ========================================================================
// Swap bytes
// ========================================================================

// -----------------------------------------------------------------------------
// See nifti_short_order implementation
ByteOrder GetByteOrder()
{
  union {
    unsigned char bb[2];
    short         ss;
  } fred;
  fred.bb[0] = 1;
  fred.bb[1] = 0;
  return (fred.ss == 1 ? LittleEndian : BigEndian);
}

// ------------------------------------------------------------------------
void swap16(char *a, char *b, long n)
{
  char c;
  for (int i = 0; i < n * 2; i += 2) {
    c = a[i];
    a[i] = b[i+1];
    b[i+1] = c;
  }
}

// ------------------------------------------------------------------------
void swap32(char *a, char *b, long n)
{
  char c;
  for (int i = 0; i < n * 4; i += 4) {
    c = a[i];
    a[i] = b[i+3];
    b[i+3] = c;
    c = a[i+1];
    a[i+1] = b[i+2];
    b[i+2] = c;
  }
}

// ------------------------------------------------------------------------
void swap64(char *a, char *b, long n)
{
  char c;
  for (int i = 0; i < n * 8; i += 8) {
    c = a[i];
    a[i] = b[i+7];
    b[i+7] = c;
    c = a[i+1];
    a[i+1] = b[i+6];
    b[i+6] = c;
    c = a[i+2];
    a[i+2] = b[i+5];
    b[i+5] = c;
    c = a[i+3];
    a[i+3] = b[i+4];
    b[i+4] = c;
  }
}


} // namespace mirtk

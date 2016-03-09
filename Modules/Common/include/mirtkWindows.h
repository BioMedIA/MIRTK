/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
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

#ifndef MIRTK_Windows_H
#define MIRTK_Windows_H


namespace mirtk {


// Windows is missing M_PI constants
#define M_PI 3.14159265358979323846
// Disable min/max macros/functions of windows.h
#define NOMINMAX
// Windows specific header file
#include <windows.h>
#include <stdio.h>
#include <stdarg.h>
// some code uses uint instead of unsigned int
typedef unsigned int uint;
// snprintf part of C++11 but not supported by all VS compilers
#define snprintf c99_snprintf
// copysign function
#if defined(_MSC_VER) && _MSC_VER < 1800
#  define copysign _copysign
#endif

// ----------------------------------------------------------------------------
inline int c99_vsnprintf(char* str, size_t size, const char* format, va_list ap)
{
  int count = -1;
  if (size  !=  0) count = _vsnprintf_s(str, size, _TRUNCATE, format, ap);
  if (count == -1) count = _vscprintf(format, ap);
  return count;
}

// ----------------------------------------------------------------------------
inline int c99_snprintf(char* str, size_t size, const char* format, ...)
{
  int count;
  va_list ap;
  va_start(ap, format);
  count = c99_vsnprintf(str, size, format, ap);
  va_end(ap);
  return count;
}


} // namespace mirtk

#endif // MIRTK_Windows_H

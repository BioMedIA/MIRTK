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

#ifndef MIRTK_Memory_H
#define MIRTK_Memory_H


// =============================================================================
// C/C++ library functions
// =============================================================================

#include <cstring>
#include <utility>
#include <memory>

namespace mirtk {


template <class T>
using UniquePtr = std::unique_ptr<T>;

template <class T>
using SharedPtr = std::shared_ptr<T>;

template <class T>
using WeakPtr = std::weak_ptr<T>;

template <class T>
SharedPtr<T> NewShared()
{
  return std::make_shared<T>();
}

template <class T, class... Args>
SharedPtr<T> NewShared(Args&&... args)
{
  return std::make_shared<T>(args...);
}

using std::memset;
using std::memcpy;
using std::memmove;
using std::memcmp;
using std::swap;

/// Byte order of each word in memory
enum ByteOrder
{
  UnknownByteOrder,
  LittleEndian,
  BigEndian
};

/// Get byte order of this system
ByteOrder GetByteOrder();

/// Swap bytes of a single word
void swap16(char *, char *, long);

/// Swap bytes of two word
void swap32(char *, char *, long);

/// Swap bytes of four word
void swap64(char *, char *, long);

/// Returns the peak (maximum so far) resident set size (physical
/// memory use) measured in bytes, or zero if the value cannot be
/// determined on this OS.
size_t GetPeakRSS();

/// Returns the current resident set size (physical memory use) measured
/// in bytes, or zero if the value cannot be determined on this OS.
size_t GetCurrentRSS();


} // namespace mirtk

// =============================================================================
// Allocate/Deallocate N-D arrays, Delete
// =============================================================================

#include "mirtk/Allocate.h"
#include "mirtk/Deallocate.h"


#endif // MIRTK_Memory_H

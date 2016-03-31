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

#ifndef MIRTK_Memory_H
#define MIRTK_Memory_H


// =============================================================================
// C/C++ library functions
// =============================================================================

#include <cstring>
#include <utility>
#include <memory>

namespace mirtk {


using std::unique_ptr;
using std::shared_ptr;

using std::memset;
using std::memcpy;
using std::memmove;
using std::memcmp;
using std::swap;

void swap16(char *, char *, long);
void swap32(char *, char *, long);
void swap64(char *, char *, long);


} // namespace mirtk

// =============================================================================
// Allocate/Deallocate N-D arrays, Delete
// =============================================================================

#include "mirtk/Allocate.h"
#include "mirtk/Deallocate.h"


#endif // MIRTK_Memory_H

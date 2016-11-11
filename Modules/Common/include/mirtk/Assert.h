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

#ifndef MIRTK_Assert_H
#define MIRTK_Assert_H

#include "mirtk/Stream.h"

// -----------------------------------------------------------------------------
#ifndef NDEBUG
#  define mirtkAssert(condition, message)                                      \
     do {                                                                      \
       if (!(condition)) {                                                     \
         cerr << "Assertion `" #condition "` failed in " << __FILE__           \
                   << " line " << __LINE__ << ": " << message << endl;         \
         exit(EXIT_FAILURE);                                                   \
       }                                                                       \
     } while (false)
#else
#  define mirtkAssert(condition, message) do { } while (false)
#endif


#endif // MIRTK_Assert_H

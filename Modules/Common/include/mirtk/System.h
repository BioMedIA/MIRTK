/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#ifndef MIRTK_System_H
#define MIRTK_System_H

#include "mirtk/String.h"


namespace mirtk {


/// Get current date in the format "%d %b %Y"
/// @sa std::put_time
string GetDate();

/// Get current time in the format "%H:%M:%S %Z"
/// @sa std::put_time
string GetTime();

/// Get current date and time in the format "%c %Z"
/// @sa std::put_time
string GetDateTime();

/// Get name of user executing this program
string GetUser();


} // namespace mirtk

#endif // MIRTK_System_H

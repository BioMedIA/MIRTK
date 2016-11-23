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

#ifndef MIRTK_Stream_H
#define MIRTK_Stream_H

// Import common C/C++ library stream functions/types into mirtk namespace
// such that they can be used within this namespace without std:: prefix

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

namespace mirtk {

using std::ios;
using std::ios_base;
using std::getline;

// Input streams
using std::istream;
using std::ifstream;
using std::istringstream;

// Output streams
using std::ostream;
using std::ofstream;
using std::ostringstream;
using std::cin;
using std::cout;
using std::cerr;
using std::flush;

// Bidirectional stream
using std::iostream;
using std::stringstream;

// Output manipulators
using std::streamsize;
using std::setw;
using std::setprecision;
using std::left;
using std::right;
using std::fixed;
using std::scientific;
using std::endl;


} // namespace mirtk

#endif // MIRTK_Stream_H

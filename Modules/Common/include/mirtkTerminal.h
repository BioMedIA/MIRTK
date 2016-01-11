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

#ifndef MIRTK_Terminal_H
#define MIRTK_Terminal_H

#include <mirtkConfig.h> // WINDOWS
#include <mirtkStream.h>

#include <cstdio>
#if !WINDOWS
#  include <unistd.h>
#endif


namespace mirtk {


// =============================================================================
// Global terminal options
// =============================================================================

/// Whether to use color output to STDOUT
extern bool stdout_color;

// =============================================================================
// Command help
// =============================================================================

/// Check if given option is a terminal option
bool IsTerminalOption(const char *);

/// Parse terminal option
void ParseTerminalOption(int &, int &, char *[]);

/// Print terminal command-line options
void PrintTerminalOptions(ostream &);

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
inline bool StdOutIsRedirected()
{
#if WINDOWS
  return false;
#else
 return !isatty(fileno(stdout));
#endif
}

// =============================================================================
// Terminal colors
// =============================================================================

// -----------------------------------------------------------------------------
/// Color escape sequences
extern const char *xreset;

extern const char *xblack;
extern const char *xred;
extern const char *xgreen;
extern const char *xyellow;
extern const char *xblue;
extern const char *xmagenta;
extern const char *xcyan;
extern const char *xwhite;

extern const char *xbrightblack;
extern const char *xbrightred;
extern const char *xbrightgreen;
extern const char *xbrightyellow;
extern const char *xbrightblue;
extern const char *xbrightmagenta;
extern const char *xbrightcyan;
extern const char *xbrightwhite;

extern const char *xboldblack;
extern const char *xboldred;
extern const char *xboldgreen;
extern const char *xboldyellow;
extern const char *xboldblue;
extern const char *xboldmagenta;
extern const char *xboldcyan;
extern const char *xboldwhite;
extern const char *xboldbrightblack;
extern const char *xboldbrightred;
extern const char *xboldbrightgreen;
extern const char *xboldbrightyellow;
extern const char *xboldbrightblue;
extern const char *xboldbrightmagenta;
extern const char *xboldbrightcyan;
extern const char *xboldbrightwhite;


} // namespace mirtk

#endif // MIRTK_Terminal_H

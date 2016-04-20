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

#include "mirtk/CommonExport.h"

#include "mirtk/Config.h" // WINDOWS
#include "mirtk/Stream.h"

#include <cstdio>
#ifndef WINDOWS
#  include <unistd.h>
#endif


namespace mirtk {


// =============================================================================
// Global terminal options
// =============================================================================

/// Whether to use color output to STDOUT
MIRTK_Common_EXPORT extern bool stdout_color;

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
#ifdef WINDOWS
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
MIRTK_Common_EXPORT extern const char *xreset;

MIRTK_Common_EXPORT extern const char *xblack;
MIRTK_Common_EXPORT extern const char *xred;
MIRTK_Common_EXPORT extern const char *xgreen;
MIRTK_Common_EXPORT extern const char *xyellow;
MIRTK_Common_EXPORT extern const char *xblue;
MIRTK_Common_EXPORT extern const char *xmagenta;
MIRTK_Common_EXPORT extern const char *xcyan;
MIRTK_Common_EXPORT extern const char *xwhite;

MIRTK_Common_EXPORT extern const char *xbrightblack;
MIRTK_Common_EXPORT extern const char *xbrightred;
MIRTK_Common_EXPORT extern const char *xbrightgreen;
MIRTK_Common_EXPORT extern const char *xbrightyellow;
MIRTK_Common_EXPORT extern const char *xbrightblue;
MIRTK_Common_EXPORT extern const char *xbrightmagenta;
MIRTK_Common_EXPORT extern const char *xbrightcyan;
MIRTK_Common_EXPORT extern const char *xbrightwhite;

MIRTK_Common_EXPORT extern const char *xboldblack;
MIRTK_Common_EXPORT extern const char *xboldred;
MIRTK_Common_EXPORT extern const char *xboldgreen;
MIRTK_Common_EXPORT extern const char *xboldyellow;
MIRTK_Common_EXPORT extern const char *xboldblue;
MIRTK_Common_EXPORT extern const char *xboldmagenta;
MIRTK_Common_EXPORT extern const char *xboldcyan;
MIRTK_Common_EXPORT extern const char *xboldwhite;
MIRTK_Common_EXPORT extern const char *xboldbrightblack;
MIRTK_Common_EXPORT extern const char *xboldbrightred;
MIRTK_Common_EXPORT extern const char *xboldbrightgreen;
MIRTK_Common_EXPORT extern const char *xboldbrightyellow;
MIRTK_Common_EXPORT extern const char *xboldbrightblue;
MIRTK_Common_EXPORT extern const char *xboldbrightmagenta;
MIRTK_Common_EXPORT extern const char *xboldbrightcyan;
MIRTK_Common_EXPORT extern const char *xboldbrightwhite;


} // namespace mirtk

#endif // MIRTK_Terminal_H

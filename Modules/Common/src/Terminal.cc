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

#include "mirtk/Config.h" // WINDOWS
#include "mirtk/Terminal.h"
#include "mirtk/Options.h"
#include "mirtk/Stream.h"


namespace mirtk {


/// Whether to use color output to STDOUT
MIRTK_Common_EXPORT bool stdout_color = false;

// =============================================================================
// Command help
// =============================================================================

// -----------------------------------------------------------------------------
bool IsTerminalOption(const char *arg)
{
  _option = NULL;
  if      (strcmp(arg, "-color")   == 0) _option = "-color";
  else if (strcmp(arg, "-nocolor") == 0) _option = "-nocolor";
  return (_option != NULL);
}

// -----------------------------------------------------------------------------
void ParseTerminalOption(int &OPTIDX, int &argc, char *argv[])
{
  if (OPTION("-color")) {
    stdout_color = true;
  } else if (OPTION("-nocolor")) {
    stdout_color = false;
  }
}

// -----------------------------------------------------------------------------
void PrintTerminalOptions(ostream &out)
{
  out << endl;
  out << "Terminal options:" << endl;
  out << "  -[no]color                   Enable/disable colored output. (default: " << (stdout_color ? "on" : "off") << ")" << endl;
}

// =============================================================================
// Terminal colors
// =============================================================================

MIRTK_Common_EXPORT const char *xreset       = "\x1b[0m";
MIRTK_Common_EXPORT const char *xblack       = "\x1b[30m";
MIRTK_Common_EXPORT const char *xred         = "\x1b[31m";
MIRTK_Common_EXPORT const char *xgreen       = "\x1b[32m";
MIRTK_Common_EXPORT const char *xyellow      = "\x1b[33m";
MIRTK_Common_EXPORT const char *xblue        = "\x1b[34m";
MIRTK_Common_EXPORT const char *xmagenta     = "\x1b[35m";
MIRTK_Common_EXPORT const char *xcyan        = "\x1b[36m";
MIRTK_Common_EXPORT const char *xwhite       = "\x1b[37m";
MIRTK_Common_EXPORT const char *xboldblack   = "\x1b[30;1m";
MIRTK_Common_EXPORT const char *xboldred     = "\x1b[31;1m";
MIRTK_Common_EXPORT const char *xboldgreen   = "\x1b[32;1m";
MIRTK_Common_EXPORT const char *xboldyellow  = "\x1b[33;1m";
MIRTK_Common_EXPORT const char *xboldblue    = "\x1b[34;1m";
MIRTK_Common_EXPORT const char *xboldmagenta = "\x1b[35;1m";
MIRTK_Common_EXPORT const char *xboldcyan    = "\x1b[36;1m";
MIRTK_Common_EXPORT const char *xboldwhite   = "\x1b[37;1m";

MIRTK_Common_EXPORT const char *xbrightblack       = "\x1b[90m";
MIRTK_Common_EXPORT const char *xbrightred         = "\x1b[91m";
MIRTK_Common_EXPORT const char *xbrightgreen       = "\x1b[92m";
MIRTK_Common_EXPORT const char *xbrightyellow      = "\x1b[93m";
MIRTK_Common_EXPORT const char *xbrightblue        = "\x1b[94m";
MIRTK_Common_EXPORT const char *xbrightmagenta     = "\x1b[95m";
MIRTK_Common_EXPORT const char *xbrightcyan        = "\x1b[96m";
MIRTK_Common_EXPORT const char *xbrightwhite       = "\x1b[97m";
MIRTK_Common_EXPORT const char *xboldbrightblack   = "\x1b[90;1m";
MIRTK_Common_EXPORT const char *xboldbrightred     = "\x1b[91;1m";
MIRTK_Common_EXPORT const char *xboldbrightgreen   = "\x1b[92;1m";
MIRTK_Common_EXPORT const char *xboldbrightyellow  = "\x1b[93;1m";
MIRTK_Common_EXPORT const char *xboldbrightblue    = "\x1b[94;1m";
MIRTK_Common_EXPORT const char *xboldbrightmagenta = "\x1b[95;1m";
MIRTK_Common_EXPORT const char *xboldbrightcyan    = "\x1b[96;1m";
MIRTK_Common_EXPORT const char *xboldbrightwhite   = "\x1b[97;1m";

// =============================================================================
// Static initializer
// =============================================================================

namespace {

const char *COLOR_ENABLED_TERMINALS[] = {
  "xterm",
  "xterm-256color",
  NULL // List terminated by NULL pointer
};

/// Static initializer for stdout_color
struct SetDefaultColorMode
{
  SetDefaultColorMode()
  {
#ifdef WINDOWS
    char *term;
    size_t len;
    errno_t err = _dupenv_s(&term, &len, "TERM");
    if (err != 0) term = nullptr;
#else
    const char *term = getenv("TERM");
#endif
    if (term == nullptr) {
      stdout_color = false;
      return;
    }
    const char **color_term = COLOR_ENABLED_TERMINALS;
    while (*color_term) {
      if (strcmp(term, *color_term) == 0) break;
      ++color_term;
    }
#ifdef WINDOWS
    free(term);
#endif
    stdout_color = (*color_term && !StdOutIsRedirected());
  }
};
static SetDefaultColorMode color_mode_initializer;

}


} // namespace mirtk

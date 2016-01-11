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

#include <mirtkTerminal.h>
#include <mirtkOptions.h>
#include <mirtkStream.h>


namespace mirtk {


/// Whether to use color output to STDOUT
bool stdout_color = false;

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

const char *xreset       = "\x1b[0m";
const char *xblack       = "\x1b[30m";
const char *xred         = "\x1b[31m";
const char *xgreen       = "\x1b[32m";
const char *xyellow      = "\x1b[33m";
const char *xblue        = "\x1b[34m";
const char *xmagenta     = "\x1b[35m";
const char *xcyan        = "\x1b[36m";
const char *xwhite       = "\x1b[37m";
const char *xboldblack   = "\x1b[30;1m";
const char *xboldred     = "\x1b[31;1m";
const char *xboldgreen   = "\x1b[32;1m";
const char *xboldyellow  = "\x1b[33;1m";
const char *xboldblue    = "\x1b[34;1m";
const char *xboldmagenta = "\x1b[35;1m";
const char *xboldcyan    = "\x1b[36;1m";
const char *xboldwhite   = "\x1b[37;1m";

const char *xbrightblack       = "\x1b[90m";
const char *xbrightred         = "\x1b[91m";
const char *xbrightgreen       = "\x1b[92m";
const char *xbrightyellow      = "\x1b[93m";
const char *xbrightblue        = "\x1b[94m";
const char *xbrightmagenta     = "\x1b[95m";
const char *xbrightcyan        = "\x1b[96m";
const char *xbrightwhite       = "\x1b[97m";
const char *xboldbrightblack   = "\x1b[90;1m";
const char *xboldbrightred     = "\x1b[91;1m";
const char *xboldbrightgreen   = "\x1b[92;1m";
const char *xboldbrightyellow  = "\x1b[93;1m";
const char *xboldbrightblue    = "\x1b[94;1m";
const char *xboldbrightmagenta = "\x1b[95;1m";
const char *xboldbrightcyan    = "\x1b[96;1m";
const char *xboldbrightwhite   = "\x1b[97;1m";

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
    const char *term = getenv("TERM");
    if (!term) term = "dumb";
    const char **color_term = COLOR_ENABLED_TERMINALS;
    while (*color_term) {
      if (strcmp(term, *color_term) == 0) break;
      ++color_term;
    }
    stdout_color = (*color_term && !StdOutIsRedirected());
  }
};
static SetDefaultColorMode color_mode_initializer;

}


} // namespace mirtk

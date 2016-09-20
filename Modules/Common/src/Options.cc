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

#include "mirtk/Options.h"

#include "mirtk/String.h"
#include "mirtk/Stream.h"
#include "mirtk/Version.h"   // PrintVersion/PrintRevision
#include "mirtk/Terminal.h"  // PrintTerminalOptions
#include "mirtk/Parallel.h"  // PrintParallelOptions
#include "mirtk/Profiling.h" // PrintProfilingOptions

#include "mirtk/CommonExport.h"


namespace mirtk {


// =============================================================================
// Global standard options
// =============================================================================

// Default: No verbose output messages
MIRTK_Common_EXPORT int verbose = 0;

// Default: No debug output
MIRTK_Common_EXPORT int debug = 0;

// =============================================================================
// Help
// =============================================================================

MIRTK_Common_EXPORT int         _posargc                = 0;
MIRTK_Common_EXPORT int         _numposarg              = -1;
MIRTK_Common_EXPORT bool        _discard_parsed_posargs = false;
MIRTK_Common_EXPORT bool        _discard_parsed_options = false;
MIRTK_Common_EXPORT const char *_option                 = NULL;

// -----------------------------------------------------------------------------
bool IsStandardOption(const char *arg)
{
  _option = NULL;
  if      (strcmp(arg, "-v")        == 0) _option = "-v";
  else if (strcmp(arg, "-verbose")  == 0) _option = "-verbose";
  else if (strcmp(arg, "-debug")    == 0) _option = "-debug";
  else if (strcmp(arg, "-revision") == 0) _option = "-revision";
  else if (strcmp(arg, "-version")  == 0 || strcmp(arg, "--version") == 0) _option = "-version";
  return (_option != NULL);
}

// -----------------------------------------------------------------------------
void ParseStandardOption(int &OPTIDX, int &argc, char *argv[])
{
  if (OPTION("-v") || OPTION("-verbose")) {
    if (HAS_ARGUMENT) verbose  = atoi(ARGUMENT);
    else              verbose += 1;
  } else if (OPTION("-debug")) {
    if (HAS_ARGUMENT) debug  = atoi(ARGUMENT);
    else              debug += 1;
  } else if (OPTION("-version") || OPTION("--version")) {
    if (HAS_ARGUMENT) {
      if (!FromString(ARGUMENT, version) || version > current_version) {
        cerr << "Invalid [-]-version argument" << endl;
        exit(1);
      }
    } else {
      PrintVersion(cout, EXECNAME);
      exit(0);
    }
  } else if (OPTION("-revision")) {
    PrintRevision(cout);
    exit(0);
  }
}

// -----------------------------------------------------------------------------
void PrintStandardOptions(ostream &out)
{
  out << endl;
  out << "Standard options:" << endl;
  out << "  -v, -verbose [n]             Increase/Set verbosity of output messages. (default: " << verbose << ")" << endl;
  out << "  -debug [level]               Increase/Set debug level for output of intermediate results. (default: " << debug << ")" << endl;
  out << "  -[-]version [major.minor]    Print version and exit or set version to emulate." << endl;
  out << "  -revision                    Print revision (or version) number only and exit." << endl;
  out << "  -h, -[-]help                 Print help and exit." << endl;
}

// -----------------------------------------------------------------------------
void PrintCommonOptions(ostream &out)
{
  PrintStandardOptions (out);
  PrintTerminalOptions (out);
  PrintParallelOptions (out);
  PrintProfilingOptions(out);
}

// =============================================================================
// Private functions used by command-line parsing macros
// =============================================================================

// -----------------------------------------------------------------------------
inline bool _IsOptionName(const char *arg)
{
  return arg[0] == '-' && (arg[1] == '-' || !IsNumber(arg));
}

// -----------------------------------------------------------------------------
int _GetNumberOfPositionalArguments(int argc, char *argv[])
{
  int n = 1;
  while (n < argc && !_IsOptionName(argv[n])) n++;
  return n - 1;
}

// -----------------------------------------------------------------------------
void _DiscardArgument(int &i, int &argc, char *argv[])
{
  for (int j = i; j < argc; j++) argv[j] = argv[j+1];
  argv[argc--] = NULL;
  i--;
}

// -----------------------------------------------------------------------------
bool _IsOption(int &i, int &argc, char *argv[], const char *opt)
{
  if ((opt == NULL && _IsOptionName(argv[i])) || (opt != NULL && strcmp(argv[i], opt) == 0)) {
    if (_discard_parsed_options) _DiscardArgument(i, argc, argv);
    _option = argv[i];
    return true;
  } else {
    _option = NULL;
    return false;
  }
}

// -----------------------------------------------------------------------------
bool _IsArgument(int i, int &argc, char *argv[])
{
  return (i+1 < argc && !_IsOptionName(argv[i+1]));
}

// -----------------------------------------------------------------------------
char *_GetPositionalArgument(int i, int &argc, char *argv[])
{
  if (i >= argc) {
    cerr << "Error: Not all required arguments specified!" << endl;
    exit(1);
  }
  if (_posargc < i) _posargc = i;
  char *arg = argv[i];
  return arg;
}

// -----------------------------------------------------------------------------
char *_GetOptionArgument(int &i, int &argc, char *argv[])
{
  if (_option) {
    i++;
    if (i >= argc) {
      cerr << "Error: Option";
      if (_option) cerr << " " << _option;
      cerr << " requires more arguments than specified!" << endl;
      exit(1);
    }
  } else {
    if (i >= argc) {
      cerr << "Error: Not all required arguments specified!" << endl;
      exit(1);
    }
    if (_posargc < i) _posargc = i;
  }
  char *arg = argv[i];
  if (_discard_parsed_options) _DiscardArgument(i, argc, argv);
  return arg;
}


} // namespace mirtk

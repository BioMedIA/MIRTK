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

#include "mirtk/Profiling.h"

#include "mirtk/Options.h"
#include "mirtk/Stream.h"
#include "mirtk/Math.h"


namespace mirtk {


// =============================================================================
// Global profiling options
// =============================================================================

// Default: No profiling of execution time
int debug_time = 0;

// Default: Use seconds for CPU time and milliseconds for GPU time
TimeUnit debug_time_unit = TIME_IN_DEFAULT_UNIT;

// =============================================================================
// Command help
// =============================================================================

// -----------------------------------------------------------------------------
bool IsProfilingOption(const char *arg)
{
  _option = NULL;
  if      (strcmp(arg, "-profile")      == 0) _option = "-profile";
  else if (strcmp(arg, "-profile-unit") == 0) _option = "-profile-unit";
  return (_option != NULL);
}

// -----------------------------------------------------------------------------
void ParseProfilingOption(int &OPTIDX, int &argc, char *argv[])
{
  if (OPTION("-profile")) {
    if (HAS_ARGUMENT) debug_time  = atoi(ARGUMENT);
    else              debug_time += 1;
  } else if (OPTION("-profile-unit")) {
    const char *arg = ARGUMENT;
    if      (strcmp(arg, "msecs") == 0) debug_time_unit = TIME_IN_MILLISECONDS;
    else if (strcmp(arg, "secs")  == 0) debug_time_unit = TIME_IN_SECONDS;
    else {
      cerr << "Error: Invalid argument for option -profile-unit!" << endl;
      exit(1);
    }
  }
}

// -----------------------------------------------------------------------------
#ifdef USE_TIMING
void PrintProfilingOptions(ostream &out)
{
  out << endl;
  out << "Profiling options:" << endl;
  out << "  -profile [n]                 Increase/Set verbosity of time measurements. (default: 0)" << endl;
#ifdef USE_CUDA
  out << "  -profile-unit <msecs|secs>   Unit of time measurements. (default: secs [CPU], msecs [GPU])" << endl;
#else
  out << "  -profile-unit <msecs|secs>   Unit of time measurements. (default: secs)" << endl;
#endif
}
#else
void PrintProfilingOptions(ostream &) {}
#endif

// =============================================================================
// Profiling output
// =============================================================================

// -----------------------------------------------------------------------------
void PrintElapsedTime(const char *section, double t, TimeUnit unit)
{
  const unsigned int section_width = 64;
  char section_buffer[section_width+1];
  char unit_buffer   [8];
  snprintf(section_buffer, section_width+1, "%s:", section);
  if (debug_time_unit == TIME_IN_DEFAULT_UNIT) {
    if (unit == TIME_IN_MILLISECONDS) {
      if (t >= 1.e+3) {
        t   *= 1.e-3;
        unit = TIME_IN_SECONDS;
      }
    } else if (unit == TIME_IN_SECONDS) {
      if (t < 1.e-1) {
        t   *= 1.e+3;
        unit = TIME_IN_MILLISECONDS;
      }
    }
  } else {
    if (debug_time_unit == TIME_IN_MILLISECONDS) {
      if (unit == TIME_IN_SECONDS)      t *= 1.e+3;
    } else if (debug_time_unit == TIME_IN_SECONDS) {
      if (unit == TIME_IN_MILLISECONDS) t *= 1.e-3;
    }
    unit = debug_time_unit;
  }
  if      (unit == TIME_IN_MILLISECONDS) snprintf(unit_buffer, 6, "msec");
  else if (unit == TIME_IN_SECONDS)      snprintf(unit_buffer, 6, "sec");
  else                                   snprintf(unit_buffer, 6, "undef");
  printf("Time for %-*s %10.3f %s\n", section_width, section_buffer, t, unit_buffer);
  fflush(stdout);
}

// -----------------------------------------------------------------------------
string ElapsedTimeToString(double t, TimeUnit unit, TimeFormat fmt, int w, char c, bool left)
{
  string str;
  if (fmt == TIME_FORMAT_UNITS) {
    if (unit == TIME_IN_DEFAULT_UNIT) {
      if (debug_time_unit == TIME_IN_MILLISECONDS) {
        if (unit == TIME_IN_SECONDS)      t *= 1.e+3;
      } else if (debug_time_unit == TIME_IN_SECONDS) {
        if (unit == TIME_IN_MILLISECONDS) t *= 1.e-3;
      }
      unit = debug_time_unit;
    }
    str = ToString(t, w, c, left);
    if      (unit == TIME_IN_MILLISECONDS) str += " msec";
    else if (unit == TIME_IN_SECONDS)      str += " sec";
    else                                   str += " undef";
  } else {
    int h, m, s;
    if (unit == TIME_IN_DEFAULT_UNIT) unit = debug_time_unit;
    if (unit == TIME_IN_MILLISECONDS) t *= 1.e-3;
    s = iround(t);
    m = s / 60;
    s = s % 60;
    if (fmt == TIME_FORMAT_MIN_SEC) {
      str  = ToString(m, w, c, left);
      str += " min ";
      str += ToString(s, w, c, left);
      str += " sec";
    } else {
      h = m / 60;
      m = m % 60;
      if (fmt == TIME_FORMAT_H_MIN_SEC) {
        str  = ToString(h, w, c, left);
        str += " h ";
        str += ToString(m, w, c, left);
        str += " min ";
        str += ToString(s, w, c, left);
        str += " sec";
      } else {
        str  = ToString(h, 2, '0');
        str += ":";
        str += ToString(m, 2, '0');
        str += ":";
        str += ToString(s, 2, '0');
      }
    }
  }
  return str;
}


} // namespace mirtk

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

#ifndef MIRTK_Profiling_H
#define MIRTK_Profiling_H

#include "mirtk/CommonExport.h"

#include "mirtk/Stream.h"
#include "mirtk/String.h"

#include <ctime>
#ifdef HAVE_TBB
// TBB includes windows header which defines min/max macros otherwise
#  ifndef NOMINMAX
#    define NOMINMAX
#    define MIRTK_UNDEF_NOMINMAX
#  endif
#  include <tbb/tick_count.h>
#  ifdef MIRTK_UNDEF_NOMINMAX
#    undef MIRTK_UNDEF_NOMINMAX
#    undef NOMINMAX
#  endif
#endif


namespace mirtk {


// -----------------------------------------------------------------------------
enum TimeUnit
{
  TIME_IN_DEFAULT_UNIT,
  TIME_IN_MILLISECONDS,
  TIME_IN_SECONDS
};

enum TimeFormat
{
  TIME_FORMAT_UNITS,       ///< Print elapsed time using time units
  TIME_FORMAT_HHMMSS,      ///< Print elapsed time with format "HH:MM:SS"
  TIME_FORMAT_H_MIN_SEC,   ///< Print elapsed time with format "[H h] [M min] [S sec]"
  TIME_FORMAT_MIN_SEC      ///< Print elapsed time with format "[M min] [S sec]"
};

// =============================================================================
// Global profiling options
// =============================================================================

/// Enable/disable profiling of execution time.
///
/// Should be set in the main function before any processing starts, e.g.,
/// depending on a command-line flag (for example -v -verbose).
/// If less or equal to zero, no timing measurements are printed to screen.
/// Otherwise, whether a timing measure is output or not depends on the
/// set debugging level.
MIRTK_Common_EXPORT extern int debug_time;

/// Time unit to use for output of time measurements.
MIRTK_Common_EXPORT extern TimeUnit debug_time_unit;

// =============================================================================
// Command help
// =============================================================================

/// Check if given option is a profiling option
bool IsProfilingOption(const char *);

/// Parse profiling option
void ParseProfilingOption(int &, int &, char *[]);

/// Print profiling command-line options
void PrintProfilingOptions(ostream &);

// =============================================================================
// CPU Profiling
// =============================================================================

/// Print elapsed time for profiled section
void PrintElapsedTime(const char *, double, TimeUnit = TIME_IN_SECONDS);

/// Convert elapsed time given in the specified units to string of given format
///
/// \param[in] t      Elapsed time in specified \p units.
/// \param[in] units  Units of time measurement.
/// \param[in] fmt    Time string format.
/// \param[in] w      Width of time field when \p fmt is not TIME_FORMAT_HHMMSS.
/// \param[in] c      Character used to fill time field.
/// \param[in] left   Whether to print time value left justified.
string ElapsedTimeToString(double t, TimeUnit units = TIME_IN_SECONDS,
                           TimeFormat fmt = TIME_FORMAT_HHMMSS,
                           int w = 0, char c = ' ', bool left = false);

// -----------------------------------------------------------------------------
/// Start measurement of execution time of current code block
///
/// @code
/// {
///   MIRTK_START_TIMING();
///   // do some work here
///   MIRTK_END_TIMING("example section");
/// }
/// @endcode
///
/// @sa MIRTK_END_TIMING
#ifdef MIRTK_WITH_PROFILING
#  ifdef HAVE_TBB
#    define MIRTK_START_TIMING()   tbb::tick_count t_start = tbb::tick_count::now()
#  else
#    define MIRTK_START_TIMING()   clock_t t_start = clock()
#  endif
#else
#  define MIRTK_START_TIMING()   do {} while (false)
#endif

// -----------------------------------------------------------------------------
/// Reset measurement of starting execution time of current code block
///
/// @code
/// {
///   MIRTK_START_TIMING();
///   // do some work here
///   MIRTK_END_TIMING("first part");
///   MIRTK_RESET_TIMING();
///   // do some work here
///   MIRTK_END_TIMING("second part");
/// }
/// @endcode
///
/// @sa MIRTK_END_TIMING
#ifdef MIRTK_WITH_PROFILING
#  ifdef HAVE_TBB
#    define MIRTK_RESET_TIMING()   t_start = tbb::tick_count::now()
#  else
#    define MIRTK_RESET_TIMING()   t_start = clock()
#  endif
#else
#  define MIRTK_RESET_TIMING()   do {} while (false)
#endif

// -----------------------------------------------------------------------------
/// End measurement of execution time of current code block.
///
/// @code
/// {
///   MIRTK_START_TIMING();
///   // do some work here
///   MIRTK_END_TIMING("example section");
/// }
/// @endcode
///
/// @note Whether or not the execution time is actually being measured and
///       printed to screen is decided at compile time depending on the
///       MIRTK_WITH_PROFILING flag.
///
/// @sa MIRTK_START_TIMING
#ifdef MIRTK_WITH_PROFILING
#  ifdef HAVE_TBB
#    define MIRTK_END_TIMING(section)                                          \
       do {                                                                    \
         ostringstream oss;                                                    \
         oss << section;                                                       \
         PrintElapsedTime(oss.str().c_str(),                                   \
                          (tbb::tick_count::now() - t_start).seconds());       \
       } while (false)
#  else
#    define MIRTK_END_TIMING(section)                                          \
       do {                                                                    \
         ostringstream oss;                                                    \
         oss << section;                                                       \
         PrintElapsedTime(oss.str().c_str(),                                   \
                          static_cast<double>(clock() - t_start)               \
                        / static_cast<double>(CLOCKS_PER_SEC));                \
       } while (false)
#  endif
#else
#  define MIRTK_END_TIMING(section)   do {} while (false)
#endif

// -----------------------------------------------------------------------------
/// End measurement of execution time of current code block.
///
/// In the following example, the timing measurement are only printed if
/// the global execution time debugging level is greater or equal the
/// specified debugging level. Hence, if debug_time is 1, only the time of
/// the following example section is printed. If debug_time is greater than 1,
/// also the times needed for each execution of inner block are printed together
/// with the summarizing total execution time of the example section. Otherwise,
/// if debug_time is less than 1, no time measurements are printed at all.
/// @code
/// {
///   MIRTK_START_TIMING();
///   // do some work here
///   {
///     MIRTK_START_TIMING();
///     // do some part of the work here
///     MIRTK_DEBUG_TIMING(2, "example part");
///   }
///   // possibly combine results and maybe do some clean up work here
///   MIRTK_DEBUG_TIMING(1, "example section");
/// }
/// @endcode
///
/// @note Whether or not the execution time is actually being measured and
///       printed to screen is decided at runtime depending on the global
///       variable debug_time.
///
/// @sa MIRTK_START_TIMING
#ifdef MIRTK_WITH_PROFILING
#  ifdef HAVE_TBB
#    define MIRTK_DEBUG_TIMING(level, section)                                 \
       do {                                                                    \
         if (debug_time >= level) {                                            \
           ostringstream oss;                                                  \
           oss << section;                                                     \
           PrintElapsedTime(oss.str().c_str(),                                 \
                            (tbb::tick_count::now() - t_start).seconds());     \
         }                                                                     \
       } while (false)
#  else
#    define MIRTK_DEBUG_TIMING(level, section)                                 \
       do {                                                                    \
         if (debug_time >= level) {                                            \
           ostringstream oss;                                                  \
           oss << section;                                                     \
           PrintElapsedTime(oss.str().c_str(),                                 \
                            static_cast<double>(clock() - t_start)             \
                          / static_cast<double>(CLOCKS_PER_SEC));              \
         }                                                                     \
       } while (false)
#  endif
#else
#  define MIRTK_DEBUG_TIMING(level, section)   do {} while (false)
#endif

// =============================================================================
// GPU Profiling
// =============================================================================

// -----------------------------------------------------------------------------
/// Start measurement of execution time of current code block.
///
/// @code
/// {
///   MIRTKCU_START_TIMING();
///   // launch CUDA kernel here
///   MIRTKCU_END_TIMING("example section");
/// }
/// @endcode
///
/// @sa IRTKCU_END_TIMING
#ifdef MIRTK_WITH_PROFILING
#  define MIRTKCU_START_TIMING()                                               \
           cudaEvent_t e_start, e_stop;                                        \
           CudaSafeCall( cudaEventCreate(&e_start) );                          \
           CudaSafeCall( cudaEventCreate(&e_stop) );                           \
           CudaSafeCall( cudaEventRecord(e_start, 0) )
#else
#  define MIRTKCU_START_TIMING()   do {} while (false)
#endif

// -----------------------------------------------------------------------------
/// Reset measurement of starting execution time of current code block
///
/// @code
/// {
///   MIRTKCU_START_TIMING();
///   // do some work here
///   MIRTKCU_END_TIMING("first part");
///   MIRTKCU_RESET_TIMING();
///   // do some work here
///   MIRTKCU_END_TIMING("second part");
/// }
/// @endcode
///
/// @sa MIRTKCU_END_TIMING
#ifdef MIRTK_WITH_PROFILING
#  define MIRTKCU_RESET_TIMING()   CudaSafeCall( cudaEventRecord(e_start, 0) )
#else
#  define MIRTKCU_RESET_TIMING()   do {} while (false)
#endif

// -----------------------------------------------------------------------------
/// Print interim measurement of execution time of current code block.
///
/// @code
/// {
///   MIRTKCU_START_TIMING();
///   // launch CUDA kernel here
///   MIRTKCU_INTERIM_TIMING("example kernel");
///   CudaSafeCall( cudaDeviceSynchronize() );
///   MIRTKCU_END_TIMING("example section");
/// }
/// @endcode
///
/// @note Whether or not the execution time is actually being measured and
///       printed to screen is decided at compile time depending on the
///       MIRTK_WITH_PROFILING flag.
///
/// @sa MIRTKCU_START_TIMING
#ifdef MIRTK_WITH_PROFILING
#  define MIRTKCU_INTERIM_TIMING(section)                                      \
     do {                                                                      \
       float t_elapsed;                                                        \
       CudaSafeCall( cudaEventRecord(e_stop, 0) );                             \
       CudaSafeCall( cudaEventSynchronize(e_stop) );                           \
       CudaSafeCall( cudaEventElapsedTime(&t_elapsed, e_start, e_stop) );      \
       ostringstream oss;                                                      \
       oss << section << " [interim]";                                         \
       PrintElapsedTime(oss.str().c_str(), t_elapsed, TIME_IN_MILLISECONDS);   \
     } while (false)
#else
#  define MIRTKCU_INTERIM_TIMING(section)   do {} while (false)
#endif

// -----------------------------------------------------------------------------
/// End measurement of execution time of current code block.
///
/// @code
/// {
///   MIRTKCU_START_TIMING();
///   // launch CUDA kernel here
///   MIRTKCU_END_TIMING("example section");
/// }
/// @endcode
///
/// @note Whether or not the execution time is actually being measured and
///       printed to screen is decided at compile time depending on the
///       MIRTK_WITH_PROFILING flag.
///
/// @sa MIRTKCU_START_TIMING
#ifdef MIRTK_WITH_PROFILING
#  define MIRTKCU_END_TIMING(section)                                          \
     do {                                                                      \
       float t_elapsed;                                                        \
       CudaSafeCall( cudaEventRecord(e_stop, 0) );                             \
       CudaSafeCall( cudaEventSynchronize(e_stop) );                           \
       CudaSafeCall( cudaEventElapsedTime(&t_elapsed, e_start, e_stop) );      \
       CudaSafeCall( cudaEventDestroy(e_start) );                              \
       CudaSafeCall( cudaEventDestroy(e_stop) );                               \
       ostringstream oss;                                                      \
       oss << section << " [GPU]";                                             \
       PrintElapsedTime(oss.str().c_str(), t_elapsed, TIME_IN_MILLISECONDS);   \
     } while (false)
#else
#  define MIRTKCU_END_TIMING(section)   do {} while (false)
#endif

// -----------------------------------------------------------------------------
/// Print interim measurement of execution time of current code block.
///
/// In the following example, the timing measurement are only printed if
/// the global execution time debugging level is greater or equal the
/// specified debugging level. Hence, if debug_time is 1, only the time of
/// the following example section is printed. If debug_time is greater than 1,
/// also the times needed for each execution of inner block are printed together
/// with the summarizing total execution time of the example section. Otherwise,
/// if debug_time is less than 1, no time measurements are printed at all.
/// @code
/// {
///   MIRTKCU_START_TIMING();
///   // launch CUDA kernel here
///   MIRTKCU_DEBUG_INTERIM_TIMING(1, "kernel name");
/// }
/// @endcode
///
/// @note Whether or not the execution time is actually being measured and
///       printed to screen is decided at runtime depending on the global
///       variable debug_time.
///
/// @sa MIRTKCU_START_TIMING
#ifdef MIRTK_WITH_PROFILING
#  define MIRTKCU_DEBUG_INTERIM_TIMING(level, section)                         \
     if (debug_time >= level) IRTKCU_INTERIM_TIMING(section)
#else
#  define MIRTKCU_DEBUG_INTERIM_TIMING(level, section)   do {} while (false)
#endif

// ---------------------------------------------------------------------------
/// End measurement of execution time of current code block.
///
/// In the following example, the timing measurement are only printed if
/// the global execution time debugging level is greater or equal the
/// specified debugging level. Hence, if debug_time is 1, only the time of
/// the following example section is printed. If debug_time is greater than 1,
/// also the times needed for each execution of inner block are printed together
/// with the summarizing total execution time of the example section. Otherwise,
/// if debug_time is less than 1, no time measurements are printed at all.
/// @code
/// {
///   MIRTKCU_START_TIMING();
///   // launch CUDA kernel here
///   MIRTKCU_DEBUG_TIMING(1, "kernel name");
/// }
/// @endcode
///
/// @note Whether or not the execution time is actually being measured and
///       printed to screen is decided at runtime depending on the global
///       variable debug_time.
///
/// @sa MIRTKCU_START_TIMING
#ifdef MIRTK_WITH_PROFILING
#  define MIRTKCU_DEBUG_TIMING(level, section)                                 \
     if (debug_time >= level) MIRTKCU_END_TIMING(section);                     \
     else do {                                                                 \
       CudaSafeCall( cudaEventDestroy(e_start) );                              \
       CudaSafeCall( cudaEventDestroy(e_stop) );                               \
     } while (false)
#else
#  define MIRTKCU_DEBUG_TIMING(level, section)   do {} while (false)
#endif


} // namespace mirtk

#endif // MIRTK_Profiling_H

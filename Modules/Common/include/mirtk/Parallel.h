/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2017 Andreas Schuh
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

#ifndef MIRTK_Parallel_H
#define MIRTK_Parallel_H

#include "mirtk/CommonExport.h"

#include "mirtk/Stream.h"
#include "mirtk/Memory.h"

#ifndef MIRTK_COMMON_WITH_TBB_MALLOC
#  define MIRTK_COMMON_WITH_TBB_MALLOC 0
#endif

#ifdef HAVE_TBB
// TBB includes windows header which defines min/max macros otherwise
#  ifndef NOMINMAX
#    define NOMINMAX
#    define MIRTK_UNDEF_NOMINMAX
#  endif
#  include <tbb/task_scheduler_init.h>
#  include <tbb/blocked_range.h>
#  include <tbb/blocked_range2d.h>
#  include <tbb/blocked_range3d.h>
#  include <tbb/parallel_for.h>
#  include <tbb/parallel_reduce.h>
#  include <tbb/concurrent_queue.h>
#  if MIRTK_COMMON_WITH_TBB_MALLOC
#    include <tbb/scalable_allocator.h>
#    include <tbb/cache_aligned_allocator.h>
#  endif
#  include <tbb/mutex.h>
#  ifdef MIRTK_UNDEF_NOMINMAX
#    undef MIRTK_UNDEF_NOMINMAX
#    undef NOMINMAX
#  endif
#endif


namespace mirtk {


// =============================================================================
// Global parallelization options
// =============================================================================

/// Enable/disable GPU acceleration
MIRTK_Common_EXPORT extern bool use_gpu;

/// Debugging level of GPU code
MIRTK_Common_EXPORT extern int debug_gpu;

/// Debugging level of TBB code
MIRTK_Common_EXPORT extern int tbb_debug;

// =============================================================================
// Command help
// =============================================================================

/// Check if given option is a parallelization option
bool IsParallelOption(const char *);

/// Parse parallelization option
void ParseParallelOption(int &, int &, char *[]);

/// Print parallelization command-line options
void PrintParallelOptions(ostream &);

// =============================================================================
// Multi-threading support using Intel's TBB
// =============================================================================

// -----------------------------------------------------------------------------
// If TBB is available and WITH_TBB is set to ON, use TBB to execute
// any parallelizable code concurrently
//
// Attention: DO NOT define TBB_DEPRECATED by default or before including the
//            other TBB header files, in particular parallel_for. The deprecated
//            behavior of parallel_for is to not choose the chunk size (grainsize)
//            automatically!
//
// http://software.intel.com/sites/products/documentation/doclib/tbb_sa/help/tbb_userguide/Automatic_Chunking.htm
#ifdef HAVE_TBB


// Import used TBB types into mirtk namespace
using tbb::task_scheduler_init;
using tbb::blocked_range;
using tbb::blocked_range2d;
using tbb::blocked_range3d;
using tbb::parallel_for;
using tbb::parallel_reduce;
using tbb::concurrent_queue;
using tbb::mutex;
using tbb::split;

#if MIRTK_COMMON_WITH_TBB_MALLOC
using tbb::scalable_allocator;
using tbb::cache_aligned_allocator;
#endif

// A task scheduler is created/terminated automatically by TBB since
// version 2.2. It is recommended by Intel not to instantiate any task
// scheduler manually. However, in order to support the -threads option
// which can be used to limit the number of threads, a global task scheduler
// instance is created and the -threads argument passed on to its initialize
// method by ParseParallelOption. There should be no task scheduler created/
// terminated in any of the MIRTK library functions and classes.
MIRTK_Common_EXPORT extern UniquePtr<task_scheduler_init> tbb_scheduler;


// -----------------------------------------------------------------------------
// Otherwise, use dummy implementations of TBB classes/functions which allows
// developers to write parallelizable code as if TBB was available and yet
// executes the code serially due to the lack of TBB (or WITH_TBB set to OFF).
// This avoids code duplication and unnecessary conditional code compilation.
#else // HAVE_TBB

template <class T>
using scalable_allocator = std::allocator<T>;

template <class T>
using cache_aligned_allocator = std::allocator<T>;

/// Dummy type used to distinguish split constructor from copy constructor
struct split {};

/// Helper for initialization of task scheduler
class task_scheduler_init
{
public:
  task_scheduler_init(int) {}
  void terminate() {}
};

/// One-dimensional range
template <typename T>
class blocked_range
{
  T _lbound;
  T _ubound;
public:
  blocked_range(T l, T u)         : _lbound(l), _ubound(u) {}
  blocked_range(T l, T u, size_t) : _lbound(l), _ubound(u) {}
  T begin() const { return _lbound; }
  T end()   const { return _ubound; }
};

/// Two-dimensional range
template <typename T>
class blocked_range2d
{
  blocked_range<T> _rows;
  blocked_range<T> _cols;

public:

  blocked_range2d(T rl, T ru,
                  T cl, T cu)
  :
    _rows (rl, ru),
    _cols (cl, cu)
  {
  }

  blocked_range2d(T rl, T ru, size_t,
                  T cl, T cu, size_t)
  :
    _rows (rl, ru),
    _cols (cl, cu)
  {
  }

  const blocked_range<T> &rows() const { return _rows; }
  const blocked_range<T> &cols() const { return _cols; }
};

/// Three-dimensional range
template <typename T>
class blocked_range3d
{
  blocked_range<T> _pages;
  blocked_range<T> _rows;
  blocked_range<T> _cols;

public:

  blocked_range3d(T pl, T pu,
                  T rl, T ru,
                  T cl, T cu)
  :
    _pages(pl, pu),
    _rows (rl, ru),
    _cols (cl, cu)
  {
  }

  blocked_range3d(T pl, T pu, size_t,
                  T rl, T ru, size_t,
                  T cl, T cu, size_t)
  :
    _pages(pl, pu),
    _rows (rl, ru),
    _cols (cl, cu)
  {
  }

  const blocked_range<T> &pages() const { return _pages; }
  const blocked_range<T> &rows() const { return _rows; }
  const blocked_range<T> &cols() const { return _cols; }
};

/// parallel_for dummy template function which executes the body serially
template <class Range, class Body>
void parallel_for(const Range &range, const Body &body) {
  body(range);
}

/// parallel_reduce dummy template function which executes the body serially
template <class Range, class Body>
void parallel_reduce(const Range &range, Body &body) {
  body(range);
}


#endif // HAVE_TBB


} // namespace mirtk

#endif // MIRTK_Parallel_H

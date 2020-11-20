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

#ifdef HAVE_TBB
#  include <tbb/tbb_stddef.h>
#endif

#include "mirtk/Parallel.h"

#include "mirtk/String.h"
#include "mirtk/Stream.h"
#include "mirtk/Options.h"

namespace mirtk {


// =============================================================================
// Global parallelization options
// =============================================================================

// Default: Disable GPU acceleration
bool use_gpu = false;

// Default: No debugging of GPU code
int debug_gpu = 0;

// Default: No debugging of TBB code
int tbb_debug = 0;

#ifdef HAVE_TBB
UniquePtr<global_control> tbb_global_control;
#endif

// =============================================================================
// Command help
// =============================================================================

// -----------------------------------------------------------------------------
bool IsParallelOption(const char *arg)
{
  _option = NULL;
  if      (strcmp(arg, "-cpu")     == 0) _option = "-cpu";
  else if (strcmp(arg, "-gpu")     == 0) _option = "-gpu";
  else if (strcmp(arg, "-threads") == 0) _option = "-threads";
  return (_option != NULL);
}

// -----------------------------------------------------------------------------
void ParseParallelOption(int &OPTIDX, int &argc, char *argv[])
{
  if (OPTION("-cpu")) {
    use_gpu = false;
  } else if (OPTION("-gpu")) {
#ifdef USE_CUDA
    use_gpu = true;
#else    
    cerr << "WARNING: Program compiled without GPU support using CUDA." << endl;
    use_gpu = false;
#endif
  } else if (OPTION("-threads")) {
    int no_threads = 0; // parse argument even if unused
    if (!FromString(ARGUMENT, no_threads)) {
      cerr << "Invalid -threads argument, must be an integer!" << endl;
      exit(1);
    }
#ifdef HAVE_TBB
    if (no_threads < 0) {
      tbb_global_control.reset(nullptr);
    } else {
      // Before TBB 4.4 Update 3, value of 0 and 1 not allowed for max_allowed_parallelism.
      #if TBB_INTERFACE_VERSION < 9003
        if (no_threads == 0) no_threads = 1;
      #endif
      // +1 for master thread, no_threads workers; -threads 0 --> serial execution
      tbb_global_control.reset(new global_control(global_control::max_allowed_parallelism, no_threads + 1));
    }
#endif
  }
}

// -----------------------------------------------------------------------------
#if defined(HAVE_TBB) || defined(USE_CUDA)
void PrintParallelOptions(ostream &out)
{
  out << endl;
  out << "Parallelization options:" << endl;
#ifdef HAVE_TBB
  out << "  -threads <n>                 Use maximal <n> threads for parallel execution. (default: automatic)" << endl;
#endif
#ifdef USE_CUDA
  out << "  -gpu                         Enable  GPU acceleration."   << (use_gpu ? " (default)" : "") << endl;
  out << "  -cpu                         Disable GPU acceleration."   << (use_gpu ? "" : " (default)") << endl;
#endif
}
#else
void PrintParallelOptions(ostream &) {}
#endif // HAVE_TBB || USE_CUDA


} // namespace mirtk

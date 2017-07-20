/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016-2017 Imperial College London
 * Copyright 2016-2017 Andreas Schuh
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

#ifndef MIRTK_Exception_H
#define MIRTK_Exception_H

#include "mirtk/String.h"
#include "mirtk/Stream.h"


namespace mirtk {


/// Enumeration of error types / exit codes
enum ErrorType
{
  ERR_None            = 0,
  ERR_Unknown         = 1,
  ERR_LogicError      = 2,
  ERR_InvalidArgument = 3,
  ERR_RuntimeError    = 4,
  ERR_IOError         = 5,
  ERR_Memory          = 6,
  ERR_NotImplemented  = 42
};

/// Convert error type to string
template <>
inline string ToString(const ErrorType &value, int w, char c, bool left)
{
  const char *str;
  switch (value) {
    case ERR_None:            { str = "NoError";         } break;
    case ERR_Unknown:         { str = "UnknownError";    } break;
    case ERR_LogicError:      { str = "LogicError";      } break;
    case ERR_InvalidArgument: { str = "InvalidArgument"; } break;
    case ERR_RuntimeError:    { str = "RuntimeError";    } break;
    case ERR_IOError:         { str = "IOError";         } break;
    case ERR_Memory:          { str = "MemoryError";     } break;
    case ERR_NotImplemented:  { str = "NotImplemented";  } break;
  }
  return ToString(str, w, c, left);
}

/// Raise error in function
///
/// The current implementation prints the error message to STDERR and terminates
/// the program with exit code 1. In future releases, when all library code has
/// been rewritten to use this function, a suitable runtime exception may be
/// thrown instead.
///
/// \param[in] err  Type of error/exception.
/// \param[in] func Name of function which is throwing the error (i.e., __func__).
/// \param[in] args Error message. The given arguments are converted to strings
///                 using the ToString template function. These strings are then
///                 concatenated to produce the complete error message.
template <typename... Args>
void Throw(ErrorType err, const char *func, Args... args)
{
  Print(cerr, err, ": ", func, ": ", args...);
  cerr << endl;
  exit(static_cast<int>(err));
}


} // namespace mirtk

#endif // MIRTK_Exception_H

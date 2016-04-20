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

/**
 * \file  mirtk/Umfpack.h
 * \brief Interface to UMFPACK routines.
 *
 * \attention Include this header in internal files such as .cc translation units only!
 */

#ifndef MIRTK_Umfpack_H
#define MIRTK_Umfpack_H

#include <umfpack.h>


namespace mirtk {


// -----------------------------------------------------------------------------
/// Translate UMFPACK status code to message string
const char *umfpack_status_message(int status)
{
  switch (status) {
    case UMFPACK_OK:                            return "Success";
#ifdef UMFPACK_WARNING_singular_matrix
    case UMFPACK_WARNING_singular_matrix:       return "Matrix is singular";
#endif
#ifdef UMFPACK_WARNING_determinant_underflow
    case UMFPACK_WARNING_determinant_underflow: return "Determinant underflow";
#endif
#ifdef UMFPACK_WARNING_determinant_overflow
    case UMFPACK_WARNING_determinant_overflow:  return "Determinant overflow";
#endif
#ifdef UMFPACK_ERROR_out_of_memory
    case UMFPACK_ERROR_out_of_memory:           return "Out of memory";
#endif
#ifdef UMFPACK_ERROR_invalid_Numeric_object
    case UMFPACK_ERROR_invalid_Numeric_object:  return "Invalid Numeric object";
#endif
#ifdef UMFPACK_ERROR_invalid_Symbolic_object
    case UMFPACK_ERROR_invalid_Symbolic_object: return "Invalid Symbolic object";
#endif
#ifdef UMFPACK_ERROR_argument_missing
    case UMFPACK_ERROR_argument_missing:        return "Argument missing";
#endif
#ifdef UMFPACK_ERROR_n_nonpositive
    case UMFPACK_ERROR_n_nonpositive:           return "N non-positive";
#endif
#ifdef UMFPACK_ERROR_invalid_matrix
    case UMFPACK_ERROR_invalid_matrix:          return "Invalid matrix";
#endif
#ifdef UMFPACK_ERROR_different_pattern
    case UMFPACK_ERROR_different_pattern:       return "Different pattern";
#endif
#ifdef UMFPACK_ERROR_invalid_system
    case UMFPACK_ERROR_invalid_system:          return "Invalid system";
#endif
#ifdef UMFPACK_ERROR_invalid_permutation
    case UMFPACK_ERROR_invalid_permutation:     return "Invalid permutation";
#endif
#ifdef UMFPACK_ERROR_internal_error
    case UMFPACK_ERROR_internal_error:          return "Internal error";
#endif
#ifdef UMFPACK_ERROR_file_IO
    case UMFPACK_ERROR_file_IO:                 return "File I/O error";
#endif
#ifdef UMFPACK_ERROR_ordering_failed
    case UMFPACK_ERROR_ordering_failed:         return "Ordering failed";
#endif
  }
  return "Unknown status code";
}


} // namespace mirtk

#endif // MIRTK_Umfpack_H_

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

#ifndef MIRTK_Deallocate_H
#define MIRTK_Deallocate_H

namespace mirtk {


/// Delete object
template <typename Type>
inline void Delete(Type *&p)
{
  delete p;
  p = NULL;
}

/// Deallocate 1D array
template <typename Type>
inline void Deallocate(Type *&p)
{
  delete[] p;
  p = NULL;
}

/// Deallocate 2D array stored in contiguous memory block
///
/// \param[in] matrix Previously allocated array or \c NULL.
/// \param[in] data   Contiguous memory used by this array, but managed
///                   separately. If not \c NULL, only the pointers in the
///                   array are deallocated, but not the \c data memory itself.
///                   Otherwise, also the contiguous data memory block is freed.
template <class Type>
inline void Deallocate(Type **&matrix, void *data = NULL)
{
  if (matrix) {
    if (matrix[0] != data) {
      delete[] matrix[0];
    }
    delete[] matrix;
    matrix = NULL;
  }
}

/// Deallocate 3D array stored in contiguous memory block
///
/// \param[in] matrix Previously allocated array or \c NULL.
/// \param[in] data   Contiguous memory used by this array, but managed
///                   separately. If not \c NULL, only the pointers in the
///                   array are deallocated, but not the \c data memory itself.
///                   Otherwise, also the contiguous data memory block is freed.
template <class Type>
inline void Deallocate(Type ***&matrix, void *data = NULL)
{
  if (matrix) {
    if (matrix[0][0] != data) {
      delete[] matrix[0][0];
    }
    delete[] matrix[0];
    delete[] matrix;
    matrix = NULL;
  }
}

/// Deallocate 4D array stored in contiguous memory block
///
/// \param[in] matrix Previously allocated array or \c NULL.
/// \param[in] data   Contiguous memory used by this array, but managed
///                   separately. If not \c NULL, only the pointers in the
///                   array are deallocated, but not the \c data memory itself.
///                   Otherwise, also the contiguous data memory block is freed.
template <class Type>
inline void Deallocate(Type ****&matrix, void *data = NULL)
{
  if (matrix) {
    if (matrix[0][0][0] != data) {
      delete[] matrix[0][0][0];
    }
    delete[] matrix[0][0];
    delete[] matrix[0];
    delete[] matrix;
    matrix = NULL;
  }
}


} // namespace mirtk

#endif // MIRTK_Deallocate_H

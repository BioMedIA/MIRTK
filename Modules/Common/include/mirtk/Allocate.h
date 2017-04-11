/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2017 Imperial College London
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

#ifndef MIRTK_Allocate_H
#define MIRTK_Allocate_H

#include "mirtk/Exception.h"


namespace mirtk {


// Usage:
// \code
// // Allocate an array of 10 pointers to a grey image and initialize them to 0
// GreyImage **images1;
// PAllocate(images1, 10);
// GreyImage **images2 = PAllocate<irtkGreyImage>(10);
// // Allocate an array of 5 integer values without initialization
// int *values1;
// Allocate(values1, 6);
// // Allocate a multi-dimensional array using existing 1D memory
// int **values2 = Allocate<int>(2, 3, values1);
// // Free memory again
// Deallocate(images1);
// Deallocate(images2);
// // Deallocate multi-dimensional array either as follows
// Deallocate(values2, values1);
// Deallocate(values1);
// // or only the multi-dimensional array, but no separate deallocation of the
// // 1D memory which is as well freed already by the following Deallocate call
// // Deallocate(values2);
// \endcode


using std::nothrow;


// =============================================================================
// 1D array
// =============================================================================

// -----------------------------------------------------------------------------
/// Allocate 1D array
template <class Type>
inline void Allocate(Type *&matrix, int n)
{
  // Set pointer to nullptr if memory size is not positive
  if (n <= 0) {
    matrix = nullptr;
    return;
  }
  // Allocate data memory
  if ((matrix = new (nothrow) Type[n]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", n * sizeof(Type), " bytes");
  }
}

// -----------------------------------------------------------------------------
/// Allocate 1D array
template <class Type>
inline Type *Allocate(int n)
{
  Type *matrix;
  Allocate(matrix, n);
  return matrix;
}

// -----------------------------------------------------------------------------
/// Allocate 1D array and initialize it
template <class Type>
inline void CAllocate(Type *&matrix, int n, const Type &init = Type())
{
  // Allocate data memory
  Allocate(matrix, n);
  // Initialize data memory
  if (matrix) {
    for (int i = 0; i < n; i++) matrix[i] = init;
  }
}

// -----------------------------------------------------------------------------
/// Allocate 1D array and initialize it
template <class Type>
inline Type *CAllocate(int n, const Type *init = nullptr)
{
  Type *matrix;
  const Type default_value = Type();
  if (!init) init = &default_value;
  CAllocate(matrix, n, *init);
  return matrix;
}

// -----------------------------------------------------------------------------
/// Allocate 1D array of pointers initialized to nullptr
template <typename Type>
inline void PAllocate(Type **&matrix, int n)
{
  // Allocate data memory
  Allocate(matrix, n);
  // Initialize data memory
  if (matrix) memset(matrix, 0, n * sizeof(Type *));
}

// -----------------------------------------------------------------------------
/// Allocate 1D array of pointers initialized to nullptr
template <typename Type>
inline Type **PAllocate(int n)
{
  Type **matrix;
  PAllocate(matrix, n);
  return matrix;
}

// =============================================================================
// 2D array
// =============================================================================

// -----------------------------------------------------------------------------
/// Allocate 2D array stored in contiguous memory block
template <class Type>
inline void CAllocate(Type **&matrix, int x, int y, const Type &init = Type())
{
  const size_t n = static_cast<size_t>(x) * static_cast<size_t>(y);
  // Set pointer to nullptr if memory size is not positive
  if (n <= 0u) {
    matrix = nullptr;
    return;
  }
  // Allocate pointers
  if ((matrix = new (nothrow) Type *[y]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", y * sizeof(Type *), " bytes");
  }
  // Allocate data memory
  if ((matrix[0] = new (nothrow) Type[n]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", n * sizeof(Type), " bytes");
  }
  // Initialize data memory
  for (size_t i = 0; i < n; ++i) matrix[0][i] = init;
  // Initialize pointers
  for (int j = 1; j < y; ++j) {
    matrix[j] = matrix[j-1] + x;
  }
}

// -----------------------------------------------------------------------------
/// Allocate 2D array stored in contiguous memory block and initialize it
template <typename Type>
inline Type **CAllocate(int x, int y, const Type *init = nullptr)
{
  Type **matrix;
  const Type default_value = Type();
  if (!init) init = &default_value;
  CAllocate(matrix, x, y, *init);
  return matrix;
}

// -----------------------------------------------------------------------------
/// Allocate 2D array stored in contiguous memory block
template <class Type>
inline void Allocate(Type **&matrix, int x, int y, Type *data = nullptr)
{
  const size_t n = static_cast<size_t>(x) * static_cast<size_t>(y);
  // Set pointer to nullptr if memory size is not positive
  if (n <= 0u) {
    matrix = nullptr;
    return;
  }
  // Allocate pointers
  if ((matrix = new (nothrow) Type *[y]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", y * sizeof(Type *), " bytes");
  }
  // Allocate data memory
  if (data) matrix[0] = data;
  else if ((matrix[0] = new (nothrow) Type[n]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", n * sizeof(Type), " bytes");
  }
  // Initialize pointers
  for (int j = 1; j < y; ++j) {
    matrix[j] = matrix[j-1] + x;
  }
}

// -----------------------------------------------------------------------------
/// Allocate 2D array stored in contiguous memory block
template <typename Type>
inline Type **Allocate(int x, int y, Type *data = nullptr)
{
  Type **matrix;
  Allocate(matrix, x, y, data);
  return matrix;
}

// -----------------------------------------------------------------------------
/// Reshape 2D array stored in contiguous memory block
template <class Type>
inline Type **Reshape(Type **matrix, int x, int y)
{
  // Keep data memory
  Type * const data = matrix[0];
  // Deallocate old pointers
  delete[] matrix;
  // Allocate new pointers
  if ((matrix = new (nothrow) Type *[y]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", y * sizeof(Type *), " bytes");
  }
  // Restore data memory
  matrix[0] = data;
  // Initialize new pointers
  for (int j = 1; j < y; ++j) {
    matrix[j] = matrix[j-1] + x;
  }
  return matrix;
}

// =============================================================================
// 3D array
// =============================================================================

// -----------------------------------------------------------------------------
/// Allocate 3D array stored in contiguous memory block
template <class Type>
inline void CAllocate(Type ***&matrix, int x, int y, int z, const Type &init = Type())
{
  const size_t n = static_cast<size_t>(x)
                 * static_cast<size_t>(y)
                 * static_cast<size_t>(z);
  // Set pointer to nullptr if memory size is not positive
  if (n <= 0u) {
    matrix = nullptr;
    return;
  }
  // Allocate pointers
  if ((matrix = new (nothrow) Type **[z]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", z * sizeof(Type **), " bytes");
  }
  if ((matrix[0] = new (nothrow) Type *[z*y]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", z * y * sizeof(Type *), " bytes");
  }
  // Allocate data memory
  if ((matrix[0][0] = new (nothrow) Type[n]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", n * sizeof(Type), " bytes");
  }
  // Initialize data memory
  for (size_t i = 0; i < n; ++i) matrix[0][0][i] = init;
  // Initialize pointers
  for (int k = 0; k < z; ++k) {
    if (k > 0) matrix[k] = matrix[k-1] + y;
    for (int j = 0; j < y; ++j) {
      matrix[k][j] = matrix[0][0] + (k*y + j) * x;
    }
  }
}

// -----------------------------------------------------------------------------
/// Allocate 3D array stored in contiguous memory block and initialize it
template <typename Type>
inline Type ***CAllocate(int x, int y, int z, const Type *init = nullptr)
{
  Type ***matrix;
  const Type default_value = Type();
  if (!init) init = &default_value;
  CAllocate(matrix, x, y, z, *init);
  return matrix;
}

// -----------------------------------------------------------------------------
/// Allocate 3D array stored in contiguous memory block
template <class Type>
inline void Allocate(Type ***&matrix, int x, int y, int z, Type *data = nullptr)
{
  const size_t n = static_cast<size_t>(x)
                 * static_cast<size_t>(y)
                 * static_cast<size_t>(z);
  // Set pointer to nullptr if memory size is not positive
  if (n <= 0u) {
    matrix = nullptr;
    return;
  }
  // Allocate pointers
  if ((matrix = new (nothrow) Type **[z]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", z * sizeof(Type **), " bytes");
  }
  if ((matrix[0] = new (nothrow) Type *[z*y]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", z * y * sizeof(Type *), " bytes");
  }
  // Allocate data memory
  if (data) matrix[0][0] = data;
  else if ((matrix[0][0] = new (nothrow) Type[n]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", n * sizeof(Type), " bytes");
  }
  // Initialize pointers
  for (int k = 0; k < z; ++k) {
    if (k > 0) matrix[k] = matrix[k-1] + y;
    for (int j = 0; j < y; ++j) {
      matrix[k][j] = matrix[0][0] + (k * y + j) * x;
    }
  }
}

// -----------------------------------------------------------------------------
/// Allocate 3D array stored in contiguous memory block
template <typename Type>
inline Type ***Allocate(int x, int y, int z, Type *data = nullptr)
{
  Type ***matrix;
  Allocate(matrix, x, y, z, data);
  return matrix;
}

// -----------------------------------------------------------------------------
/// Reshape 3D array stored in contiguous memory block
template <class Type>
inline Type ***Reshape(Type ***matrix, int x, int y, int z)
{
  // Keep data memory
  Type * const data = matrix[0][0];
  // Deallocate old pointers
  delete[] matrix[0];
  delete[] matrix;
  // Allocate new pointers
  if ((matrix = new (nothrow) Type **[z]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", z * sizeof(Type **), " bytes");
  }
  if ((matrix[0] = new (nothrow) Type *[z*y]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", z * y * sizeof(Type *), " bytes");
  }
  // Restore data memory
  matrix[0][0] = data;
  // Initialize new pointers
  for (int k = 0; k < z; ++k) {
    if (k > 0) matrix[k] = matrix[k-1] + y;
    for (int j = 0; j < y; ++j) {
      matrix[k][j] = matrix[0][0] + (k * y + j) * x;
    }
  }
  return matrix;
}

// =============================================================================
// 4D array
// =============================================================================

// -----------------------------------------------------------------------------
/// Allocate 4D array stored in contiguous memory block and initialize it
template <class Type>
inline void CAllocate(Type ****&matrix, int x, int y, int z, int t, const Type &init = Type())
{
  const size_t n = static_cast<size_t>(x)
                 * static_cast<size_t>(y)
                 * static_cast<size_t>(z)
                 * static_cast<size_t>(t);
  // Set pointer to nullptr if memory size is not positive
  if (n <= 0u) {
    matrix = nullptr;
    return;
  }
  // Allocate pointers
  if ((matrix = new (nothrow) Type ***[t]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", t * sizeof(Type ***), " bytes");
  }
  if ((matrix[0] = new (nothrow) Type **[t*z]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", t * z * sizeof(Type **), " bytes");
  }
  if ((matrix[0][0] = new (nothrow) Type *[t*z*y]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", t * z * y * sizeof(Type *), " bytes");
  }
  // Allocate data memory
  if ((matrix[0][0][0] = new (nothrow) Type[n]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", n * sizeof(Type), " bytes");
  }
  // Initialize data memory
  for (size_t i = 0; i < n; ++i) matrix[0][0][0][i] = init;
  // Initialize pointers
  for (int l = 0; l < t; ++l) {
    if (l > 0) matrix[l] = matrix[l-1] + z;
    for (int k = 0; k < z; ++k) {
      matrix[l][k] = matrix[0][0] + (l * z + k) * y;
      for (int j = 0; j < y; ++j) {
        matrix[l][k][j] = matrix[0][0][0] + ((l * z + k) * y + j) * x;
      }
    }
  }
}

// -----------------------------------------------------------------------------
/// Allocate 4D array stored in contiguous memory block and initialize it
template <typename Type>
inline Type ****CAllocate(int x, int y, int z, int t, const Type *init = nullptr)
{
  Type ****matrix;
  const Type default_value = Type();
  if (!init) init = &default_value;
  CAllocate(matrix, x, y, z, t, *init);
  return matrix;
}

// -----------------------------------------------------------------------------
/// Allocate 4D array stored in contiguous memory block
template <class Type>
inline void Allocate(Type ****&matrix, int x, int y, int z, int t, Type *data = nullptr)
{
  const size_t n = static_cast<size_t>(x)
                 * static_cast<size_t>(y)
                 * static_cast<size_t>(z)
                 * static_cast<size_t>(t);
  // Set pointer to nullptr if memory size is not positive
  if (n <= 0u) {
    matrix = nullptr;
    return;
  }
  // Allocate pointers
  if ((matrix = new (nothrow) Type ***[t]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", t * sizeof(Type ***), " bytes");
  }
  if ((matrix[0] = new (nothrow) Type ** [t*z]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", t * z * sizeof(Type **), " bytes");
  }
  if ((matrix[0][0] = new (nothrow) Type * [t*z*y]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", t * z * y * sizeof(Type *), " bytes");
  }
  // Allocate data memory
  if (data) matrix[0][0][0] = data;
  else if ((matrix[0][0][0] = new (nothrow) Type[n]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", n * sizeof(Type), " bytes");
  }
  // Initialize pointers
  for (int l = 0; l < t; ++l) {
    if (l > 0) matrix[l] = matrix[l-1] + z;
    for (int k = 0; k < z; ++k) {
      matrix[l][k] = matrix[0][0] + (l * z + k) * y;
      for (int j = 0; j < y; ++j) {
        matrix[l][k][j] = matrix[0][0][0] + ((l * z + k) * y + j) * x;
      }
    }
  }
}

// -----------------------------------------------------------------------------
/// Allocate 4D array stored in contiguous memory block
template <typename Type>
inline Type ****Allocate(int x, int y, int z, int t, Type *data = nullptr)
{
  Type ****matrix;
  Allocate(matrix, x, y, z, t, data);
  return matrix;
}

// -----------------------------------------------------------------------------
/// Reshape 4D array stored in contiguous memory block
template <class Type>
inline Type ****Reshape(Type ****matrix, int x, int y, int z, int t)
{
  // Keep data memory
  Type * const data = matrix[0][0][0];
  // Deallocate old pointers
  delete[] matrix[0][0];
  delete[] matrix[0];
  delete[] matrix;
  // Allocate new pointers
  if ((matrix = new (nothrow) Type ***[t]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", t * sizeof(Type ***), " bytes");
  }
  if ((matrix[0] = new (nothrow) Type **[t*z]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", t * z * sizeof(Type **), " bytes");
  }
  if ((matrix[0][0] = new (nothrow) Type *[t*z*y]) == nullptr) {
    Throw(ERR_Memory, __FUNCTION__, "Failed to allocate ", t * z * y * sizeof(Type *), " bytes");
  }
  // Restore data memory
  matrix[0][0][0] = data;
  // Initialize new pointers
  for (int l = 0; l < t; ++l) {
    if (l > 0) matrix[l] = matrix[l-1] + z;
    for (int k = 0; k < z; ++k) {
      matrix[l][k] = matrix[0][0] + (l * z + k) * y;
      for (int j = 0; j < y; ++j) {
        matrix[l][k][j] = matrix[0][0][0] + ((l * z + k) * y + j) * x;
      }
    }
  }
  return matrix;
}


} // namespace mirtk

#endif // MIRTK_Allocate_H

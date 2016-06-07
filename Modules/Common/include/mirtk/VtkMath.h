/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#ifndef MIRTK_VtkMath_H
#define MIRTK_VtkMath_H

#include "vtkConfigure.h" // VTK version macros 


// See http://www.paraview.org/Bug/view.php?id=14164
#if VTK_MAJOR_VERSION == 6 && VTK_MINOR_VERSION == 0 && VTK_PATCH_VERSION == 0
#  define isnan    ::std::isnan
#  define isinf    ::std::isinf
#  define isfinite ::std::isfinite
#endif

#include <vtkMath.h> // DO NOT use double quotes which would cause endless recursive
                     // include of this file when filesystem is not case sensitive!

#if VTK_MAJOR_VERSION == 6 && VTK_MINOR_VERSION == 0 && VTK_PATCH_VERSION == 0
#  undef isnan
#  undef isinf
#  undef isfinite
#endif

namespace mirtk { namespace vtkmath {


// -----------------------------------------------------------------------------
/// Subtract two 2D vectors
inline void Subtract2D(const double a[2], const double b[2], double c[2])
{
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
}

// -----------------------------------------------------------------------------
/// Compute squared distance between 2D points
inline double Distance2BetweenPoints2D(const double a[2], const double b[2])
{
  const double dx = b[0] - a[0];
  const double dy = b[1] - a[1];
  return dx * dx + dy * dy;
}


} } // namespace mirtk::vtkmath


#endif // MIRTK_VtkMath_H

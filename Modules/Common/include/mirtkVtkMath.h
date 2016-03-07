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


// See http://www.paraview.org/Bug/view.php?id=14164
#if VTK_MAJOR_VERSION == 6 && VTK_MINOR_VERSION == 0 && VTK_PATCH_VERSION == 0
#  define isnan    ::std::isnan
#  define isinf    ::std::isinf
#  define isfinite ::std::isfinite
#endif

#include <vtkMath.h>

#if VTK_MAJOR_VERSION == 6 && VTK_MINOR_VERSION == 0 && VTK_PATCH_VERSION == 0
#  undef isnan
#  undef isinf
#  undef isfinite
#endif


#endif // MIRTK_VtkMath_H

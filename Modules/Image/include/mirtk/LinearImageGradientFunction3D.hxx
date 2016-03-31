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

#ifndef MIRTK_LinearImageGradientFunction3D_HXX
#define MIRTK_LinearImageGradientFunction3D_HXX

#include "mirtk/LinearImageGradientFunction3D.h"
#include "mirtk/LinearImageGradientFunction.hxx"


namespace mirtk {


// -----------------------------------------------------------------------------
template <class TImage>
GenericLinearImageGradientFunction3D<TImage>
::GenericLinearImageGradientFunction3D()
{
  this->NumberOfDimensions(3);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearImageGradientFunction3D<TImage>::GradientType
GenericLinearImageGradientFunction3D<TImage>
::GetInside(double x, double y, double z, double t) const
{
  return this->GetInside3D(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearImageGradientFunction3D<TImage>::GradientType
GenericLinearImageGradientFunction3D<TImage>
::GetWithPaddingInside(double x, double y, double z, double t) const
{
  return this->GetWithPaddingInside3D(x, y, z, t);
}


} // namespace mirtk

#endif // MIRTK_LinearImageGradientFunction3D_HXX

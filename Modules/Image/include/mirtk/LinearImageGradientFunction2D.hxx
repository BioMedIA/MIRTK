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

#ifndef MIRTK_LinearImageGradientFunction2D_HXX
#define MIRTK_LinearImageGradientFunction2D_HXX

#include "mirtk/LinearImageGradientFunction2D.h"
#include "mirtk/LinearImageGradientFunction.hxx"


namespace mirtk {


// -----------------------------------------------------------------------------
template <class TImage>
GenericLinearImageGradientFunction2D<TImage>
::GenericLinearImageGradientFunction2D()
{
  this->NumberOfDimensions(2);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearImageGradientFunction2D<TImage>::GradientType
GenericLinearImageGradientFunction2D<TImage>
::GetInside(double x, double y, double z, double t) const
{
  return this->GetInside2D(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearImageGradientFunction2D<TImage>::GradientType
GenericLinearImageGradientFunction2D<TImage>
::GetWithPaddingInside(double x, double y, double z, double t) const
{
  return this->GetWithPaddingInside2D(x, y, z, t);
}


} // namespace mirtk

#endif // MIRTK_LinearImageGradientFunction2D_HXX

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

#include "mirtk/IntensityCrossCorrelation.h"

#include "mirtk/Math.h"
#include "mirtk/ObjectFactory.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(IntensityCrossCorrelation);


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
IntensityCrossCorrelation
::IntensityCrossCorrelation(const char *name)
:
  ImageSimilarity(name)
{
}

// -----------------------------------------------------------------------------
IntensityCrossCorrelation
::IntensityCrossCorrelation(const IntensityCrossCorrelation &other)
:
  ImageSimilarity(other)
{
}

// -----------------------------------------------------------------------------
IntensityCrossCorrelation::~IntensityCrossCorrelation()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double IntensityCrossCorrelation::Evaluate()
{
  VoxelType *tgt = Target()->Data();
  VoxelType *src = Source()->Data();

  double x = .0, y = .0, xy = .0, x2 = .0, y2 = .0;
  int    n =  0;
  for (int idx = 0; idx < NumberOfVoxels(); ++idx, ++tgt, ++src) {
    if (IsForeground(idx)) {
      x  += *tgt;
      y  += *src;
      xy += (*tgt) * (*src);
      x2 += (*tgt) * (*tgt);
      y2 += (*tgt) * (*tgt);
      ++n;
    }
  }

  if (n == 0) return .0;
  return (xy - (x * y) / n) / (sqrt(x2 - x * x / n) * sqrt(y2 - y *y / n));
}

// -----------------------------------------------------------------------------
bool IntensityCrossCorrelation
::NonParametricGradient(const RegisteredImage *image, GradientImageType *gradient)
{
  cerr << "IntensityCrossCorrelation::NonParametricGradient: Not implemented" << endl;
  exit(1);

  // Apply chain rule to obtain gradient w.r.t y = T(x)
  MultiplyByImageGradient(image, gradient);
  return true;
}


} // namespace mirtk

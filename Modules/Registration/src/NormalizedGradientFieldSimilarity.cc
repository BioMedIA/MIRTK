/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Stefan Pszczolkowski Parraguez, Andreas Schuh
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

#include "mirtk/NormalizedGradientFieldSimilarity.h"

#include "mirtk/Math.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
NormalizedGradientFieldSimilarity
::NormalizedGradientFieldSimilarity(const char *name, double weight)
:
  GradientFieldSimilarity(name, weight),
  _TargetNoise(1.0),
  _SourceNoise(1.0),
  _TargetNormalization (100.0),
  _SourceNormalization (100.0)
{
}

// -----------------------------------------------------------------------------
NormalizedGradientFieldSimilarity
::NormalizedGradientFieldSimilarity(const NormalizedGradientFieldSimilarity &other)
:
  GradientFieldSimilarity(other),
  _TargetNoise(other._TargetNoise),
  _SourceNoise(other._SourceNoise),
  _TargetNormalization(other._TargetNormalization),
  _SourceNormalization(other._SourceNormalization)
{
}

// -----------------------------------------------------------------------------
NormalizedGradientFieldSimilarity::~NormalizedGradientFieldSimilarity()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool NormalizedGradientFieldSimilarity::SetWithPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Target noise") == 0) {
    return FromString(value, _TargetNoise) && _TargetNoise > 0;
  } else if (strcmp(param, "Source noise") == 0) {
    return FromString(value, _SourceNoise) && _SourceNoise > 0;
  }
  return GradientFieldSimilarity::SetWithPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList NormalizedGradientFieldSimilarity::Parameter() const
{
  ParameterList params = GradientFieldSimilarity::Parameter();
  Insert(params, "Target noise", ToString(_TargetNoise));
  Insert(params, "Source noise", ToString(_SourceNoise));
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double NormalizedGradientFieldSimilarity::NormalizationFactor(RegisteredImage *image, double noise) const
{
  double norm = .0;
  int    n    = 0;

  VoxelType *dx = image->Data(0, 0, 0, 1);
  VoxelType *dy = image->Data(0, 0, 0, 2);
  VoxelType *dz = image->Data(0, 0, 0, 3);

  for (int idx = 0; idx < NumberOfVoxels(); ++idx, ++dx, ++dy, ++dz) {
    if (IsForeground(idx)) {
      norm += sqrt((*dx) * (*dx) + (*dy) * (*dy) + (*dz) * (*dz));
      ++n;
    }
  }

  return ((n > 0) ? (norm * noise / n) : .0);
}

// -----------------------------------------------------------------------------
void NormalizedGradientFieldSimilarity::Update(bool hessian)
{
  const bool initial_update = _InitialUpdate;
  GradientFieldSimilarity::Update(hessian); // sets _InitialUpdate = false
  if (initial_update || _Target->Transformation()) {
    _TargetNormalization = NormalizationFactor(_Target, _TargetNoise);
    if (_TargetNormalization == 0) {
      cerr << "NormalizedGradientFieldSimilarity::Update: Target image has no structure!" << endl;
      exit(1);
    }
  }
  if (initial_update || _Source->Transformation()) {
    _SourceNormalization = NormalizationFactor(_Source, _SourceNoise);
    if (_SourceNormalization == 0) {
      cerr << "NormalizedGradientFieldSimilarity::Update: Source image has no structure!" << endl;
      exit(1);
    }
  }
}


} // namespace mirtk

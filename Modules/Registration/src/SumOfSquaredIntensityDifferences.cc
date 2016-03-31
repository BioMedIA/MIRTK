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

#include "mirtk/SumOfSquaredIntensityDifferences.h"

#include "mirtk/Math.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/ObjectFactory.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(SumOfSquaredIntensityDifferences);


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace SumOfSquaredIntensityDifferencesUtils {


// -----------------------------------------------------------------------------
// Types
typedef SumOfSquaredIntensityDifferences::VoxelType    VoxelType;
typedef SumOfSquaredIntensityDifferences::GradientType GradientType;

// -----------------------------------------------------------------------------
/// Sum the squared intensity differences
struct EvaluateSumOfSquaredDifferences : public VoxelReduction
{
  SumOfSquaredIntensityDifferences *_Sim;
  double                            _Sum;
  int                               _Cnt;

  EvaluateSumOfSquaredDifferences(SumOfSquaredIntensityDifferences *sim)
  :
    _Sim(sim), _Sum(.0), _Cnt(0)
  {}

  EvaluateSumOfSquaredDifferences(const EvaluateSumOfSquaredDifferences &rhs)
  :
    _Sim(rhs._Sim), _Sum(rhs._Sum), _Cnt(rhs._Cnt)
  {}

  void split(const EvaluateSumOfSquaredDifferences &lhs)
  {
    _Sum = .0;
    _Cnt =  0;
  }

  void join(const EvaluateSumOfSquaredDifferences &rhs)
  {
    _Sum += rhs._Sum;
    _Cnt += rhs._Cnt;
  }

  void operator ()(int i, int j, int k, int, VoxelType *t, VoxelType *s)
  {
    if (_Sim->IsForeground(i, j, k)) {
      _Sum += static_cast<double>(pow(*t - *s, 2));
      ++_Cnt;
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate gradient of sum of squared intensity differences
struct EvaluateSumOfSquaredDifferencesGradient : public VoxelFunction
{
  SumOfSquaredIntensityDifferences *_Sim;

  EvaluateSumOfSquaredDifferencesGradient(SumOfSquaredIntensityDifferences *sim)
  :
    _Sim(sim)
  {}

  void operator ()(int i, int j, int k, int, const VoxelType *t, const VoxelType *s, GradientType *g)
  {
    if (_Sim->IsForeground(i, j, k)) {
      *g = -2.0 * static_cast<double>(*t - *s);
    } else {
      *g = .0;
    }
  }
};


} // namespace SumOfSquaredIntensityDifferencesUtils
using namespace SumOfSquaredIntensityDifferencesUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
SumOfSquaredIntensityDifferences
::SumOfSquaredIntensityDifferences(const char *name)
:
  ImageSimilarity(name),
  _MaxSqDiff(1.0), _Value(.0), _N(0)
{
}

// -----------------------------------------------------------------------------
SumOfSquaredIntensityDifferences
::SumOfSquaredIntensityDifferences(const SumOfSquaredIntensityDifferences &other)
:
  ImageSimilarity(other),
  _MaxSqDiff(other._MaxSqDiff), _Value(other._Value), _N(other._N)
{
}

// -----------------------------------------------------------------------------
SumOfSquaredIntensityDifferences::~SumOfSquaredIntensityDifferences()
{
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void SumOfSquaredIntensityDifferences::Initialize()
{
  // Initialize base class
  ImageSimilarity::Initialize();

  // Determine maximum possible intensity difference
  double tmin, tmax, smin, smax;
  Target()->InputImage()->GetMinMaxAsDouble(&tmin, &tmax);
  Source()->InputImage()->GetMinMaxAsDouble(&smin, &smax);
  _MaxSqDiff = (tmin - smax) * (tmin - smax) > (tmax - smin) * (tmax - smin) ?
               (tmin - smax) * (tmin - smax) : (tmax - smin) * (tmax - smin);
  if (_MaxSqDiff == .0) _MaxSqDiff = 1.0;
  _Value = .0, _N =  0;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void SumOfSquaredIntensityDifferences::Update(bool gradient)
{
  // Upate base class and moving image(s)
  ImageSimilarity::Update(gradient);
  // Evaluate sum of squared differences over all voxels
  EvaluateSumOfSquaredDifferences ssd(this);
  ParallelForEachVoxel(_Domain, _Target, _Source, ssd);
  _Value = ssd._Sum, _N = ssd._Cnt;
}

// -----------------------------------------------------------------------------
void SumOfSquaredIntensityDifferences::Exclude(const blocked_range3d<int> &region)
{
  EvaluateSumOfSquaredDifferences ssd(this);
  ParallelForEachVoxel(region, _Target, _Source, ssd);
  _Value -= ssd._Sum, _N -= ssd._Cnt;
}

// -----------------------------------------------------------------------------
void SumOfSquaredIntensityDifferences::Include(const blocked_range3d<int> &region)
{
  EvaluateSumOfSquaredDifferences ssd(this);
  ParallelForEachVoxel(region, _Target, _Source, ssd);
  _Value += ssd._Sum, _N += ssd._Cnt;
}

// -----------------------------------------------------------------------------
double SumOfSquaredIntensityDifferences::Evaluate()
{
  return ((_N > 0) ? (_Value / (_N * _MaxSqDiff)) : .0);
}

// -----------------------------------------------------------------------------
bool SumOfSquaredIntensityDifferences
::NonParametricGradient(const RegisteredImage *image, GradientImageType *gradient)
{
  // Compute gradient of similarity w.r.t given moving image
  EvaluateSumOfSquaredDifferencesGradient eval(this);
  ParallelForEachVoxel(_Domain, (image == Target() ? Source() : Target()), image, gradient, eval);
  if (_N > 0) *gradient /= _N * _MaxSqDiff;

  // Apply chain rule to obtain gradient w.r.t y = T(x)
  MultiplyByImageGradient(image, gradient);
  return true;
}


} // namespace mirtk

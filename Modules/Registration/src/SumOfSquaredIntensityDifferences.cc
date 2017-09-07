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
  double                            _Scale;

  EvaluateSumOfSquaredDifferencesGradient(SumOfSquaredIntensityDifferences *sim, double scale = 1.)
  :
    _Sim(sim), _Scale(2. * scale)
  {}

  void operator ()(int i, int j, int k, int, const VoxelType *t, const VoxelType *s, GradientType *g)
  {
    if (_Sim->IsForeground(i, j, k)) {
      *g = -_Scale * static_cast<double>(*t - *s);
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
  _MinTargetIntensity(.0), _MaxTargetIntensity(1.0),
  _MinSourceIntensity(.0), _MaxSourceIntensity(1.0),
  _MaxSqDiff(1.0), _SumSqDiff(.0),
  _NumberOfForegroundVoxels(0)
{
}

// -----------------------------------------------------------------------------
void SumOfSquaredIntensityDifferences
::CopyAttributes(const SumOfSquaredIntensityDifferences &other)
{
  _MinTargetIntensity       = other._MinTargetIntensity;
  _MaxTargetIntensity       = other._MaxTargetIntensity;
  _MinSourceIntensity       = other._MinSourceIntensity;
  _MaxSourceIntensity       = other._MaxSourceIntensity;
  _MaxSqDiff                = other._MaxSqDiff;
  _SumSqDiff                = other._SumSqDiff;
  _NumberOfForegroundVoxels = other._NumberOfForegroundVoxels;
}

// -----------------------------------------------------------------------------
SumOfSquaredIntensityDifferences
::SumOfSquaredIntensityDifferences(const SumOfSquaredIntensityDifferences &other)
:
  ImageSimilarity(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
SumOfSquaredIntensityDifferences &SumOfSquaredIntensityDifferences
::operator =(const SumOfSquaredIntensityDifferences &other)
{
  if (this != &other) {
    ImageSimilarity::operator =(other);
    CopyAttributes(other);
  }
  return *this;
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

  // Determine maximum intensities and maximum possible intensity difference
  Target()->InputImage()->GetMinMaxAsDouble(_MinTargetIntensity, _MaxTargetIntensity);
  Source()->InputImage()->GetMinMaxAsDouble(_MinSourceIntensity, _MaxSourceIntensity);
  _MaxSqDiff = pow(max(abs(_MaxSourceIntensity - _MinTargetIntensity),
                       abs(_MaxTargetIntensity - _MinSourceIntensity)), 2);
  if (_MaxSqDiff == .0) _MaxSqDiff = 1.0;

  // Initialize sum of squared differences
  _SumSqDiff = .0;
  _NumberOfForegroundVoxels =  0;
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
  _SumSqDiff = ssd._Sum, _NumberOfForegroundVoxels = ssd._Cnt;
}

// -----------------------------------------------------------------------------
void SumOfSquaredIntensityDifferences::Exclude(const blocked_range3d<int> &region)
{
  EvaluateSumOfSquaredDifferences ssd(this);
  ParallelForEachVoxel(region, _Target, _Source, ssd);
  _SumSqDiff -= ssd._Sum, _NumberOfForegroundVoxels -= ssd._Cnt;
}

// -----------------------------------------------------------------------------
void SumOfSquaredIntensityDifferences::Include(const blocked_range3d<int> &region)
{
  EvaluateSumOfSquaredDifferences ssd(this);
  ParallelForEachVoxel(region, _Target, _Source, ssd);
  _SumSqDiff += ssd._Sum, _NumberOfForegroundVoxels += ssd._Cnt;
}

// -----------------------------------------------------------------------------
double SumOfSquaredIntensityDifferences::Evaluate()
{
  if (_NumberOfForegroundVoxels == 0) return .0;
  return _SumSqDiff / (_NumberOfForegroundVoxels * _MaxSqDiff);
}

// -----------------------------------------------------------------------------
bool SumOfSquaredIntensityDifferences
::NonParametricGradient(const RegisteredImage *image, GradientImageType *gradient)
{
  // Gradient is zero when input has no foreground
  if (_NumberOfForegroundVoxels <= 0) {
    *gradient = 0.;
    return true;
  }

  // Normalization factor
  double norm = 1. / (_NumberOfForegroundVoxels * _MaxSqDiff);

  // Compute gradient of similarity w.r.t given moving image
  EvaluateSumOfSquaredDifferencesGradient eval(this, norm);
  ParallelForEachVoxel(_Domain, (image == Target() ? Source() : Target()), image, gradient, eval);

  // Apply chain rule to obtain gradient w.r.t y = T(x)
  MultiplyByImageGradient(image, gradient);
  return true;
}


} // namespace mirtk

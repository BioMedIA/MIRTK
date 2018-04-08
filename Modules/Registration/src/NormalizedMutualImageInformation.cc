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

#include "mirtk/NormalizedMutualImageInformation.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/BSpline.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/ObjectFactory.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(NormalizedMutualImageInformation);


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace NormalizedMutualImageInformationUtils {


// -----------------------------------------------------------------------------
template <class BinType>
inline double ValToRange(const Histogram1D<BinType> &hist, double val)
{
  if (val < hist.Min()) return .0;
  if (val > hist.Max()) return static_cast<double>(hist.NumberOfBins() - 1);
  return hist.ValToRange(val);
}

// -----------------------------------------------------------------------------
class CalculateGradient : public VoxelFunction
{
  typedef NormalizedMutualImageInformation::BinType BinType;

  const NormalizedMutualImageInformation *_This;
  const Histogram2D<BinType>             &_LogJointHistogram;
  const Histogram1D<BinType>             &_LogMarginalXHistogram;
  const Histogram1D<BinType>             &_LogMarginalYHistogram;
  double                                  _NormalizedMutualInformation;
  double                                  _NormalizedJointEntropy;

  /// Round value to fixed precision to avoid floating point differences of
  /// NMI gradient in a symmetric (SVFFD) registration when target and source
  /// images are exchanged in order to obtain inverse consistent results.
  ///
  /// Moreover, noted that optimization would perform more iterations with this
  /// fixed precision NMI derivative than without, i.e., yielding higher final NMI.
  inline double RoundDerivativeValue(double x) const
  {
    return 1e-9 * static_cast<double>(static_cast<long long>(1e+9 * static_cast<double>(x)));
  }

public:

  CalculateGradient(const NormalizedMutualImageInformation *_this,
                    const Histogram2D<BinType> &logJointHistogram,
                    const Histogram1D<BinType> &logMarginalXHistogram,
                    const Histogram1D<BinType> &logMarginalYHistogram,
                    double je_norm, double nmi)
  :
    _This(_this),
    _LogJointHistogram(logJointHistogram),
    _LogMarginalXHistogram(logMarginalXHistogram),
    _LogMarginalYHistogram(logMarginalYHistogram),
    _NormalizedMutualInformation(nmi),
    _NormalizedJointEntropy(je_norm)
  {}

  template <class TImage, class TIntensity, class TGradient>
  void operator ()(const TImage &, int idx, const TIntensity *tgt, const TIntensity *src, TGradient *deriv)
  {
    if (_This->IsForeground(idx)) {
      const int    target_nbins = _LogMarginalXHistogram.NumberOfBins();
      const int    source_nbins = _LogMarginalYHistogram.NumberOfBins();
      const double target_value = ValToRange(_LogMarginalXHistogram, *tgt);
      const double source_value = ValToRange(_LogMarginalYHistogram, *src);

      int t1 = static_cast<int>(     target_value ) - 1;
      int t2 = static_cast<int>(ceil(target_value)) + 1;
      int s1 = static_cast<int>(     source_value ) - 1;
      int s2 = static_cast<int>(ceil(source_value)) + 1;

      if (t1 <  0           ) t1 = 0;
      if (t2 >= target_nbins) t2 = target_nbins - 1;
      if (s1 <  0           ) s1 = 0;
      if (s2 >= source_nbins) s2 = source_nbins - 1;

      double jointEntropyGrad  = 0.;
      double targetEntropyGrad = 0.;
      double sourceEntropyGrad = 0.;
      double w;

      for (int t = t1; t <= t2; ++t)
      for (int s = s1; s <= s2; ++s) {
        w = BSpline<double>::B  (static_cast<double>(t) - target_value) *
            BSpline<double>::B_I(static_cast<double>(s) - source_value);
        jointEntropyGrad  += w * static_cast<double>(_LogJointHistogram(t, s));
        targetEntropyGrad += w * static_cast<double>(_LogMarginalXHistogram(t));
        sourceEntropyGrad += w * static_cast<double>(_LogMarginalYHistogram(s));
      }

      (*deriv) = (targetEntropyGrad + sourceEntropyGrad - _NormalizedMutualInformation * jointEntropyGrad) / _NormalizedJointEntropy;
      (*deriv) = RoundDerivativeValue(*deriv);
    }
  }
};


} // namespace NormalizedMutualImageInformationUtils
using namespace NormalizedMutualImageInformationUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void NormalizedMutualImageInformation
::CopyAttributes(const NormalizedMutualImageInformation &other)
{
  _TargetHistogram = other._TargetHistogram;
  _SourceHistogram = other._SourceHistogram;
  _TargetEntropy   = other._TargetEntropy;
  _SourceEntropy   = other._SourceEntropy;
  _JointEntropy    = other._JointEntropy;
}

// -----------------------------------------------------------------------------
NormalizedMutualImageInformation
::NormalizedMutualImageInformation(const char *name)
:
  HistogramImageSimilarity(name),
  _TargetHistogram(0),
  _SourceHistogram(0),
  _TargetEntropy(NaN),
  _SourceEntropy(NaN),
  _JointEntropy(NaN)
{
}

// -----------------------------------------------------------------------------
NormalizedMutualImageInformation
::NormalizedMutualImageInformation(const NormalizedMutualImageInformation &other)
:
  HistogramImageSimilarity(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
NormalizedMutualImageInformation &NormalizedMutualImageInformation
::operator =(const NormalizedMutualImageInformation &other)
{
  if (this != &other) {
    HistogramImageSimilarity::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
NormalizedMutualImageInformation::~NormalizedMutualImageInformation()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void NormalizedMutualImageInformation::Update(bool gradient)
{
  // Update base class and moving image(s)
  HistogramImageSimilarity::Update(gradient);

  // Update marginal histograms
  _Histogram.HistogramX(_TargetHistogram);
  _Histogram.HistogramY(_SourceHistogram);

  _TargetEntropy = _TargetHistogram.Entropy();
  _SourceEntropy = _SourceHistogram.Entropy();
  _JointEntropy  = _Histogram.JointEntropy();
}

// -----------------------------------------------------------------------------
double NormalizedMutualImageInformation::Evaluate()
{
  return 2.0 - (_TargetEntropy + _SourceEntropy) / _JointEntropy;
}

// -----------------------------------------------------------------------------
double NormalizedMutualImageInformation::RawValue(double value) const
{
  return 2.0 - HistogramImageSimilarity::RawValue(value);
}

// -----------------------------------------------------------------------------
// This code is based on an idea from Marc Modat for computing the NMI
// derivative as suggested in his niftyreg package
bool NormalizedMutualImageInformation
::NonParametricGradient(const RegisteredImage *image, GradientImageType *gradient)
{
  MIRTK_START_TIMING();

  RegisteredImage     *fixed;
  Histogram2D<BinType> jhist(0, 0);
  Histogram1D<BinType> xhist(0);
  Histogram1D<BinType> yhist(0);

  // Swap target and source if similarity derived w.r.t transformed "target"
  if (image == Target()) {
    fixed = Source();
    jhist = _Histogram.Transposed();
    xhist = _SourceHistogram;
    yhist = _TargetHistogram;
  } else {
    fixed = Target();
    jhist = _Histogram;
    xhist = _TargetHistogram;
    yhist = _SourceHistogram;
  }

  // Log transform histograms
  jhist.Log();
  xhist.Log();
  yhist.Log();

  // Normalized mutual information value
  double nmi = (_TargetEntropy + _SourceEntropy) / _JointEntropy;

  // Joint entropy value times normalisation factor, including negative sign
  double je_norm = - _JointEntropy * _Histogram.NumberOfSamples();

  // Evaluate similarity gradient w.r.t given transformed image
  CalculateGradient eval(this, jhist, xhist, yhist, je_norm, nmi);
  memset(gradient->Data(), 0, _NumberOfVoxels * sizeof(GradientType));
  ParallelForEachVoxel(fixed, image, gradient, eval);

  // Apply chain rule to obtain gradient w.r.t y = T(x)
  MultiplyByImageGradient(image, gradient);

  MIRTK_DEBUG_TIMING(2, "similarity gradient computation (NMI)");
  return true;
}


} // namespace

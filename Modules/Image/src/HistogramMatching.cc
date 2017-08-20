/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
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

#include "mirtk/HistogramMatching.h"

#include "mirtk/Histogram1D.h"
#include "mirtk/Profiling.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
HistogramMatching<VoxelType>::HistogramMatching()
:
  _NumberOfBins(512),
  _NumberOfSteps(10),
  _CutOff(.01)
{
}

// =============================================================================
// Execution
// =============================================================================

// ---------------------------------------------------------------------------
template <class VoxelType>
void HistogramMatching<VoxelType>::Initialize()
{
  // Check reference image
  if (_Reference == nullptr) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Reference image required");
  }

  // Initialize base class
  ImageToImage<VoxelType>::Initialize();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HistogramMatching<VoxelType>::Run()
{
  MIRTK_START_TIMING();

  this->Initialize();

  const ImageType * const input  = this->Input();
  const ImageType * const target = this->Reference();
  ImageType       * const output = this->Output();

  // Compute source image histogram
  Histogram1D<int> src_hist;
  double src_min, src_max;
  input->GetMinMaxAsDouble(&src_min, &src_max);
  src_hist.Min(src_min);
  src_hist.Max(src_max);
  src_hist.PutNumberOfBins(_NumberOfBins);
  for (int vox = 0; vox < input->NumberOfVoxels(); ++vox) {
    if (input->IsForeground(vox)) {
      src_hist.AddSample(input->GetAsDouble(vox));
    }
  }

  // Compute target image histogram
  Histogram1D<int> tgt_hist;
  double tgt_min, tgt_max;
  _Reference->GetMinMaxAsDouble(&tgt_min, &tgt_max);
  tgt_hist.Min(tgt_min);
  tgt_hist.Max(tgt_max);
  tgt_hist.PutNumberOfBins(_NumberOfBins);
  for (int vox = 0; vox < target->NumberOfVoxels(); ++vox) {
    if (target->IsForeground(vox)) {
      tgt_hist.AddSample(target->GetAsDouble(vox));
    }
  }

  // Compute linear maps
  double pct  = _CutOff / 100.;
  double step = (1. - 2. * _CutOff / 100.) / _NumberOfSteps;

  Array<double> src_levels, tgt_levels;
  src_levels.reserve(_NumberOfSteps + 3);
  tgt_levels.reserve(_NumberOfSteps + 3);

  src_levels.push_back(src_min);
  tgt_levels.push_back(tgt_min);
  for (int i = 0; i <= _NumberOfSteps; ++i, pct += step) {
    double src_level = src_hist.CDFToVal(pct);
    double tgt_level = tgt_hist.CDFToVal(pct);
    if (tgt_level - tgt_levels.back() < 1e-6 * (tgt_max - tgt_min)) {
      cerr << this->NameOfClass() << ": " << pct * 100. << " percentile collapses in reference, skipping" << endl;
    } else {
      src_levels.push_back(src_level);
      tgt_levels.push_back(tgt_level);
    }
  }
  src_levels.push_back(src_max);
  tgt_levels.push_back(tgt_max);

  // Map input intensities
  double value, bg = -inf;
  if (target->HasBackgroundValue()) {
    bg = target->GetBackgroundValueAsDouble();
  } else {
    bg = tgt_min - 1.;
  }
  for (int bin, vox = 0; vox < input->NumberOfVoxels(); ++vox) {
    if (input->IsForeground(vox)) {
      value = input->GetAsDouble(vox);
      for (bin = 0; bin < static_cast<int>(src_levels.size()); ++bin) {
        if (value <= src_levels[bin]) break;
      }
      if (bin == 0) {
        value = tgt_levels.front();
      } else if (bin >= static_cast<int>(src_levels.size() - 1)) {
        value = tgt_levels.back();
      } else {
        value  = (value - src_levels[bin - 1]) / (src_levels[bin] - src_levels[bin - 1]);
        value *= tgt_levels[bin] - tgt_levels[bin - 1];
        value += tgt_levels[bin-1];
      }
    } else {
      value = bg;
    }
    output->PutAsDouble(vox, value);
  }
  if (input->HasBackgroundValue()) {
    output->PutBackgroundValueAsDouble(bg);
  }
  if (input->HasMask()) {
    output->PutMask(new BinaryImage(*input->GetMask()), true);
  }

  this->Finalize();

  MIRTK_DEBUG_TIMING(5, this->NameOfClass());
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class HistogramMatching<float>;
template class HistogramMatching<double>;


} // namespace mirtk

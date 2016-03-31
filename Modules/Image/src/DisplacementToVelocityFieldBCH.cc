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

#include "mirtk/DisplacementToVelocityFieldBCH.h"

#include "mirtk/Deallocate.h"
#include "mirtk/ImageFunction.h"
#include "mirtk/GaussianBlurring2D.h"
#include "mirtk/LieBracketImageFilter.h"
#include "mirtk/VelocityToDisplacementFieldSS.h"
#include "mirtk/NaryVoxelFunction.h"


namespace mirtk {


// ===========================================================================
// Construction/Destruction
// ===========================================================================

// ---------------------------------------------------------------------------
template <class VoxelType>
DisplacementToVelocityFieldBCH<VoxelType>::DisplacementToVelocityFieldBCH()
:
  _ExponentialFilter(new VelocityToDisplacementFieldSS<VoxelType>),
  _CustomExponentialFilter(false),
  _NumberOfIterations(8),
  _NumberOfTerms(3),
  _UseJacobian(false),
  _SmoothVelocities(false),
  _dv(NULL), _l1(NULL), _l2(NULL), _l3(NULL), _l4(NULL)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
DisplacementToVelocityFieldBCH<VoxelType>::~DisplacementToVelocityFieldBCH()
{
  if (!_CustomExponentialFilter) delete _ExponentialFilter;
}

// ===========================================================================
// Filter implementation
// ===========================================================================

// ----------------------------------------------------------------------------
template <class VoxelType>
void DisplacementToVelocityFieldBCH<VoxelType>
::ExponentialFilter(VelocityToDisplacementField<VoxelType> *filter)
{
  if (!_CustomExponentialFilter) delete _ExponentialFilter;
  _ExponentialFilter       = filter;
  _CustomExponentialFilter = (_ExponentialFilter != NULL);
}

// ----------------------------------------------------------------------------
template <class VoxelType>
void DisplacementToVelocityFieldBCH<VoxelType>::UpperIntegrationLimit(double t)
{
  _ExponentialFilter->UpperIntegrationLimit(t);
}

// ----------------------------------------------------------------------------
template <class VoxelType>
double DisplacementToVelocityFieldBCH<VoxelType>::UpperIntegrationLimit()
{
  return _ExponentialFilter->UpperIntegrationLimit();
}

// ----------------------------------------------------------------------------
template <class VoxelType>
void DisplacementToVelocityFieldBCH<VoxelType>::NumberOfSteps(int n)
{
  _ExponentialFilter->NumberOfSteps(n);
}

// ----------------------------------------------------------------------------
template <class VoxelType>
int DisplacementToVelocityFieldBCH<VoxelType>::NumberOfSteps()
{
  return _ExponentialFilter->NumberOfSteps();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void DisplacementToVelocityFieldBCH<VoxelType>::Initialize()
{
  // Initialize base class
  DisplacementToVelocityField<VoxelType>::Initialize();

  // Get attributes of input vector field
  const ImageAttributes &grid = this->Input()->Attributes();

  // Check input
  if ((grid._z <= 1 && grid._t != 2 && grid._t != 3) || (grid._z > 1 && grid._t != 3)) {
    cerr << this->NameOfClass() << "::Initialize: Input must be a valid vector field" << endl;
    exit(1);
  }

  // Check parameters
  if (_NumberOfTerms < 2 || _NumberOfTerms > 6) {
    cerr << this->NameOfClass() << "::Initialize: Number of BCH terms must be 2-6" << endl;
    exit(1);
  }

  // Ensure that exponential filter is set
  if (_ExponentialFilter == NULL) {
    // Developer may by mistake passed a NULL  pointer to SetExponentialFilter,
    // so let them know about this mistake rather than instantiating a default
    // exponential filter here. The default filter is created in the constructor.
    cerr << this->NameOfClass() << "::Initialize: No filter for exponentiation of velocity field set" << endl;
    exit(1);
  }

  // Allocate intermediate images
  switch (_NumberOfTerms) {
    // Attention: Must be in descending order and without break statements!
    case 6: _l4 = new ImageType(grid);
    case 5: _l3 = new ImageType(grid);
    case 4: _l2 = new ImageType(grid);
    case 3: _l1 = new ImageType(grid);
    case 2: _dv = new ImageType(grid);
  };

  // Initialize output velocity field
  // v_0 = 0
  // v_1 = exp(-v_0) ° phi = phi
  this->Output()->CopyFrom(this->Input()->Data());

  // Initialize exponential filter
  _ExponentialFilter->Input(0, this->Output());
  _ExponentialFilter->Input(1, this->Input ());
  _ExponentialFilter->Output(_dv);
  _ExponentialFilter->ComputeInverse(true);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void DisplacementToVelocityFieldBCH<VoxelType>::Finalize()
{
  // Deallocate intermediate images
  Delete(_dv);
  Delete(_l1);
  Delete(_l2);
  Delete(_l3);
  Delete(_l4);
  // Finalize base class
  DisplacementToVelocityField<VoxelType>::Finalize();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void DisplacementToVelocityFieldBCH<VoxelType>::Run()
{
  MIRTK_START_TIMING();

  // Do the initial set up
  this->Initialize();

  // Get pointer to output (AFTER initialization of base class!)
  ImageType *v = this->Output();

  // Iteratively update velocity field using Baker-Campbell-Hausdorff formula
  for (int n = 0; n < _NumberOfIterations; n++) {
    // Compute dv = exp(-v) ° phi
    _ExponentialFilter->Run();
    // Smooth to stabilize computation
    // (see Insight Journal article of Vercauteren et al. at http://hdl.handle.net/10380/3060)
    if (_SmoothVelocities) {
      GaussianBlurring<VoxelType> blur(2.0 * v->GetXSize(), 2.0 * v->GetYSize(),
                                       v->Z() > 1 ? 2.0 * v->GetZSize() : .0);
      //blur.Input (v);
      //blur.Output(v);
      //blur.Run();
      blur.Input (_dv);
      blur.Output(_dv);
      blur.Run();
    }
    // Calculate required Lie brackets
    if (_NumberOfTerms > 2) liebracket(_l1,  v,  _dv, _UseJacobian); //          [v, dv]
    if (_NumberOfTerms > 3) liebracket(_l2,  v,  _l1, _UseJacobian); // [v,      [v, dv]]
    if (_NumberOfTerms > 4) liebracket(_l3, _dv, _l1, _UseJacobian); // [dv,     [v, dv]]
    if (_NumberOfTerms > 5) liebracket(_l4, _dv, _l2, _UseJacobian); // [dv, [v, [v, dv]]]
    // Compute update using the BCH formula
    NaryVoxelFunction::EvaluateBCHFormula bch;
    if      (_NumberOfTerms == 2) ParallelForEachScalar(v, _dv,                     v, bch);
    else if (_NumberOfTerms == 3) ParallelForEachScalar(v, _dv, _l1,                v, bch);
    else if (_NumberOfTerms == 4) ParallelForEachScalar(v, _dv, _l1, _l2,           v, bch);
    else if (_NumberOfTerms == 5) ParallelForEachScalar(v, _dv, _l1, _l2, _l3,      v, bch);
    else if (_NumberOfTerms == 6) ParallelForEachScalar(v, _dv, _l1, _l2, _l3, _l4, v, bch);
  }

  // Do the final cleaning up
  this->Finalize();

  MIRTK_DEBUG_TIMING(1, "DisplacementToVelocityFieldBCH");
}

// ===========================================================================
// Explicit template instantiations
// ===========================================================================

template class DisplacementToVelocityFieldBCH<float>;
template class DisplacementToVelocityFieldBCH<double>;


} // namespace mirtk

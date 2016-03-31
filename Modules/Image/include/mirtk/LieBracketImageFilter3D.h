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

#ifndef MIRTK_LieBracketImageFilter3D_H
#define MIRTK_LieBracketImageFilter3D_H


#include "mirtk/LieBracketImageFilter.h"

#include "mirtk/Matrix.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"


namespace mirtk {


/**
 * Image filter for computation of Lie bracket of two 3D vector fields.
 */

template <class TVoxel>
class LieBracketImageFilter3D : public LieBracketImageFilter<TVoxel>
{
  mirtkImageFilterMacro(LieBracketImageFilter3D, TVoxel);

public:

  /// Type of superclass of this concrete image filter
  typedef LieBracketImageFilter<TVoxel> Superclass;

protected:

  /// World to image matrix (excluding translation)
  Matrix _MatW2I;

  /// Compute 1st order derivatives of given vector field
  void Jacobian(Matrix &, const ImageType &, int, int, int);

  /// Initialize filter
  virtual void Initialize();

public:

  /// Constructor
  LieBracketImageFilter3D();

  /// Destructor
  virtual ~LieBracketImageFilter3D();

  // Import other overloads
  using Superclass::Output;

  /// Set output
  virtual void Output(ImageType *);

  /// Run filter on every voxel
  virtual void Run();

  /// Run filter on single voxel
  virtual void Run(double [3], int, int, int);

  /// Run filter on single voxel and component
  virtual double Run(int, int, int, int);

};

///////////////////////////////////////////////////////////////////////////////
// Inline definitions
///////////////////////////////////////////////////////////////////////////////

// --------------------------------------------------------------------------
template <class VoxelType>
LieBracketImageFilter3D<VoxelType>::LieBracketImageFilter3D()
:
  _MatW2I(3, 3)
{
  _MatW2I.Ident();
}

// --------------------------------------------------------------------------
template <class VoxelType>
LieBracketImageFilter3D<VoxelType>::~LieBracketImageFilter3D()
{
}

// --------------------------------------------------------------------------
template <class VoxelType>
void LieBracketImageFilter3D<VoxelType>::Output(ImageType *output)
{
  Superclass::Output(output);
  _MatW2I = output->GetWorldToImageMatrix()(0, 0, 3, 3);
}

// --------------------------------------------------------------------------
template <class VoxelType>
void LieBracketImageFilter3D<VoxelType>::Initialize()
{
  // Initialize base class
  Superclass::Initialize();

  // Ensure that input vector fields are 3D
  // Note that base class already ensures that inputs and output have same attributes
  // Accept also a 2D slice of a 3D vector field as input
  if (this->GetInput(0)->T() != 3) {
    cerr << this->NameOfClass() << "::Initialize: Input images are no 3D vector fields" << endl;
    exit(1);
  }
}

// ---------------------------------------------------------------------------
// Note: Using NN extrapolation at boundary
template <class VoxelType>
void LieBracketImageFilter3D<VoxelType>
::Jacobian(Matrix &jac, const ImageType &v, int i, int j, int k)
{
  int a, b;
  // Finite difference in x dimension
  if (i <= 0) {
    a = i;
    b = i + 1;
  } else if (i >= v.GetX()-1) {
    a = i - 1;
    b = i;
  } else {
    a = i - 1;
    b = i + 1;
  }
  jac(0, 0) = 0.5 * (v(b, j, k, 0) - v(a, j, k, 0));
  jac(1, 0) = 0.5 * (v(b, j, k, 1) - v(a, j, k, 1));
  jac(2, 0) = 0.5 * (v(b, j, k, 2) - v(a, j, k, 2));
  // Finite difference in y dimension
  if (j <= 0) {
    a = j;
    b = j + 1;
  } else if (j >= v.GetY()-1) {
    a = j - 1;
    b = j;
  } else {
    a = j - 1;
    b = j + 1;
  }
  jac(0, 1) = 0.5 * (v(i, b, k, 0) - v(i, a, k, 0));
  jac(1, 1) = 0.5 * (v(i, b, k, 1) - v(i, a, k, 1));
  jac(2, 1) = 0.5 * (v(i, b, k, 2) - v(i, a, k, 2));
  // Finite difference in z dimension
  if (k <= 0) {
    a = k;
    b = k + 1;
  } else if(k >= v.GetZ()-1) {
    a = k - 1;
    b = k;
  } else {
    a = k - 1;
    b = k + 1;
  }
  jac(0, 2) = 0.5 * (v(i, j, b, 0) - v(i, j, a, 0));
  jac(1, 2) = 0.5 * (v(i, j, b, 1) - v(i, j, a, 1));
  jac(2, 2) = 0.5 * (v(i, j, b, 2) - v(i, j, a, 2));
  // Project derivative from image to world space
  jac *= _MatW2I;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
double LieBracketImageFilter3D<VoxelType>::Run(int i, int j, int k, int t)
{
  Matrix lJ(3, 3), rJ(3, 3);

  const ImageType &lv = *this->GetInput(0);
  const ImageType &rv = *this->GetInput(1);

  const VoxelType &lx = lv(i, j, k, 0);
  const VoxelType &ly = lv(i, j, k, 1);
  const VoxelType &lz = lv(i, j, k, 2);
  const VoxelType &rx = rv(i, j, k, 0);
  const VoxelType &ry = rv(i, j, k, 1);
  const VoxelType &rz = rv(i, j, k, 2);

  Jacobian(lJ, lv, i, j, k);
  Jacobian(rJ, rv, i, j, k);

  return (lJ(t, 0) * rx - lx * rJ(t, 0)) +
         (lJ(t, 1) * ry - ly * rJ(t, 1)) +
         (lJ(t, 2) * rz - lz * rJ(t, 2));
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void LieBracketImageFilter3D<VoxelType>::Run(double vec[3], int i, int j, int k)
{
  Matrix lJ(3, 3), rJ(3, 3);

  const ImageType &lv = *this->GetInput(0);
  const ImageType &rv = *this->GetInput(1);

  const VoxelType &lx = lv(i, j, k, 0);
  const VoxelType &ly = lv(i, j, k, 1);
  const VoxelType &lz = lv(i, j, k, 2);
  const VoxelType &rx = rv(i, j, k, 0);
  const VoxelType &ry = rv(i, j, k, 1);
  const VoxelType &rz = rv(i, j, k, 2);

  Jacobian(lJ, lv, i, j, k);
  Jacobian(rJ, rv, i, j, k);

  vec[0] = (lJ(0, 0) * rx - lx * rJ(0, 0)) +
           (lJ(0, 1) * ry - ly * rJ(0, 1)) +
           (lJ(0, 2) * rz - lz * rJ(0, 2));

  vec[1] = (lJ(1, 0) * rx - lx * rJ(1, 0)) +
           (lJ(1, 1) * ry - ly * rJ(1, 1)) +
           (lJ(1, 2) * rz - lz * rJ(1, 2));

  vec[2] = (lJ(2, 0) * rx - lx * rJ(2, 0)) +
           (lJ(2, 1) * ry - ly * rJ(2, 1)) +
           (lJ(2, 2) * rz - lz * rJ(2, 2));
}

namespace LieBracketImageFilter3DUtils {

// ---------------------------------------------------------------------------
/// Parallelizable body of LieBracketImageFilter3D::Run.
template <class VoxelType>
class Run
{
private:

  typedef LieBracketImageFilter3D<VoxelType> FilterType;
  typedef typename FilterType::ImageType     ImageType;

  FilterType *_LieBracketFilter; ///< Lie bracket filter
  ImageType  *_Output;           ///< Output vector field

public:

  /// Constructor
  Run(FilterType *filter, ImageType *output)
  :
    _LieBracketFilter(filter),
    _Output(output)
  {}

  /// Run Lie bracket filter for each voxel in specified range
  void operator ()(const blocked_range3d<int> &r) const
  {
    double vec[3];
    for (int k = r.pages().begin(); k != r.pages().end(); ++k)
    for (int j = r.rows ().begin(); j != r.rows ().end(); ++j)
    for (int i = r.cols ().begin(); i != r.cols ().end(); ++i) {
      _LieBracketFilter->Run(vec, i, j, k);
      _Output->PutAsDouble(i, j, k, 0, vec[0]);
      _Output->PutAsDouble(i, j, k, 1, vec[1]);
      _Output->PutAsDouble(i, j, k, 2, vec[2]);
    }
  }
};

} // namespace LieBracketImageFilter3DUtils

// ---------------------------------------------------------------------------
template <class VoxelType>
void LieBracketImageFilter3D<VoxelType>::Run()
{
  blocked_range3d<int> voxels(0, this->GetInput(0)->Z(),
                              0, this->GetInput(0)->Y(),
                              0, this->GetInput(0)->X());
  LieBracketImageFilter3DUtils::Run<VoxelType> run(this, this->Output());
  MIRTK_START_TIMING();
  this->Initialize();
  parallel_for(voxels, run);
  this->Finalize();
  MIRTK_DEBUG_TIMING(2, "LieBracketImageFilter");
}


} // namespace mirtk

#endif // MIRTK_LieBracketImageFilter3D_H

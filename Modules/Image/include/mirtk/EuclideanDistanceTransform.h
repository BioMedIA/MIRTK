/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
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

#ifndef MIRTK_EuclideanDistanceTransform_H
#define MIRTK_EuclideanDistanceTransform_H

#include "mirtk/ImageToImage.h"


namespace mirtk {


template <class TVoxel>
class EuclideanDistanceTransform : public ImageToImage<TVoxel>
{
  mirtkInPlaceImageFilterMacro(EuclideanDistanceTransform, TVoxel);

public:

  /// 2D or 3D distance transform
  enum Mode { DT_2D, DT_3D };

protected:

  /// 2D or 3D distance transform
  Mode _distanceTransformMode;

  /// Calculate the Vornoi diagram
  int edtVornoiEDT(long *, long);

  /// Calculate 2D distance transform
  void edtComputeEDT_2D(char *, long *, long, long);

  /// Calculate 3D distance transform
  void edtComputeEDT_3D(char *, long *, long, long, long);

  /// Calculate the Vornoi diagram for anisotripic voxel sizes
  int edtVornoiEDT_anisotropic(VoxelType *, long, double);

  /// Calculate 2D distance transform for anisotripic voxel sizes
  void edtComputeEDT_2D_anisotropic(const VoxelType *, VoxelType *, long, long, double, double);

  /// Calculate 3D distance transform for anisotripic voxel sizes
  void edtComputeEDT_3D_anisotropic(const VoxelType *, VoxelType *, long, long, long, double, double, double);

public:

  /// Default constructor
  EuclideanDistanceTransform(Mode = DT_3D);

  /// Destructor (empty).
  ~EuclideanDistanceTransform() {};

  // Run distance transform
  virtual void Run();

  // Get Radial
  virtual void Radial();

  // Get Radial+Thickness
  virtual void TRadial();
};


} // namespace mirtk

#endif // MIRTK_EuclideanDistanceTransform_H

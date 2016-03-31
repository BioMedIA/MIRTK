/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2011, 2016 Imperial College London
 * Copyright 2011       Paul Aljabar
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

#include "mirtk/CityBlockDistanceTransform.h"

#include "mirtk/Array.h"
#include "mirtk/BaseImage.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
CityBlockDistanceTransform<VoxelType>::CityBlockDistanceTransform()
{
  _flipType = FlipNone;
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void CityBlockDistanceTransform<VoxelType>::Initialize2D()
{
  // May need to flip input...
  GenericImage<VoxelType> *input  = const_cast<GenericImage<VoxelType> *>(this->Input());
  GenericImage<VoxelType> *output = this->Output();

  // Make an array with an extra layer in each direction.
  ImageAttributes attr = input->Attributes();

  // One of the dimensions should be a singleton.
  if (attr._x != 1 && attr._y != 1 && attr._z != 1){
  	cerr << "CityBlockDistanceTransform<VoxelType>::Initialize2D: One spatial dimension should be a singleton." << endl;
  	exit(1);
  }

  // Ensure the singleton spatial dimension is the z-dimension.

  // Default.
  _flipType = FlipNone;
 
  if (attr._x == 1) {
  	input ->FlipXZ(0);
  	output->FlipXZ(0);
  	_flipType = FlipXZ;
  } else if (attr._y == 1) {
  	input ->FlipYZ(0);
  	output->FlipYZ(0);
  	_flipType = FlipYZ;
  }

  // Get the attributes again as they may have changed.
  attr = input->Attributes();

  // Increment the planar dimensions.
  attr._x += 2;
  attr._y += 2;

  _data.Initialize(attr);

  const int nx = attr._x;
  const int ny = attr._y;
  const int nt = attr._t;

  for (int l = 0; l < nt; ++l) {

    // Outer rows and columns are duplicates of outer rows and
    // columns of input image.
    for (int j = 1; j < ny - 1; ++j) {
      _data.Put(0   , j, 0, l, input->Get(0   , j-1, 0, l) > 0 ? 1 : 0);
      _data.Put(nx-1, j, 0, l, input->Get(nx-3, j-1, 0, l) > 0 ? 1 : 0);
    }
    for (int i = 1; i < nx - 1; ++i) {
      _data.Put(i, 0   , 0, l, input->Get(i-1, 0   , 0, l) > 0 ? 1 : 0);
      _data.Put(i, ny-1, 0, l, input->Get(i-1, ny-3, 0, l) > 0 ? 1 : 0);
    }

    // Copy original image into interior.
    for (int j = 1; j < ny - 1; ++j)
    for (int i = 1; i < nx - 1; ++i) {
      _data.Put(i, j, 0, l, input->Get(i-1, j-1, 0, l) > 0 ? 1 : 0);
    }
  }

  // Only need to use the first four offsets in the 2D case.
  _offsets.Initialize(nx, ny, CONNECTIVITY_6);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void CityBlockDistanceTransform<VoxelType>::Initialize3D()
{
  const GenericImage<VoxelType> *input  = this->Input();

  // Make an array with an extra layer in each direction.
  ImageAttributes attr = input->Attributes();
  attr._x += 2;
  attr._y += 2;
  attr._z += 2;

  _data.Initialize(attr);

  const int nx = _data.X();
  const int ny = _data.Y();
  const int nz = _data.Z();
  const int nt = _data.T();

  for (int l = 0; l < nt; ++l) {

    // Outer layers are duplicates of outer layer of input image.
    for (int j = 1; j < ny - 1; ++j)
    for (int i = 1; i < nx - 1; ++i) {
      _data.Put(i, j, 0   , l, input->Get(i-1, j-1, 0   , l) > 0 ? 1 : 0);
      _data.Put(i, j, nz-1, l, input->Get(i-1, j-1, nz-3, l) > 0 ? 1 : 0);
  	}

    for (int k = 1; k < nz - 1; ++k)
    for (int i = 1; i < nx - 1; ++i) {
      _data.Put(i, 0   , k, l, input->Get(i-1, 0   , k-1, l) > 0 ? 1 : 0);
      _data.Put(i, ny-1, k, l, input->Get(i-1, ny-3, k-1, l) > 0 ? 1 : 0);
    }

    for (int k = 1; k < nz - 1; ++k)
    for (int j = 1; j < ny - 1; ++j) {
      _data.Put(0   , j, k, l, input->Get(0   , j-1, k-1, l) > 0 ? 1 : 0);
      _data.Put(nx-1, j, k, l, input->Get(nx-3, j-1, k-1, l) > 0 ? 1 : 0);
    }

    // Copy original image into interior.
    for (int k = 1; k < nz - 1; ++k)
    for (int j = 1; j < ny - 1; ++j)
    for (int i = 1; i < nx - 1; ++i) {
      _data.Put(i, j, k, l, input->Get(i-1, j-1, k-1, l) > 0 ? 1 : 0);
    }
  }

  // Offsets for searching in the neighbourhood of a voxel.
  _offsets.Initialize(nx, ny, CONNECTIVITY_6);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void CityBlockDistanceTransform<VoxelType>::Initialize()
{
  ImageToImage<VoxelType>::Initialize();

  const int nx = this->Input()->X();
  const int ny = this->Input()->Y();
  const int nz = this->Input()->Z();

  if (nx == 1 || ny == 1 || nz == 1) {
    this->Initialize2D();
  } else {
    this->Initialize3D();
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void CityBlockDistanceTransform<VoxelType>::Run2D()
{
  VoxelType  val;
  int        borderVoxelCount = 0, objectVoxelCount = 0;
  GreyPixel *ptr2current, *ptr2offset;

  GenericImage<VoxelType> * const output = this->Output();

  const int nx = _data.X();
  const int ny = _data.Y();
  const int nz = _data.Z();
  const int nt = _data.T();

  if (nz != 1) {
    cerr << "CityBlockDistanceTransform<VoxelType>::Run2D(): Expect one slice only" << endl;
    exit(1);
  }

  // Storage for recording the indices of border voxels.
  Array<int> borderIndices(nx * ny);

  for (int l = 0; l < nt; ++l) {
    do {

      // Increment the distance for all current object voxels.
      objectVoxelCount = 0;

      for (int j = 1; j < ny - 1; ++j)
      for (int i = 1; i < nx - 1; ++i) {
        if (_data(i, j, 0, l) > 0) {
          ++objectVoxelCount;
          val = output->Get(i-1, j-1, 0, l);
          output->Put(i-1, j-1, 0, l, 1 + val);
        }
      }

      if (objectVoxelCount == 0){
        cerr << "irtkCityBlockDistanceTransform<VoxelType>::Run2D() : No object voxels." << endl;
        break;
      }

      // Remove the border from the current object(s).
      borderVoxelCount = 0;
      GreyPixel * const ptr2start = _data.Data(0, 0, 0, l);
      for (int j = 1; j < ny - 1; ++j)
      for (int i = 1; i < nx - 1; ++i) {
        if (_data(i, j, 0, l) > 0) {
          // Object voxel. Check face neighbourhood for background.
          ptr2current = _data.Data(i, j, 0, l);
          // 2D case: Only need to check first 4 neighbours in offset list.
          for (int m = 0; m < 4; ++m) {
            ptr2offset = ptr2current + _offsets(m);
            if (*ptr2offset < 1) {
              // Border voxel.
              borderIndices[borderVoxelCount] = static_cast<int>(ptr2current - ptr2start);
              ++borderVoxelCount;
              break;
            }
          }
        }
      }

      // Remove the current border voxels.
      for (int n = 0; n < borderVoxelCount; ++n){
	      ptr2current = ptr2start + borderIndices[n];
	      (*ptr2current) = 0;
      }

      // Update count of object voxels.
      objectVoxelCount -= borderVoxelCount;

    } while (objectVoxelCount > 0);
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void CityBlockDistanceTransform<VoxelType>::Run3D()
{
  VoxelType  val;
  int        borderVoxelCount = 0, objectVoxelCount = 0;
  GreyPixel *ptr2current, *ptr2offset;

  GenericImage<VoxelType> * const output = this->Output();

  const int nx = _data.X();
  const int ny = _data.Y();
  const int nz = _data.Z();
  const int nt = _data.T();

  // Storage for recording the indices of border voxels.
  Array<int> borderIndices(nx * ny * nz);

  for (int l = 0; l < nt; ++l) {

    do {

      // Increment the distance for all current object voxels.
      objectVoxelCount = 0;
      for (int k = 1; k < nz - 1; ++k)
      for (int j = 1; j < ny - 1; ++j)
      for (int i = 1; i < nx - 1; ++i) {
        if (_data(i, j, k, l) > 0) {
          ++objectVoxelCount;
          val = output->Get(i-1, j-1, k-1, l);
          output->Put(i-1, j-1, k-1, l, 1 + val);
        }
      }

      if (objectVoxelCount == 0){
        cerr << "irtkCityBlockDistanceTransform<VoxelType>::Run3D() : No object voxels." << endl;
        break;
      }


      // Remove the border from the current object.
      borderVoxelCount = 0;
      GreyPixel * const ptr2start = _data.Data(0, 0, 0, l);
      for (int k = 1; k < nz - 1; ++k)
      for (int j = 1; j < ny - 1; ++j)
      for (int i = 1; i < nx - 1; ++i) {
        if (_data(i, j, k, l) > 0) {
          // Object voxel. Check face neighbourhood for background (6-neighbourhood).
          ptr2current = _data.Data(i, j, k, l);
          for (int m = 0; m < 6; ++m) {
            ptr2offset = ptr2current + _offsets(m);
            if (*ptr2offset < 1) {
              // Border voxel.
              borderIndices[borderVoxelCount] = static_cast<int>(ptr2current - ptr2start);
              ++borderVoxelCount;
              break;
            }
          }
        }
      }

      // Remove the border voxels.
      for (int n = 0; n < borderVoxelCount; ++n){
        ptr2current = ptr2start + borderIndices[n];
        (*ptr2current) = 0;
      }

      // Update count of object voxels.
      objectVoxelCount -= borderVoxelCount;

    } while (objectVoxelCount > 0);

  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void CityBlockDistanceTransform<VoxelType>::Run()
{
  this->Initialize();

  const int nx = this->Input()->X();
  const int ny = this->Input()->Y();
  const int nz = this->Input()->Z();

  if (nx == 1 || ny == 1 || nz == 1) {
    this->Run2D();
  } else {
    this->Run3D();
  }

  this->Finalize();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void CityBlockDistanceTransform<VoxelType>::Finalize()
{
  // If image is 2D, may need to undo a dimension flip.
  if (_data.Z() == 1){
    GenericImage<VoxelType> *input  = const_cast<GenericImage<VoxelType> *>(this->Input());
    GenericImage<VoxelType> *output = this->Output();
    if (_flipType == FlipXZ) {
      input ->FlipXZ(0);
      output->FlipXZ(0);
    } else if (_flipType == FlipYZ) {
      input ->FlipYZ(0);
      output->FlipYZ(0);
    }
  }

  _data.Clear();
  _flipType = FlipNone;

  ImageToImage<VoxelType>::Finalize();
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class CityBlockDistanceTransform<RealPixel>;
template class CityBlockDistanceTransform<GreyPixel>;


} // namespace mirtk

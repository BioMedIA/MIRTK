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

// TODO: Compare transformed image to expected result
#include "gtest/gtest.h"

#include "mirtk/RegisteredImage.h"

#include "mirtk/Profiling.h"
#include "mirtk/GenericImage.h"
#include "mirtk/RigidTransformation.h"
#include "mirtk/BSplineFreeFormTransformation3D.h"
#include "mirtk/MultiLevelFreeFormTransformation.h"

namespace mirtk {


// ===========================================================================
// Helper
// ===========================================================================

// ---------------------------------------------------------------------------
/// Fill test image such that each voxel's value corresponds to its index
template <class TVoxel>
void fill_test_image(GenericImage<TVoxel> &image)
{
  const int numvox = image.GetNumberOfVoxels();
  if (static_cast<double>(numvox - 1) > voxel_limits<TVoxel>::max()) {
    cerr << "fill_test_image: Overflow!" << endl;
    exit(1);
  }
  for (int idx = 0; idx < image.GetNumberOfVoxels(); ++idx)
  {
    image.Put(idx, static_cast<TVoxel>(idx));
  }
}

// ===========================================================================
// Tests
// ===========================================================================

// ---------------------------------------------------------------------------
TEST(RegisteredImage, IdentityTransformation)
{

}

// ---------------------------------------------------------------------------
TEST(RegisteredImage, GlobalAndLocalTransformation)
{
  // Allocate images
  ImageAttributes attr(64, 64, 32);
  RegisteredImage::InputImageType image(attr);
  RegisteredImage source;
  // Initialize untransformed source image such that voxels can be distinguished
  fill_test_image(image);
  // Prepare global transformation which translates image one voxel in y
  RigidTransformation global;
  global.PutTranslationY(attr._dy);
  // Prepare FFD which translates image two voxels in x
  BSplineFreeFormTransformation3D local(attr, attr._dx, attr._dy, attr._dz);
  for (int xdof = 0; xdof < local.NumberOfCPs(); ++xdof) {
    local.Put(xdof, 2.0 * attr._dx);
  }
  // Initialize moving image
  MultiLevelFreeFormTransformation mffd(global);
  mffd.PushLocalTransformation(&local);
  //BinaryImage mask(attr); mask = true; source.PutMask(&mask);
  source.InputImage(&image);
  source.Transformation(&mffd);
  source.Initialize(attr);
  // Update moving image
  {
    MIRTK_START_TIMING();
    source.Update(true, false, false, true);
    MIRTK_END_TIMING("RegisteredImage::Update(false)");
  }
  {
    MIRTK_START_TIMING();
    source.Update(true, true, false, true);
    MIRTK_END_TIMING("RegisteredImage::Update(true)");
  }
  mffd.PopLocalTransformation();
}


} // namespace mirtk

// ===========================================================================
// Main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

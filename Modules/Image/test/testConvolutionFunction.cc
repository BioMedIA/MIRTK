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

#include "gtest/gtest.h"

#include "mirtk/GenericImage.h"
#include "mirtk/ConvolutionFunction.h"
#include "mirtk/ScalarGaussian.h"
#include "mirtk/ScalarFunctionToImage.h"

using namespace mirtk;
using namespace mirtk::ConvolutionFunction;

static const char *test_file   = NULL;
static const char *result_file = NULL;
static double      sigma       = 1.0;
static double      bgvalue     = 0.0;

// ===========================================================================
// Kernels
// ===========================================================================

GenericImage<RealPixel> Gaussian(double sigma, double dx)
{
  sigma = sigma / dx;
  ScalarGaussian                   gaussian(sigma, 1, 1, 0, 0, 0);
  GenericImage<RealPixel>          kernel(2 * int(round(4.0 * sigma)) + 1, 1, 1);
  ScalarFunctionToImage<RealPixel> generator;
  generator.Input (&gaussian);
  generator.Output(&kernel);
  generator.Run();
  RealPixel sum = .0;
  for (int i = 0; i < kernel.X(); ++i) sum += kernel(i, 0, 0, 0);
  for (int i = 0; i < kernel.X(); ++i) kernel(i, 0, 0, 0) /= sum;
  return kernel;
}

// ===========================================================================
// Tests
// ===========================================================================

// ---------------------------------------------------------------------------
TEST(ConvolutionFunction, ConvolveMirroredForeground)
{
  GenericImage<RealPixel> image (test_file);
  GenericImage<RealPixel> output(image.GetImageAttributes());
  GenericImage<RealPixel> xkernel = Gaussian(sigma, image.GetXSize());
  GenericImage<RealPixel> ykernel = Gaussian(sigma, image.GetYSize());
  GenericImage<RealPixel> zkernel = Gaussian(sigma, image.GetZSize());
  image.PutBackgroundValueAsDouble(bgvalue, true);
  ConvolveMirroredForegroundInX<RealPixel> xconv(&image, &xkernel);
  ParallelForEachVoxel(image.GetImageAttributes(), image, output, xconv);
  ConvolveMirroredForegroundInY<RealPixel> yconv(&image, &ykernel);
  ParallelForEachVoxel(image.GetImageAttributes(), output, image, yconv);
  ConvolveMirroredForegroundInZ<RealPixel> zconv(&image, &zkernel);
  ParallelForEachVoxel(image.GetImageAttributes(), image, output, zconv);
  output.Write(result_file);
}

// ===========================================================================
// Main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  if (argc < 3 || argc > 5) {
    cerr << "usage: " << argv[0] << " <input> <output> [sigma] [background]" << endl;
    exit(1);
  }
  test_file   = argv[1];
  result_file = argv[2];
  if (argc > 3) sigma   = atof(argv[3]);
  if (argc > 4) bgvalue = atof(argv[4]);
  return RUN_ALL_TESTS();
}

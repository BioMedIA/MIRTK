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

// TODO: Compare downsampled images to expected result instead of just writing files.

#include "gtest/gtest.h"

#include "mirtk/GenericImage.h"
#include "mirtk/Downsampling.h"
#include "mirtk/GaussianBlurring.h"
#include "mirtk/Resampling.h"
#include "mirtk/ScalarGaussian.h"
//#include "mirtk/ScalarGaussianDx.h"
#include "mirtk/ScalarFunctionToImage.h"
#include "mirtk/ConvolutionFunction.h"
#include "mirtk/InterpolateImageFunction.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace mirtk;
using namespace mirtk::ConvolutionFunction;

static const char *test_file   = NULL;
static const char *result_file = NULL;
static double      sigma       = 0.0;
static double      bgvalue     = 0.0;
static int         factor      = 2;

// ===========================================================================
// Tests
// ===========================================================================

// ---------------------------------------------------------------------------
TEST(Downsampling, GaussianPyramid1)
{
  char buffer[128];
  const int               levels = 4;
  GenericImage<RealPixel> pyramid[levels];
  pyramid[0].Read(test_file);
  pyramid[0].PutBackgroundValueAsDouble(bgvalue);
  for (int l = 1; l < levels; ++l) {
    GaussianBlurring<RealPixel> blur(.7 * pyramid[l-1].GetXSize(),
                                     .7 * pyramid[l-1].GetYSize(),
                                     .7 * pyramid[l-1].GetZSize());
    blur.Input (&pyramid[l-1]);
    blur.Output(&pyramid[l]);
    blur.Run();
    Resampling<RealPixel> resample(pyramid[l-1].GetX() / 2,
                                   pyramid[l-1].GetY() / 2,
                                   pyramid[l-1].GetZ() / 2);
    resample.Interpolator(InterpolateImageFunction::New(Interpolation_Linear, &pyramid[l]));
    resample.Input (&pyramid[l]);
    resample.Output(&pyramid[l]);
    resample.Run();
    snprintf(buffer, 128, "pyramid1_level_%d", l);
    delete resample.Interpolator();
    pyramid[l].Write(buffer);
  }
}

// ---------------------------------------------------------------------------
//TEST(Downsampling, GaussianPyramid2)
//{
//  char buffer[128];
//  const int               levels = 4;
//  GenericImage<RealPixel> pyramid[levels];
//  pyramid[0].Read(test_file);
//  pyramid[0].PutBackgroundValueAsDouble(bgvalue, true);
//  ScalarGaussian   gaussian (0.7, 1.0, 1.0, .0, .0, .0);
//  ScalarGaussianDx gaussianD(0.7, 1.0, 1.0, .0, .0, .0);
//  GenericImage<RealPixel> kernel (2*round(4.*0.7)+1, 1, 1);
//  GenericImage<RealPixel> kernelD(2*round(4.*0.7)+1, 1, 1);
//  ScalarFunctionToImage<RealPixel> gaussianSource;
//  gaussianSource.Input (&gaussian);
//  gaussianSource.Output(&kernel);
//  gaussianSource.Run();
//  gaussianSource.Input (&gaussianD);
//  gaussianSource.Output(&kernelD);
//  gaussianSource.Run();
//  RealPixel sum = .0;
//  for (int i = 0; i < kernel.X(); ++i) sum += kernel(i, 0, 0);
//  for (int i = 0; i < kernel.X(); ++i) kernel(i, 0, 0) /= sum;
//  for (int l = 1; l < levels; ++l) {
//    ImageAttributes attr = pyramid[l-1].GetImageAttributes();
//    GenericImage<RealPixel> tmp;
//    tmp = pyramid[l] = pyramid[l-1];
//    ConvolveExtendedForegroundInX<RealPixel> convx(&pyramid[l-1], &kernel);
//    ParallelForEachVoxel(attr, pyramid[l-1], pyramid[l], convx);
//    ConvolveExtendedForegroundInY<RealPixel> convy(&pyramid[l-1], &kernel);
//    ParallelForEachVoxel(attr, pyramid[l], tmp, convy);
//    ConvolveExtendedForegroundInZ<RealPixel> convz(&pyramid[l-1], &kernel);
//    ParallelForEachVoxel(attr, tmp, pyramid[l], convz);
//    Resampling<RealPixel> resample(pyramid[l-1].X() / 2,
//                                   pyramid[l-1].Y() / 2,
//                                   pyramid[l-1].Z() / 2);
//    resample.Interpolator(InterpolateImageFunction::New(Interpolation_Linear, &pyramid[l]));
//    resample.Input (&pyramid[l]);
//    resample.Output(&pyramid[l]);
//    resample.Run();
//    snprintf(buffer, 128, "pyramid2_level_%d", l);
//    delete resample.Interpolator();
//    pyramid[l].Write(buffer);
//  }
//}

/*
// ---------------------------------------------------------------------------
TEST(Downsampling, Run)
{
  GreyImage      image(test_file);
  ScalarGaussian gaussian(sigma, 1.0, 1.0, .0, .0, .0);
  image.PutBackgroundValueAsDouble(bgvalue, true);
  Downsampling<GreyPixel> downsampler(factor);
  double x = -0.5, y = -0.5, z = -0.5;
  image.ImageToWorld(x, y, z);
  cout << "Start of image region before: (" << x << ", " << y << ", " << z << ")" << endl;
  downsampler.Kernel(&gaussian, round(4.0 * sigma));
  downsampler.Input (&image);
  downsampler.Output(&image);
  downsampler.Run();
  x = -0.5, y = -0.5, z = -0.5;
  image.ImageToWorld(x, y, z);
  cout << "Start of image region after:  (" << x << ", " << y << ", " << z << ")" << endl;
  image.Write(result_file);
}
*/

// ===========================================================================
// Main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  if (argc < 3 || argc > 6) {
    cerr << "usage: " << argv[0] << " <input> <output> [factor] [sigma] [background]" << endl;
    exit(1);
  }
  test_file   = argv[1];
  result_file = argv[2];
  if (argc > 3) factor  = atoi(argv[3]);
  if (argc > 4) sigma   = atof(argv[4]);
  if (argc > 5) bgvalue = atof(argv[5]);
  return RUN_ALL_TESTS();
}

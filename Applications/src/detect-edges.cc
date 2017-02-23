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

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/IOConfig.h"

#include "mirtk/GenericImage.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/Resampling.h"
#include "mirtk/GaussianBlurring.h"
#include "mirtk/GaussianBlurringWithPadding.h"
#include "mirtk/GradientImageFilter.h"
#include "mirtk/ConvolutionFunction.h"

#include "mirtk/LinearInterpolateImageFunction.hxx"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <input> <output> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Convolves the input image with an edge detection operator." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input    Input greyscale image." << endl;
  cout << "  output   Output edge map." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -sigma <float>   Standard deviation of Gaussian blurring filter. (default: 0)" << endl;
  cout << "  -central         Use central differences. (default)" << endl;
  cout << "  -sobel           Use Sobel operator." << endl;
  cout << "  -prewitt         Use Prewitt operator." << endl;
  cout << "  -laplace         Use Laplace operator." << endl;
  cout << "  -differential    Use differential edge detection." << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
enum EdgeOperator {
  GRADIENT,
  SOBEL,
  PREWITT,
  LAPLACE,
  DIFFERENTIAL
};

// -----------------------------------------------------------------------------
/// Resample image to isotropic voxel size
void MakeIsotropic(RealImage &image)
{
  const double ds = min(min(image.GetXSize(), image.GetYSize()), image.GetZSize());
  GenericLinearInterpolateImageFunction<ByteImage> interp;
  Resampling<RealPixel> resampling(ds, ds, ds);
  resampling.Input(&image);
  resampling.Output(&image);
  resampling.Interpolator(&interp);
  resampling.Run();
}

// -----------------------------------------------------------------------------
RealImage ComputeCentralDifferences(const RealImage &image)
{
  RealImage output;
  GradientImageFilter<RealPixel> gradient;
  gradient.UseOrientation(false);
  gradient.UseVoxelSize(false);
  gradient.Input (&image);
  gradient.Output(&output);
  gradient.Run();
  return output;
}

// -----------------------------------------------------------------------------
/// Voxel function for parallel convolution of image with discrete Laplace operator
struct LaplaceOperator : public VoxelFunction
{
  const RealImage *_Image;
  LaplaceOperator(const RealImage *image = NULL)
  :
    _Image(image)
  {}

  void operator ()(int i, int j, int k, int, RealPixel *out) const
  {
    const RealImage &f = *_Image;
    if (0 < i && i < f.X()-1 && 0 < j && j < f.Y()-1 && 0 < k && k < f.Z()-1) {
      *out = .25 * (  f(i-1, j, k) + f(i+1, j, k)
                    + f(i, j-1, k) + f(i, j+1, k)
                    + f(i, j, k-1) + f(i, j, k+1) - 6.0 * f(i, j, k));
    }
  }
};

// -----------------------------------------------------------------------------
/// Convolve intensity image by the finite difference Laplacian kernel
RealImage ApplyLaplaceOperator(const RealImage &image)
{
  RealImage output(image);
  ParallelForEachVoxel(LaplaceOperator(&image), output.Attributes(), output);
  return output;
}

// -----------------------------------------------------------------------------
/// http://en.wikipedia.org/wiki/Edge_detection#Differential_edge_detection
struct DifferentialEdgeFunction : public VoxelFunction
{
  const RealImage *_Image;
  DifferentialEdgeFunction(const RealImage *image = NULL)
  :
    _Image(image)
  {}

  void operator ()(int i, int j, int k, int, RealPixel *out) const
  {
    const RealImage &f = *_Image;
    if (0 < i && i < f.X()-1 && 0 < j && j < f.Y()-1 && 0 < k && k < f.Z()-1) {
      // First order derivatives
      double fx = .5 * (f(i+1, j, k) - f(i-1, j, k));
      double fy = .5 * (f(i, j+1, k) - f(i, j-1, k));
      double fz = .5 * (f(i, j, k+1) - f(i, j, k-1));
      // Second order derivatives
      double fijk = f(i, j, k);
      double fxx = .25 * (f(i+1, j,   k  ) - fijk             - fijk             + f(i-1, j,   k  ));
      double fxy = .25 * (f(i+1, j+1, k  ) - f(i-1, j+1, k  ) - f(i+1, j-1, k  ) + f(i-1, j-1, k  ));
      double fxz = .25 * (f(i+1, j,   k+1) - f(i-1, j,   k+1) - f(i+1, j,   k-1) + f(i-1, j,   k-1));
      double fyy = .25 * (f(i,   j+1, k  ) - fijk             - fijk             + f(i,   j-1, k  ));
      double fyz = .25 * (f(i,   j+1, k+1) - f(i,   j-1, k+1) - f(i,   j+1, k-1) + f(i,   j-1, k-1));
      double fzz = .25 * (f(i,   j,   k+1) - fijk             - fijk             + f(i,   j,   k-1));
      // Second order derivative in the direction of the image gradient
      *out =        fx * fx * fxx + fy * fy * fyy + fz * fz * fzz
           + 2.0 * (fx * fy * fxy + fx * fz * fxz + fy * fz * fyz);
//      if (fabs(*out) > 1000.) *out = copysign(1000.0, *out);
    }
  }
};

// -----------------------------------------------------------------------------
/// http://en.wikipedia.org/wiki/Edge_detection#Differential_edge_detection
RealImage ComputeDifferentialEdgeFunction(const RealImage &image)
{
  RealImage output(image.Attributes());
  ParallelForEachVoxel(DifferentialEdgeFunction(&image), output.Attributes(), output);
  return output;
}

// -----------------------------------------------------------------------------
RealImage ApplySobelOperator(const RealImage &image)
{
  typedef ConvolutionFunction::ConvolveInX<RealPixel> ConvX;
  typedef ConvolutionFunction::ConvolveInY<RealPixel> ConvY;
  typedef ConvolutionFunction::ConvolveInZ<RealPixel> ConvZ;

  RealPixel h[3] = { 1.0, 2.0, 1.0};
  RealPixel g[3] = {-1.0,  .0, 1.0};

  const ImageAttributes &attr = image.Attributes();
  RealImage gx(attr), gy(attr), gz(attr), tmp(attr);

  ParallelForEachVoxel(ConvX(&image, g, 3, false), attr, image, gx);
  ParallelForEachVoxel(ConvY(&gx,    h, 3, false), attr, gx,    tmp);
  ParallelForEachVoxel(ConvZ(&gx,    h, 3, false), attr, tmp,   gx);

  ParallelForEachVoxel(ConvX(&image, h, 3, false), attr, image, gy);
  ParallelForEachVoxel(ConvY(&gy,    g, 3, false), attr, gy,   tmp);
  ParallelForEachVoxel(ConvZ(&gy,    h, 3, false), attr, tmp,   gy);

  ParallelForEachVoxel(ConvX(&image, h, 3, false), attr, image, gz);
  ParallelForEachVoxel(ConvY(&gz,    h, 3, false), attr, gz,   tmp);
  ParallelForEachVoxel(ConvZ(&gz,    g, 3, false), attr, tmp,   gz);

  RealImage gm(attr);
  for (int i = 0; i < image.NumberOfVoxels(); ++i) {
    gm(i) = sqrt(gx(i) * gx(i) + gy(i) * gy(i) + gz(i) * gz(i));
  }
  return gm;
}

// -----------------------------------------------------------------------------
RealImage ApplyPrewittOperator(const RealImage &image)
{
  cerr << "Prewitt operator not implemented" << endl;
  exit(1);
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  REQUIRES_POSARGS(2);

  const char *image_name  = POSARG(1);
  const char *output_name = POSARG(2);

  double       sigma   = 0.;
  RealPixel    padding = NaN;
  EdgeOperator op      = GRADIENT;

  for (ALL_OPTIONS) {
    if      (OPTION("-sigma")) sigma = atof(ARGUMENT);
    else if (OPTION("-gradient") || OPTION("-central")) op = GRADIENT;
    else if (OPTION("-sobel"))        op = SOBEL;
    else if (OPTION("-prewitt"))      op = PREWITT;
    else if (OPTION("-laplace"))      op = LAPLACE;
    else if (OPTION("-differential")) op = DIFFERENTIAL;
    else if (OPTION("-padding")) PARSE_ARGUMENT(padding);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read image
  InitializeIOLibrary();
  RealImage image(image_name);

  // Blur image
  if (sigma > .0) {
    if (IsNaN(padding)) {
      GaussianBlurring<RealPixel> blur(sigma);
      blur.Input (&image);
      blur.Output(&image);
      blur.Run();
    } else {
      GaussianBlurringWithPadding<RealPixel> blur(sigma, padding);
      blur.Input (&image);
      blur.Output(&image);
      blur.Run();
    }
  }

  // Convolve image
  RealImage output(image.Attributes());
  switch (op) {
    case GRADIENT:     output = ComputeCentralDifferences(image);       break;
    case SOBEL:        output = ApplySobelOperator(image);              break;
    case PREWITT:      output = ApplyPrewittOperator(image);            break;
    case LAPLACE:      output = ApplyLaplaceOperator(image);            break;
    case DIFFERENTIAL: output = ComputeDifferentialEdgeFunction(image); break;
  }

  // Discard
  if (!IsNaN(padding)) {
    int n = 0;
    bool outside;
    const int rl = (image.T() > 1 ? 1 : 0);
    const int rk = (image.Z() > 1 ? 1 : 0);
    const int rj = (image.Y() > 1 ? 1 : 0);
    const int ri = (image.X() > 1 ? 1 : 0);
    for (int l = 0; l < image.T(); ++l)
    for (int k = 0; k < image.Z(); ++k)
    for (int j = 0; j < image.Y(); ++j)
    for (int i = 0; i < image.X(); ++i) {
      outside = (image(i, j, k, l) <= padding);
      if (!outside) {
        for (int nl = l-rl; nl <= l+rl; ++nl)
        for (int nk = k-rk; nk <= k+rk; ++nk)
        for (int nj = j-rj; nj <= j+rj; ++nj)
        for (int ni = i-ri; ni <= i+ri; ++ni) {
          if (image.IsInside(ni, nj, nk, nl) && image(ni, nj, nk, nl) <= padding) {
            outside = true;
            nl = l + 2 * rl;
            nk = k + 2 * rk;
            nj = j + 2 * rj;
            break;
          }
        }
      }
      if (outside) {
        n++;
        output(i, j, k, l) = 0.;
      }
    }
    cout << "No. of padded outside values = " << n << endl;
  }

  // Write edge map
  output.Write(output_name);
  return 0;
}

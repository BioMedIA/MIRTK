/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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
#include "mirtk/InterpolateImageFunction.h"
#include "mirtk/Resampling.h"
#include "mirtk/ResamplingWithPadding.h"
#include "mirtk/GaussianInterpolateImageFunction.h"
#include "mirtk/GaussianInterpolateImageFunction2D.h"
#include "mirtk/ConstExtrapolateImageFunction.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  const     Indent option_desc_indent(27, 1);
  const int max_line_length = 80;
  int       line_length;
  string    mode;

  cout << endl;
  cout << "Usage: " << name << " <input> <output> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Resamples an image on a lattice with specified voxel size." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input    Input image." << endl;
  cout << "  output   Resampled output image." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -imsize <nx> <ny> <nz>   New image size in number of voxels (default: input image size)" << endl;
  cout << "  -size <dx> <dy> <dz>     New voxel size in mm (default: 1 1 1)" << endl;
  cout << "  -isotropic <m>           Resample to isotropic voxel size: m*min(dx, min(dy, dz))" << endl;
  cout << "  -interp <mode>           Interpolation mode (case insensitive):" << endl;
  line_length = max_line_length;
  for (int i = 0; i < Interpolation_Last; ++i) {
    if (i > 0) {
      cout << ", ";
      line_length += 2;
    }
    mode = ToString(static_cast<InterpolationMode>(i));
    line_length += static_cast<int>(mode.length());
    if (line_length >= max_line_length) {
      cout << endl << option_desc_indent;
      line_length = option_desc_indent.Spaces() + static_cast<int>(mode.length());
    }
    cout << mode;
    if (i == Interpolation_Linear) cout << " (default)";
  }
  cout << endl;
  cout << "  -extrap <mode>           Extrapolation mode (case insensitive):" << endl;
  line_length = max_line_length;
  for (int i = 0; i < Extrapolation_Last; ++i) {
    if (i > 0) {
      cout << ", ";
      line_length += 2;
    }
    mode = ToString(static_cast<ExtrapolationMode>(i));
    line_length += static_cast<int>(mode.length());
    if (line_length >= max_line_length) {
      cout << endl << option_desc_indent;
      line_length = option_desc_indent.Spaces() + static_cast<int>(mode.length());
    }
    cout << mode;
    if (i == Extrapolation_Const) cout << " (default)";
  }
  cout << endl;
  cout << "  -outside <value>         Constant outside value extrapolation. (default: 0)" << endl;
  cout << "  -sigma <value>           Sigma value of Gaussian interpolator. (default: 1)" << endl;
  cout << "  -padding <value>         Background padding. (default: none)" << endl;
  cout << "                           Only linear interpolation for padding!" << endl;
  cout << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
template <class TVoxel>
bool Resample(BaseImage *image, InterpolateImageFunction *interp,
              int nx, int ny, int nz, double dx, double dy, double dz)
{
  GenericImage<TVoxel> *im = dynamic_cast<GenericImage<TVoxel> *>(image);
  if (im == nullptr) return false;
  Resampling<TVoxel> resampling(nx, ny, nz, dx, dy, dz);
  resampling.Input (im);
  resampling.Output(im);
  resampling.Interpolator(interp);
  resampling.Run();
  return true;
}

// -----------------------------------------------------------------------------
template <class TVoxel>
bool Resample(BaseImage *image, InterpolateImageFunction *interp,
              int nx, int ny, int nz, double dx, double dy, double dz,
              double padding_value)
{
  GenericImage<TVoxel> *im = dynamic_cast<GenericImage<TVoxel> *>(image);
  if (im == nullptr) return false;
  ResamplingWithPadding<TVoxel> resampling(nx, ny, nz, dx, dy, dz, static_cast<TVoxel>(padding_value));
  resampling.Input (im);
  resampling.Output(im);
  resampling.Interpolator(interp);
  resampling.Run();
  return true;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // Parse arguments
  EXPECTS_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  // Read image
  InitializeIOLibrary();
  UniquePtr<BaseImage> image(BaseImage::New(input_name));

  // Parse optional arguments
  InterpolationMode interpolation_mode  = Interpolation_NN;
  ExtrapolationMode extrapolation_mode  = Extrapolation_Default;
  double            interpolation_sigma = 1.0; // Interpolation_Gaussian parameter
  double            outside_value       = 0.0; // Extrapolation_Const parameter
  double            isotropic           = .0;
  double            padding_value       = NaN;
  double            dx                  = image->XSize();
  double            dy                  = image->YSize();
  double            dz                  = image->ZSize();

  int nx, ny, nz;
  nx = ny = nz = 0;

  for (ALL_OPTIONS) {
    if (OPTION("-imsize")) {
      PARSE_ARGUMENT(nx);
      PARSE_ARGUMENT(ny);
      PARSE_ARGUMENT(nz);
    }
    else if (OPTION("-size")) {
      PARSE_ARGUMENT(dx);
      PARSE_ARGUMENT(dy);
      PARSE_ARGUMENT(dz);
    }
    else if (OPTION("-isotropic")) {
      if (HAS_ARGUMENT) { 
        PARSE_ARGUMENT(isotropic);
      } else isotropic = 1.0;
    }
    else if (OPTION("-padding")) {
      PARSE_ARGUMENT(padding_value);
    }
    else if (OPTION("-outside")) {
      PARSE_ARGUMENT(outside_value);
      extrapolation_mode = Extrapolation_Const;
    }
    else if (OPTION("-sigma")) {
      PARSE_ARGUMENT(interpolation_sigma);
    }
    else if (OPTION("-interp")) {
      PARSE_ARGUMENT(interpolation_mode);
    }
    else if (OPTION("-extrap")) {
      PARSE_ARGUMENT(extrapolation_mode);
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Create interpolator
  UniquePtr<InterpolateImageFunction> interpolator;
  interpolator.reset(InterpolateImageFunction::New(interpolation_mode, extrapolation_mode, image.get()));
  interpolator->DefaultValue(outside_value);

  // TODO: The actual type of the templated GenericGaussianInterpolateImageFunction's is not known.
  //       Add virtual parameter setters which take string arguments to InterpolateImageFunction.
  if (strcmp(interpolator->NameOfClass(), "GaussianInterpolateImageFunction") == 0) {
    reinterpret_cast<GaussianInterpolateImageFunction *>(interpolator.get())->Sigma(interpolation_sigma);
  }
  if (strcmp(interpolator->NameOfClass(), "GaussianInterpolateImageFunction2D") == 0) {
    reinterpret_cast<GaussianInterpolateImageFunction2D *>(interpolator.get())->Sigma(interpolation_sigma);
  }
  if (interpolator->Extrapolator()) {
    if (strcmp(interpolator->Extrapolator()->NameOfClass(), "ConstExtrapolateImageFunction") == 0) {
      reinterpret_cast<ConstExtrapolateImageFunction *>(interpolator->Extrapolator())->DefaultValue(outside_value);
    }
  }

  // Make voxel size isotropic
  if (isotropic > .0) {
    image->GetPixelSize(&dx, &dy, &dz);
    dx = dy = dz = min(min(dx, dy), dz) * isotropic;
    if (verbose) cout << "Resampling image to isotropic voxel size: " << dx << " [mm]";
  } else {
    if (verbose) cout << "Resampling ... ";
  }
  if (verbose) cout.flush();

  // Resample image
  if (IsNaN(padding_value)) {
    switch (image->GetDataType()) {
      case MIRTK_VOXEL_UNSIGNED_CHAR: {
        Resample<unsigned char>(image.get(), interpolator.get(), nx, ny, nz, dx, dy, dz);
      } break;
      case MIRTK_VOXEL_SHORT: {
        Resample<short>(image.get(), interpolator.get(), nx, ny, nz, dx, dy, dz);
      } break;
      case MIRTK_VOXEL_UNSIGNED_SHORT: {
        Resample<unsigned short>(image.get(), interpolator.get(), nx, ny, nz, dx, dy, dz);
      } break;
      case MIRTK_VOXEL_INT: {
        Resample<int>(image.get(), interpolator.get(), nx, ny, nz, dx, dy, dz);
      } break;
      case MIRTK_VOXEL_FLOAT: {
        Resample<float>(image.get(), interpolator.get(), nx, ny, nz, dx, dy, dz);
      } break;
      case MIRTK_VOXEL_DOUBLE: {
        Resample<double>(image.get(), interpolator.get(), nx, ny, nz, dx, dy, dz);
      } break;
      default:
        FatalError("Unsupported voxel type: " << image->GetDataType());
    }
  } else {
    switch (image->GetDataType()) {
      case MIRTK_VOXEL_UNSIGNED_CHAR: {
        Resample<unsigned char>(image.get(), interpolator.get(), nx, ny, nz, dx, dy, dz, padding_value);
      } break;
      case MIRTK_VOXEL_SHORT: {
        Resample<short>(image.get(), interpolator.get(), nx, ny, nz, dx, dy, dz, padding_value);
      } break;
      case MIRTK_VOXEL_UNSIGNED_SHORT: {
        Resample<unsigned short>(image.get(), interpolator.get(), nx, ny, nz, dx, dy, dz, padding_value);
      } break;
      case MIRTK_VOXEL_INT: {
        Resample<int>(image.get(), interpolator.get(), nx, ny, nz, dx, dy, dz, padding_value);
      } break;
      case MIRTK_VOXEL_FLOAT: {
        Resample<float>(image.get(), interpolator.get(), nx, ny, nz, dx, dy, dz, padding_value);
      } break;
      case MIRTK_VOXEL_DOUBLE: {
        Resample<double>(image.get(), interpolator.get(), nx, ny, nz, dx, dy, dz, padding_value);
      } break;
      default:
        FatalError("Unsupported voxel type: " << image->GetDataType());
    }
  }
  if (verbose) cout << " done" << endl;

  // Write image
  image->Write(output_name);

  return 0;
}

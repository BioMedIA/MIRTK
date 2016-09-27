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

#include "mirtk/ImageReader.h"
#include "mirtk/IOConfig.h"
#include "mirtk/GaussianBlurring.h"
#include "mirtk/GaussianBlurring4D.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <input> <output> [<sigma>] [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  This program blurs an input image using a Gaussian filter with a standard" << endl;
  cout << "  deviation equal to the specified sigma value. The dimensions in which the" << endl;
  cout << "  filter recursively applied can be specified by the respective options below." << endl;
  cout << "  Note that more than one of these options can be specified and that each can" << endl;
  cout << "  also be given more than once. The blurring of the respective dimensions" << endl;
  cout << "  is then performed in the order of the given blurring options." << endl;
  cout << endl;
  cout << "  The output image is by default saved as short integers if the" << endl;
  cout << "  input image was of integral voxel type, otherwise as float image." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input    Input image. Multiple channels/frames are smoothed independently." << endl;
  cout << "  output   Blurred output image." << endl;
  cout << "  sigma    Standard deviation of Gaussian smoothing kernel. When negative, its" << endl;
  cout << "           absolute value is multiplied by the average voxel size. (default: 1)" << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -3D      Blur image in all spatial directions. (default)" << endl;
  cout << "  -4D      Blur image in all dimensions." << endl;
  cout << "  -x       Blur image in x dimension." << endl;
  cout << "  -y       Blur image in y dimension." << endl;
  cout << "  -z       Blur image in z dimension." << endl;
  cout << "  -short   Set data type of output to short integers." << endl;
  cout << "  -float   Set data type of output to floating point." << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
void ParseSigma(const string &arg, double &sigma)
{
  string value;
  string units = ValueUnits(arg, &value, "signed");
  if (!FromString(value, sigma)) {
    FatalError("Invalid sigma value argument: " << value);
  }
  if (units != "signed") {
    if (sigma < 0.) {
      FatalError("Sigma value cannot be negative when units are specified");
    }
    if (units == "vox") {
      sigma = -sigma;
    } else if (units != "mm") {
      FatalError("Invalid sigma value units: " << units);
    }
  }
}

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // Parse positional arguments
  REQUIRES_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);
  double      sigma       = 1.0;

  if (NUM_POSARGS == 3) {
    ParseSigma(string(POSARG(3)), sigma);
  } else if (NUM_POSARGS > 3) {
    PrintHelp(EXECNAME);
    cout << endl;
    FatalError("Invalid number of positional arguments");
  }

  // Determine input image data type
  InitializeIOLibrary();
  UniquePtr<ImageReader> reader(ImageReader::New(input_name));
  int data_type = reader->DataType();

  // Parse and discard options
  DISCARD_PARSED_OPTIONS();

  for (ALL_OPTIONS) {
    if      (OPTION("-short")) data_type = MIRTK_VOXEL_SHORT;
    else if (OPTION("-float")) data_type = MIRTK_VOXEL_FLOAT;
    else HANDLE_COMMON_OPTION();
  }

  if (sigma == 0.) {
    FatalError("Sigma value cannot be zero!");
  }

  // Read image
  RealImage input(input_name);

  // Parse blurring options
  bool default_blurring = true;
  for (ALL_OPTIONS) {
    if (OPTION("-3D")) {
      default_blurring = false;
      GaussianBlurring<RealPixel> blur(
        sigma < 0. ? abs(sigma) * (input.XSize() + input.YSize() + input.ZSize()) / 3. : sigma
      );
      blur.Input (&input);
      blur.Output(&input);
      blur.Run();
    } else if (OPTION("-4D")) {
      default_blurring = false;
      GaussianBlurring<RealPixel> blur(
        sigma < 0. ? abs(sigma) * (input.XSize() + input.YSize() + input.ZSize() + input.TSize()) / 4. : sigma
      );
      blur.Input (&input);
      blur.Output(&input);
      blur.Run();
    } else if (OPTION("-x") || OPTION("-X")) {
      default_blurring = false;
      GaussianBlurring<RealPixel> blur(sigma < 0. ? abs(sigma) * input.XSize() : sigma);
      blur.Input (&input);
      blur.Output(&input);
      blur.RunX();
    } else if (OPTION("-y") || OPTION("-Y")) {
      default_blurring = false;
      GaussianBlurring<RealPixel> blur(sigma < 0. ? abs(sigma) * input.YSize() : sigma);
      blur.Input (&input);
      blur.Output(&input);
      blur.RunY();
    } else if (OPTION("-z") || OPTION("-Z")) {
      default_blurring = false;
      GaussianBlurring<RealPixel> blur(sigma < 0. ? abs(sigma) * input.ZSize() : sigma);
      blur.Input (&input);
      blur.Output(&input);
      blur.RunZ();
    } else if (OPTION("-sigma")) {
      // change sigma for next blur operation
      ParseSigma(string(ARGUMENT), sigma);
    } else HANDLE_UNKNOWN_OPTION();
  }

  // Default blurring if no blurring option given
  if (default_blurring) {
    GaussianBlurring<RealPixel> blur(
      sigma < 0. ? abs(sigma) * (input.XSize() + input.YSize() + input.ZSize()) / 3. : sigma
    );
    blur.Input (&input);
    blur.Output(&input);
    blur.Run();
  }

  // Write image
  switch (data_type) {
    case MIRTK_VOXEL_FLOAT:
    case MIRTK_VOXEL_DOUBLE: {
      input.Write(output_name);
    } break;
    default:
      GreyImage output(input);
      output.Write(output_name);
  }

  return 0;
}

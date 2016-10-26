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
#include "mirtk/GaussianPyramidFilter.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <input> <output> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Downsamples an image using an iterative Gaussian pyramid filter.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Input intensity imaeg.\n";
  cout << "  output   Input image downsampled n times.\n";
  cout << "\n";
  cout << "Options:\n";
  cout << "  -iterations <n>   Number of iterations. (default: 1)\n";
  PrintStandardOptions(cout);
  cout << "\n";
  cout.flush();
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
BaseImage *Downsample(BaseImage *input)
{
  GaussianPyramidFilter<VoxelType> filter(0, 1);
  filter.Input (dynamic_cast<GenericImage<VoxelType> *>(input));
  filter.Output(new GenericImage<VoxelType>());
  filter.Run();
  return filter.Output();
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  REQUIRES_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  // Default parameters
  int iter = 1;

  // Parse optional arguments
  for (ALL_OPTIONS) {
    if (OPTION("-iterations") || OPTION("-iter")) iter = atoi(ARGUMENT);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read input image
  if (verbose) cout << "Reading image ... "; cout.flush();
  InitializeIOLibrary();
  UniquePtr<BaseImage> input(BaseImage::New(input_name));
  if (verbose) cout << "done." << endl;

  // Downsample image
  if (verbose) cout << "Downsampling image ... ";
  UniquePtr<BaseImage> output;
  for (int i = 0; i < iter; ++i) {
    if (i > 0) input.swap(output);
    switch (input->GetScalarType()) {
      case MIRTK_VOXEL_UNSIGNED_CHAR:  { output.reset(Downsample<unsigned char >(input.get())); break; }
      case MIRTK_VOXEL_SHORT:          { output.reset(Downsample<short         >(input.get())); break; }
      case MIRTK_VOXEL_UNSIGNED_SHORT: { output.reset(Downsample<unsigned short>(input.get())); break; }
      case MIRTK_VOXEL_INT:            { output.reset(Downsample<int           >(input.get())); break; }
      case MIRTK_VOXEL_UNSIGNED_INT:   { output.reset(Downsample<unsigned int  >(input.get())); break; }
      case MIRTK_VOXEL_FLOAT:          { output.reset(Downsample<float         >(input.get())); break; }
      case MIRTK_VOXEL_DOUBLE:         { output.reset(Downsample<double        >(input.get())); break; }
      default: {
        FatalError("Unsupported scalar type: " << input->GetScalarType());
        exit(1);
      }
    }
  }
  if (verbose) cout << "done." << endl;

  // Write output image
  output->Write(output_name);
  return 0;
}

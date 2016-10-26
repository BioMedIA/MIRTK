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
#include "mirtk/Erosion.h"

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
  cout << "  Erodes an input image by replacing a voxel's value by the minimum\n";
  cout << "  of the values of its neighboring voxels.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Input intensity/segmentation image.\n";
  cout << "  output   Eroded output image.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -iterations <n>     Number of iterations. (default: 1)\n";
  cout << "  -connectivity <n>   Type of voxel connectivity (4, 6, 18, or 26). (default: 26)\n";
  PrintStandardOptions(cout);
  cout << "\n";
  cout.flush();
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  EXPECTS_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  int              iterations   = 1;
  ConnectivityType connectivity = CONNECTIVITY_26;

  for (ALL_OPTIONS) {
    if (OPTION("-iterations") || OPTION("-iter")) {
      PARSE_ARGUMENT(iterations);
    }
    else if (OPTION("-connectivity") || OPTION("-neighbors") || OPTION("-number-of-neighbors")) {
      PARSE_ARGUMENT(connectivity);
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  InitializeIOLibrary();
  UniquePtr<BaseImage> image(BaseImage::New(input_name));

  if (verbose) cout << "Dilating ... ", cout.flush();
  switch (image->GetDataType()) {
    case MIRTK_VOXEL_BINARY:  Erode<BinaryPixel>(image.get(), iterations, connectivity); break;
    case MIRTK_VOXEL_GREY:    Erode<GreyPixel  >(image.get(), iterations, connectivity); break;
    case MIRTK_VOXEL_REAL:    Erode<RealPixel  >(image.get(), iterations, connectivity); break;
    default: {
      RealImage other(*image);
      Erode<RealPixel>(&other, iterations, connectivity);
      *image = other;
    } break;
  }
  if (verbose) cout << "done" << endl;

  image->Write(output_name);

  return 0;
}


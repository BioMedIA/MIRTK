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

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: convert <input> <output> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Converts an image from one voxel type to another.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -char|uchar|short|ushort|float|double   Output voxel type.\n";
  cout << "  -rescale <min> <max>                    Output minimum and maximum intensity.\n";
  PrintStandardOptions(cout);
  cout << "\n";
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

  int    voxel_type = MIRTK_VOXEL_UNKNOWN;
  double min_value  = numeric_limits<double>::quiet_NaN();
  double max_value  = numeric_limits<double>::quiet_NaN();

  for (ALL_OPTIONS) {
    if      (OPTION("-char"))   voxel_type = MIRTK_VOXEL_CHAR;
    else if (OPTION("-uchar"))  voxel_type = MIRTK_VOXEL_UNSIGNED_CHAR;
    else if (OPTION("-short"))  voxel_type = MIRTK_VOXEL_SHORT;
    else if (OPTION("-ushort")) voxel_type = MIRTK_VOXEL_UNSIGNED_SHORT;
    else if (OPTION("-float"))  voxel_type = MIRTK_VOXEL_FLOAT;
    else if (OPTION("-double")) voxel_type = MIRTK_VOXEL_DOUBLE;
    else if (OPTION("-binary")) voxel_type = MIRTK_VOXEL_BINARY;
    else if (OPTION("-byte"))   voxel_type = MIRTK_VOXEL_BYTE;
    else if (OPTION("-grey"))   voxel_type = MIRTK_VOXEL_GREY;
    else if (OPTION("-real"))   voxel_type = MIRTK_VOXEL_REAL;
    else if (OPTION("-rescale")) {
      PARSE_ARGUMENT(min_value);
      PARSE_ARGUMENT(max_value);
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read image
  InitializeIOLibrary();
  UniquePtr<BaseImage> input(BaseImage::New(input_name));
  if (voxel_type == MIRTK_VOXEL_UNKNOWN) {
    voxel_type = input->GetDataType();
  }

  // Scale image
  if (!IsNaN(min_value) || !IsNaN(max_value)) {
    if (IsNaN(min_value) || IsNaN(max_value)) {
      double vmin, vmax;
      input->GetMinMaxAsDouble(vmin, vmax);
      if (IsNaN(min_value)) min_value = vmin;
      if (IsNaN(max_value)) max_value = vmax;
    }
    if (min_value >= max_value) swap(min_value, max_value);
    if (voxel_type != MIRTK_VOXEL_FLOAT && voxel_type != MIRTK_VOXEL_DOUBLE) {
      input.reset(new GenericImage<double>(*input));
    }
    input->PutMinMaxAsDouble(min_value, max_value);
  }

  // Convert image
  UniquePtr<BaseImage> output(BaseImage::New(voxel_type));
  *output = *input;

  // Write image
  output->Write(output_name);

  return 0;
}

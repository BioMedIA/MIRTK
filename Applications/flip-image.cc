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

#include <mirtkCommon.h>
#include <mirtkOptions.h>

#include <mirtkIOConfig.h>
#include <mirtkBaseImage.h>

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
  cout << "  Swaps the two image dimensions at a time, in the order of the input\n";
  cout << "  options. Both, the image data and the coordinates of the image origin\n";
  cout << "  are swapped each time an input option is processed.\n";
  cout << "\n";
  cout << "Options:\n";
  cout << "  -xy, -yx   Swap x and y dimension.\n";
  cout << "  -xz, -zx   Swap x and z dimension.\n";
  cout << "  -xt, -tx   Swap x and t dimension.\n";
  cout << "  -yz, -zy   Swap y and z dimension.\n";
  cout << "  -yt, -ty   Swap y and t dimension.\n";
  cout << "  -zt, -tz   Swap z and t dimension.\n";
  PrintStandardOptions(cout);
  cout << "\n";
  cout.flush();
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  EXPECTS_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  InitializeIOLibrary();
  unique_ptr<BaseImage> image(BaseImage::New(input_name));

  for (ALL_OPTIONS) {
    if      (OPTION("-xy") || OPTION("-yx")) image->FlipXY(true);
    else if (OPTION("-xz") || OPTION("-zx")) image->FlipXZ(true);
    else if (OPTION("-xt") || OPTION("-tx")) image->FlipXT(true);
    else if (OPTION("-yz") || OPTION("-zy")) image->FlipYZ(true);
    else if (OPTION("-yt") || OPTION("-ty")) image->FlipYT(true);
    else if (OPTION("-zt") || OPTION("-tz")) image->FlipZT(true);
    else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
  }

  image->Write(output_name);

  return 0;
}

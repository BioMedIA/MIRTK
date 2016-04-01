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
  cout << endl;
  cout << "Usage: " << name << " <input> <output> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Applies a sequence of one or more spatial reflections to an image." << endl;
  cout << "  The reflection options are processed in the order given and changing" << endl;
  cout << "  the order can change the result." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -x         Reflect x axis." << endl;
  cout << "  -y         Reflect y axis." << endl;
  cout << "  -z         Reflect z axis." << endl;
  cout << "  -xy, -yx   Swap x and y axes." << endl;
  cout << "  -xz, -zx   Swap x and z axes." << endl;
  cout << "  -yz, -zy   Swap y and z axes." << endl;
  PrintStandardOptions(cout);
  cout << endl;
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

  InitializeIOLibrary();
  unique_ptr<BaseImage> image(BaseImage::New(input_name));

  for (ALL_OPTIONS) {
    if      (OPTION("-x")) image->ReflectX();
    else if (OPTION("-y")) image->ReflectY();
    else if (OPTION("-z")) image->ReflectZ();
    else if (OPTION("-xy") || OPTION("-yx")) image->FlipXY(false);
    else if (OPTION("-xz") || OPTION("-zx")) image->FlipXZ(false);
    else if (OPTION("-yz") || OPTION("-zy")) image->FlipYZ(false);
    else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
  }

  image->Write(output_name);

  return 0;
}

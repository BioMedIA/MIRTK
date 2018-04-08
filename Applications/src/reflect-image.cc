/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2017 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
 * Copyright 2013-2017 Andreas Schuh
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

// TODO: This command is obsolete; use flip-image instead.
//       The "mirtk.subprocess.path" function can map "reflect-image"
//       requests to "flip-image" command executions to maintain
//       backwards compatibility. This would include setting
//       "-axes off" for "-xy"... options that were removed here.

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
  cout << "  Applies a  one or more spatial reflections along an image axis.\n";
  cout << "  A more generic tool that can also be used to swap two axes is\n";
  cout << "  the flip-image command.\n";
  cout << "\n";
  cout << "Options:\n";
  cout << "  -x   Reflect x axis.\n";
  cout << "  -y   Reflect y axis.\n";
  cout << "  -z   Reflect z axis.\n";
  cout << "  -t   Reflect t axis.\n";
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
  UniquePtr<BaseImage> image(BaseImage::New(input_name));

  bool modify_axes = false;
  bool reflect_x   = false;
  bool reflect_y   = false;
  bool reflect_z   = false;
  bool reflect_t   = false;

  for (ALL_OPTIONS) {
    if (OPTION("-axes")) {
      modify_axes = true;
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(modify_axes);
      }
    }
    else if (OPTION("-noaxes")) {
      modify_axes = false;
    }
    else if (OPTION("-x")) reflect_x = true;
    else if (OPTION("-y")) reflect_y = true;
    else if (OPTION("-z")) reflect_z = true;
    else if (OPTION("-t")) reflect_t = true;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (reflect_x) image->ReflectX(modify_axes);
  if (reflect_y) image->ReflectY(modify_axes);
  if (reflect_z) image->ReflectZ(modify_axes);
  if (reflect_t) image->ReflectT(modify_axes);

  image->Write(output_name);

  return 0;
}

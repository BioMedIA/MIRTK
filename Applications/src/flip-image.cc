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
#include "mirtk/BaseImage.h"

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
  cout << "  are swapped each time an input option is processed. To swap the image\n";
  cout << "  axes instead of the origin, use the :option:`-axes` option before the\n";
  cout << "  swap options. This option should be used to reorder the image dimensions\n";
  cout << "  without changing the world coordinates of the voxels. When swapping of\n";
  cout << "  the coordinate axes is enabled, the coordinates of the image origin are\n";
  cout << "  kept the same, i.e., :option:`-origin` is ignored. The default behavior\n";
  cout << "  is to swap the image data and the coordinates of the image origin. This\n";
  cout << "  may in many cases not have the desired effect, but has been this way\n";
  cout << "  already for some time.\n";
  cout << "\n";
  cout << "Options:\n";
  cout << "  -axes [on|off]     Enable/disable swapping of the coordinate axes. (default: off)\n";
  cout << "  -noaxes            Disable swapping of the coordinate axes.\n";
  cout << "  -origin [on|off]   Enable/disable swapping of the origin coordinates. (default: on)\n";
  cout << "  -noorigin          Disable swapping of the origin coordinates.\n";
  cout << "  -xy, -yx           Swap x and y dimension.\n";
  cout << "  -xz, -zx           Swap x and z dimension.\n";
  cout << "  -xt, -tx           Swap x and t dimension.\n";
  cout << "  -yz, -zy           Swap y and z dimension.\n";
  cout << "  -yt, -ty           Swap y and t dimension.\n";
  cout << "  -zt, -tz           Swap z and t dimension.\n";
  PrintStandardOptions(cout);
  cout << endl;
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
  UniquePtr<BaseImage> image(BaseImage::New(input_name));
  bool modify_origin = true;
  bool modify_axes   = false;

  for (ALL_OPTIONS) {
    if (OPTION("-xy") || OPTION("-yx")) {
      if (modify_axes) {
        image->SwapXY(true);
      } else {
        image->FlipXY(modify_origin);
      }
    }
    else if (OPTION("-xz") || OPTION("-zx")) {
      if (modify_axes) {
        image->SwapXZ(true);
      } else {
        image->FlipXZ(modify_origin);
      }
    }
    else if (OPTION("-xt") || OPTION("-tx")) {
      if (modify_axes) {
        image->SwapXT(true);
      } else {
        image->FlipXT(modify_origin);
      }
    }
    else if (OPTION("-yz") || OPTION("-zy")) {
      if (modify_axes) {
        image->SwapYZ(true);
      } else {
        image->FlipYZ(modify_origin);
      }
    }
    else if (OPTION("-yt") || OPTION("-ty")) {
      if (modify_axes) {
        image->SwapYT(true);
      } else {
        image->FlipYT(modify_origin);
      }
    }
    else if (OPTION("-zt") || OPTION("-tz")) {
      if (modify_axes) {
        image->SwapZT(true);
      } else {
        image->FlipZT(modify_origin);
      }
    }
    else if (OPTION("-axes")) {
      modify_axes = true;
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(modify_axes);
      }
    }
    else if (OPTION("-noaxes")) {
      modify_axes = false;
    }
    else if (OPTION("-origin")) {
      modify_origin = true;
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(modify_origin);
      }
    }
    else if (OPTION("-noorigin")) {
      modify_origin = false;
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  image->Write(output_name);

  return 0;
}

/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2017 Imperial College London
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
#include "mirtk/GenericImage.h"
#include "mirtk/HistogramMatching.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// ------------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <target> <source> <output> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Matches the intensity distribution of the source image to match the\n";
  cout << "  distribution of the target image using a piecewise linear function [1].\n";
  cout << "\n";
  cout << "  [1] Nyul, Udupa, and Zhang, \"New variants of a method of MRI scale standardization\",\n";
  cout << "      IEEE TMI 19(2), pp. 143-150, 2000, http://dx.doi.org/10.1109/42.836373.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  target   Reference image.\n";
  cout << "  source   Input image whose histogram should match the reference histogram.\n";
  cout << "  output   Output image with histogram matching the reference histogram.\n";
  cout << "\n";
  cout << "Options:\n";
  cout << "  -Tp <value>\n";
  cout << "    Target padding value. (default: NaN)\n";
  cout << "  -Sp <value>\n";
  cout << "    Source padding value. (default: NaN)\n";
  cout << "  -dtype short|int|float|double\n";
  cout << "    Data type of output image. (default: input data type)\n";
  PrintCommonOptions(cout);
  cout << endl;
}


// =============================================================================
// Main
// =============================================================================

// ------------------------------------------------------------------------------
int main(int argc, char **argv)
{
  REQUIRES_POSARGS(3);
  InitializeIOLibrary();

  const char *target_name = POSARG(1);
  const char *source_name = POSARG(2);
  const char *output_name = POSARG(3);

  double source_padding = NaN;
  double target_padding = NaN;
  ImageDataType dtype = MIRTK_VOXEL_UNKNOWN;

  for (ALL_OPTIONS) {
    if      (OPTION("-Tp")) PARSE_ARGUMENT(target_padding);
    else if (OPTION("-Sp")) PARSE_ARGUMENT(source_padding);
    else if (OPTION("-dtype")) PARSE_ARGUMENT(dtype);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (verbose) cout << "Reading target image ... ", cout.flush();
  RealImage target(target_name);
  target.PutBackgroundValueAsDouble(target_padding, true);
  if (verbose) cout << "done" << endl;

  if (verbose) cout << "Reading source image ... ", cout.flush();
  RealImage source;
  {
    UniquePtr<BaseImage> tmp(BaseImage::New(source_name));
    if (dtype == MIRTK_VOXEL_UNKNOWN) {
      dtype = static_cast<ImageDataType>(tmp->GetDataType());
    }
    source = *tmp;
  }
  source.PutBackgroundValueAsDouble(source_padding, true);
  if (verbose) cout << "done" << endl;

  if (verbose) cout << "Matching histogram...";
  HistogramMatching<RealPixel> match;
  match.Input(&source);
  match.Reference(&target);
  match.Output(&source);
  match.Run();
  if (verbose) cout << " done" << endl;

  if (verbose) cout << "Writing output image...";
  if (dtype != source.GetDataType()) {
    UniquePtr<BaseImage> output(BaseImage::New(dtype));
    *output = source;
    output->Write(output_name);
  } else {
    source.Write(output_name);
  }
  if (verbose) cout << " done" << endl;

  return 0;
}

/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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


// ===========================================================================
// Help
// ===========================================================================

// ---------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <image> <weights> <pbmap> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Computes the probability of a voxel to belong to given segment." << endl;
  PrintCommonOptions(cout);
  cout << endl;
}


// ===========================================================================
// Main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc, char **argv)
{
  EXPECTS_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  InitializeIOLibrary();
  RealImage image(input_name);
  const int nv = image.NumberOfVoxels();

  Array<double> mean;
  Array<double> var;

  for (ALL_OPTIONS) {
    if (OPTION("-pbmap")) {
      RealImage pbmap(ARGUMENT);
      RealPixel v, w, vsum = 0., wsum = 0.;
      for (int idx = 0; idx < nv; ++idx) {
        v = image(idx);
        w = pbmap(idx);
        vsum += w * v;
        wsum += w;
      }
      RealPixel mu = vsum / wsum;
      vsum = 0.;
      for (int idx = 0; idx < nv; ++idx) {
        v = image(idx) - mu;
        w = pbmap(idx);
        vsum += w * v * v;
      }
      mean.push_back(mu);
      var .push_back(vsum / wsum);
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  char   label;
  double prob, max_prob;

  GenericImage<char> labels(image.Attributes());
  for (int idx = 0; idx < nv; ++idx) {
    max_prob = 0.;
    label    = 0;
    for (size_t i = 0; i < mean.size(); ++i) {
      prob = exp(-.5 * pow(image(idx) - mean[i], 2) / var[i]);
      if (prob > max_prob) {
        max_prob = prob;
        label    = i+1;
      }
    }
    labels(idx) = label;
  }

  labels.Write(output_name);
  return 0;
}

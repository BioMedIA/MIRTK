/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
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

#include "mirtk/Options.h"
#include "mirtk/IOConfig.h"

#include "mirtk/BaseImage.h"
#include "mirtk/GenericImage.h"
#include "mirtk/GaussianBlurring.h"

using namespace mirtk;


void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <cine> <ED output> <ES output>" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  This program detects ED and ES phases from a cine image sets." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  cine        Input cine image sets." << endl;
  cout << "  ED output   ED phase image." << endl;
  cout << "  ES output   ES phase image." << endl;
  cout << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

int main(int argc, char* argv[] )
{
  EXPECTS_POSARGS(1);
  const char *cine_name = POSARG(1);
  const char *ed_name = nullptr;
  const char *es_name = nullptr;
  int esphase, frames;
  short cine_max, cine_min, cinedis;
  double *similarity, *smoothsimilarity, dif;
  // Check command line
  InitializeIOLibrary();

  // Parse source and target images
  for (ALL_OPTIONS) {
    if (OPTION("--output-ed") ed_name = ARGUMENT;
    else if (OPTION("--output-es") es_name = ARGUMENT;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (!ed_name && !es_name) {
    FatalError("At least one of --output-ed and --output-es required!");
  }

  // Create images
  GreyImage cine(cine_name);
  GreyImage blurred(cine.GetImageAttributes());

  GaussianBlurring<GreyPixel> gaussianBlurring(2);
  gaussianBlurring.Input (&cine);
  gaussianBlurring.Output(&blurred);
  gaussianBlurring.Run();

  ImageAttributes attr = cine.GetImageAttributes();
  Array<double> similarity(attr._t, 0.0);
  Array<double> smoothsimilarity(attr._t, 0.0);
  frames = attr._t;
  attr._t = 1;

  // Create similarity
  blurred.GetMinMax(&cine_min,&cine_max);
  cinedis = cine_max - cine_min;
  for (int i = 0; i < frames; ++i) {
    similarity[i] = 0;
    smoothsimilarity[i] = 0;
  }
  // Evaluate similarity
  for (int t = 0; t < frames; ++t) {
    for (int k = 0; k < cine.GetZ(); ++k)
    for (int j = 0; j < cine.GetY(); ++j)
    for (int i = 0; i < cine.GetX(); ++i) {
      dif = (blurred.GetAsDouble(i,j,k,t) - blurred.GetAsDouble(i,j,k,0)) / cinedis;
      similarity[t] += dif*dif;
    }
    if (verbose) {
      cout << "similarity : " << similarity[t] << endl;
    }
  }
  for (i = 0; i < frames; i++) {
    similarity[i] = sqrt(similarity[i]);
  }
  // Smooth similarity
  for (i = 1; i < frames - 1; i++) {
    smoothsimilarity[i] = (similarity[i-1] + similarity[i] + similarity[i+1]) / 3;
  }
  // Find min similarity
  dif = 0;
  for (i = 0; i < frames; i++) {
    if (dif < smoothsimilarity[i]) {
      dif = smoothsimilarity[i];
      esphase = i;
    }
  }
  if (verbose) {
    cout << "ES phase is: " << esphase << endl;
  }

  // Output
  if (ed_name) {
    cine.GetFrame(0).Write(ed_name);
  }
  if (es_name) {
    cine.GetFrame(es_phase).Write(es_name);
  }
  delete []smoothsimilarity;
  delete []similarity;
}

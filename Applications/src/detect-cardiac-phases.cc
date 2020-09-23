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
  const char *cine_name = NULL;
  const char *out_ED_name = NULL;
  const char *out_ES_name = NULL;
  int esphase, frames;
  short cine_max, cine_min, cinedis;
  double *similarity, *smoothsimilarity, dif;
  // Check command line
  REQUIRES_POSARGS(3);
  InitializeIOLibrary();

  // Parse source and target images
  cine_name = POSARG(1);
  out_ED_name = POSARG(2);
  out_ES_name = POSARG(3);

  // Create images
  GreyImage cine;
  cine.Read(cine_name);
  GreyImage blured;
  blured.Initialize(cine.GetImageAttributes());

  GaussianBlurring<GreyPixel> gaussianBlurring(2);
  gaussianBlurring.Input (&cine);
  gaussianBlurring.Output(&blured);
  gaussianBlurring.Run();

  ImageAttributes atr = cine.GetImageAttributes();
  similarity = new double[atr._t];
  smoothsimilarity = new double[atr._t];
  frames = atr._t;
  atr._t = 1;

  GreyImage out_ED(atr);
  GreyImage out_ES(atr);
  out_ED = cine.GetFrame(0);
  // Create similarity
  blured.GetMinMax(&cine_min,&cine_max);
  cinedis = cine_max - cine_min;
  for(int i = 0; i < frames; ++i){
      similarity[i] = 0;
      smoothsimilarity[i] = 0;
  }
  // Evaluate similarity
  for (int t = 0; t < frames; ++t) {
      for (int k = 0; k < cine.GetZ(); ++k) {
          for (int j = 0; j< cine.GetY(); ++j) {
              for (int i = 0; i<cine.GetX(); ++i) {
                 dif = (blured.GetAsDouble(i,j,k,t) - blured.GetAsDouble(i,j,k,0))/cinedis;
                 similarity[t] += dif*dif;
              }
          }
      }
      cout << "similarity : " << similarity[t] << endl;
  }
  for(i = 0; i < frames; i++)
      similarity[i] = sqrt(similarity[i]);
  // Smooth similarity
  for(i = 1; i < frames - 1; i++)
      smoothsimilarity[i] = (similarity[i-1] + similarity[i] + similarity[i+1])/3;
  // Find min similarity
  dif = 0;
  for(i = 0; i < frames; i++){
      if(dif < smoothsimilarity[i]){
          dif = smoothsimilarity[i];
          esphase = i;
      }
  }

  cout << "ES phase is: " << esphase << endl;

  out_ES = cine.GetFrame(esphase);
  // Output
  out_ED.Write(out_ED_name);
  out_ES.Write(out_ES_name);
  delete []smoothsimilarity;
  delete []similarity;
}
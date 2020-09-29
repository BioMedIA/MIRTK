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

#include "mirtk/ImageReader.h"
#include "mirtk/IOConfig.h"
#include "mirtk/GenericImage.h"

using namespace mirtk;


void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <image> <output> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  This program contrasts an image based on its mean pixel values and standard deviation." << endl;
  cout << "  output = clamp(input, min=max(mean - a * std, min_value), max=min(mean + a * std, max_value))." << endl;
  cout << "  By default, a = 3, min_value = MIN_GREY, max_value = MAX_GREY " << endl;

  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input    Input image." << endl;
  cout << "  output   Contrasted output image." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -std-factor <value>    Double. Standard deviation multiplication factor" << endl;
  cout << "  -min-value  <value>    Double. Min pixel value" << endl;
  cout << "  -max-value  <value>    Double. Max pixel value" << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

int main(int argc, char *argv[])
{
  // Parse positional arguments
  EXPECTS_POSARGS(2);

  // Parse and discard options
  DISCARD_PARSED_OPTIONS();
  double min_value = MIN_GREY;
  double max_value = MAX_GREY;
  double std_factor = 3.0;
  for (ALL_OPTIONS) {
    if (OPTION("-min-value")) {
        PARSE_ARGUMENT(min_value);
    } else if (OPTION("-max-value")){
        PARSE_ARGUMENT(max_value);
    } else if (OPTION("-std-factor")){
        PARSE_ARGUMENT(std_factor);
    }
    else {
        HANDLE_COMMON_OR_UNKNOWN_OPTION();
    }
  }

  InitializeIOLibrary();
  RealImage image(POSARG(1));
  char *output_name = POSARG(2);

  double average = image.Mean();
  RealPixel std = image.GetSD();

  for(int t = 0; t < image.GetT(); t++){
    for(int z = 0; z < image.GetZ(); z++){
      for(int y = 0; y < image.GetY(); y++){
        for(int x = 0; x < image.GetX(); x++){
          if(image.GetAsDouble(x, y, z, t) > average + std_factor * std){
            image.PutAsDouble(x, y, z, t, average + std_factor * std);
          }else if(image.GetAsDouble(x, y, z, t) < average - std_factor * std){
            image.PutAsDouble(x, y, z, t, average - std_factor * std);
          }
          if(image.GetAsDouble(x, y, z, t) < min_value)
            image.PutAsDouble(x, y, z, t, min_value);
          if(image.GetAsDouble(x, y, z, t) > max_value)
            image.PutAsDouble(x, y, z, t, max_value);
        }
      }
    }
  }

  image.Write(output_name);
}
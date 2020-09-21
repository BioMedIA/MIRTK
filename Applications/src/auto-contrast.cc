/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$
Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/ImageReader.h"
#include "mirtk/IOConfig.h"
#include "mirtk/GenericImage.h"

using namespace mirtk;

char *output_name = NULL;


void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <image> <output> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  This program auto contrast an image." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input    Input image." << endl;
  cout << "  output   Contrasted output image." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -Tp <value>        Target padding value all pixel < value = value" << endl;
  cout << "  -Tp2 <maxvalue>    Target padding maximum value all pixel > maxvalue = maxvalue" << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

int main(int argc, char *argv[])
{
  // Parse positional arguments
  REQUIRES_POSARGS(2);
  // Determine input image data type
  InitializeIOLibrary();

  int x, y, z, t, ok;
  double average,padding1,padding2;
  RealPixel std;
  RealImage *image;
  RealImage *output;

  padding1 = MIN_GREY;
  padding2 = MAX_GREY;

  image = new RealImage(POSARG(1));
  output = new RealImage(*image);
  output_name = POSARG(2);
  // Parse and discard options
  DISCARD_PARSED_OPTIONS();

  for (ALL_OPTIONS) {
    if (OPTION("-Tp")) {
        PARSE_ARGUMENT(padding1);
    } else if (OPTION("-Tp2")){
        PARSE_ARGUMENT(padding2);
    } else {
        HANDLE_COMMON_OR_UNKNOWN_OPTION();
    }
  }

  average = image->GetAverage();
  std = image->GetSD();

  for(t = 0; t < output->GetT(); t++){
      for(z = 0; z < output->GetZ(); z++){
          for(y = 0; y < output->GetY(); y++){
              for(x = 0; x < output->GetX(); x++){
                  if(output->GetAsDouble(x,y,z,t) > average + 3.0*std){
                      output->PutAsDouble(x,y,z,t,average + 3.0*std);
                  }else if(output->GetAsDouble(x,y,z,t) < average - 3.0*std){
                      output->PutAsDouble(x,y,z,t,average - 3.0*std);
                  }
                  if(output->GetAsDouble(x,y,z,t) < padding1)
                      output->PutAsDouble(x,y,z,t,padding1);
                  if(output->GetAsDouble(x,y,z,t) > padding2)
                      output->PutAsDouble(x,y,z,t,padding2);
              }
          }
      }
  }

  output->Write(output_name);

  delete image;
  delete output;
}
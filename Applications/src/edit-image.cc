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

#include "mirtk/BaseImage.h"
#include "mirtk/IOConfig.h"

#ifdef HAVE_MIRTK_Transformation
  #include "mirtk/Transformation.h"
  #include "mirtk/HomogeneousTransformation.h"
  #include "mirtk/MultiLevelTransformation.h"
#endif

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
  cout << "  Modifies the attributes of an image stored in the header." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input    Input image." << endl;
  cout << "  output   Output image." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -size <dx> <dy> <dz>   Voxel size   (in mm)" << endl;
  cout << "  -tsize <dt>            Voxel size   (in ms)" << endl;
  cout << "  -origin <x> <y> <z>    Image spatial origin (in mm)" << endl;
  cout << "  -torigin <t>           Image temporal origin (in ms)" << endl;
  cout << "  -orientation <x1> <x2> <x3>  <y1> <y2> <y3>  <z1> <z2> <z3>" << endl;
  cout << "                         Image orientation." << endl;
  cout << endl;
  cout << "  -copy-size <image>                      Copy voxel size." << endl;
  cout << "  -copy-origin <image>                    Copy origin." << endl;
  cout << "  -copy-orientation <image>               Copy orientation." << endl;
  cout << "  -copy-origin-orientation <image>        Alias for :option:`-copy-origin` :option:`-copy-orientation`." << endl;
  cout << "  -copy-origin-orientation-size <image>   Alias for :option:`-copy-origin` :option:`-copy-orientation` :option:`-copy-size`." << endl;
  cout << endl;
  cout << "  -reset           Set orientation, origin, and affine transformation matrix to default." << endl;
  cout << "  -reset-dof       Set affine transformation matrix to default." << endl;
#ifdef HAVE_MIRTK_Transformation
  cout << "  -dofin <file>    Apply transformation to axis, spacing and origin information" << endl;
  cout << "                   in the header. Note that any shearing that is present is" << endl;
  cout << "                   stored as additional affine transformation (c.f. -putdof)." << endl;
  cout << "  -putdof <file>   Store affine transformation in image header (NIfTI only)." << endl;
#endif // HAVE_MIRTK_Transformation
  cout << endl;
  cout << "  -swapxy   Swap the x and y axis vectors." << endl;
  cout << "  -swapxz   Swap the x and z axis vectors." << endl;
  cout << "  -swapyz   Swap the y and z axis vectors." << endl;
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

  double origin[4];
  double xsize, ysize, zsize, tsize;
  double xaxis[3], yaxis[3], zaxis[3];

  for (ALL_OPTIONS) {
    if (OPTION("-size")) {
      PARSE_ARGUMENT(xsize);
      PARSE_ARGUMENT(ysize);
      PARSE_ARGUMENT(zsize);
      image->PutPixelSize(xsize, ysize, zsize);
    }
    else if (OPTION("-tsize")) {
      PARSE_ARGUMENT(tsize);
      image->PutTSize(tsize);
    }
    else if (OPTION("-origin")) {
      PARSE_ARGUMENT(origin[0]);
      PARSE_ARGUMENT(origin[1]);
      PARSE_ARGUMENT(origin[2]);
      image->PutOrigin(origin[0], origin[1], origin[2]);
    }
    else if (OPTION("-torigin")) {
      PARSE_ARGUMENT(origin[3]);
      image->PutTOrigin(origin[3]);
    }
    else if (OPTION("-orientation")) {
      PARSE_ARGUMENT(xaxis[0]);
      PARSE_ARGUMENT(xaxis[1]);
      PARSE_ARGUMENT(xaxis[2]);
      PARSE_ARGUMENT(yaxis[0]);
      PARSE_ARGUMENT(yaxis[1]);
      PARSE_ARGUMENT(yaxis[2]);
      PARSE_ARGUMENT(zaxis[0]);
      PARSE_ARGUMENT(zaxis[1]);
      PARSE_ARGUMENT(zaxis[2]);
      image->PutOrientation(xaxis, yaxis, zaxis);
    }
    else if (OPTION("-copy-size")) {
      GreyImage target(ARGUMENT);
      target.GetPixelSize(&xsize, &ysize, &zsize);
      image->PutPixelSize(xsize, ysize, zsize);
    }
    else if (OPTION("-copy-origin")) {
      GreyImage target(ARGUMENT);
      target.GetOrigin(origin[0], origin[1], origin[2], origin[3]);
      image->PutOrigin(origin[0], origin[1], origin[2], origin[3]);
    }
    else if (OPTION("-copy-orientation")) {
      GreyImage target(ARGUMENT);
      target.GetOrientation(xaxis, yaxis, zaxis);
      image->PutOrientation(xaxis, yaxis, zaxis);
    }
    else if (OPTION("-copy-origin-orientation")) {
      GreyImage target(ARGUMENT);
      target.GetOrigin(origin[0], origin[1], origin[2], origin[3]);
      image->PutOrigin(origin[0], origin[1], origin[2], origin[3]);
      target.GetOrientation(xaxis, yaxis, zaxis);
      image->PutOrientation(xaxis, yaxis, zaxis);
    }
    else if (OPTION("-copy-origin-orientation-size")) {
      GreyImage target(ARGUMENT);
      target.GetPixelSize(&xsize, &ysize, &zsize);
      image->PutPixelSize(xsize, ysize, zsize);
      target.GetOrientation(xaxis, yaxis, zaxis);
      image->PutOrientation(xaxis, yaxis, zaxis);
      target.GetOrigin(origin[0], origin[1], origin[2], origin[3]);
      image->PutOrigin(origin[0], origin[1], origin[2], origin[3]);
    }
    #ifdef HAVE_MIRTK_Transformation
    else if (OPTION("-dofin") || OPTION("-dofin_i") || OPTION("-putdof") || OPTION("-putdof_i")) {
      Matrix m;
      UniquePtr<Transformation> dof(Transformation::New(ARGUMENT));
      HomogeneousTransformation *lin  = dynamic_cast<HomogeneousTransformation *>(dof.get());
      MultiLevelTransformation  *mffd = dynamic_cast<MultiLevelTransformation  *>(dof.get());
      if      (lin)  m = lin->GetMatrix();
      else if (mffd) m = mffd->GetGlobalTransformation()->GetMatrix();
      else {
        FatalError("Option -dofin requires a rigid, similarity, affine or MFFD transformation as argument");
      }
      if      (strcmp(OPTNAME, "-dofin"    ) == 0) image->PutAffineMatrix(m,          true);
      else if (strcmp(OPTNAME, "-dofin_i"  ) == 0) image->PutAffineMatrix(m.Invert(), true);
      else if (strcmp(OPTNAME, "-putdof"   ) == 0) image->PutAffineMatrix(m,          false);
      else if (strcmp(OPTNAME, "-putdof_i" ) == 0) image->PutAffineMatrix(m.Invert(), false);
    }
    #endif // HAVE_MIRTK_Transformation
    else if (OPTION("-reset")) {
      ImageAttributes defaultAttr;
      image->PutOrientation(defaultAttr._xaxis,defaultAttr._yaxis,defaultAttr._zaxis);
      image->PutOrigin(defaultAttr._xorigin,defaultAttr._yorigin,defaultAttr._zorigin);
      image->ResetAffineMatrix();
    }
    else if (OPTION("-reset-dof")) {
      image->ResetAffineMatrix();
    }
    else if (OPTION("-swapxy") || OPTION("-swapyx")) {
      image->GetOrientation(xaxis, yaxis, zaxis);
      image->PutOrientation(yaxis, xaxis, zaxis);
    }
    else if (OPTION("-swapxz") || OPTION("-swapzx")) {
      image->GetOrientation(xaxis, yaxis, zaxis);
      image->PutOrientation(zaxis, yaxis, xaxis);
    }
    else if (OPTION("-swapyz") || OPTION("-swapzy")) {
      image->GetOrientation(xaxis, yaxis, zaxis);
      image->PutOrientation(xaxis, zaxis, yaxis);
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (verbose) image->Print();
  image->Write(output_name);

  return 0;
}

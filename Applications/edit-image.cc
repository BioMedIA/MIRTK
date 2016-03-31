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

#include <mirtkCommon.h>
#include <mirtkOptions.h>

#include <mirtkBaseImage.h>
#include <mirtkIOConfig.h>

#ifdef HAVE_MIRTK_Transformation
  #include <mirtkTransformation.h>
  #include <mirtkHomogeneousTransformation.h>
  #include <mirtkMultiLevelTransformation.h>
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
  unique_ptr<BaseImage> image(BaseImage::New(input_name));

  double origin[4];
  double xsize, ysize, zsize;
  double xaxis[3], yaxis[3], zaxis[3];

  for (ALL_OPTIONS) {
    if (OPTION("-size")) {
      xsize = atof(ARGUMENT);
      ysize = atof(ARGUMENT);
      zsize = atof(ARGUMENT);
      image->PutPixelSize(xsize, ysize, zsize);
    }
    else if (OPTION("-tsize")) {
      image->PutTSize(atof(ARGUMENT));
    }
    else if (OPTION("-origin")) {
      origin[0] = atof(argv[1]);
      origin[1] = atof(argv[1]);
      origin[2] = atof(argv[1]);
      image->PutOrigin(origin[0], origin[1], origin[2]);
    }
    else if (OPTION("-torigin")) {
      origin[3] = atof(argv[1]);
      image->GetOrigin(origin[0], origin[1], origin[2]);
      image->PutOrigin(origin[0], origin[1], origin[2], origin[3]);
    }
    else if (OPTION("-orientation")) {
      xaxis[0] = atof(argv[1]);
      xaxis[1] = atof(argv[1]);
      xaxis[2] = atof(argv[1]);
      yaxis[0] = atof(argv[1]);
      yaxis[1] = atof(argv[1]);
      yaxis[2] = atof(argv[1]);
      zaxis[0] = atof(argv[1]);
      zaxis[1] = atof(argv[1]);
      zaxis[2] = atof(argv[1]);
      image->PutOrientation(xaxis, yaxis, zaxis);
    }
    else if (OPTION("-copy-size")) {
      GreyImage target(ARGUMENT);
      target.GetPixelSize(&xsize, &ysize, &zsize);
      image->PutPixelSize(xsize, ysize, zsize);
    }
    else if (OPTION("-copy-origin")) {
      GreyImage target(ARGUMENT);
      target.GetOrigin(origin[0], origin[1], origin[2]);
      image->PutOrigin(origin[0], origin[1], origin[2]);
    }
    else if (OPTION("-copy-orientation")) {
      GreyImage target(ARGUMENT);
      target.GetOrientation(xaxis, yaxis, zaxis);
      image->PutOrientation(xaxis, yaxis, zaxis);
    }
    else if (OPTION("-copy-origin-orientation")) {
      GreyImage target(ARGUMENT);
      target.GetOrigin(origin[0], origin[1], origin[2]);
      image->PutOrigin(origin[0], origin[1], origin[2]);
      target.GetOrientation(xaxis, yaxis, zaxis);
      image->PutOrientation(xaxis, yaxis, zaxis);
    }
    else if (OPTION("-copy-origin-orientation-size")) {
      GreyImage target(ARGUMENT);
      target.GetPixelSize(&xsize, &ysize, &zsize);
      image->PutPixelSize(xsize, ysize, zsize);
      target.GetOrientation(xaxis, yaxis, zaxis);
      image->PutOrientation(xaxis, yaxis, zaxis);
      target.GetOrigin(origin[0], origin[1], origin[2]);
      image->PutOrigin(origin[0], origin[1], origin[2]);
    }
    #ifdef HAVE_MIRTK_Transformation
    else if (OPTION("-dofin") || OPTION("-dofin_i") || OPTION("-putdof") || OPTION("-putdof_i")) {
      Matrix m;
      unique_ptr<Transformation> dof(Transformation::New(ARGUMENT));
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
    else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
  }

  if (verbose) image->Print();
  image->Write(output_name);

  return 0;
}

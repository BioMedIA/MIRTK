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
  #include "mirtk/Matrix.h"
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
  cout << "\n";
  cout << "Usage: " << name << " <input> <output> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Modifies the attributes of an image stored in the header.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Input image.\n";
  cout << "  output   Output image.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -spacing, -voxel-size, -size <dx> [<dy> [<dz> [<dt>]]]\n";
  cout << "      Voxel size (dx, dy, dz in mm, dt in ms)\n";
  cout << "  -dx <dx>\n";
  cout << "      Spatial voxel size in x dimension.\n";
  cout << "  -dy <dy>\n";
  cout << "      Spatial voxel size in y dimension.\n";
  cout << "  -dz <dz>\n";
  cout << "      Spatial voxel size in z dimension.\n";
  cout << "  -dt, tsize <dt>\n";
  cout << "      Temporal voxel size (in ms)\n";
  cout << "  -origin <x> <y> [<z> [<t>]]\n";
  cout << "      Image spatial origin (in mm)\n";
  cout << "  -torigin <t>\n";
  cout << "      Temporal image origin (in ms)\n";
  cout << "  -orientation <x1> <x2> <x3>  <y1> <y2> <y3>  <z1> <z2> <z3>\n";
  cout << "      Image orientation. The axes direction vectors are normalized to unit length.\n";
  cout << "\n";
  cout << "  -copy-spacing, -copy-voxel-size, -copy-size <image>\n";
  cout << "      Copy voxel size.\n";
  cout << "  -copy-origin <image>\n";
  cout << "      Copy origin.\n";
  cout << "  -copy-orientation <image>\n";
  cout << "      Copy orientation.\n";
  cout << "  -copy-origin-orientation <image>\n";
  cout << "      Alias for :option:`-copy-origin` :option:`-copy-orientation`.\n";
  cout << "  -copy-origin-orientation-spacing, -copy-origin-orientation-voxel-size, -copy-origin-orientation-size <image>\n";
  cout << "      Alias for :option:`-copy-origin` :option:`-copy-orientation` :option:`-copy-spacing`.\n";
  cout << "  -reset\n";
  cout << "      Set orientation, origin, and affine transformation matrix to default.\n";
  cout << "  -reset-dof\n";
  cout << "      Set affine transformation matrix to default.\n";
#ifdef HAVE_MIRTK_Transformation
  cout << "  -dofin <file>\n";
  cout << "      Apply transformation to axis, spacing and origin information\n";
  cout << "      in the header. Note that any shearing that is present is\n";
  cout << "      stored as additional affine transformation (c.f. :option:`-putdof`).\n";
  cout << "  -dofin_i <fil>\n";
  cout << "      Same as :option:`-dofin` but using the inverse transformation.\n";
  cout << "  -putdof <file>\n";
  cout << "      Store affine transformation in image header (NIfTI only).\n";
  cout << "  -putdof_i <fil>\n";
  cout << "      Same as :option:`-putdof` but using the inverse transformation.\n";
  cout << "  -dofout <file>\n";
  cout << "      Save transformation which maps the world coordinates of the\n";
  cout << "      output image to the world coordinates of the input image.\n";
  cout << "      Applying this output transformation to the output image\n";
  cout << "      using the :option:`-dofin` option restores the previous\n";
  cout << "      image origin, orientation, and voxel size.\n";
#endif // HAVE_MIRTK_Transformation
  cout << "\n";
  cout << "  -swapxy\n";
  cout << "      Swap the x and y axis vectors.\n";
  cout << "  -swapxz\n";
  cout << "      Swap the x and z axis vectors.\n";
  cout << "  -swapyz\n";
  cout << "      Swap the y and z axis vectors.\n";
  PrintStandardOptions(cout);
  cout << endl;
}

// -----------------------------------------------------------------------------
double Normalize(double v[3])
{
  double l = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] /= l, v[1] /= l, v[2] /= l;
  return l;
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

  #ifdef HAVE_MIRTK_Transformation
    ImageAttributes input_attr = image->Attributes();
    input_attr._smat.Ident();
    const Matrix input_mat = input_attr.GetImageToWorldMatrix();
    const char *dofout_name = nullptr;
  #endif

  for (ALL_OPTIONS) {
    if (OPTION("-spacing") || OPTION("-voxel-size") || OPTION("-size")) {
      image->GetPixelSize(xsize, ysize, zsize, tsize);
      PARSE_ARGUMENT(xsize);
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(ysize);
        if (HAS_ARGUMENT) PARSE_ARGUMENT(zsize);
      } else {
        ysize = zsize = xsize;
      }
      if (HAS_ARGUMENT) PARSE_ARGUMENT(tsize);
      image->PutPixelSize(xsize, ysize, zsize, tsize);
    }
    else if (OPTION("-dx")) {
      image->GetPixelSize(xsize, ysize, zsize);
      PARSE_ARGUMENT(xsize);
      image->PutPixelSize(xsize, ysize, zsize);
    }
    else if (OPTION("-dy")) {
      image->GetPixelSize(xsize, ysize, zsize);
      PARSE_ARGUMENT(ysize);
      image->PutPixelSize(xsize, ysize, zsize);
    }
    else if (OPTION("-dz")) {
      image->GetPixelSize(xsize, ysize, zsize);
      PARSE_ARGUMENT(zsize);
      image->PutPixelSize(xsize, ysize, zsize);
    }
    else if (OPTION("-dt") || OPTION("-tsize")) {
      PARSE_ARGUMENT(tsize);
      image->PutTSize(tsize);
    }
    else if (OPTION("-origin")) {
      image->GetOrigin(origin[0], origin[1], origin[2], origin[3]);
      PARSE_ARGUMENT(origin[0]);
      PARSE_ARGUMENT(origin[1]);
      if (HAS_ARGUMENT) PARSE_ARGUMENT(origin[2]);
      if (HAS_ARGUMENT) PARSE_ARGUMENT(origin[3]);
      image->PutOrigin(origin[0], origin[1], origin[2], origin[3]);
    }
    else if (OPTION("-torigin")) {
      PARSE_ARGUMENT(origin[3]);
      image->PutTOrigin(origin[3]);
    }
    else if (OPTION("-orientation")) {
      PARSE_ARGUMENT(xaxis[0]);
      PARSE_ARGUMENT(xaxis[1]);
      PARSE_ARGUMENT(xaxis[2]);
      Normalize(xaxis);
      PARSE_ARGUMENT(yaxis[0]);
      PARSE_ARGUMENT(yaxis[1]);
      PARSE_ARGUMENT(yaxis[2]);
      Normalize(yaxis);
      PARSE_ARGUMENT(zaxis[0]);
      PARSE_ARGUMENT(zaxis[1]);
      PARSE_ARGUMENT(zaxis[2]);
      Normalize(zaxis);
      image->PutOrientation(xaxis, yaxis, zaxis);
    }
    else if (OPTION("-copy-spacing") || OPTION("-copy-voxel-size") || OPTION("-copy-size")) {
      GreyImage target(ARGUMENT);
      target.GetPixelSize(&xsize, &ysize, &zsize, &tsize);
      image->PutPixelSize(xsize, ysize, zsize, tsize);
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
    else if (OPTION("-copy-origin-orientation-spacing") ||
             OPTION("-copy-origin-orientation-voxel-size") ||
             OPTION("-copy-origin-orientation-size")) {
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
    else if (OPTION("-dofout")) {
      dofout_name = ARGUMENT;
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

  #ifdef HAVE_MIRTK_Transformation
    if (dofout_name) {
      ImageAttributes output_attr = image->Attributes();
      output_attr._smat.Ident();
      const Matrix output_mat = output_attr.GetImageToWorldMatrix();
      AffineTransformation dofout;
      dofout.PutMatrix(input_mat * output_mat.Inverse());
      dofout.Write(dofout_name);
    }
  #endif

  if (verbose) image->Print();
  image->Write(output_name);

  return 0;
}

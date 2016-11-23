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
#include "mirtk/ImageReader.h"
#include "mirtk/IOConfig.h"
#include "mirtk/Transformation.h"
#include "mirtk/HomogeneousTransformation.h"
#include "mirtk/RigidTransformation.h"
#include "mirtk/ImageTransformation.h"
#include "mirtk/InterpolateImageFunction.h"

using namespace mirtk;


// ===========================================================================
// Help
// ===========================================================================

void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <source> <output> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Applies a transformation to an input image. Each voxel center of" << endl;
  cout << "  the target image is mapped by the given transformation to the space of" << endl;
  cout << "  the source image. The output intensity for the target voxel is the" << endl;
  cout << "  source image intensity interpolated at the mapped point and cast to" << endl;
  cout << "  the output data type." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  source   Source image." << endl;
  cout << "  output   Transformed source image." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -dof <file>      Affine source header transformation. (default: none)" << endl;
  cout << "  -dof_i <file>    Inverse affine source header transformation. (default: none)" << endl;
  cout << "  -dofin <file>    Transformation or 'Id'/'Identity'. (default: Id)" << endl;
  cout << "  -interp <mode>   Interpolation mode, e.g., \"NN\". (default: Linear)" << endl;
  cout << "  -target <file>   Target image. (default: source)" << endl;
  cout << "  -type, -dtype, -datatype <type>" << endl;
  cout << "      Data type of output image. (default: data type of source)" << endl;
  cout << "  -Tp <value>      Target padding value. (default: none)" << endl;
  cout << "  -Sp <value>      Source padding value. (default: none)" << endl;
  cout << "  -Tt <value>      Time point of target image. (default: torigin)" << endl;
  cout << "  -St <value>      Time point of source image. (default: torigin)" << endl;
  cout << "  -invert          Invert transformation. (default: off)" << endl;
  cout << "  -2d              Project transformation to 2D, ignoring mapped z coordinate. (default: off)" << endl;
  cout << endl;
  cout << "Interpolation modes:" << endl;
  for (int i = 0; i < Interpolation_Last; ++i) {
    cout << "  " << ToString(static_cast<InterpolationMode>(i)) << endl;
  }
  PrintCommonOptions(cout);
  cout << endl;
}

// ===========================================================================
// Main
// ===========================================================================

int main(int argc, char **argv)
{
  // Parse positional arguments
  EXPECTS_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  // Read image
  InitializeIOLibrary();
  UniquePtr<ImageReader> input_reader(ImageReader::New(input_name));
  UniquePtr<BaseImage>   source(input_reader->Run());

  // Parse optional arguments
  const char       *dof_name      = nullptr;
  bool              dof_invert    = false;
  const char       *dofin_name    = nullptr;
  const char       *target_name   = nullptr;
  InterpolationMode interpolation = Interpolation_NN;
  ImageDataType     dtype         = MIRTK_VOXEL_UNKNOWN;

  double target_t = NaN;
  double source_t = NaN;

  double source_padding = 0;
  double target_padding = NaN;
  bool  invert          = false;
  bool  twod            = false;

  for (ALL_OPTIONS) {
    if      (OPTION("-dof")   ) dof_name       = ARGUMENT, dof_invert = false;
    else if (OPTION("-dof_i") ) dof_name       = ARGUMENT, dof_invert = true;
    else if (OPTION("-dofin") ) dofin_name     = ARGUMENT;
    else if (OPTION("-target")) target_name    = ARGUMENT;
    else if (OPTION("-Tp")    ) PARSE_ARGUMENT(target_padding);
    else if (OPTION("-Sp")    ) PARSE_ARGUMENT(source_padding);
    else if (OPTION("-Tt")    ) PARSE_ARGUMENT(target_t);
    else if (OPTION("-St")    ) PARSE_ARGUMENT(source_t);
    else if (OPTION("-interp")) PARSE_ARGUMENT(interpolation);
    else HANDLE_BOOL_OPTION(invert);
    else HANDLE_BOOLEAN_OPTION("2d", twod);
    // backwards compatibility options
    else if (OPTION("-nn")     ) interpolation  = Interpolation_NN;
    else if (OPTION("-linear") ) interpolation  = Interpolation_Linear;
    else if (OPTION("-bspline")) interpolation  = Interpolation_BSpline;
    else if (OPTION("-cspline")) interpolation  = Interpolation_CSpline;
    else if (OPTION("-sinc")   ) interpolation  = Interpolation_Sinc;
    else if (OPTION("-sbased") ) interpolation  = Interpolation_SBased;
    else if (OPTION("-type") || OPTION("-dtype") || OPTION("-datatype")) {
      PARSE_ARGUMENT(dtype);
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }
  if (dtype == MIRTK_VOXEL_UNKNOWN) {
    dtype = static_cast<ImageDataType>(source->GetDataType());
  }

  // Apply affine header transformation
  if (dof_name) {
    UniquePtr<Transformation> t(Transformation::New(dof_name));
    HomogeneousTransformation *lin = dynamic_cast<HomogeneousTransformation *>(t.get());
    if (lin) {
      Matrix mat = lin->GetMatrix();
      if (dof_invert) mat.Invert();
      source->PutAffineMatrix(mat, true);
    } else {
      FatalError("Header transformation (-dof|-dof_i) can only be affine");
    }
  }

  // Read target image
  UniquePtr<BaseImage> target;
  if (target_name) {
    UniquePtr<ImageReader> target_reader(ImageReader::New(target_name));
    target.reset(target_reader->Run());
  }

  // Instantiate interpolator
  UniquePtr<InterpolateImageFunction> interpolator(InterpolateImageFunction::New(interpolation));

  // Initialize output image
  // Note: Always use floating point for intermediate interpolated image values!
  if (target) {
    if (target->T() != source->T() || (target->GetDataType() != MIRTK_VOXEL_FLOAT && target->GetDataType() != MIRTK_VOXEL_DOUBLE)) {
      UniquePtr<BaseImage> tmp(new RealImage(target->Attributes(), source->T()));
      tmp->PutTSize(source->GetTSize());
      for (int l = 0; l < tmp->T(); ++l)
      for (int k = 0; k < tmp->Z(); ++k)
      for (int j = 0; j < tmp->Y(); ++j)
      for (int i = 0; i < tmp->X(); ++i) {
        tmp->PutAsDouble(i, j, k, l, target->GetAsDouble(i, j, k, 0));
      }
      target.reset(tmp.release());
    }
  } else {
    target.reset(new RealImage(source->Attributes()));
    if (!IsNaN(target_padding)) {
      const int nvox = source->NumberOfVoxels();
      for (int vox = 0; vox < nvox; ++vox) {
        target->PutAsDouble(vox, source->GetAsDouble(vox));
      }
    }
  }
  if (IsNaN(target_padding)) target_padding = -inf;

  // Set temporal offset
  if (!IsNaN(target_t)) target->PutTOrigin(target_t);
  if (!IsNaN(source_t)) source->PutTOrigin(source_t);

  // Instantiate image transformation
  UniquePtr<Transformation> transformation;
  if (dofin_name == NULL || strcmp(dofin_name, "identity") == 0
                         || strcmp(dofin_name, "Identity") == 0
                         || strcmp(dofin_name, "Id")       == 0) {
    // Create identity transformation
    transformation.reset(new RigidTransformation());
  } else {
    // Read transformation
    transformation.reset(Transformation::New(dofin_name));
  }

  // Create image transformation filter
  ImageTransformation imagetransformation;
  imagetransformation.Input(source.get());
  imagetransformation.Transformation(transformation.get());
  imagetransformation.Output(target.get());
  imagetransformation.TargetPaddingValue(target_padding);
  imagetransformation.SourcePaddingValue(source_padding);
  imagetransformation.Interpolator(interpolator.get());
  imagetransformation.Invert(invert);
  imagetransformation.TwoD(twod);

  // Transform source image
  imagetransformation.Run();
  if (imagetransformation.NumberOfSingularPoints() > 0) {
    ostringstream msg;
    msg << "Transformation is non-invertible at "
        << imagetransformation.NumberOfSingularPoints() << " point";
    if (imagetransformation.NumberOfSingularPoints() > 1) msg << 's';
    Warning(msg.str());
  }

  // Convert to output data type
  if (target->GetDataType() != dtype) {
    UniquePtr<BaseImage> tmp(BaseImage::New(dtype));
    *tmp = *target;
    target.reset(tmp.release());
  }

  // Write the transformed image
  target->Write(output_name);

  return 0;
}

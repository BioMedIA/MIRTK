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
#include "mirtk/ResamplingWithPadding.h"

using namespace mirtk;


// ===========================================================================
// Help
// ===========================================================================

void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <source> <output> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Applies a transformation to an input image. Each voxel center of\n";
  cout << "  the target image is mapped by the given transformation to the space of\n";
  cout << "  the source image. The output intensity for the target voxel is the\n";
  cout << "  source image intensity interpolated at the mapped point and cast to\n";
  cout << "  the output data type.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  source   Source image.\n";
  cout << "  output   Transformed source image.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -dofin <file>\n";
  cout << "      Transformation or 'Id'/'Identity'. (default: Id)\n";
  cout << "  -invert [on|off], -noinvert\n";
  cout << "      Enable/disable inversion of :option:`-dofin` transformation. (default: off)\n";
  cout << "  -interp, -interpolation <mode>\n";
  cout << "      Interpolation mode, e.g., \"NN\". (default: Linear)\n";
  cout << "  -target <file>\n";
  cout << "      Target image. (default: source)\n";
  cout << "  -target-affdof <file>\n";
  cout << "      Affine target header transformation. (default: none)\n";
  cout << "  -target-invdof <file>\n";
  cout << "      Inverse affine target header transformation. (default: none)\n";
  cout << "  -source-affdof <file>\n";
  cout << "      Affine source header transformation. (default: none)\n";
  cout << "  -source-invdof <file>\n";
  cout << "      Inverse affine source header transformation. (default: none)\n";
  cout << "  -apply-affdof\n";
  cout << "      Apply affine header transformation to output image header.\n";
  cout << "      When this option is not specified, the output image attributes\n";
  cout << "      are identical to the :option:`-target` image. When this option\n";
  cout << "      is given, the :option:`-target-affdof` or :option:`-target-invdof`,\n";
  cout << "      respectively, is applied to the output image attributes.\n";
  cout << "  -spacing, -voxel-size <dx> [<dy> [<dz> [<dt>]]]\n";
  cout << "      Voxel size of output image. (default: :option:`-target` spacing)\n";
  cout << "  -type, -dtype, -datatype <type>\n";
  cout << "      Data type of output image. (default: data type of source)\n";
  cout << "  -Tp, -target-padding <value>\n";
  cout << "      Target padding value. (default: none)\n";
  cout << "  -Sp, -source-padding <value>\n";
  cout << "      Source padding value. (default: none)\n";
  cout << "  -padding <value>\n";
  cout << "      Set both :option:`-target-padding` and :option:`-source-padding` to same value.\n";
  cout << "  -Tt, -target-time <value>\n";
  cout << "      Time point of target image. (default: torigin)\n";
  cout << "  -St, -source-time <value>\n";
  cout << "      Time point of source image. (default: torigin)\n";
  cout << "  -2d [on|off], -no2d\n";
  cout << "      Project transformed points to 2D, i.e., ignore mapped z coordinate. (default: off)\n";
  cout << "  -3d\n";
  cout << "      Alias for :option:`-2d off`.\n";
  cout << "\n";
  cout << "Interpolation modes:\n";
  for (int i = 0; i < Interpolation_Last; ++i) {
    InterpolationMode mode = static_cast<InterpolationMode>(i);
    if (mode == InterpolationWithoutPadding(mode)) {
      cout << "  " << ToString(static_cast<InterpolationMode>(i));
      if (mode != Interpolation_Default) cout << " [with padding]";
      cout << "\n";
    }
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
  const char       *srcdof_name   = nullptr;
  bool              srcdof_invert = false;
  const char       *dofin_name    = nullptr;
  const char       *target_name   = nullptr;
  const char       *tgtdof_name   = nullptr;
  bool              tgtdof_invert = false;
  bool              affdof_apply  = false;
  InterpolationMode interpolation = Interpolation_NN;
  ImageDataType     dtype         = MIRTK_VOXEL_UNKNOWN;
  double            spacing[3]    = {0., 0., 0.};

  double target_t = NaN;
  double source_t = NaN;

  double source_padding = 0;
  double target_padding = NaN;
  bool   invert         = false;
  bool   twod           = false;

  for (ALL_OPTIONS) {
    if (OPTION("-dofin") ) {
      dofin_name = ARGUMENT;
    }
    else if (OPTION("-target")) {
      target_name = ARGUMENT;
    }
    else if (OPTION("-target-affdof")) {
      tgtdof_name   = ARGUMENT;
      tgtdof_invert = false;
    }
    else if (OPTION("-target-invdof")) {
      tgtdof_name   = ARGUMENT;
      tgtdof_invert = true;
    }
    else if (OPTION("-source-affdof") || OPTION("-dof")) {
      srcdof_name   = ARGUMENT;
      srcdof_invert = false;
    }
    else if (OPTION("-source-invdof") || OPTION("-dof_i")) {
      srcdof_name   = ARGUMENT;
      srcdof_invert = true;
    }
    else if (OPTION("-apply-affdof")) {
      affdof_apply = true;
    }
    else if (OPTION("-spacing") || OPTION("-voxel-size")) {
      PARSE_ARGUMENT(spacing[0]);
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(spacing[1]);
        if (HAS_ARGUMENT) {
          PARSE_ARGUMENT(spacing[2]);
        }
      } else {
        spacing[1] = spacing[2] = spacing[0];
      }
    }
    else if (OPTION("-padding")) {
      PARSE_ARGUMENT(target_padding);
      source_padding = target_padding;
    }
    else if (OPTION("-Tp") || OPTION("-target-padding")) {
      PARSE_ARGUMENT(target_padding);
    }
    else if (OPTION("-Sp") || OPTION("-source-padding")) {
      PARSE_ARGUMENT(source_padding);
    }
    else if (OPTION("-Tt") || OPTION("-target-time")) {
      PARSE_ARGUMENT(target_t);
    }
    else if (OPTION("-St") || OPTION("-source-time")) {
      PARSE_ARGUMENT(source_t);
    }
    else if (OPTION("-interp") || OPTION("-interpolation")) {
      PARSE_ARGUMENT(interpolation);
    }
    else HANDLE_BOOL_OPTION(invert);
    else HANDLE_BOOLEAN_OPTION("2d", twod);
    else if (OPTION("-3d")) twod = false;
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

  // Instantiate interpolator
  UniquePtr<InterpolateImageFunction> interpolator(InterpolateImageFunction::New(interpolation));

  // Initialize output image
  // Note: Always use floating point for intermediate interpolated image values!
  UniquePtr<RealImage> target;
  if (target_name) {
    UniquePtr<ImageReader> reader(ImageReader::New(target_name));
    UniquePtr<BaseImage> image(reader->Run());
    if (image->T() == source->T()) {
      target.reset(dynamic_cast<RealImage *>(image.get()));
    }
    if (target) {
      image.release();
    } else {
      target.reset(new RealImage(image->Attributes(), source->T()));
      target->PutTSize(source->GetTSize());
      for (int l = 0; l < target->T(); ++l)
      for (int k = 0; k < target->Z(); ++k)
      for (int j = 0; j < target->Y(); ++j)
      for (int i = 0; i < target->X(); ++i) {
        target->PutAsDouble(i, j, k, l, image->GetAsDouble(i, j, k, 0));
      }
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

  // Resample to desired output spacing
  if (spacing[0] > 0. || spacing[1] > 0. || spacing[2] > 0.) {
    double dx, dy, dz;
    target->GetPixelSize(dx, dy, dz);
    if (!fequal(dx, spacing[0]) || !fequal(dy, spacing[1]) || !fequal(dz, spacing[2])) {
      if (spacing[0] > 0.) dx = spacing[0];
      if (spacing[1] > 0.) dy = spacing[1];
      if (spacing[2] > 0.) dz = spacing[2];
      if (IsInf(target_padding)) {
        Resampling<RealPixel> resampler(dx, dy, dz);
        resampler.Input(target.get());
        resampler.Output(target.get());
        resampler.Interpolator(interpolator.get());
        resampler.Run();
      } else {
        ResamplingWithPadding<RealPixel> resampler(dx, dy, dz, target_padding);
        resampler.Input(target.get());
        resampler.Output(target.get());
        resampler.Run();
      }
    }
  }

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

  // Set affine header transformation
  if (srcdof_name) {
    UniquePtr<Transformation> t(Transformation::New(srcdof_name));
    HomogeneousTransformation *lin = dynamic_cast<HomogeneousTransformation *>(t.get());
    if (lin) {
      Matrix mat = lin->GetMatrix();
      if (srcdof_invert) mat.Invert();
      source->PutAffineMatrix(mat, false);
    } else {
      FatalError("Source header transformation must be affine");
    }
  }
  if (tgtdof_name) {
    UniquePtr<Transformation> t(Transformation::New(tgtdof_name));
    HomogeneousTransformation *lin = dynamic_cast<HomogeneousTransformation *>(t.get());
    if (lin) {
      Matrix mat = lin->GetMatrix();
      if (tgtdof_invert) mat.Invert();
      target->PutAffineMatrix(mat, affdof_apply);
    } else {
      FatalError("Target header transformation must be affine");
    }
  } else if (!target_name && srcdof_name) {
    target->PutAffineMatrix(source->GetAffineMatrix(), affdof_apply);
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

  // Reset affine header transformation
  if (!affdof_apply && (tgtdof_name || (!target_name && srcdof_name))) {
    target->ResetAffineMatrix();
  }

  // Write the transformed image
  if (target->GetDataType() != dtype) {
    UniquePtr<BaseImage> output(BaseImage::New(dtype));
    *output = *target;
    output->Write(output_name);
  } else {
    target->Write(output_name);
  }

  return 0;
}

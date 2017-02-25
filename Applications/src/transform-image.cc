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
#include "mirtk/HashImage.h"
#include "mirtk/ImageReader.h"
#include "mirtk/IOConfig.h"
#include "mirtk/Transformation.h"
#include "mirtk/HomogeneousTransformation.h"
#include "mirtk/RigidTransformation.h"
#include "mirtk/ImageTransformation.h"
#include "mirtk/InterpolateImageFunction.h"
#include "mirtk/ResamplingWithPadding.h"
#include "mirtk/NearestNeighborInterpolateImageFunction.h"

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
  cout << "  -interpolation, -interp <mode>\n";
  cout << "      Interpolation mode, e.g., \"NN\". (default: Linear)\n";
  cout << "  -labels [<n>...|all]\n";
  cout << "      Transform the specified segmentation labels. When no arguments or \"all\" is\n";
  cout << "      given, all input labels are transformed. When :option:`-interpolation` mode\n";
  cout << "      is nearest neighbor (\"NN\"), the input label image is directly resampled\n";
  cout << "      by assigning the output the nearest label value. Otherwise, the default,\n";
  cout << "      the binary mask corresponding to each label is resampled using the selected\n";
  cout << "      interpolation mode (e.g., linear) and the resulting fuzzy segmentations\n";
  cout << "      converted back into a hard segmentation by assigning the label with the\n";
  cout << "      highest interpolated value (pseudo-probability). This results in smoother\n";
  cout << "      transformed hard segmentations with reduced NN interpolation artifacts.\n";
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
// Auxiliaries
// ===========================================================================

typedef GenericImage<double> CoordMap;
typedef HashImage<short> FuzzySeg;

// ---------------------------------------------------------------------------
struct ComposeImageToWorldWithDisplacementMap
{
  const Matrix   *_ImageToWorld;
  const CoordMap *_Displacement;
  const Matrix   *_WorldToImage;
  CoordMap       *_CompositeMap;

  void operator ()(const blocked_range3d<int> &re) const
  {
    double u, v, w, x, y, z;
    for (int k = re.pages().begin(); k < re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j < re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i < re.cols ().end(); ++i) {
      u = i, v = j, w = k;
      Transform(*_ImageToWorld, u, v, w, x, y, z);
      x += _Displacement->Get(i, j, k, 0);
      y += _Displacement->Get(i, j, k, 1);
      z += _Displacement->Get(i, j, k, 2);
      Transform(*_WorldToImage, x, y, z, u, v, w);
      _CompositeMap->Put(i, j, k, 0, u);
      _CompositeMap->Put(i, j, k, 1, v);
      _CompositeMap->Put(i, j, k, 2, w);
    }
  }
};

// ---------------------------------------------------------------------------
/// Evaluate map of target voxel indices to continuous source voxel indices
int EvaluateTargetToSourceMap(CoordMap &map, const BaseImage *source,
                              const Transformation *dof, bool invert = false)
{
  map = 0.;
  const double t  = source->GetTOrigin();
  const double t0 = map.GetTOrigin();
  int nsingular = 0;
  if (invert) {
    nsingular = dof->InverseDisplacement(map, t, t0);
  } else {
    dof->Displacement(map, t, t0);
  }
  MIRTK_START_TIMING();
  const Matrix i2w = map.GetImageToWorldMatrix();
  const Matrix w2i = source->GetWorldToImageMatrix();
  ComposeImageToWorldWithDisplacementMap compose;
  compose._ImageToWorld = &i2w;
  compose._Displacement = &map;
  compose._WorldToImage = &w2i;
  compose._CompositeMap = &map;
  parallel_for(blocked_range3d<int>(0, map.Z(), 0, map.Y(), 0, map.X()), compose);
  MIRTK_DEBUG_TIMING(1, "composing maps");
  return nsingular;
}

// ---------------------------------------------------------------------------
BinaryImage Binarize(const GreyImage &seg, GreyPixel label)
{
  BinaryImage mask(seg);
  const int nvox = seg.NumberOfVoxels();
  for (int idx = 0; idx < nvox; ++idx) {
    mask(idx) = (seg(idx) == label ? 1 : 0);
  }
  return mask;
}

// ---------------------------------------------------------------------------
struct ResampleLabelBody
{
  const InterpolateImageFunction *_LabelMask;
  const CoordMap                 *_CoordMap;
  FuzzySeg                        _Output;

  ResampleLabelBody(const CoordMap *map, const InterpolateImageFunction *interp)
  :
    _CoordMap(map),
    _LabelMask(interp),
    _Output(map->Attributes(), 1)
  {
    _Output.DefaultValue(0);
  }

  ResampleLabelBody(ResampleLabelBody &other, split)
  :
    _CoordMap(other._CoordMap),
    _LabelMask(other._LabelMask),
    _Output(other._Output.Attributes())
  {
    _Output.DefaultValue(0);
  }

  void join(ResampleLabelBody &other)
  {
    _Output.Data().insert(other._Output.Begin(), other._Output.End());
  }

  void operator ()(const blocked_range<int> &slices)
  {
    double u, v, w, prob;
    for (int k = slices.begin(); k < slices.end(); ++k)
    for (int j = 0; j < _Output.Y(); ++j)
    for (int i = 0; i < _Output.X(); ++i) {
      u = _CoordMap->Get(i, j, k, 0);
      v = _CoordMap->Get(i, j, k, 1);
      w = _CoordMap->Get(i, j, k, 2);
      prob = _LabelMask->Evaluate(u, v, w);
      if (prob > 0.) {
        prob = clamp(prob, 0., 1.) * voxel_limits<FuzzySeg::VoxelType>::max();
        _Output.Put(i, j, k, static_cast<FuzzySeg::VoxelType>(prob));
      }
    }
  }
};

// ---------------------------------------------------------------------------
FuzzySeg ResampleLabel(const GreyImage &seg, GreyPixel label,
                       const CoordMap &map, InterpolationMode mode)
{
  MIRTK_START_TIMING();
  const BinaryImage mask = Binarize(seg, label);
  UniquePtr<InterpolateImageFunction> interp(InterpolateImageFunction::New(mode));
  interp->DefaultValue(0.);
  interp->Input(&mask);
  interp->Initialize();
  ResampleLabelBody resample(&map, interp.get());
#ifdef HAVE_TBB
  int n = 0;
  const int nvox = seg.NumberOfSpatialVoxels();
  for (int idx = 0; idx < nvox; ++idx) {
    if (seg(idx) == label) ++n;
  }
  if (n > iceil(.75 * nvox)) {
    resample(blocked_range<int>(0, map.Z()));
  } else {
    parallel_reduce(blocked_range<int>(0, map.Z(), map.Z() / 8), resample);
  }
#else
  resample(blocked_range<int>(0, map.Z()));
#endif
  MIRTK_DEBUG_TIMING(1, "label resampling");
  return resample._Output;
}

// ---------------------------------------------------------------------------
template <class T>
Array<T> ToArray(const OrderedSet<T> &values)
{
  size_t i = 0;
  Array<T> array(values.size());
  for (auto value : values) {
    array[i++] = value;
  }
  return array;
}

// ---------------------------------------------------------------------------
struct MakeHardSegmentationBody
{
  const GreyImage        *_Source;
  const CoordMap         *_CoordMap;
  const Array<GreyPixel> *_Labels;
  const Array<FuzzySeg>  *_Probs;
  FuzzySeg::VoxelType     _MinProb;
  GreyImage              *_Output;

  void operator ()(const blocked_range3d<int> &re) const
  {
    int u, v, w;
    Array<GreyPixel> candidates;
    candidates.reserve(_Labels->size());
    FuzzySeg::VoxelType prob, max_prob;
    for (int k = re.pages().begin(); k < re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j < re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i < re.cols ().end(); ++i) {
      candidates.clear();
      max_prob = 0;
      for (size_t l = 0; l < _Labels->size(); ++l) {
        prob = (*_Probs)[l].Get(i, j, k);
        if (prob > _MinProb && prob >= max_prob) {
          max_prob = prob;
          candidates.push_back((*_Labels)[l]);
        }
      }
      if (candidates.size() == 1) {
        _Output->Put(i, j, k, candidates[0]);
      } else {
        u = iround(_CoordMap->Get(i, j, k, 0));
        v = iround(_CoordMap->Get(i, j, k, 1));
        w = iround(_CoordMap->Get(i, j, k, 2));
        if (_Source->IsInside(u, v, w)) {
          _Output->Put(i, j, k, _Source->Get(u, v, w));
        }
      }
    }
  }
};

// ---------------------------------------------------------------------------
GreyImage *MakeHardSegmentation(const GreyImage        &seg,
                                const Array<FuzzySeg>  &probs,
                                const Array<GreyPixel> &labels,
                                const CoordMap         &map,
                                FuzzySeg::VoxelType eps = 0)
{
  MIRTK_START_TIMING();
  UniquePtr<GreyImage> output(new GreyImage(map.Attributes(), 1));
  MakeHardSegmentationBody body;
  body._Source   = &seg;
  body._CoordMap = &map;
  body._Labels   = &labels;
  body._Probs    = &probs;
  body._MinProb  = eps;
  body._Output   = output.get();
  parallel_for(blocked_range3d<int>(0, map.Z(), 0, map.Y(), 0, map.X()), body);
  MIRTK_DEBUG_TIMING(1, "making hard segmentation");
  return output.release();
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

  // Parse optional arguments
  const char       *srcdof_name   = nullptr;
  bool              srcdof_invert = false;
  const char       *dofin_name    = nullptr;
  const char       *target_name   = nullptr;
  const char       *tgtdof_name   = nullptr;
  bool              tgtdof_invert = false;
  bool              affdof_apply  = false;
  InterpolationMode interpolation = Interpolation_Linear;
  ImageDataType     dtype         = MIRTK_VOXEL_UNKNOWN;
  double            spacing[3]    = {0., 0., 0.};
  bool              all_labels    = false;
  OrderedSet<GreyPixel> labels;

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
    else if (OPTION("-labels") || OPTION("-label")) {
      GreyPixel label;
      do {
        string arg = ToLower(ARGUMENT);
        if (arg == "all") {
          all_labels = true;
        } else if (FromString(arg, label)) {
          labels.insert(label);
        } else {
          FatalError("Invalid -labels option argument: " << arg);
        }
      } while (HAS_ARGUMENT);
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

  // Initialize I/O library
  InitializeIOLibrary();
  UniquePtr<BaseImage> output;

  // Instantiate image transformation
  int nsingular = 0;
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

  if (labels.empty() && !all_labels) {

    // Read source image
    UniquePtr<BaseImage> source(BaseImage::New(input_name));

    // Instantiate image interpolator
    UniquePtr<InterpolateImageFunction> interpolator;
    interpolator.reset(InterpolateImageFunction::New(interpolation));

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

    // Transform source intensity image
    ImageTransformation imagetransformation;
    imagetransformation.Input(source.get());
    imagetransformation.Transformation(transformation.get());
    imagetransformation.Output(target.get());
    imagetransformation.TargetPaddingValue(target_padding);
    imagetransformation.SourcePaddingValue(source_padding);
    imagetransformation.Interpolator(interpolator.get());
    imagetransformation.Invert(invert);
    imagetransformation.TwoD(twod);
    imagetransformation.Run();
    nsingular = imagetransformation.NumberOfSingularPoints();

    output.reset(target.release());

  } else {

    // Read input segmentation
    GreyImage source(input_name);
    const int num_src_vox = source.NumberOfSpatialVoxels();
    const GreyPixel bg_label = static_cast<GreyPixel>(source_padding);
    if (all_labels || labels.empty()) {
      GreyPixel label;
      labels.clear();
      for (int idx = 0; idx < num_src_vox; ++idx) {
        label = source(idx);
        labels.insert(source(idx));
      }
    } else {
      for (int idx = 0; idx < num_src_vox; ++idx) {
        if (labels.find(source(idx)) == labels.end()) {
          source(idx) = bg_label;
        }
      }
      labels.insert(bg_label);
    }

    // Output segmentation attributes
    ImageAttributes attr = source.Attributes();
    if (target_name) {
      BinaryImage temp(target_name);
      attr = temp.Attributes();
    }

    // Set desired output voxel size
    if (spacing[0] > .0) {
      attr._x  = iceil(attr._x * attr._dx / spacing[0]);
      attr._dx = spacing[0];
    }
    if (spacing[1] > .0) {
      attr._y  = iceil(attr._y * attr._dy / spacing[1]);
      attr._dy = spacing[1];
    }
    if (spacing[2] > .0) {
      attr._z  = iceil(attr._z * attr._dz / spacing[2]);
      attr._dz = spacing[2];
    }

    // Set temporal offset
    if (!IsNaN(target_t)) attr._torigin = target_t;
    if (!IsNaN(source_t)) source.PutTOrigin(source_t);

    // Set affine header transformation
    if (srcdof_name) {
      UniquePtr<Transformation> t(Transformation::New(srcdof_name));
      HomogeneousTransformation *lin = dynamic_cast<HomogeneousTransformation *>(t.get());
      if (lin) {
        Matrix mat = lin->GetMatrix();
        if (srcdof_invert) mat.Invert();
        source.PutAffineMatrix(mat, false);
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
        attr.PutAffineMatrix(mat, affdof_apply);
      } else {
        FatalError("Target header transformation must be affine");
      }
    } else if (!target_name && srcdof_name) {
      attr.PutAffineMatrix(source.GetAffineMatrix(), affdof_apply);
    }

    // Transform source segmentation using NN interpolation
    if (interpolation == Interpolation_NN) {

      if (verbose) {
        cout << "Resampling label image using NN interpolation...";
        cout.flush();
      }
      output.reset(new GreyImage(attr));
      GenericNearestNeighborInterpolateImageFunction<GreyImage> nn;
      ImageTransformation imagetransformation;
      imagetransformation.Input(&source);
      imagetransformation.Transformation(transformation.get());
      imagetransformation.Output(output.get());
      imagetransformation.SourcePaddingValue(source_padding);
      imagetransformation.Interpolator(&nn);
      imagetransformation.Invert(invert);
      imagetransformation.TwoD(twod);
      imagetransformation.Run();
      nsingular = imagetransformation.NumberOfSingularPoints();
      if (verbose) cout << " done" << endl;

    // Transform source labels individually and obtain transformed
    // hard segmentation from individual fuzzy segmentation images
    } else {

      if (verbose) {
        cout << "Evaluating target to source index map...";
        cout.flush();
      }
      CoordMap map(attr, 3);
      nsingular = EvaluateTargetToSourceMap(map, &source, transformation.get(), invert);
      if (verbose) cout << " done" << endl;

      Array<FuzzySeg> probs;
      probs.reserve(labels.size());
      for (auto label : labels) {
        if (verbose) {
          cout << "Resampling label " << label << "...";
          cout.flush();
        }
        probs.push_back(ResampleLabel(source, label, map, interpolation));
        if (verbose) cout << " done" << endl;
      }

      if (verbose) {
        cout << "Making hard segmentation...";
        cout.flush();
      }
      FuzzySeg::VoxelType eps = .5 * voxel_limits<FuzzySeg::VoxelType>::max();
      output.reset(MakeHardSegmentation(source, probs, ToArray(labels), map, eps));
      if (verbose) cout << " done" << endl;

    }
  }

  // Report number of singular points
  if (nsingular > 0) {
    ostringstream msg;
    msg << "Transformation is non-invertible at " << nsingular << " point";
    if (nsingular > 1) msg << 's';
    Warning(msg.str());
  }

  // Reset affine header transformation
  if (!affdof_apply && (tgtdof_name || (!target_name && srcdof_name))) {
    output->ResetAffineMatrix();
  }

  // Write the transformed image
  if (dtype != MIRTK_VOXEL_UNKNOWN && output->GetDataType() != dtype) {
    UniquePtr<BaseImage> image(BaseImage::New(dtype));
    *image = *output;
    image->Write(output_name);
  } else {
    output->Write(output_name);
  }

  return 0;
}

/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
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

#include "mirtk/IOConfig.h"

#include "mirtk/Matrix.h"
#include "mirtk/GenericImage.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/InterpolateImageFunction.h"
#include "mirtk/DataStatistics.h"

#ifdef HAVE_MIRTK_Transformation
#  include "mirtk/Transformation.h"
#  include "mirtk/HomogeneousTransformation.h"
#  include "mirtk/FluidFreeFormTransformation.h"
#endif

using namespace mirtk;
using namespace mirtk::data::statistic;


// =============================================================================
// Help
// =============================================================================

/// Print program usage information
void PrintHelp(const char* name)
{
  cout << "\n";
#ifdef HAVE_MIRTK_Transformation
  cout << "Usage: " << name << " <output> -images <images.lst> [options]\n";
  cout << "       " << name << " <output> -image <image1> [-dof <dof1>...] -image <image2> [-dof <dof2>...]... [options]\n";
  cout << "       " << name << " <output> <image1> [<dof1>...] <image2> [<dof2>...]... [options]\n";
#else
  cout << "Usage: " << name << " <output> -images <images.lst> [options]\n";
  cout << "       " << name << " <output> -image <image1> -image <image2> ... [options]\n";
  cout << "       " << name << " <output> <image1> <image2>... [options]\n";
#endif
  cout << "       " << name << " <output> <sequence> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Computes voxel-wise average image of the given (transformed) input images.\n";
#ifdef HAVE_MIRTK_Transformation
  cout << "\n";
  cout << "  The input images are optionally transformed by one or more image-to-average space\n";
  cout << "  transformation (:option:`-dof`) or its inverse (:option:`-dof_i`), respectively.\n";
  cout << "  When more than one transformation is given, these are concatenated via funtion\n";
  cout << "  composition.\n";
#endif
  cout << "\n";
  cout << "Required arguments:\n";
  cout << "  <output>\n";
  cout << "      Voxel-wise average image.\n";
  cout << "\n";
  cout << "Input options:\n";
  cout << "  -image <file>\n";
  cout << "      A single input image.\n";
  cout << "  -images <file>\n";
  cout << "      Text file with N lines containing the file name of each input image,\n";
  cout << "      optionally a transformation file name (see :option:`-dof` and :option:`-dof_i`),\n";
  cout << "      optionally followed by a weight for a weighted output average (default weight is 1).\n";
  cout << "      The first line of the text file must specify the common base directory\n";
  cout << "      of all relative image and transformation file paths occurring on the\n";
  cout << "      subsequent N lines. A path starting with './' must be relative to the\n";
  cout << "      directory containing the input text file itself.\n";
  cout << "  -delim, -delimiter <c>\n";
  cout << "      Delimiter used in :option:`-images` file.\n";
  cout << "      (default: ',' for .csv, '\\t' for .tsv, and ' ' otherwise)\n";
#ifdef HAVE_MIRTK_Transformation
  cout << "  -invert\n";
  cout << "      Invert transformations specified in :option:`-images` file.\n";
  cout << "  -dof <file>\n";
  cout << "      Specifies a transformation to be applied to the preceeding image. (default: none)\n";
  cout << "  -dof_i <file>\n";
  cout << "      Specifies a transformation whose inverse is to be applied to the\n";
  cout << "      preceeding image similar to :option:`-dof`. (default: none)\n";
#endif
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -reference, -target <file>|max-size|max-space|average\n";
  cout << "      When <file> name specified, use the attributes of this reference image for average image.\n";
  cout << "      When \"max-size\" is specified, the input image with the most number of voxels is used as\n";
  cout << "      reference, whereas \"max-space\" selects the image which occupies the largest amount world space.\n";
  cout << "      By default, an \"average\" voxel size and image orientation is computed from the input images\n";
  cout << "      and the image size adjusted such that the image bounding boxes of all (affinely transformed)\n";
  cout << "      input images are contained within the average image space. (default: average)\n";
  cout << "  -voxel-size, -spacing, -size <dx> [<dy> [<dz>]]\n";
  cout << "      Voxel size of intensity average image. (default: :option:`-reference`)\n";
  cout << "  -padding <value>\n";
  cout << "      Input padding and output background value. No input padding if not specified.\n";
  cout << "      Output background value zero by default or minimum average intensity minus 1. (default: 0)\n";
  cout << "  -normalize, -normalization <mode>\n";
  cout << "      Input intensity normalization:\n";
  cout << "      - ``none``:    Use input intensities unmodified. (default)\n";
  cout << "      - ``mean``:    Divide by mean foreground value.\n";
  cout << "      - ``median``:  Divide by median foreground value.\n";
  cout << "      - ``z-score``: Subtract mean and divide by standard deviation.\n";
  cout << "      - ``unit``:    Rescale input intensities to [0, 1].\n";
  cout << "      - ``dist``:    Rescale input intensities to average mean and standard deviation of input images.\n";
  cout << "  -rescale, -rescaling <mode>|<min> <max>\n";
  cout << "      Linear rescaling of average intensity values:\n";
  cout << "      - ``none``: No rescaling of averaged intensities. (default)\n";
  cout << "      - ``unit``: Rescale average intensities to [0, 1].\n";
  cout << "      - ``dist``: Rescale average intensities to average mean and standard deviation of input images.\n";
  cout << "      - ``<min> <max>``: Rescale average intensities to specified output range.\n";
  cout << "  -margin <n>\n";
  cout << "      Crop/pad average image and ensure a margin of <n> voxels at each boundary. (default: -1/off)\n";
  cout << "  -interpolation, -interp <mode>\n";
  cout << "      Interpolation mode, e.g., NN, Linear, BSpline, Cubic, Sinc. (default: Linear)\n";
  cout << "  -label <value>\n";
  cout << "      Segmentation label of which to create an average probability map.\n";
  cout << "  -datatype, -dtype, -type char|uchar|short|float|double\n";
  cout << "      Data type of output image. The intermediate average image always has floating point data type.\n";
  cout << "      When this option is given, this average is cast to the respective output type before writing\n";
  cout << "      the result to the output image file. (default: float)\n";
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Types
// =============================================================================

// Type of input images
typedef RealPixel                                   InputType;
typedef GenericImage<InputType>                     InputImage;
typedef GenericInterpolateImageFunction<InputImage> InputImageFunction;

// Type of average image
typedef GenericImage<float>  AverageImage;

// Type of intensity image normalization
enum NormalizationMode
{
  Normalization_None,         ///< Use input intensity values unmodified
  Normalization_Mean,         ///< Divide input image values by mean intensity
  Normalization_Median,       ///< Divide input image values by median intensity
  Normalization_ZScore,       ///< Subtract mean intensity and divide by standard deviation
  Normalization_MeanStDev,    ///< Rescale input to average mean and standard deviation
  Normalization_UnitRange     ///< Rescale input intensities to [0, 1]
};

// Type of average image rescaling functions
enum RescalingMode
{
  Rescaling_None,        ///< No rescaling of average (normalized) input intensities
  Rescaling_MeanStDev,   ///< Rescale output to input mean and standard deviation
  Rescaling_Range        ///< Rescale output image to specified output range
};

// =============================================================================
// Auxiliary functions
// =============================================================================

/// Read image file names and times from input list file
int read_image_list_file(const char *image_list_name,
                         Array<string> &names,
                         Array<Array<string> > &dofs,
                         Array<double> &weights,
                         const char *delim = " ")
{
  Array<string> dof_names;
  string   base_dir = Directory(image_list_name); // Base directory for relative paths
  ifstream iff(image_list_name);                  // List input file stream
  string   line;                                  // Input line
  double   weight;                                // Image weight
  int      l = 0;                                 // Line number

  names  .clear();
  dofs   .clear();
  weights.clear();

  // Read base directory for image files
  if (!getline(iff, line)) {
    cerr << "Error: Cannot parse list file " << image_list_name << endl;
    return 0;
  }
  const bool discard_empty = true;
  const bool handle_quotes = true;
  auto columns = Split(line, delim, 0, discard_empty, handle_quotes);
  if (columns[0].empty()) columns[0] = ".";
  if (!base_dir.empty() && columns[0].front() != PATHSEP) {
    if (columns[0] != ".") base_dir += PATHSEP + columns[0];
  } else {
    base_dir = columns[0];
  }
  l++;

  while (getline(iff, line)) {
    l++;
    // Ignore comment lines starting with # character
    if (line[0] == '#') continue;
    // Split line at delimiter
    columns = Split(line, delim, 0, discard_empty, handle_quotes);
    if (columns.empty()) continue;
    // Insert columns into output containers
    if (!base_dir.empty() && columns[0].front() != PATHSEP) {
      names.push_back(base_dir + PATHSEP + columns[0]);
    } else {
      names.push_back(columns[0]);
    }
    auto ndofs = columns.size() - 1;
    if (FromString(columns.back(), weight)) {
      weights.push_back(weight);
      ndofs -= 1;
    } else {
      weights.push_back(1.0);
    }
    dof_names.clear();
    if (ndofs > 0) {
      dof_names.reserve(ndofs);
      for (size_t i = 1; i <= ndofs; ++i) {
        if (!columns[i].empty()) {
          if (!base_dir.empty() && columns[i].front() != PATHSEP) {
            dof_names.push_back(base_dir + PATHSEP + columns[i]);
          } else {
            dof_names.push_back(columns[i]);
          }
        }
      }
    }
    dofs.push_back(dof_names);
  }

  return static_cast<int>(names.size());
}

// -----------------------------------------------------------------------------
#ifdef HAVE_MIRTK_Transformation
void Read(const string &image_name, const Array<string> &imdof_names, const Array<bool> &imdof_invert,
          InputImage &image, Array<UniquePtr<Transformation> > &dofs, Array<bool> &invert)
{
  dofs.clear();
  invert.clear();
  image.Read(image_name.c_str());
  if (!imdof_names.empty()) {
    size_t i;
    for (i = 0; i < imdof_names.size(); ++i) {
      if (imdof_names[i].empty()) continue;
      dofs.push_back(UniquePtr<Transformation>(Transformation::New(imdof_names[i].c_str())));
      invert.push_back(imdof_invert[i]);
    }
    if (!dofs.empty()) {
      HomogeneousTransformation *lin;
      MultiLevelTransformation *mffd;
      FluidFreeFormTransformation *fluid;
      for (i = 0; i < dofs.size(); ++i) {
        lin = dynamic_cast<HomogeneousTransformation *>(dofs[i].get());
        if (!lin) {
          mffd = dynamic_cast<MultiLevelTransformation *>(dofs[i].get());
          if (!mffd || mffd->NumberOfLevels() > 0) break;
          lin = mffd->GetGlobalTransformation();
        }
        Matrix mat = lin->GetMatrix();
        if (invert[i]) mat.Invert();
        image.PutAffineMatrix(mat, true);
      }
      if (i > 0) {
        dofs  .erase(dofs  .begin(), dofs  .begin() + i);
        invert.erase(invert.begin(), invert.begin() + i);
      }
      if (!dofs.empty()) {
        fluid = dynamic_cast<FluidFreeFormTransformation *>(dofs.front().get());
        if (fluid && !invert.front()) {
          Matrix mat = fluid->GetGlobalTransformation()->GetMatrix();
          image.PutAffineMatrix(mat, true);
          fluid->GetGlobalTransformation()->Reset();
        }
      }
    }
  }
}
#endif // HAVE_MIRTK_Transformation

// -----------------------------------------------------------------------------
UniquePtr<bool[]> ForegroundMaskArray(const BaseImage *image)
{
  const int num = image->NumberOfSpatialVoxels();
  UniquePtr<bool[]> mask(new bool[num]);
  for (int idx = 0; idx < num; ++idx) {
    mask[idx] = (!IsNaN(image->GetAsDouble(idx)) && image->IsForeground(idx));
  }
  return mask;
}

// -----------------------------------------------------------------------------
struct AddVoxelValueToAverage : public VoxelFunction
{
  AverageImage       *_Average;
  InputImageFunction *_Image;
  double              _Weight;
  int                 _Label;
  double              _Scale;
  double              _Intercept;
  double              _Min;
  double              _Max;

  #ifdef HAVE_MIRTK_Transformation
    const Array<UniquePtr<Transformation> > *_Transformation;
    const Array<bool>                       *_Invert;
    const GenericImage<double>              *_Displacement;
  #endif

  template <class T>
  void operator()(int i, int j, int k, int, T *avg)
  {
    double x = i, y = j, z = k;
    _Average->ImageToWorld(x, y, z);
    #ifdef HAVE_MIRTK_Transformation
      if (_Displacement) {
        x += _Displacement->Get(i, j, k, 0);
        y += _Displacement->Get(i, j, k, 1);
        z += _Displacement->Get(i, j, k, 2);
      } else if (_Transformation) {
        // Note: Input transformation is from image to average!
        const int ndofs = static_cast<int>(_Transformation->size());
        for (int n = ndofs - 1; n >= 0; --n) {
          if ((*_Invert)[n]) (*_Transformation)[n]->Transform(x, y, z);
          else               (*_Transformation)[n]->Inverse  (x, y, z);
        }
      }
    #endif // HAVE_MIRTK_Transformation
    _Image->Input()->WorldToImage(x, y, z);
    if (_Label > 0) {
      if (static_cast<int>(_Image->Evaluate(x, y, z)) == _Label) {
        *avg += static_cast<T>(_Weight);
      }
    } else {
      double value = _Image->Evaluate(x, y, z);
      if (!IsNaN(value) && value != _Image->DefaultValue()) {
        value = _Scale * clamp(value, _Min, _Max) + _Intercept;
        if (IsNaN(*avg)) {
          *avg = static_cast<T>(_Weight * value);
        } else {
          *avg += static_cast<T>(_Weight * value);
        }
      }
    }
  }
};

// -----------------------------------------------------------------------------
bool Add(AverageImage &average, InputImage &image,
         #ifdef HAVE_MIRTK_Transformation
           const Array<UniquePtr<Transformation> > &dof,
           const Array<bool>                       &inv,
         #endif // HAVE_MIRTK_Transformation
         double weight, int label,
         InterpolationMode interpolation = Interpolation_Linear,
         NormalizationMode normalization = Normalization_None,
         double avg_mean = NaN, double avg_sigma = NaN)
{
  AddVoxelValueToAverage add;
  // Determine min/max intensity used to clamp interpolated values
  const double bg = image.GetBackgroundValueAsDouble();
  const InputType * const data = image.Data();
  const int num = image.NumberOfSpatialVoxels();
  UniquePtr<bool[]> mask = ForegroundMaskArray(&image);
  Extrema::Calculate(add._Min, add._Max, num, data, mask.get());
  if (IsNaN(add._Min) || IsNaN(add._Max)) return false;
  // Compute normalization parameters
  add._Scale = 1.;
  add._Intercept = 0.;
  if (normalization == Normalization_Mean) {
    double mean = Mean::Calculate(num, data, mask.get());
    if (!fequal(mean, 0.)) {
      add._Scale = 1. / mean;
    }
  } else if (normalization == Normalization_Median) {
    double median = Median::Calculate(num, data, mask.get());
    if (!fequal(median, 0.)) {
      add._Scale = 1. / median;
    }
  } else if (normalization == Normalization_ZScore) {
    double mean, sigma;
    NormalDistribution::Calculate(mean, sigma, num, data, mask.get());
    if (fequal(sigma, 0.)) return false;
    add._Scale = 1. / sigma;
    add._Intercept = - mean / sigma;
    double zscore_min = (add._Min - mean) / sigma;
    double zscore_max = (add._Max - mean) / sigma;
    double output_min = ((IsNaN(bg) || bg <= 0.) ? 0. : bg + .01);
    double output_max = output_min + (zscore_max - zscore_min);
    double zscore_mul = (output_max - output_min) / (zscore_max - zscore_min);
    double zscore_add = output_min - zscore_mul * zscore_min;
    add._Scale     = zscore_mul * add._Scale;
    add._Intercept = zscore_mul * add._Intercept + zscore_add;
  } else if (normalization == Normalization_MeanStDev) {
    double mean, sigma;
    NormalDistribution::Calculate(mean, sigma, num, data, mask.get());
    if (fequal(avg_sigma, 0.) || fequal(sigma, 0.)) return false;
    add._Scale = avg_sigma / sigma;
    add._Intercept = avg_mean - add._Scale / sigma;
  } else if (normalization == Normalization_UnitRange) {
    double range = add._Max - add._Min;
    if (fequal(range, 0.)) return false;
    add._Scale = 1. / (add._Max - add._Min);
    add._Intercept = - add._Scale * add._Min;
  }
  // Mask array no longer needed
  mask.reset();
  // Set interpolation function
  if (label > 0) interpolation = Interpolation_NN;
  UniquePtr<InputImageFunction> interp;
  interp.reset(InputImageFunction::New(interpolation, &image));
  interp->DefaultValue(bg);
  interp->Initialize();
  // Precompute displacement vector field
  #ifdef HAVE_MIRTK_Transformation
    bool cache = false;
    GenericImage<double> disp;
    const int ndofs = static_cast<int>(dof.size());
    for (int n = ndofs - 1; n >= 0; --n) {
      if (dof[n]->RequiresCachingOfDisplacements()) {
        cache = true;
        break;
      }
    }
    if (cache) {
      const double t  = average.GetTOrigin();
      const double t0 = -1.;
      disp.Initialize(average.Attributes(), 3);
      // Note: Input transformation is from image to average!
      for (int n = ndofs - 1; n >= 0; --n) {
        // Attention: Must specify both t and t0 to call the overloaded
        //            [Inverse]Displacement function which does not reset
        //            the current displacement vectors to zero before!
        if (inv[n]) dof[n]->Displacement       (disp, t, t0);
        else        dof[n]->InverseDisplacement(disp, t, t0);
      }
    }
    add._Transformation = &dof;
    add._Invert         = &inv;
    add._Displacement   = (cache ? &disp : nullptr);
  #endif // HAVE_MIRTK_Transformation
  // Transform, normalize, and add image values
  add._Average = &average;
  add._Weight  = weight;
  add._Label   = label;
  add._Image   = interp.get();
  ParallelForEachVoxel(average.Attributes(), average, add);
  return true;
}

// -----------------------------------------------------------------------------
/// Voxel function for parallel convolution of 2D image with discrete Laplace operator
class LaplaceOperator2D : public VoxelFunction
{
  const AverageImage *_Image;
  const double2       _NeighborWeight;
  const double        _CenterWeight;

public:

  LaplaceOperator2D(const AverageImage *image)
  :
    _Image(image),
    _NeighborWeight(make_double2(_Image->XSize() * _Image->XSize(),
                                 _Image->YSize() * _Image->YSize())),
    _CenterWeight(2. * _NeighborWeight.x + 2. * _NeighborWeight.y)
  {}

  void operator ()(int i, int j, int, int, AverageImage::VoxelType *out) const
  {
    const AverageImage &f = *_Image;
    if (f.IsInsideForeground(i, j)) {
      double invalue  = f(i, j);
      double outvalue = 0.;
      if (f.IsInsideForeground(i-1, j)) outvalue += f(i-1, j) * _NeighborWeight.x;
      else                              outvalue += invalue   * _NeighborWeight.x;
      if (f.IsInsideForeground(i+1, j)) outvalue += f(i+1, j) * _NeighborWeight.x;
      else                              outvalue += invalue   * _NeighborWeight.x;
      if (f.IsInsideForeground(i, j-1)) outvalue += f(i, j-1) * _NeighborWeight.y;
      else                              outvalue += invalue   * _NeighborWeight.y;
      if (f.IsInsideForeground(i, j+1)) outvalue += f(i, j+1) * _NeighborWeight.y;
      else                              outvalue += invalue   * _NeighborWeight.y;
      *out = voxel_cast<AverageImage::VoxelType>(outvalue - invalue * _CenterWeight);
    }
  }
};

// -----------------------------------------------------------------------------
/// Voxel function for parallel convolution of 3D image with discrete Laplace operator
class LaplaceOperator3D : public VoxelFunction
{
  const AverageImage *_Image;
  const double3       _NeighborWeight;
  const double        _CenterWeight;

public:

  LaplaceOperator3D(const AverageImage *image)
  :
    _Image(image),
    _NeighborWeight(make_double3(_Image->XSize() * _Image->XSize(),
                                 _Image->YSize() * _Image->YSize(),
                                 _Image->ZSize() * _Image->ZSize())),
    _CenterWeight(2. * _NeighborWeight.x + 2. * _NeighborWeight.y + 2. * _NeighborWeight.z)
  {}

  void operator ()(int i, int j, int k, int, AverageImage::VoxelType *out) const
  {
    const AverageImage &f = *_Image;
    if (f.IsInsideForeground(i, j, k)) {
      double invalue  = f(i, j, k);
      double outvalue = 0.;
      if (f.IsInsideForeground(i-1, j, k)) outvalue += f(i-1, j, k) * _NeighborWeight.x;
      else                                 outvalue += invalue      * _NeighborWeight.x;
      if (f.IsInsideForeground(i+1, j, k)) outvalue += f(i+1, j, k) * _NeighborWeight.x;
      else                                 outvalue += invalue      * _NeighborWeight.x;
      if (f.IsInsideForeground(i, j-1, k)) outvalue += f(i, j-1, k) * _NeighborWeight.y;
      else                                 outvalue += invalue      * _NeighborWeight.y;
      if (f.IsInsideForeground(i, j+1, k)) outvalue += f(i, j+1, k) * _NeighborWeight.y;
      else                                 outvalue += invalue      * _NeighborWeight.y;
      if (f.IsInsideForeground(i, j, k-1)) outvalue += f(i, j, k-1) * _NeighborWeight.z;
      else                                 outvalue += invalue      * _NeighborWeight.z;
      if (f.IsInsideForeground(i, j, k+1)) outvalue += f(i, j, k+1) * _NeighborWeight.z;
      else                                 outvalue += invalue      * _NeighborWeight.z;
      *out = voxel_cast<AverageImage::VoxelType>(outvalue - invalue * _CenterWeight);
    }
  }
};

// -----------------------------------------------------------------------------
/// Convolve intensity image by the finite difference Laplacian kernel
AverageImage ApplyLaplaceOperator(const AverageImage &image)
{
  AverageImage output(image);
  output.ClearBackgroundValue();
  if (output.Z() == 1 || output.ZSize() == 0.) {
    ParallelForEachVoxel(LaplaceOperator2D(&image), output.Attributes(), output);
  } else {
    ParallelForEachVoxel(LaplaceOperator3D(&image), output.Attributes(), output);
  }
  return output;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  InputImage    sequence, image;
  Array<string> image_name;
  Array<double> image_weight;
  int           nimages;

  #ifdef HAVE_MIRTK_Transformation
    Array<UniquePtr<Transformation> > dofs;
    Array<Array<string> >             imdof_name;
    Array<Array<bool> >               imdof_invert;
    Array<bool>                       invert;
  #endif // HAVE_MIRTK_Transformation

  InitializeIOLibrary();

  // ---------------------------------------------------------------------------
  // Command-line parsing

  // Required positional argument(s)
  REQUIRES_POSARGS(1);
  const char *output_name = POSARG(1);

  // Optional positional argument(s)
  int nposarg;
  for (nposarg = NUM_POSARGS + 1; nposarg < argc; ++nposarg) {
    bool img = (strcmp(argv[nposarg], "-image") == 0);
    bool inv = (strcmp(argv[nposarg], "-dof_i") == 0);
    bool dof = inv || (strcmp(argv[nposarg], "-dof")   == 0);
    if (!img && !dof) {
      if (argv[nposarg][0] == '-') break;
      size_t l = strlen(argv[nposarg]);
      dof = ((l > 4u && strcmp(&argv[nposarg][l-4u], ".dof")    == 0) ||
             (l > 7u && strcmp(&argv[nposarg][l-7u], ".dof.gz") == 0));
    } else {
      ++nposarg;
      if (nposarg == argc) {
        PrintHelp(EXECNAME);
        cout << endl;
        FatalError("Option " << argv[nposarg-1] << " requires an argument");
      }
    }
    if (dof) {
      #ifdef HAVE_MIRTK_Transformation
        imdof_name  .resize(image_name.size());
        imdof_invert.resize(image_name.size());
        imdof_name  .back().push_back(argv[nposarg]);
        imdof_invert.back().push_back(inv);
      #else // HAVE_MIRTK_Transformation
        FatalError("Cannot apply an image transformation, rebuild with Transformation module enabled!");
      #endif // HAVE_MIRTK_Transformation
    } else {
      image_name.push_back(argv[nposarg]);
    }
  }
  --nposarg;

  nimages = static_cast<int>(image_name.size());
  image_weight.resize(nimages, 1.0);
  #ifdef HAVE_MIRTK_Transformation
    imdof_name  .resize(nimages);
    imdof_invert.resize(nimages);
  #endif // HAVE_MIRTK_Transformation

  // Parse arguments
  const char        *image_list_name  = nullptr;
  const char        *image_list_delim = nullptr;
  string             reference_name;
  double             padding          = numeric_limits<double>::quiet_NaN();
  InterpolationMode  interpolation    = Interpolation_Linear;
  NormalizationMode  normalization    = Normalization_None;
  RescalingMode      rescaling        = Rescaling_None;
  double             output_min       = -inf;
  double             output_max       = +inf;
  bool               invert_dofs      = false;
  int                label            = -1;
  ImageDataType      dtype            = MIRTK_VOXEL_FLOAT;
  int                margin           = -1;
  bool               sharpen          = false;
  double             dx = .0, dy = .0, dz = .0;

  for (ARGUMENTS_AFTER(nposarg)) {
    if (OPTION("-images")) {
      image_list_name = ARGUMENT;
    }
    else if (OPTION("-delim") || OPTION("-delimiter")) {
      image_list_delim = ARGUMENT;
    }
    else if (OPTION("-reference") || OPTION("-target")) {
      const string farg = ARGUMENT;
      const string larg = ToLower(farg);
      if (larg == "maxsize" || larg == "max-size") {
        reference_name = "max-size";
      } else if (larg == "maxvol" || larg == "maxvolume" || larg == "max-vol" || larg == "max-volume" ||
                 larg == "maxarea" || larg == "max-area" || larg == "maxspace" || larg == "max-space") {
        reference_name = "max-space";
      } else {
        reference_name = farg;
      }
    }
    else if (OPTION("-spacing") || OPTION("-voxel-size") || OPTION("-size")) {
      PARSE_ARGUMENT(dx);
      dy = dz = dx;
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(dy);
        if (HAS_ARGUMENT) PARSE_ARGUMENT(dz);
        else              dz = .0;
      }
    }
    else if (OPTION("-normalization") || OPTION("-normalize")) {
      if (HAS_ARGUMENT) {
        const string arg = ToLower(ARGUMENT);
        bool bval;
        if (FromString(arg, bval)) {
          normalization = (bval ? Normalization_Mean : Normalization_None);
        } else if (arg == "none") {
          normalization = Normalization_None;
        } else if (arg == "mean") {
          normalization = Normalization_Mean;
        } else if (arg == "median") {
          normalization = Normalization_Median;
        } else if (arg == "zscore" || arg == "z-score") {
          normalization = Normalization_ZScore;
        } else if (arg == "unit") {
          normalization = Normalization_UnitRange;
        } else if (arg == "dist") {
          normalization = Normalization_MeanStDev;
        } else {
          FatalError("Invalid -normalization mode: " << arg);
        }
      } else {
        normalization = Normalization_Mean;
      }
    }
    else if (OPTION("-rescaling") || OPTION("-rescale")) {
      if (HAS_ARGUMENT) {
        const string arg = ToLower(ARGUMENT);
        if (HAS_ARGUMENT) {
          if (!FromString(arg, output_min)) {
            FatalError("Option -resaling must either have one string argument or two floating point numbers as arguments");
          }
          PARSE_ARGUMENT(output_max);
          rescaling = Rescaling_Range;
        } else {
          bool bval;
          if (FromString(arg, bval)) {
            rescaling = (bval ? Rescaling_MeanStDev : Rescaling_None);
          } else {
            const string larg = ToLower(arg);
            if (larg == "none") {
              rescaling = Rescaling_None;
            } else if (arg == "unit") {
              output_min = 0.;
              output_max = 1.;
              rescaling = Rescaling_Range;
            } else if (arg == "dist") {
              rescaling = Rescaling_MeanStDev;
            } else {
              FatalError("Invalid -rescaling mode: " << arg);
            }
          }
        }
      } else {
        rescaling = Rescaling_MeanStDev;
      }
    }
    else HANDLE_BOOL_OPTION(sharpen);
    else HANDLE_BOOLEAN_OPTION("invert", invert_dofs);
    else if (OPTION("-padding")) PARSE_ARGUMENT(padding);
    else if (OPTION("-margin")) PARSE_ARGUMENT(margin);
    else if (OPTION("-label")) PARSE_ARGUMENT(label);
    else if (OPTION("-interpolation") || OPTION("-interp")) PARSE_ARGUMENT(interpolation);
    else if (OPTION("-datatype") || OPTION("-dtype") || OPTION("-type")) {
      PARSE_ARGUMENT(dtype);
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (label > 0 && IsNaN(padding)) padding = .0;

  // ---------------------------------------------------------------------------
  // Collect (further) input image meta-data...

  // ...from files specified as positional arguments
  if (image_name.size() == 1) {
    if (verbose) cout << "Reading sequence ... ", cout.flush();
    sequence.Read(image_name.front().c_str());
    nimages = sequence.T();
    if (nimages < 2) {
      if (verbose) cout << endl;
      FatalError("Input sequence contains only one temporal frame!");
    }
    if (verbose) cout << " done\n" << endl;
    image_weight.resize(nimages, 1.0);
  }

  // ...from files specified in list file
  if (image_list_name) {
    Array<string>         names;
    Array<Array<string> > imdofs;
    Array<double>         weights;
    if (image_list_delim == nullptr) {
      const string ext = Extension(image_list_name);
      if      (ext == ".csv") image_list_delim = ",";
      else if (ext == ".tsv") image_list_delim = "\t";
      else                    image_list_delim = " ";
    }
    const int n = read_image_list_file(image_list_name, names, imdofs, weights, image_list_delim);
    if (n == 0) {
      if (verbose) cout << endl;
      FatalError("Failed to parse input list file " << image_list_name << "!");
    }
    if (verbose) cout << "Found " << n << " images in list file " << image_list_name << endl;
    nimages += n;
    image_name  .insert(image_name  .end(), names  .begin(), names  .end());
    image_weight.insert(image_weight.end(), weights.begin(), weights.end());
    imdof_name  .insert(imdof_name  .end(), imdofs .begin(), imdofs .end());
    imdof_invert.resize(nimages);
    for (int i = nimages - n; i < nimages; ++i) {
      imdof_invert[i].resize(imdof_name[i].size(), invert_dofs);
    }
  }

  // Check that any input image is given
  if (nimages < 1) {
    PrintHelp(EXECNAME);
    FatalError("No input image(s) specified!");
  }

  // Normalize weights
  double wsum = 0.;
  for (int n = 0; n < nimages; ++n) {
    wsum += image_weight[n];
  }
  if (wsum == 0.) {
    FatalError("Sum of image weights is zero!");
  }
  for (int n = 0; n < nimages; ++n) {
    image_weight[n] /= wsum;
  }

  // ---------------------------------------------------------------------------
  // Collect input properties if needed
  double min_input = -inf, max_input = +inf;
  double sum_of_means  = 0.;
  double sum_of_sigmas = 0.;
  Array<ImageAttributes> attrs;
  bool calc_input_distribution = false;
  bool collect_attrs           = false;
  if (reference_name.empty() || reference_name == "max-size" || reference_name == "max-space") {
    collect_attrs = sequence.IsEmpty();
  }
  if (normalization == Normalization_MeanStDev || rescaling == Rescaling_MeanStDev) {
    calc_input_distribution = true;
  }
  if (collect_attrs || calc_input_distribution) {
    if (verbose) {
      cout << "Determine";
      if (collect_attrs) {
        if (reference_name == "max-size" || reference_name == "max-space") {
          cout << " image with " << reference_name;
        } else {
          cout << " common field-of-view";
        }
      }
      if (calc_input_distribution) {
        if (collect_attrs) cout << ", and";
        cout << " average intensity distribution";
      }
      cout << "...";
      cout.flush();
    }
    for (int n = 0; n < nimages; ++n) {
      #ifdef HAVE_MIRTK_Transformation
        Read(image_name[n], imdof_name[n], imdof_invert[n], image, dofs, invert);
        dofs.clear(), invert.clear(); // unused here
      #else // HAVE_MIRTK_Transformation
        image.Read(image_name[n].c_str());
      #endif // HAVE_MIRTK_Transformation
      if (image.T() > 1) {
        if (verbose) cout << " failed" << endl;
        FatalError("Image " << (n + 1) << " has four dimensions!");
      }
      if (collect_attrs) {
        attrs.push_back(OrthogonalFieldOfView(image.Attributes()));
      }
      image.PutBackgroundValueAsDouble(padding);
      if (calc_input_distribution) {
        double mean, sigma, min_value, max_value;
        const int num = image.NumberOfSpatialVoxels();
        UniquePtr<bool[]> mask = ForegroundMaskArray(&image);
        NormalDistribution::Calculate(mean, sigma, num, image.Data(), mask.get());
        Extrema::Calculate(min_value, max_value, num, image.Data(), mask.get());
        if (min_value < min_input) min_input = min_value;
        if (max_value > max_input) max_input = max_value;
        sum_of_means  += mean;
        sum_of_sigmas += sigma;
      }
    }
    if (verbose) cout << " done" << endl;
  }

  const double avg_mean  = sum_of_means  / nimages;
  const double avg_sigma = sum_of_sigmas / nimages;

  // ---------------------------------------------------------------------------
  // Average image attributes (field-of-view)
  ImageAttributes fov;
  if (!reference_name.empty() && reference_name != "max-size" && reference_name != "max-space") {
    if (verbose) cout << "Using attributes of reference image" << endl;
    BinaryImage reference(reference_name.c_str());
    fov = reference.Attributes();
  } else if (sequence.IsEmpty()) {
    if (reference_name == "max-size") {
      int vox = 0;
      for (const auto &attr : attrs) {
        if (attr.NumberOfSpatialPoints() > vox) {
          vox = attr.NumberOfSpatialPoints();
          fov = attr;
        }
      }
    } else if (reference_name == "max-space") {
      double vol = 0.;
      for (const auto &attr : attrs) {
        if (attr.Space() > vol) {
          vol = attr.Space();
          fov = attr;
        }
      }
    } else {
      fov = OverallFieldOfView(attrs);
    }
  } else {
    fov = sequence.Attributes();
    fov._t = 1, fov._dt *= nimages;
  }

  // Set desired output voxel size
  if (dx > .0) {
    fov._x  = iceil(fov._x * fov._dx / dx);
    fov._dx = dx;
  }
  if (dy > .0) {
    fov._y  = iceil(fov._y * fov._dy / dy);
    fov._dy = dy;
  }
  if (dz > .0) {
    fov._z  = iceil(fov._z * fov._dz / dz);
    fov._dz = dz;
  }

  // ---------------------------------------------------------------------------
  // Compute average image
  AverageImage average(fov);
  average = NaN;
  average.PutBackgroundValueAsDouble(NaN);
  for (int n = 0; n < nimages; ++n) {
    if (sequence.IsEmpty()) {
      if (verbose) {
        cout << "Add image " << setw(3) << (n+1) << " out of " << nimages << "... ";
        cout.flush();
      }
      #ifdef HAVE_MIRTK_Transformation
        Read(image_name[n], imdof_name[n], imdof_invert[n], image, dofs, invert);
      #else // HAVE_MIRTK_Transformation
        image.Read(image_name[n].c_str());
      #endif // HAVE_MIRTK_Transformation
    } else {
      if (verbose) {
        cout << "Add frame " << setw(3) << (n+1) << " out of " << nimages << "... ";
        cout.flush();
      }
      sequence.GetFrame(image, n);
      #ifdef HAVE_MIRTK_Transformation
        dofs.clear();
        invert.clear();
      #endif
    }
    image.PutBackgroundValueAsDouble(padding, true);

    bool ok;
    #ifdef HAVE_MIRTK_Transformation
      ok = Add(average, image, dofs, invert, image_weight[n], label, interpolation, normalization, avg_mean, avg_sigma);
      dofs.clear(), invert.clear();
    #else
      ok = Add(average, image, image_weight[n], label, interpolation, normalization, avg_mean, avg_sigma);
    #endif
    if (ok) {
      if (verbose) cout << " done" << endl;
    } else {
      if (verbose) cout << " skip" << endl;
      if (normalization == Normalization_UnitRange || normalization == Normalization_ZScore || normalization == Normalization_MeanStDev) {
        Warning("Input image " << image_name[n] << " either contains no or constant foreground values!");
      } else {
        Warning("Input image " << image_name[n] << " contains no foreground values!");
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Replace NaN values by suitable background value
  AverageImage::VoxelType min_value, max_value;
  average.GetMinMax(min_value, max_value);
  AverageImage::VoxelType bg = min_value - 1.;
  if (bg > 0.) bg = 0.;
  for (int idx = 0; idx < average.NumberOfVoxels(); ++idx) {
    if (IsNaN(average(idx))) {
      average(idx) = bg;
    }
  }
  average.PutBackgroundValueAsDouble(bg);

  // ---------------------------------------------------------------------------
  // Crop/pad average image
  if (margin >= 0) {
    if (verbose) {
      cout << "Cropping average image adding margin of " << margin << " background layers...";
      cout.flush();
    }
    average.CropPad(margin);
    if (verbose) cout << " done" << endl;
  }

  // ---------------------------------------------------------------------------
  // Sharpen average image (cf. itkLaplacianSharpeningImageFilter, ANTs AverageImages)
  if (sharpen) {
    // TODO: Implement this as reusable MIRTK image filter
    if (verbose) {
      cout << "Sharpening average image using Laplacian operator...";
      cout.flush();
    }
    const int num = average.NumberOfSpatialVoxels();
    UniquePtr<bool[]> mask = ForegroundMaskArray(&average);
    typedef AverageImage::VoxelType Real;
    AverageImage laplace = ApplyLaplaceOperator(average);
    Real mean = Mean::Calculate(num, average.Data(), mask.get());
    Real min_value, max_value;
    average.GetMinMax(min_value, max_value);
    laplace.PutMinMax(min_value, max_value);
    if (debug) {
      average.Write("debug_average_intensity.nii.gz");
      laplace.Write("debug_average_laplacian.nii.gz");
    }
    for (int idx = 0; idx < num; ++idx) {
      if (mask[idx]) average(idx) -= laplace(idx);
    }
    if (debug) {
      average.Write("debug_average_sharpened.nii.gz");
    }
    Real shift = mean - Mean::Calculate(num, average.Data(), mask.get());
    for (int idx = 0; idx < num; ++idx) {
      if (mask[idx]) {
        average(idx) = clamp(average(idx) + shift, min_value, max_value);
      }
    }
    if (verbose) cout << " done" << endl;
  }

  // ---------------------------------------------------------------------------
  // Replace background value by specified padding value
  if (!IsNaN(padding)) {
    AverageImage::VoxelType bg;
    bg = voxel_cast<AverageImage::VoxelType>(padding);
    for (int idx = 0; idx < average.NumberOfVoxels(); ++idx) {
      if (average.IsBackground(idx)) average(idx) = bg;
    }
    average.PutBackgroundValueAsDouble(padding);
  }

  // ---------------------------------------------------------------------------
  // Rescale to input mean and standard deviation
  if (rescaling != Rescaling_None) {
    double scale = 1., shift = 0.;
    const int num = average.NumberOfSpatialVoxels();
    UniquePtr<bool[]> mask = ForegroundMaskArray(&average);
    if (rescaling == Rescaling_MeanStDev) {
      if (verbose) {
        cout << "Rescaling average to input mean and standard deviation...";
        cout.flush();
      }
      double mean, sigma;
      NormalDistribution::Calculate(mean, sigma, num, average.Data(), mask.get());
      if (!fequal(avg_sigma, 0.) && !fequal(sigma, 0.)) {
        scale = avg_sigma / sigma;
      }
      shift = avg_mean - scale * mean;
    } else if (rescaling == Rescaling_Range) {
      if (verbose) {
        cout << "Rescaling average to output range [" << output_min << ", " << output_max << "]...";
        cout.flush();
      }
      double min_value, max_value;
      Extrema::Calculate(min_value, max_value, num, average.Data(), mask.get());
      if (!fequal(max_value, min_value)) {
        scale = (output_max - output_min) / (max_value - min_value);
      }
      shift = output_min - scale * min_value;
    }
    double min_value = +inf;
    if (!fequal(scale, 1.) || !fequal(shift, 0.)) {
      double value;
      for (int idx = 0; idx < num; ++idx) {
        if (mask[idx]) {
          value = clamp(scale * double(average(idx)) + shift, min_input, max_input);
          if (value < min_value) min_value = value;
          average(idx) = voxel_cast<AverageImage::VoxelType>(value);
        }
      }
    }
    if (verbose) {
      cout << " done" << endl;
    }
    double bg = average.GetBackgroundValueAsDouble();
    if (min_value < bg) {
      bg = min(0., min_value - 1.);
      AverageImage::VoxelType padding = static_cast<AverageImage::VoxelType>(bg);
      if (verbose) {
        cout << "Input -padding value larger than rescaled minimum, pad output with: " << padding << endl;
      }
      for (int idx = 0; idx < num; ++idx) {
        if (!mask[idx]) average(idx) = padding;
      }
      average.PutBackgroundValueAsDouble(bg);
    }
  }

  // ---------------------------------------------------------------------------
  // Write average image
  if (average.GetDataType() != dtype) {
    UniquePtr<BaseImage> output(BaseImage::New(dtype));
    *output = average;
    output->Write(output_name);
  } else {
    average.Write(output_name);
  }

  return 0;
}

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

#ifdef HAVE_MIRTK_Transformation
#  include "mirtk/Transformation.h"
#  include "mirtk/HomogeneousTransformation.h"
#  include "mirtk/FluidFreeFormTransformation.h"
#endif

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

/// Print program usage information
void PrintHelp(const char* name)
{
  cout << endl;
#ifdef HAVE_MIRTK_Transformation
  cout << "Usage: " << name << " <output> -images <images.lst> [options]" << endl;
  cout << "       " << name << " <output> -image <image1> [-dof <dof1>...] -image <image2> [-dof <dof2>...]... [options]" << endl;
  cout << "       " << name << " <output> <image1> [<dof1>...] <image2> [<dof2>...]... [options]" << endl;
#else
  cout << "Usage: " << name << " <output> -images <images.lst> [options]" << endl;
  cout << "       " << name << " <output> -image <image1> -image <image2> ... [options]" << endl;
  cout << "       " << name << " <output> <image1> <image2>... [options]" << endl;
#endif
  cout << "       " << name << " <output> <sequence> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Computes voxelwise average image of the given (transformed) input images." << endl;
#ifdef HAVE_MIRTK_Transformation
  cout << endl;
  cout << "  The input images are optionally transformed by an image-to-average space transformation" << endl;
  cout << "  (:option:`-dof`) or its inverse (:option:`-dof_i`)." << endl;
#endif
  cout << endl;
  cout << "Required arguments:" << endl;
  cout << "  <output>                 Voxel-wise average image." << endl;
  cout << "  -images <file>           Text file with N lines containing the file name of each input image," << endl;
  cout << "                           optionally a transformation file name (see :option:`-dof` and :option:`-dof_i`)," << endl;
  cout << "                           optionally followed by a weight for a weighted output average (default weight is 1)." << endl;
  cout << "                           The first line of the text file must specify the common base directory" << endl;
  cout << "                           of all relative image and transformation file paths occurring on the" << endl;
  cout << "                           subsequent N lines. A path starting with './' must be relative to the" << endl;
  cout << "                           directory containing the input text file itself." << endl;
  cout << "  -delim, -delimiter <c>   Delimiter used in :option:`-images` file." << endl;
  cout << "                           (default: ',' for .csv, '\\t' for .tsv, and ' ' otherwise)" << endl;
  cout << "  -invert                  Invert transformations specified in :option:`-images` file." << endl;
  cout << "  -image <file>            A single input image." << endl;
#ifdef HAVE_MIRTK_Transformation
  cout << "  -dof <file>              Specifies a transformation to be applied to the preceeding image. (default: none)" << endl;
  cout << "  -dof_i <file>            Specifies a transformation whose inverse is to be applied to the" << endl;
  cout << "                           preceeding image similar to :option:`-dof`. (default: none)" << endl;
#endif
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -reference <file>          Reference image for output image attributes.\n";
  cout << "  -size <dx> [<dy> [<dz>]]   Voxel size of intensity average image. (default: average of input images)" << endl;
  cout << "  -padding <value>           Input padding and output background value. No input padding if not specified." << endl;
  cout << "                             Output background value zero by default or minimum average intensity minus 1. (default: 0)" << endl;
  cout << "  -interp <mode>             Interpolation mode, e.g., NN, Linear, BSpline, Cubic, Sinc. (default: Linear)" << endl;
  cout << "  -label <value>             Segmentation label of which to create an average probability map." << endl;
  cout << "  -datatype, -dtype, -type char|uchar|short|float|double" << endl;
  cout << "      Data type of output image. The intermediate average image always has floating point data type." << endl;
  cout << "      When this option is given, this average is cast to the respective output type before writing" << endl;
  cout << "      the result to the output image file. (default: float)" << endl;
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
struct AddVoxelValueToAverage : public VoxelFunction
{
  AverageImage       *_Average;
  InputImageFunction *_Image;
  double              _Weight;
  int                 _Label;

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
      const T bg = static_cast<T>(_Average->GetBackgroundValueAsDouble());
      if ((IsNaN(bg) && IsNaN(*avg)) || (*avg == bg)) {
        *avg = static_cast<T>(_Weight * _Image->Evaluate(x, y, z));
      } else {
        *avg += static_cast<T>(_Weight * _Image->Evaluate(x, y, z));
      }
    }
  }
};

// -----------------------------------------------------------------------------
void Add(AverageImage &average, InputImage &image,
         #ifdef HAVE_MIRTK_Transformation
           const Array<UniquePtr<Transformation> > &dof,
           const Array<bool>                       &inv,
         #endif // HAVE_MIRTK_Transformation
         double            weight        = 1.0,
         int               label         = -1,
         InterpolationMode interpolation = Interpolation_Linear)
{
  AddVoxelValueToAverage add;
  if (label > 0) interpolation = Interpolation_NN;
  UniquePtr<InputImageFunction> interp;
  interp.reset(InputImageFunction::New(interpolation, &image));
  interp->Initialize();
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
  add._Average = &average;
  add._Weight  = weight;
  add._Label   = label;
  add._Image   = interp.get();
  ParallelForEachVoxel(average.Attributes(), average, add);
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
  const char        *reference_name   = nullptr;
  double             padding          = numeric_limits<double>::quiet_NaN();
  InterpolationMode  interpolation    = Interpolation_Linear;
  bool               invert_dofs      = false;
  int                label            = -1;
  ImageDataType      dtype            = MIRTK_VOXEL_FLOAT;
  int                margin           = -1;
  double             dx = .0, dy = .0, dz = .0;

  for (ARGUMENTS_AFTER(nposarg)) {
    if (OPTION("-images")) {
      image_list_name = ARGUMENT;
    }
    else if (OPTION("-delim") || OPTION("-delimiter")) {
      image_list_delim = ARGUMENT;
    }
    else if (OPTION("-reference")) {
      reference_name = ARGUMENT;
    }
    else if (OPTION("-size")) {
      PARSE_ARGUMENT(dx);
      dy = dz = dx;
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(dy);
        if (HAS_ARGUMENT) PARSE_ARGUMENT(dz);
        else              dz = .0;
      }
    }
    else if (OPTION("-padding"))   PARSE_ARGUMENT(padding);
    else if (OPTION("-invert"))    invert_dofs = true;
    else if (OPTION("-margin"))    PARSE_ARGUMENT(margin);
    else if (OPTION("-label"))     PARSE_ARGUMENT(label);
    else if (OPTION("-interp"))    PARSE_ARGUMENT(interpolation);
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

  // Compute sum of weights
  double wsum = .0;
  for (int n = 0; n < nimages; ++n) {
    wsum += image_weight[n];
  }

  // Compute voxel-wise average of (interpolated) input intensities on
  // common (mapped) image domain of all input images; for example,
  // used for anatomical intensity atlas construction

  // Compute field-of-view which contains all images
  if (verbose) cout << "Determine field-of-view which contains all images...", cout.flush();
  ImageAttributes fov;
  if (reference_name) {
    GreyImage reference(reference_name);
    fov = reference.Attributes();
  } else if (sequence.IsEmpty()) {
    Array<ImageAttributes> attr;
    for (int n = 0; n < nimages; ++n) {
      #ifdef HAVE_MIRTK_Transformation
        Read(image_name[n], imdof_name[n], imdof_invert[n], image, dofs, invert);
        dofs.clear(); // unused here
        invert.clear();
      #else // HAVE_MIRTK_Transformation
        image.Read(image_name[n].c_str());
      #endif // HAVE_MIRTK_Transformation
      if (image.T() > 1) {
        if (verbose) cout << " failed" << endl;
        FatalError("Image " << (n+1) << " has four dimensions!");
      }
      attr.push_back(OrthogonalFieldOfView(image.Attributes()));
    }
    fov = OverallFieldOfView(attr);
  } else {
    fov = sequence.Attributes();
    fov._t = 1, fov._dt *= nimages;
  }
  if (verbose) cout << " done" << endl;

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

  // Compute average image
  AverageImage average(fov);
  average = padding;
  average.PutBackgroundValueAsDouble(padding);
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

    #ifdef HAVE_MIRTK_Transformation
      Add(average, image, dofs, invert, image_weight[n], label, interpolation);
      dofs.clear();
      invert.clear();
    #else
      Add(average, image, image_weight[n], label, interpolation);
    #endif

    if (verbose) cout << " done" << endl;
  }
  if (wsum > .0) average /= wsum;

  // Crop/pad average image
  if (margin >= 0) {
    average.CropPad(margin);
  }

  // Replace NaN background values
  if (IsNaN(average.GetBackgroundValueAsDouble())) {
    AverageImage::VoxelType min_value, max_value;
    average.GetMinMax(min_value, max_value);
    padding = min(.0, double(min_value) - 1.0);
    for (int idx = 0; idx < average.NumberOfVoxels(); ++idx) {
      if (IsNaN(average(idx))) {
        average(idx) = voxel_cast<AverageImage::VoxelType>(padding);
      }
    }
    average.PutBackgroundValueAsDouble(padding);
  }

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

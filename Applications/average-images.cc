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

#include <mirtkCommon.h>
#include <mirtkOptions.h>

#include <mirtkIOConfig.h>

#include <mirtkMatrix.h>
#include <mirtkGenericImage.h>
#include <mirtkVoxelFunction.h>
#include <mirtkInterpolateImageFunction.h>

#ifdef HAVE_MIRTK_Transformation
#  include <mirtkTransformation.h>
#  include <mirtkHomogeneousTransformation.h>
#  include <mirtkFluidFreeFormTransformation.h>
#endif

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

/// Print program usage information
void PrintHelp(const char* name)
{
  cout << endl;
  cout << "Usage: " << name << " <output> -images <images.lst> [options]" << endl;
  cout << "       " << name << " <output> -image <image1> [-dof <dof1>] -image <image2> [-dof <dof2>]... [options]" << endl;
  cout << "       " << name << " <output> <image1> [<dof1>] <image2> [<dof2>]... [options]" << endl;
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
  cout << "                            optionally followed by a weight for a weighted output average (default weight is 1)." << endl;
  cout << "                           The first line of the text file must specify the common base directory" << endl;
  cout << "                           of all relative image and transformation file paths occurring on the" << endl;
  cout << "                           subsequent N lines. A path starting with './' must be relative to the" << endl;
  cout << "                           directory containing the input text file itself." << endl;
  cout << "  -image <file>            A single input image." << endl;
#ifdef HAVE_MIRTK_Transformation
  cout << "  -dof <file>              Specifies a transformation to be applied to the preceeding image. (default: none)" << endl;
  cout << "  -dof_i <file>            Specifies a transformation whose inverse is to be applied to the" << endl;
  cout << "                           preceeding image. (default: none)" << endl;
#endif
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -size <dx> [<dy> [<dz>]]   Voxel size of intensity average image. (default: average of input images)" << endl;
  cout << "  -padding <value>           Input padding and output background value. No input padding if not specified." << endl;
  cout << "                             Output background value zero by default or minimum average intensity minus 1. (default: 0)" << endl;
  cout << "  -interp <mode>             Interpolation mode, e.g., NN, Linear, BSpline, Cubic, Sinc. (default: Linear)" << endl;
  cout << "  -label <value>             Segmentation label of which to create an average probability map." << endl;
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

// Type of output image
typedef float                     OutputType;
typedef GenericImage<OutputType>  OutputImage;

// =============================================================================
// Auxiliary functions
// =============================================================================

/// Read image file names and times from input list file
int read_image_list_file(const char        *image_list_name,
                         Array<string>     &names,
                         Array<string>     &dofs,
                         Array<OutputType> &weights)
{
  string   base_dir = Directory(image_list_name); // Base directory containing image files
  ifstream iff(image_list_name);                  // List input file stream
  string   line;                                  // Input line
  string   part;                                  // Line part
  string   entry;                                 // Column entry
  bool     inquote = false;                       // Parsing parts of a quoted string
  int      l = 0;                                 // Line number

  names  .clear();
  dofs   .clear();
  weights.clear();

  // Read base directory for image files
  if (!getline(iff, line)) {
    cerr << "Error: Cannot parse list file " << image_list_name << endl;
    return 0;
  }

  if (!base_dir.empty() && line[0] != PATHSEP) {
    if (line != ".") base_dir += PATHSEP + line;
  } else {
    base_dir  = line;
  }
  l++;

  while (getline(iff, line)) {
    l++;
    // Ignore comment lines starting with # character
    if (line[0] == '#') continue;
    // Split line at space
    istringstream lss(line);
    int c = 0;
    while (getline(lss, part, ' ')) {
      if (inquote) entry += ' ';
      if (part.empty()) continue;
      if (part[0] == '"') {
        part.erase(0, 1);
        inquote = !inquote;
      }
      if (part[part.size()-1] == '"') {
        part.erase(part.size()-1, 1);
        inquote = !inquote;
      }
      entry += part;
      if (inquote) continue;
      // Add column entry to the proper output container
      if (c == 0) {
        if (!base_dir.empty() && entry[0] != PATHSEP) entry = base_dir + PATHSEP + entry;
        names.push_back(entry);
      } else if (c == 1 && !isdigit(entry[0])) {
        if (!base_dir.empty() && entry[0] != PATHSEP) entry = base_dir + PATHSEP + entry;
        dofs.push_back(entry);
      } else if (c == 2 || (c == 1 && isdigit(entry[0]))) {
        if (c == 1) dofs.push_back("");
        istringstream iss(entry);
        OutputType          weight;
        iss >> weight;
        if (!iss) {
          cerr << "Error: " << image_list_name << ":" << l << ": Failed to parse weight!" << endl;
          names  .clear();
          weights.clear();
          return 0;
        }
        weights.push_back(weight);
      } else {
        cerr << "Error: " << image_list_name << ":" << l << ": Too many entries!"<< endl;
        names  .clear();
        weights.clear();
        return 0;
      }
      // Parse next column
      entry.clear();
      c++;
    }
    if (c == 0) continue; // Skip empty lines
    // Add default for non-specified columns
    dofs   .resize(names.size());
    weights.resize(names.size(), 1.0);
  }

  return static_cast<int>(names.size());
}

// -----------------------------------------------------------------------------
#ifdef HAVE_MIRTK_Transformation
void Read(const string &image_name, const string    &imdof_name,
          InputImage   &image,      Transformation *&imdof,
          bool imdof_invert = false)
{
  image.Read(image_name.c_str());
  if (imdof_name.empty()) {
    imdof = NULL;
  } else {
    imdof = Transformation::New(imdof_name.c_str());
    HomogeneousTransformation   *lin   = dynamic_cast<HomogeneousTransformation   *>(imdof);
    FluidFreeFormTransformation *fluid = dynamic_cast<FluidFreeFormTransformation *>(imdof);
    if (fluid && imdof_invert) lin = fluid->GetGlobalTransformation();
    if (lin) {
      Matrix mat = lin->GetMatrix();
      if (imdof_invert) mat.Invert();
      image.PutAffineMatrix(mat, true);
      Delete(imdof);
    }
    if (fluid && lin) fluid->GetGlobalTransformation()->Reset();
  }
}
#endif // HAVE_MIRTK_Transformation

// -----------------------------------------------------------------------------
struct AddVoxelValueToAverage : public VoxelFunction
{
  OutputImage        *_Average;
  InputImageFunction *_Image;
  double              _Weight;
  int                 _Label;

  #ifdef HAVE_MIRTK_Transformation
    Transformation       *_Transformation;
    GenericImage<double> *_Displacement;
    bool                  _Invert;
  #endif

  void operator()(int i, int j, int k, int, OutputType *avg)
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
        if (_Invert) _Transformation->Transform(x, y, z);
        else         _Transformation->Inverse  (x, y, z);
      }
    #endif // HAVE_MIRTK_Transformation
    _Image->Input()->WorldToImage(x, y, z);
    if (_Label > 0) {
      if (static_cast<int>(_Image->Evaluate(x, y, z)) == _Label) {
        *avg += static_cast<OutputType>(_Weight);
      }
    } else {
      const OutputType bg = static_cast<OutputType>(_Average->GetBackgroundValueAsDouble());
      if ((IsNaN(bg) && IsNaN(*avg)) || (*avg == bg)) {
        *avg = static_cast<OutputType>(_Weight * _Image->Evaluate(x, y, z));
      } else {
        *avg += static_cast<OutputType>(_Weight * _Image->Evaluate(x, y, z));
      }
    }
  }
};

// -----------------------------------------------------------------------------
void Add(OutputImage &average, InputImage &image, int label = -1,
         #ifdef HAVE_MIRTK_Transformation
           Transformation   *transformation = NULL,
           bool              invert         = false,
         #endif // HAVE_MIRTK_Transformation
         OutputType        weight         = 1.0,
         InterpolationMode interpolation  = Interpolation_Linear)
{
  if (label > 0) interpolation = Interpolation_NN;
  AddVoxelValueToAverage add;
  #ifdef HAVE_MIRTK_Transformation
    GenericImage<double> disp;
    if (transformation && transformation->RequiresCachingOfDisplacements()) {
      disp.Initialize(average.Attributes(), 3);
      // Note: Input transformation is from image to average!
      if (invert) transformation->Displacement(disp);
      else        transformation->InverseDisplacement(disp);
    }
    add._Transformation = transformation;
    add._Displacement   = (disp.IsEmpty() ? NULL : &disp);
    add._Invert         = invert;
  #endif // HAVE_MIRTK_Transformation
  add._Average = &average;
  add._Weight  = weight;
  add._Label   = label;
  add._Image   = InputImageFunction::New(interpolation, &image);
  add._Image->Initialize();
  ParallelForEachVoxel(average.Attributes(), average, add);
  delete add._Image;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  InputImage        sequence, image;
  Array<string>     image_name;
  Array<OutputType> image_weight;
  int               nimages;

  #ifdef HAVE_MIRTK_Transformation
    Transformation *imdof = NULL;
    Array<string>   imdof_name;
    Array<bool>     imdof_invert;
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
        imdof_name  .back() = argv[nposarg];
        imdof_invert.back() = inv;
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
  const char        *image_list_name = NULL;
  const char        *reference_name  = NULL;
  double             padding         = numeric_limits<double>::quiet_NaN();
  InterpolationMode  interpolation   = Interpolation_Linear;
  bool               voxelwise       = false;
  int                label           = -1;
  int                margin          = -1;
  double             dx = .0, dy = .0, dz = .0;

  for (ARGUMENTS_AFTER(nposarg)) {
    if      (OPTION("-images" ))   image_list_name = ARGUMENT;
    else if (OPTION("-reference")) reference_name = ARGUMENT;
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
    else if (OPTION("-voxelwise")) voxelwise = true;
    else if (OPTION("-margin"))    PARSE_ARGUMENT(margin);
    else if (OPTION("-label"))     PARSE_ARGUMENT(label);
    else if (OPTION("-interp"))    PARSE_ARGUMENT(interpolation);
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
  }

  // ...from files specified in list file
  if (image_list_name) {
    Array<string> names;
    Array<string> imdofs;
    Array<OutputType>   weights;
    int n = read_image_list_file(image_list_name, names, imdofs, weights);
    if (n == 0) {
      if (verbose) cout << endl;
      FatalError("Failed to parse input list file " << image_list_name << "!");
    }
    if (verbose) cout << "Found " << n << " images in list file " << image_list_name << endl;
    nimages += n;
    image_name  .insert(image_name  .end(), names  .begin(), names  .end());
    image_weight.insert(image_weight.end(), weights.begin(), weights.end());
    imdof_name  .insert(imdof_name  .end(), imdofs .begin(), imdofs .end());
    imdof_invert.resize(nimages, false);
  }

  // Check that any input image is given
  if (nimages < 1) {
    PrintHelp(EXECNAME);
    FatalError("No input image(s) specified!");
  }

  // Compute sum of weights
  OutputType wsum = .0;
  for (int n = 0; n < nimages; ++n) wsum += image_weight[n];

  // Compute voxel-wise average of (interpolated) input intensities on
  // common (mapped) image domain of all input images; for example,
  // used for anatomical intensity atlas construction

  // Compute field-of-view which contains all images
  if (verbose) cout << "Determine field-of-view which contains all images...", cout.flush();
  ImageAttributes fov;
  if (reference_name) {
    GreyImage reference(reference_name);
    fov = reference.Attributes();
  } else {
    Array<ImageAttributes> attr;
    for (int n = 0; n < nimages; ++n) {
      #ifdef HAVE_MIRTK_Transformation
        Read(image_name[n], imdof_name[n], image, imdof, imdof_invert[n]);
        Delete(imdof);
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
  OutputImage average(fov);
  average = static_cast<OutputType>(padding);
  average.PutBackgroundValueAsDouble(padding);
  for (int n = 0; n < nimages; ++n) {
    if (sequence.IsEmpty()) {
      if (verbose) {
        cout << "Add image " << setw(3) << (n+1) << " out of " << nimages << "... ";
        cout.flush();
      }
      #ifdef HAVE_MIRTK_Transformation
        Read(image_name[n], imdof_name[n], image, imdof, imdof_invert[n]);
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
        imdof = NULL;
      #endif
    }

    Add(average, image, label, 
        #ifdef HAVE_MIRTK_Transformation
          imdof, imdof_invert[n],
        #endif
        image_weight[n], interpolation);

    #ifdef HAVE_MIRTK_Transformation
      Delete(imdof);
    #endif

    if (verbose) cout << " done" << endl;
  }
  if (wsum > .0) average /= wsum;

  // Crop/pad average image
  if (margin >= 0) average.CropPad(margin);

  // Replace NaN background values
  if (IsNaN(average.GetBackgroundValueAsDouble())) {
    OutputImage::VoxelType min_value, max_value;
    average.GetMinMax(min_value, max_value);
    padding = min(.0, double(min_value) - 1.0);
    for (int idx = 0; idx < average.NumberOfVoxels(); ++idx) {
      if (IsNaN(average(idx))) {
        average(idx) = voxel_cast<OutputImage::VoxelType>(padding);
      }
    }
    average.PutBackgroundValueAsDouble(padding);
  }

  // Write average image
  average.Write(output_name);


  return 0;
}

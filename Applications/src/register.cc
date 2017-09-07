/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2017 Andreas Schuh
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

#include "mirtk/NumericsConfig.h"
#include "mirtk/IOConfig.h"
#include "mirtk/TransformationConfig.h"
#include "mirtk/RegistrationConfig.h"
#if MIRTK_Registration_WITH_Deformable
#  include "mirtk/DeformableConfig.h"
#endif

#include "mirtk/GenericImage.h"
#include "mirtk/GenericRegistrationFilter.h"
#include "mirtk/GenericRegistrationLogger.h"
#include "mirtk/GenericRegistrationDebugger.h"

#include "mirtk/Transformation.h"
#include "mirtk/HomogeneousTransformation.h"
#include "mirtk/RigidTransformation.h"
#include "mirtk/SimilarityTransformation.h"
#include "mirtk/AffineTransformation.h"
#include "mirtk/LinearFreeFormTransformation3D.h"

#if MIRTK_Registration_WITH_PointSet
#  include "mirtk/PointSetIO.h"
#  include "vtkSmartPointer.h"
#  include "vtkPolyData.h"
#endif

using namespace mirtk;


// =============================================================================
// Version and help
// =============================================================================

/// Print command synopsis
void PrintSynopsis(const char *name)
{
  cout << "Usage: " << name << " -images <images.lst> [options]\n";
  cout << "       " << name << " -image <image1> [-dof <dof1>] -image <image2> [-dof <dof2>]... [options]\n";
  cout << "       " << name << " -image <image1> <image2>... [options]\n";
  cout << "       " << name << " <image1> <image2>... [options]\n";
  cout << "       " << name << " <image_sequence> [options]\n";
#if MIRTK_Registration_WITH_PointSet
  cout << "       " << name << " -pset <pointset1> [-dof <dof1>] -pset <pointset2> [-dof <dof2>]... [options]\n";
  cout << "       " << name << " <dataset1> <dataset2>... [options]\n";
  cout << "       " << name << " -pset <pointset1> <pointset2>... [options]\n";
  cout << "       " << name << " -image <image1> -points <pointset1> [-dof <dof1>]... [options]\n";
#endif
}

/// Print brief program usage information
void PrintUsage(const char* name)
{
  cout << "\n";
  PrintSynopsis(name);
  cout << "\n";
  cout << "Required options:\n";
  cout << "  -dofout <file>          Write transformation to specified file.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -model <name>           Transformation model(s). (default: Rigid+Affine+FFD)\n";
  cout << "  -image <file>...        Input image(s) to be registered.\n";
  cout << "  -output <file>          Write (first) transformed source image to specified file.\n";
#if MIRTK_Registration_WITH_PointSet
  cout << "  -pset <file>...         Input points, curve(s), surface(s), and/or other simplicial complex(es).\n";
#endif
  cout << "  -dof <file>             Affine transformation to be applied to the preceeding image and/or polydata.\n";
  cout << "  -dof_i <file>           Apply inverse of pre-computed affine transformation instead (cf. -dof).\n";
  cout << "  -mask <file>            Reference mask which defines the domain within which to evaluate the\n";
  cout << "                          energy function (i.e. data fidelity terms). (default: none)\n";
  cout << "  -reset-mask [yes|no]    Set value of all :option:`-mask` voxels to active (1). (default: no)\n";
  cout << "  -dofin <file>           Initial transformation estimate. (default: guess or identity)\n";
  cout << "  -par <name> <value>     Specify parameter value directly as command argument.\n";
  cout << "  -parin <file>           Read parameters from configuration file. If \"stdin\" or \"cin\",\n";
  cout << "                          the parameters are read from standard input instead. (default: none)\n";
  cout << "  -parout <file>          Write parameters to the named configuration file. (default: none)\n";
  cout << "  -v, -verbose [n]        Increase/Set verbosity of output messages. (default: " << verbose << ")\n";
  cout << "  -h -[-]help             Print complete help and exit.\n";
  cout << "\n";
}

/// Print full program usage information
void PrintHelp(const char* name)
{
  cout << "\n";
  PrintSynopsis(name);
  cout << "\n";
  cout << "Description:\n";
  cout << "  Registers a set of images, polygonal surface meshes, and/or point clouds (e.g. fiducial markers).\n";
  cout << "  The set of input images can be comprised of multiple channels (e.g., acquired with different imaging\n";
  cout << "  modalities) at different time points. For longitudinal data, the temporal origin in the NIfTI header\n";
  cout << "  identifies the time point that each input image belongs to. How all input images and polydata sets are\n";
  cout << "  registered with one another is determined by an energy function. This energy function is formulated in\n";
  cout << "  a simplified math expression using MATLAB-style indexing for the individual input files, i.e.,\n";
  cout << "\n";
  cout << "  Energy function =      SIM[Image dissimilarity](I(1), I(2:end) o T)...\n";
  cout << "                  +      PDM[Point set distance](T o P(1), P(2:end))...\n";
  cout << "                  + 1e-3 BE[Bending energy](T)...\n";
  cout << "                  +    0 LE[Linear energy](T)...\n";
  cout << "                  +    0 TP[Topology preservation](T)...\n";
  cout << "                  +    0 VP[Volume preservation](T)...\n";
  cout << "                  +    0 LogJac[LogJac penalty](T)...\n";
  cout << "                  +    0 NegJac[NegJac penalty](T)...\n";
  cout << "                  +    0 Sparsity(T)\n";
  cout << "\n";
  cout << "  where only energy terms with non-zero weights are active during the registration.\n";
  cout << "  The image dissimilarity term is only added if at least two input images are given.\n";
  cout << "  Similarly, only with at least two input polydata sets, the PDM term is added.\n";
  cout << "  These energy terms are referenced in the configuration file using their respective\n";
  cout << "  identifier in square brackets. For example, to change the weight of the bending energy\n";
  cout << "  smoothness term, add \"Bending energy weight = 0.01\" to the configuration file or\n";
  cout << "  use the -par option on the command-line. To enable volume preservation, set the\n";
  cout << "  parameter \"Volume preservation weight\" to a positive value.\n";
  cout << "\n";
  cout << "Output options:\n";
  cout << "  -dofout <file>    Write transformation to specified file.\n";
  cout << "  -output <file>    Write (first) transformed source image to specified file.\n";
  cout << "                    Given the flexibility of the energy function formulation,\n";
  cout << "                    this option may not always give the desired output and is\n";
  cout << "                    limited to standard pairwise image registration using a\n";
  cout << "                    single image dissimilarity term. Use transform-image command\n";
  cout << "                    with the :option:`-dofout` file as input otherwise.\n";
  cout << "\n";
  cout << "Input options:\n";
  cout << "  -image <file>...\n";
  cout << "      Input image to be registered. The order in which the images are\n";
  cout << "      given on the command-line corresponds to the image indices used\n";
  cout << "      in the registration energy formulation. The first input image has\n";
  cout << "      index 1 and is denoted by the symbol I(1) in the image similarity term.\n";
#if MIRTK_Registration_WITH_PointSet
  cout << "  -pset <file>...\n";
  cout << "      Input points, curve(s), surface(s), or other simplicial complex(es).\n";
  cout << "      The order in which the point set files are given on the command-line\n";
  cout << "      corresponds to the point set indices used in the registration energy\n";
  cout << "      formulation. The first input point set has index 1 and is denoted by\n";
  cout << "      the symbol P(1), C(1), or S(1), respectively, in the point set distance\n";
  cout << "      measure (PDM) term.\n";
#endif
  cout << "  -dof <file>\n";
  cout << "      Specifies a pre-computed affine transformation to be applied to\n";
  cout << "      the preceeding image and/or polydata. Note the difference to :option:`-dofin`,\n";
  cout << "      which specifies an initial estimate of the :option:`-dofout` transformation.\n";
  cout << "      (default: none)\n";
  cout << "  -dof_i <file>\n";
  cout << "      Specifies a pre-computed affine transformation whose inverse is to be\n";
  cout << "      applied to the preceeding image and/or polydata (see :option:`-dof`). (default: none)\n";
  cout << "  -images <file>\n";
  cout << "      Text file with N lines containing the file names of each\n";
  cout << "      2D/3D input image, optionally followed by a space character and the\n";
  cout << "      associated (relative) acquisition time (in ms). If no acquisition\n";
  cout << "      time is given, the temporal offset of the image header is used.\n";
  cout << "      Images specified as positional arguments are read before the images\n";
  cout << "      named in this (then optional) text file and the index of the images\n";
  cout << "      given in the text file are offset by the number of images named\n";
  cout << "      as positional arguments or using :option:`-image`, respectively.\n";
  cout << "  -psets <file>\n";
  cout << "      Equivalent to :option:`-images`, but for point sets (:option:`-pset`).\n";
  cout << "  -dofs <file>\n";
  cout << "      Read names of affine transformations of input images from\n";
  cout << "      specified text file. The affine transformations are set as\n";
  cout << "      if the images had been transformed by these before. Hence,\n";
  cout << "      these transformations are not part of the output transformation.\n";
  cout << "      The text file must contain on each line the basename of the input image\n";
  cout << "      followed by the corresponding affine transformation file. Input images\n";
  cout << "      not included are assumed to not be further pre-transformed. (default: none)\n";
  cout << "  -dofins <file>\n";
  cout << "      Read pairwise transformations from the files specified\n";
  cout << "      in the given text file which lists for each pair-wise\n";
  cout << "      transformation the corresponding target image, followed\n";
  cout << "      by the source image and the name of the transformation file.\n";
  cout << "      Note that the target and source image names must be listed\n";
  cout << "      in the :option:`-images` list file. (default: none)\n";
  cout << "  -dofin <file>\n";
  cout << "      Read initial transformation from file if :option:`-dofins` not specified.\n";
  cout << "      When no initial guess is given, and the first transformation model is a\n";
  cout << "      homogeneous transformation, a translation which aligns image foreground centers\n";
  cout << "      of mass is used. The identity mapping ('Id' or 'identity') is used as initial guess\n";
  cout << "      for deformable models unless the <file> argument is 'guess' to align the centers.\n";
  cout << "      If the given transformation cannot be used directly as starting point of the\n";
  cout << "      registration, it will be approximated by an instance of the chosen transformation\n";
  cout << "      model at the initial resolution level. The input transformation may thus be of\n";
  cout << "      different type than the output transformation of the registration.\n";
  cout << "      (default: guess or Id/identity)\n";
  cout << "  -mask <file>\n";
  cout << "      Reference mask which defines the domain within which to evaluate the\n";
  cout << "      energy function (i.e. image similarity). The registered images will\n";
  cout << "      thus be resampled within the corresponding domain of the world system.\n";
  cout << "      By default, the foreground of the target image defines this domain.\n";
  cout << "\n";
  cout << "Configuration options:\n";
  cout << "  -parin <file>\n";
  cout << "      Read parameters from configuration file. If <file> is \"stdin\" or \"cin\",\n";
  cout << "      the parameters are read from standard input instead. (default: none)\n";
  cout << "  -par <name> <value>\n";
  cout << "      Specify any parameter value usually found in :option:`-parin` file as command argument.\n";
  cout << "      The other configuration options are convenient shortcuts for commonly customized parameters.\n";
  cout << "  -model <m1>[+<m2>...]\n";
  cout << "      \"Transformation model\". Multiple models can be concatenated using a plus sign.\n";
  cout << "      (default: Rigid+Affine+FFD)\n";
  cout << "  -multi-level-model, -composition None|Sum|Fluid|LogSum\n";
  cout << "      \"Multi-level transformation model\" used to combine global affine transformation\n";
  cout << "      with the local free-form deformations at each resolution level. (default: Sum)\n";
  cout << "  -sim <value>\n";
  cout << "      Specifies concrete \"Image (dis-)similarity\" measure of SIM \"Energy function\" term.\n";
  cout << "      Most often used measures are MSE/SSD, NMI, and NCC/LNCC. (default: NMI)\n";
  cout << "  -bins <n>\n";
  cout << "      \"No. of bins\" used for NMI image similarity measure.\n";
  cout << "  -window <width> [<units> [<type>]]\n";
  cout << "      Local window used for (local) NCC/LNCC image similarity measure.\n";
  cout << "      A window <width> of zero corresponds to a global NCC measure.\n";
  cout << "      The <units> can be either \"vox\" (default) or \"mm\". The type\n";
  cout << "      of the window can be \"box\" (default), \"sigma\", \"fwhm\", or \"fwtm\",\n";
  cout << "      where the latter correspond to a Gaussian window with either the specified\n";
  cout << "      standard deviation, full width at half maximum, or full width at tenth maximum,\n";
  cout << "      respectively. The default is a box window with uniform weights for each voxel.\n";
  cout << "  -interp, -interpolation <mode>\n";
  cout << "      \"Image interpolation\" mode. (default: \"Fast linear [with padding]\")\n";
  cout << "  -extrap, -extrapolation <mode>\n";
  cout << "      \"Image extrapolation\" mode. (default: \"Default\" for used interpolation mode)\n";
  cout << "  -levels <from> [<to>]\n";
  cout << "      Image/FFD resolution levels. The <from> number corresponds to the \"No. of levels\",\n";
  cout << "      and the <to> number is the final level which is 1 by default.\n";
  cout << "      When images are given as input, the default number of levels is 4 and 1 otherwise.\n";
  cout << "  -level <n>\n";
  cout << "      Alias for :option:`-levels` <n> <n> which only performs the registration on a single level.\n";
  cout << "  -bg, -background, -padding <value>\n";
  cout << "      \"Background value\" (threshold) of input and output images (default: none)\n";
  cout << "  -ds <width>\n";
  cout << "      \"Control point spacing\" of free-form deformation on highest resolution level. (default: 4x min voxel size)\n";
  cout << "  -be <w>\n";
  cout << "      \"Bending energy weight\" of free-form deformation. (default: 0.001)\n";
  cout << "  -le <w> [<lambda>]\n";
  cout << "      \"Linear energy weight\" of free-form deformation. (default: 0)\n";
  cout << "  -tp <w>\n";
  cout << "      \"Topology preservation weight\" of free-form deformation. (default: 0)\n";
  cout << "  -vp <w>\n";
  cout << "      \"Volume preservation weight\" of free-form deformation. (default: 0)\n";
  cout << "  -lj, -jl, -log-jac, -jac <w>\n";
  cout << "      \"LogJac penalty weight\" of free-form deformation. For a classic FFD transformation model\n";
  cout << "      this penalty term is equivalent to the volume preservation term. When applied to the SVFFD\n";
  cout << "      model, however, this penalty applies to the Jacobian determinant of the velocity field. (default: 0)\n";
  cout << "  -nj, -neg-jac <w>\n";
  cout << "      \"NegJac penalty weight\" of free-form deformation. For a classic FFD transformation model\n";
  cout << "      this penalty term is equivalent to the volume preservation term. When applied to the SVFFD\n";
  cout << "      model, however, this penalty applies to the Jacobian determinant of the velocity field. (default: 0)\n";
  cout << "  -parout <file>\n";
  cout << "      Write parameters to the named configuration file. Note that after initialization of\n";
  cout << "      the registration, an interim configuration file is written. This file is overwritten\n";
  cout << "      once the registration finished with final configuration used during the course of the\n";
  cout << "      registration. (default: none)\n";
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliary functions
// =============================================================================

/// Read image/polydata file names and times from input list file
int read_input_list_file(const char *input_list_name, Array<string>& names, Array<double> &times)
{
  // Base directory of relative file paths
  const string base_dir = Directory(input_list_name);

  string   prefix;               // Common file name prefix
  string   suffix;               // Common file name suffix/extension
  ifstream iff(input_list_name); // List input file stream
  string   line;                 // Input line
  string   part;                 // Line part
  string   entry;                // Column entry
  bool     inquote = false;      // Parsing parts of a quoted string
  int      l = 0;                // Line number

  names.clear();
  times.clear();

  // Read base directory and default extensionf for image files
  while (getline(iff, line)) {
    l++;
    if (line.empty()) continue;
    // Discard leading white space characters
    line = line.substr(line.find_first_not_of(" \t"));
    // Ignore comment lines starting with # character
    if (line[0] == '#') continue;
    // Split first line at spaces
    int c = 0;
    istringstream lss(line);
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
      if (c == 0) {
        if (entry.length() > 1 && entry[0] == '.' && entry[1] == PATHSEP) {
          entry = entry.substr(2);
        }
        if (entry[0] != PATHSEP) {
          if (base_dir.empty()) prefix = "";
          else prefix = base_dir + PATHSEP;
          if (entry != ".") prefix += entry;
        } else {
          prefix = entry;
        }
      } else if (c == 1) {
        suffix = entry;
      } else {
        cerr << "Error: " << input_list_name << ":" << l << ": Too many entries!"<< endl;
        return 0;
      }
      entry.clear();
      c++;
    }
    if (c != 0) break; // Skip empty lines
  }
  if (l == 0) {
    cerr << "Error: Cannot parse list file " << input_list_name << endl;
    return 0;
  }

  // Process following lines
  while (getline(iff, line)) {
    l++;
    // Discard leading white space characters
    if (line.empty()) continue;
    line = line.substr(line.find_first_not_of(" \t"));
    // Ignore comment lines starting with # character
    if (line[0] == '#') continue;
    // Split line at spaces
    int c = 0;
    istringstream lss(line);
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
        if (!prefix.empty() && entry[0] != PATHSEP) entry = prefix + entry;
        names.push_back(entry + suffix);
      } else if (c == 1) {
        double time;
        if (!FromString(entry.c_str(), time)) {
          cerr << "Error: " << input_list_name << ":" << l << ": Failed to parse time value!" << endl;
          names.clear();
          times.clear();
          return 0;
        }
        times.push_back(time);
      } else {
        cerr << "Error: " << input_list_name << ":" << l << ": Too many entries!"<< endl;
        names.clear();
        times.clear();
        return 0;
      }
      // Parse next column
      entry.clear();
      c++;
    }
    if (0 < c && c < 2) {
      times.push_back(numeric_limits<double>::quiet_NaN()); // no time specified
    }
  }

  return static_cast<int>(names.size());
}

/// Read transformation file names and times from input list file
int read_dofin_list_file(const char *list_name, const Array<string> &images, Array<int> &targets, Array<int> &sources, Array<string> &dofs)
{
  string   base_dir = Directory(list_name); // Base directory containing transformation files
  ifstream iff(list_name);                  // List input file stream
  string   line;                            // Input line
  string   part;                            // Line part
  string   entry;                           // Column entry
  bool     inquote = false;                 // Parsing quoted string or not
  int      l = 0;                           // Line number

  targets.clear();
  sources.clear();
  dofs.clear();

  // Read base directory for transformation files
  if (!getline(iff, line)) {
    cerr << "Error: Cannot parse list file " << list_name << endl;
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
    // Split line at comma
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
      if (c == 0 || c == 1) {
        int n = 0;
        while (static_cast<unsigned int>(n) < images.size() && entry != BaseName(images[n])) n++;
        if (static_cast<unsigned int>(n) == images.size()) {
          cerr << "Error: " << list_name << ":" << l << ": Unknown "<< ((c == 0) ? "target" : "source") << " image!" << endl;
          dofs   .clear();
          targets.clear();
          sources.clear();
          return 0;
        }
        if (c == 0) targets.push_back(n);
        else        sources.push_back(n);
      } else if (c == 2) {
        if (!base_dir.empty() && entry[0] != PATHSEP) entry = base_dir + PATHSEP + entry;
        dofs.push_back(entry);
      } else {
        cerr << "Error: " << list_name << ":" << l << ": Too many comma separated entries!"<< endl;
        dofs   .clear();
        targets.clear();
        sources.clear();
        return 0;
      }
      // Parse next column
      entry.clear();
      c++;
    }
    if (c == 0) continue; // Skip empty lines
    if (c < 2) {
      cerr << "Error: Missing source image and transformation name at line " << l << endl;
    } else if (c < 3) {
      cerr << "Error: Missing transformation name at line " << l << endl;
    }
    if (c != 3) {
      dofs   .clear();
      targets.clear();
      sources.clear();
      return 0;
    }
  }

  return static_cast<int>(dofs.size());
}

// -----------------------------------------------------------------------------
bool IsIdentity(const string &name)
{
  return name.empty() || name == "identity" || name == "Identity" || name == "Id";
}

// =============================================================================
// Read input
// =============================================================================

// -----------------------------------------------------------------------------
struct ConcurrentImageReader
{
  enum Error { None, InvalidDoF };

  void operator()(const blocked_range<size_t> &re) const
  {
    HomogeneousTransformation *lin;
    for (size_t n = re.begin(); n != re.end(); ++n) {
      UniquePtr<BaseImage> &image = (*_Image)[n];
      if (image.get() == NULL) {
        image.reset(BaseImage::New(_ImageName[n].c_str()));
        if (!IsIdentity(_DoFName[n])) {
          UniquePtr<Transformation> dof(Transformation::New(_DoFName[n].c_str()));
          lin = dynamic_cast<HomogeneousTransformation *>(dof.get());
          if (lin) {
            Matrix mat = lin->GetMatrix();
            if (_DoFInvert[n]) mat.Invert();
            image->PutAffineMatrix(mat, true);
          } else {
            _Error[n] = InvalidDoF;
          }
        }
      }
      _Error[n] = None;
    }
  }

  static void Run(const Array<string>           &fname,
                  const Array<string>           &tname,
                  const Array<bool>             &tinv,
                  Array<UniquePtr<BaseImage> > &image,
                  Error                         *error)
  {
    ConcurrentImageReader body;
    body._ImageName = fname;
    body._DoFName   = tname;
    body._DoFInvert = tinv;
    body._Image     = &image;
    body._Error     = error;
    blocked_range<size_t> idx(0u, fname.size());
//    parallel_for(idx, body);
    body(idx);
  }

private:

  Array<string>                  _ImageName;
  Array<string>                  _DoFName;
  Array<bool>                    _DoFInvert;
  Array<UniquePtr<BaseImage> > *_Image;
  Error                         *_Error;
};

// -----------------------------------------------------------------------------
#if MIRTK_Registration_WITH_PointSet
struct ConcurrentPointSetReader
{
  enum Error { None, CannotOpenFile, InvalidType, EmptyPointSet };

private:

  Array<string>                         _FileName;
  Array<string>                         _DoFName;
  Array<bool>                           _DoFInvert;
  Array<vtkSmartPointer<vtkPointSet> > *_PointSet;
  Error                                *_Error;

public:

  void operator()(const blocked_range<size_t> &re) const
  {
    double p[3];
    Array<vtkSmartPointer<vtkPointSet> > &pointset = *_PointSet;
    for (size_t n = re.begin(); n != re.end(); ++n) {
      pointset[n] = ReadPointSet(_FileName[n].c_str());
      if (pointset[n]) {
        vtkPoints *points = pointset[n]->GetPoints();
        if (points->GetNumberOfPoints() == 0) {
          _Error[n] = EmptyPointSet;
        } else {
          if (!IsIdentity(_DoFName[n])) {
            UniquePtr<Transformation> t(Transformation::New(_DoFName[n].c_str()));
            for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
              points->GetPoint(i, p);
              if (_DoFInvert[n]) t->Inverse  (p[0], p[1], p[2]);
              else               t->Transform(p[0], p[1], p[2]);
              points->SetPoint(i, p);
            }
          }
          _Error[n] = None;
        }
      } else {
        _Error[n] = CannotOpenFile;
      }
    }
  }

  static void Run(const Array<string>                  &fname,
                  const Array<string>                  &tname,
                  const Array<bool>                    &tinv,
                  Array<vtkSmartPointer<vtkPointSet> > &pointset,
                  Error                                *error)
  {
    ConcurrentPointSetReader body;
    body._FileName  = fname;
    body._DoFName   = tname;
    body._DoFInvert = tinv;
    body._PointSet  = &pointset;
    body._Error     = error;
    blocked_range<size_t> idx(0, fname.size());
//    parallel_for(idx, body);
    body(idx);
  }

};
#endif // MIRTK_Registration_WITH_PointSet

// =============================================================================
// Main function
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  // Initialize libraries / object factories
  InitializeNumericsLibrary();
  InitializeIOLibrary();
  InitializeTransformationLibrary();
  InitializeRegistrationLibrary();
  #if MIRTK_Registration_WITH_Deformable
    InitializeDeformableLibrary();
  #endif

  // Default verbosity, output can be suppressed using option "-v 0"
  verbose = 1;

  // Print brief help if called without any arguments
  if (argc == 1) {
    PrintUsage(EXECNAME);
    exit(1);
  }

  // Positional arguments
  REQUIRES_POSARGS(0);

  Array<string> image_names;
  Array<string> imdof_names;
  Array<bool>   imdof_invert;

  Array<string> pset_names;
  Array<string> pdof_names;
  Array<bool>   pdof_invert;

  for (ALL_POSARGS) {
    size_t len = strlen(POSARG(ARGIDX));
    if (len > 4 && (strcmp(POSARG(ARGIDX) + len - 4u, ".vtk") == 0 ||
                    strcmp(POSARG(ARGIDX) + len - 4u, ".VTK") == 0)) {
      pset_names.push_back(POSARG(ARGIDX));
    } else {
      image_names.push_back(POSARG(ARGIDX));
    }
  }
  imdof_names .resize(image_names.size(), "identity");
  pdof_names  .resize(pset_names .size(), "identity");
  imdof_invert.resize(image_names.size(), false);
  pdof_invert .resize(pset_names .size(), false);

  // Optional arguments
  bool debug_output_level_prefix = true;
  const char *image_list_name    = nullptr;
  const char *dofin_list_name    = nullptr;
  const char *pset_list_name     = nullptr;
  const char *dofin_name         = nullptr;
  const char *dofout_name        = nullptr;
  const char *imgout_name        = nullptr;
  const char *tgtdof_name        = nullptr;
  const char *parin_name         = nullptr;
  const char *parout_name        = nullptr;
  const char *mask_name          = nullptr;
  bool        reset_mask         = false;
  ParameterList params;

  enum {
    UnknownTransformation,
    ImageTransformation,
    PointSetTransformation,
    ImageAndPointSetTransformation
  } doftype = UnknownTransformation;

  for (ALL_OPTIONS) {
    if (OPTION("-image")) {
      do {
        image_names .push_back(ARGUMENT);
        imdof_names .push_back("identity");
        imdof_invert.push_back(false);
      } while (HAS_ARGUMENT);
      if (doftype == PointSetTransformation) doftype = ImageAndPointSetTransformation;
      else                                   doftype = ImageTransformation;
    }
    else if (OPTION("-points") || OPTION("-pointset") || OPTION("-pset") ||
             OPTION("-mesh")   || OPTION("-polydata") || OPTION("-poly")) {
#if MIRTK_Registration_WITH_PointSet
      const char *arg = ARGUMENT;
      const string ext = Extension(arg);
      if (strcmp(OPTNAME, "-points") == 0 && (ext == ".txt" || ext == ".csv")) {
        pset_list_name = arg;
        doftype = UnknownTransformation;
      } else {
        for(;;) {
          pset_names.push_back(arg);
          pdof_names .push_back("identity");
          pdof_invert.push_back(false);
          if (!HAS_ARGUMENT) break;
          arg = ARGUMENT;
        }
        if (doftype == ImageTransformation) doftype = ImageAndPointSetTransformation;
        else                                doftype = PointSetTransformation;
      }
#else // MIRTK_Registration_WITH_PointSet
      cerr << "Point set input only supported when built with MIRTK PointSet module!" << endl;
      exit(1);
#endif // MIRTK_Registration_WITH_PointSet
    }
    else if (OPTION("-dof") || OPTION("-dof_i")) {
      const bool  invert   = (strcmp(OPTNAME, "-dof_i") == 0);
      const char *name = ARGUMENT;
      switch (doftype) {
        case ImageTransformation:
          imdof_names .back() = name;
          imdof_invert.back() = invert;
          break;
        case PointSetTransformation:
          pdof_names .back() = name;
          pdof_invert.back() = invert;
          break;
        case ImageAndPointSetTransformation:
          imdof_names .back() = pdof_names .back() = name;
          imdof_invert.back() = pdof_invert.back() = invert;
          break;
        default:
          cerr << EXECNAME << ": -dof option must be preceeded by -image and/or -poly/-polydata" << endl;
          exit(1);
      }
    }
    else if (OPTION("-images")) {
      image_list_name = ARGUMENT;
    }
    else if (OPTION("-pointsets") || OPTION("-psets") || OPTION("-meshes") || OPTION("-polys")) {
      pset_list_name = ARGUMENT;
    }
    else if (OPTION("-dofins")) dofin_list_name = ARGUMENT;
    else if (OPTION("-dofin" )) dofin_name      = ARGUMENT;
    else if (OPTION("-dofout")) dofout_name     = ARGUMENT;
    else if (OPTION("-output")) imgout_name     = ARGUMENT;
    else if (OPTION("-disp"))   tgtdof_name     = ARGUMENT;
    else if (OPTION("-mask"))   mask_name       = ARGUMENT;
    else HANDLE_BOOLEAN_OPTION("reset-mask", reset_mask);
    else if (OPTION("-nodebug-level-prefix")) {
      debug_output_level_prefix = false;
    }
    // Parameter
    else if (OPTION("-par")) {
      const char *param = ARGUMENT;
      const char *value = ARGUMENT;
      Insert(params, param, value);
    }
    else if (OPTION("-parin" )) {
      parin_name = ARGUMENT;
    }
    else if (OPTION("-parout")) {
      parout_name = ARGUMENT;
    }
    // Shortcuts for often used -par "<parameter> = <value>"
    else if (OPTION("-model")) {
      Insert(params, "Transformation model", ARGUMENT);
    }
    else if (OPTION("-multi-level-model") || OPTION("-composition")) {
      Insert(params, "Multi-level transformation model", ARGUMENT);
    }
    else if (OPTION("-sim")) {
      Insert(params, "Image (dis-)similarity measure", ARGUMENT);
    }
    else if (OPTION("-bins")) {
      int n;
      PARSE_ARGUMENT(n);
      Insert(params, "No. of bins", n);
    }
    else if (OPTION("-window-size") || OPTION("-window")) {
      float width;
      PARSE_ARGUMENT(width);
      string units = "vox";
      if (HAS_ARGUMENT) units = ARGUMENT;
      string type  = "box";
      if (HAS_ARGUMENT) type = ARGUMENT;
      Insert(params, string("Local window size [") + type + string("]"), ToString(width) + units);
    }
    else if (OPTION("-pdm")) {
      Insert(params, "Point set distance measure", ARGUMENT);
    }
    else if (OPTION("-level")) {
      int n;
      PARSE_ARGUMENT(n);
      Insert(params, "First level", n);
      Insert(params, "Final level", n);
    }
    else if (OPTION("-levels")) {
      int n, m = 1;
      PARSE_ARGUMENT(n);
      if (HAS_ARGUMENT) PARSE_ARGUMENT(m);
      Insert(params, "First level", n);
      Insert(params, "Final level", m);
    }
    else if (OPTION("-interp") || OPTION("-interpolation")) {
      Insert(params, "Image interpolation mode", ARGUMENT);
    }
    else if (OPTION("-extrap") || OPTION("-extrapolation")) {
      Insert(params, "Image extrapolation mode", ARGUMENT);
    }
    else if (OPTION("-ds")) {
      double ds;
      PARSE_ARGUMENT(ds);
      Insert(params, "Control point spacing", ds);
    }
    else if (OPTION("-be")) {
      double w;
      PARSE_ARGUMENT(w);
      Insert(params, "Bending energy weight", w);
    }
    else if (OPTION("-le")) {
      double w;
      PARSE_ARGUMENT(w);
      Insert(params, "Linear energy weight", w);
      if (HAS_ARGUMENT) {
        double lambda;
        PARSE_ARGUMENT(lambda);
        Insert(params, "Linear energy lambda", lambda / w);
      }
    }
    else if (OPTION("-vp")) {
      double w;
      PARSE_ARGUMENT(w);
      Insert(params, "Volume preservation weight", w);
    }
    else if (OPTION("-tp")) {
      double w;
      PARSE_ARGUMENT(w);
      Insert(params, "Topology preservation weight", w);
    }
    else if (OPTION("-lj") || OPTION("-jl") || OPTION("-logjac") || OPTION("-log-jac") || OPTION("-jac")) {
      double w;
      PARSE_ARGUMENT(w);
      Insert(params, "LogJac penalty weight", w);
    }
    else if (OPTION("-nj") || OPTION("-negjac") || OPTION("-neg-jac")) {
      double w;
      PARSE_ARGUMENT(w);
      Insert(params, "NegJac penalty weight", w);
    }
    else if (OPTION("-padding") || OPTION("-bg") || OPTION("-background")) {
      double v;
      PARSE_ARGUMENT(v);
      Insert(params, "Background value", v);
    }
    // Unknown option
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (!dofout_name && !imgout_name) {
    FatalError("Either -dofout (recommended) or -output option argument required! Use -dofout none to not write any output.");
  }

  // TODO: Read initial transformations from list file (cf. obsolete tdreg)
  if (dofin_list_name) {
    FatalError("Option -dofins currently not implemented");
  }

  // ---------------------------------------------------------------------------
  // Print version information
  if (verbose) {
    PrintVersion(cout, EXECNAME);
    if (verbose > 1) cout.put('\n');
  }

  MIRTK_START_TIMING();

  // ---------------------------------------------------------------------------
  // Read configuration
  GenericRegistrationFilter registration;

  if (verbose > 2) cout << "Reading configuration ..." << endl;
  // Set registration filter input and configuration
  bool parin_stdin = (parin_name && (strcmp(parin_name, "stdin") == 0 ||
                                     strcmp(parin_name, "STDIN") == 0 ||
                                     strcmp(parin_name, "cin")   == 0));
  // 1. Read parameters from file
  if (!parin_stdin && parin_name && !registration.Read(parin_name, verbose > 2)) {
    cerr << "Failed to read configuration from file \"" << parin_name << "\"!" << endl;
    exit(1);
  }
  // 2. Set parameters provided as command arguments
  if (!registration.Parameter(params)) {
    cerr << "Failed to parse configuration given as command arguments!" << endl;
    exit(1);
  }
  // 3. Add any parameters read from standard input stream
  if (parin_stdin) {
    if (verbose) {
      cout << "\nEnter additional parameters now (press Ctrl-D to continue):" << endl;
    }
    if (!registration.Read(cin, verbose > 2)) {
      cerr << "Failed to read configuration from standard input stream!" << endl;
      exit(1);
    }
  }
  if (verbose > 2 || (parin_stdin && verbose > 0)) {
    if (verbose > 2) cout << "\n";
    cout << "Reading configuration ... done";
    if (verbose > 2) cout << "\n";
    cout << endl;
  }

  // ---------------------------------------------------------------------------
  // Determine number of input images
  Array<double>                image_times;
  Array<UniquePtr<BaseImage> > images;

  if (image_names.size() == 1) {
    UniquePtr<BaseImage> sequence(BaseImage::New(image_names[0].c_str()));
    const int nframes = sequence->GetT();
    if (nframes < 2) {
      if (verbose > 1) cout << " failed\n" << endl;
      cerr << "Error: Input sequence contains only one temporal frame!" << endl;
      exit(1);
    }
    image_names.resize(nframes, image_names[0]);
    image_times.resize(nframes);
    images.resize(nframes);
    for (int n = 0; n < nframes; ++n) {
      BaseImage *frame = NULL;
      sequence->GetFrame(frame, n);
      images[n].reset(frame);
      image_times[n] = frame->GetTOrigin();
    }
  }

  if (image_list_name) {
    Array<string> names;
    Array<double> times;
    if (read_input_list_file(image_list_name, names, times) == 0) {
      if (verbose > 1) cout << endl;
      cerr << "Error: Failed to parse input list file " << image_list_name << "!" << endl;
      exit(1);
    }
    image_names.insert(image_names.end(), names.begin(), names.end());
    image_times.insert(image_times.end(), times.begin(), times.end());
  }

  // ---------------------------------------------------------------------------
  // Determine number of input point sets
  Array<double> pset_times;

  if (pset_list_name) {
    Array<string> names;
    Array<double> times;
    if (read_input_list_file(pset_list_name, names, times) == 0) {
      if (verbose > 1) cout << endl;
      cerr << "Error: Failed to parse input list file " << pset_list_name << "!" << endl;
      exit(1);
    }
    pset_names.insert(pset_names.end(), names.begin(), names.end());
    pset_times.insert(pset_times.end(), times.begin(), times.end());
  }

  // ---------------------------------------------------------------------------
  // Parse energy formula to determine which input files are really used
  registration.ParseEnergyFormula(static_cast<int>(image_names.size()),
                                  static_cast<int>(pset_names .size()), 1);

  const int nimages = registration.NumberOfRequiredImages();
  image_names .resize(nimages);
  image_times .resize(nimages);
  imdof_names .resize(nimages);
  imdof_invert.resize(nimages);
  images      .resize(nimages);

  const int npsets = registration.NumberOfRequiredPointSets();
  pset_names .resize(npsets);
  pset_times .resize(npsets);
  pdof_names .resize(npsets);
  pdof_invert.resize(npsets);

  if (verbose > 2) {
    if (nimages > 0) {
      cout << "Input images:\n";
      for (size_t n = 0; n < image_names.size(); ++n) {
        cout << "  t = " << setw(5) << right << image_times[n] << ": " << image_names[n] << "\n";
      }
    }
    if (npsets > 0) {
      if (nimages > 0) cout << '\n';
      cout << "Input point sets:\n";
      for (size_t n = 0; n < pset_names.size(); ++n) {
        cout << "  t = " << setw(5) << right << pset_times[n] << ": " << pset_names[n] << "\n";
      }
    }
    cout << endl;
  }

  // ---------------------------------------------------------------------------
  // Read input images
  {
    if (nimages > 0 && verbose > 1) {
      cout << "Reading images ..........";
      cout.flush();
    }

    ConcurrentImageReader::Error *errors = new ConcurrentImageReader::Error[nimages];
    ConcurrentImageReader::Run(image_names, imdof_names, imdof_invert, images, errors);

    int nerrors = 0;
    for (int n = 0; n < nimages; ++n) {
      if (errors[n] == ConcurrentImageReader::InvalidDoF) {
        ++nerrors;
        if (verbose > 1 && nerrors == 1) cout << " failed\n" << endl;
        cerr << "Error: Implicit transformation of image " << (n+1) << " can only be affine" << endl;
      }
    }
    delete[] errors;
    if (nerrors > 0) exit(1);

    if (nimages > 0 && verbose > 1) {
      cout << " done" << endl;
    }
  }

  // Set input images
  for (int n = 0; n < nimages; ++n) {
    images[n]->PutTOrigin(image_times[n]);
    #if 1  // TODO: Fix the actual issue so this is not needed!
      images[n]->PutAffineMatrix(images[n]->GetAffineMatrix(), true);
      if (!images[n]->GetAffineMatrix().IsIdentity()) {
        if (verbose > 0) cout << endl;
        Warning("Input image has shearing component in affine matrix (NIfTI sform)!"
                "\nThis may potentially result in a suboptimal output transformation!"
                "\nConsider pre-transforming the image with the given affine transformation."
                "\nThis issue has yet to be fixed properly within the Registration module.");
      }
    #endif
    registration.AddInput(images[n].get());
  }

  // ---------------------------------------------------------------------------
  // Read input point sets
#if MIRTK_Registration_WITH_PointSet
  Array<vtkSmartPointer<vtkPointSet> > psets(npsets);
  {
    if (npsets > 0 && verbose > 1) {
      cout << "Reading point sets ......";
      cout.flush();
    }

    ConcurrentPointSetReader::Error *errors = new ConcurrentPointSetReader::Error[npsets];
    ConcurrentPointSetReader::Run(pset_names, pdof_names, pdof_invert, psets, errors);

    int  nerrors = 0;
    for (int n = 0; n < npsets; ++n) {
      if (errors[n] != ConcurrentPointSetReader::None) {
        if (verbose > 1 && nerrors == 0) cout << " failed\n" << endl;
        ++nerrors;
      }
      switch (errors[n]) {
        case ConcurrentPointSetReader::CannotOpenFile: {
          cerr << "Error: Cannot open file " << pset_names[n] << "!" << endl;
        } break;
        case ConcurrentPointSetReader::InvalidType: {
          cerr << "Error: Input point set " << pset_names[n] << " must be vtkPointSet!" << endl;
        } break;
        case ConcurrentPointSetReader::EmptyPointSet: {
          cerr << "Error: Input point set " << pset_names[n] << " is empty!" << endl;
        } break;
        case ConcurrentPointSetReader::None: break;
      }
    }
    delete[] errors;
    if (nerrors > 0) exit(1);

    if (npsets > 0 && verbose > 1) {
      cout << " done" << endl;
    }
  }

  // Set input point sets
  for (int n = 0; n < npsets; ++n) {
    registration.AddInput(psets[n], pset_times[n]);
  }
#endif // MIRTK_Registration_WITH_PointSet

  // ---------------------------------------------------------------------------
  // Read target transformation
  UniquePtr<Transformation> tgtdof;
  if (tgtdof_name) {
    if (IsIdentity(tgtdof_name)) {
      tgtdof.reset(new RigidTransformation());
    } else {
      if (Transformation::CheckHeader(tgtdof_name)) {
        tgtdof.reset(Transformation::New(tgtdof_name));
      } else {
        GenericImage<double> disp(tgtdof_name);
        if (disp.T() != 3) {
          FatalError("Input -disp image must have 3 components/channels!");
        }
        tgtdof.reset(new LinearFreeFormTransformation3D(disp));
      }
    }
    registration.TargetTransformation(tgtdof.get());
  }

  // ---------------------------------------------------------------------------
  // Initialize registration
  GenericRegistrationLogger   logger;
  GenericRegistrationDebugger debugger("mirtk_");
  debugger.LevelPrefix(debug_output_level_prefix);

  logger.Verbosity(verbose - 1);
  if ((debug_time == 0 && verbose > 0) ||
      (debug_time  > 0 && verbose > 1)) {
    registration.AddObserver(logger);
  }
  if (debug) {
    registration.AddObserver(debugger);
  }

  Transformation *dofout = nullptr;
  registration.Output(&dofout);

  // Read mask which defines domain on which similarity is evaluated
  UniquePtr<BinaryImage> mask;
  if (mask_name) {
    mask.reset(new BinaryImage(mask_name));
    if (reset_mask) *mask = 1;
    registration.Domain(mask.get());
  }

  // Guess unset parameters, must be called before MakeInitialGuess and Write -parout
  registration.GuessParameter();

  // Write initial parameters to file
  //
  // Note: Will be overwritten once the registration finished successfully
  //       with the actual parameters used below.
  if (parout_name) {
    registration.Write(parout_name);
  }

  // Read initial transformation
  UniquePtr<Transformation> dofin;
  if (dofin_name) {
    if (IsIdentity(dofin_name)) {
      dofin.reset(new RigidTransformation());
    } else if (ToLower(dofin_name) == "guess") {
      dofin.reset(registration.MakeInitialGuess());
    } else {
      dofin.reset(Transformation::New(dofin_name));
    }
    registration.InitialGuess(dofin.get());
  } else {
    registration.InitialGuess(tgtdof.get());
  }

  MIRTK_DEBUG_TIMING(1, "reading input data");

  // ---------------------------------------------------------------------------
  // Run registration
  const clock_t start_cpu_time = clock();
#ifdef HAVE_TBB
  tbb::tick_count start_wall_time = tbb::tick_count::now();
#endif

  registration.Run();
  registration.DeleteObserver(logger);
  registration.DeleteObserver(debugger);

  if (verbose) {
    cout << "\n";
    double sec = static_cast<double>(clock() - start_cpu_time) / static_cast<double>(CLOCKS_PER_SEC);
    #ifdef HAVE_TBB
      cout << "CPU time is " << ElapsedTimeToString(sec, TIME_IN_SECONDS, TIME_FORMAT_H_MIN_SEC, 2) << "\n";
      sec = (tbb::tick_count::now() - start_wall_time).seconds();
    #endif
    cout << "Finished in " << ElapsedTimeToString(sec, TIME_IN_SECONDS, TIME_FORMAT_H_MIN_SEC, 2) << endl;
  }

  // Write final transformation
  if (dofout_name && strcmp(dofout_name, "none") != 0 &&
                     strcmp(dofout_name, "None") != 0 &&
                     strcmp(dofout_name, "NONE") != 0) {
    if (dofout->TypeOfClass() == TRANSFORMATION_SIMILARITY) {
      // Write affine transformation instead, because most other programs
      // cannot deal with the new similarity transformation type (yet)
      // TODO: Update other tools (e.g., rview) to handle SimilarityTransformation
      AffineTransformation aff(*static_cast<SimilarityTransformation *>(dofout));
      aff.Write(dofout_name);
    } else {
      dofout->Write(dofout_name);
    }
  }

  // Write actual parameters used to file
  if (parout_name) registration.Write(parout_name);

  // Write (first) transformed source image
  if (imgout_name) {
    int t = -1, s = -1;
    for (int n = 0; n < registration.NumberOfImages(); ++n) {
      if (t == -1 && registration.IsTargetImage(n)) t = n;
      if (s == -1 && registration.IsSourceImage(n)) s = n;
    }
    if (t == -1 || s == -1 || t == s) {
      delete dofout;
      FatalError("Sorry, could not determine which input image to resample on which target domain to write -output image. Use transform-image command instead.");
    }
    RegisteredImage output;
    RegisteredImage::InputImageType input(*images[s]);
    input.PutBackgroundValueAsDouble(registration.BackgroundValue(s));
    output.InputImage(&input);
    output.Transformation(dofout);
    output.InterpolationMode(registration.InterpolationMode(s));
    output.ExtrapolationMode(registration.ExtrapolationMode(s));
    output.PutBackgroundValueAsDouble(registration.BackgroundValue(s));
    output.Initialize(images[t]->Attributes());
    output.Update(true, false, false, true);
    output.Write(imgout_name);
  }

  // Clean up
  delete dofout;

  return 0;
}

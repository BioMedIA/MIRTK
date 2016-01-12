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

#include <mirtkNumericsConfig.h>
#include <mirtkImageIOConfig.h>
#include <mirtkTransformationConfig.h>
#include <mirtkRegistrationConfig.h>
#if MIRTK_Registration_WITH_Deformable
#  include <mirtkDeformableConfig.h>
#endif

#include <mirtkGenericImage.h>
#include <mirtkGenericRegistrationFilter.h>
#include <mirtkGenericRegistrationLogger.h>
#include <mirtkGenericRegistrationDebugger.h>

#include <mirtkTransformation.h>
#include <mirtkHomogeneousTransformation.h>
#include <mirtkRigidTransformation.h>
#include <mirtkSimilarityTransformation.h>
#include <mirtkAffineTransformation.h>

#if MIRTK_Registration_WITH_PointSet
#  include <mirtkPointSetUtils.h>
#  include <vtkSmartPointer.h>
#  include <vtkPolyData.h>
#  include <vtkOBJReader.h>
#  include <vtkPLYReader.h>
#  include <vtkSTLReader.h>
#  include <vtkDataSetReader.h>
#  include <vtkXMLPolyDataReader.h>
#  include <vtkTriangleFilter.h>
#  include <vtkDataSetSurfaceFilter.h>
#endif

using namespace mirtk;


// =============================================================================
// Version and help
// =============================================================================

/// Print command synopsis
void PrintSynopsis(const char *name)
{
  cout << "Usage: " << name << " -images <images.lst> [options]" << endl;
  cout << "       " << name << " -image <image1> [-dof <dof1>] -image <image2> [-dof <dof2>]... [options]" << endl;
  cout << "       " << name << " -image <image1> <image2>... [options]" << endl;
  cout << "       " << name << " <image1> <image2>... [options]" << endl;
  cout << "       " << name << " <image_sequence> [options]" << endl;
#if MIRTK_Registration_WITH_PointSet
  cout << "       " << name << " -pset <pointset1> [-dof <dof1>] -pset <pointset2> [-dof <dof2>]... [options]" << endl;
  cout << "       " << name << " <dataset1> <dataset2>... [options]" << endl;
  cout << "       " << name << " -pset <pointset1> <pointset2>... [options]" << endl;
  cout << "       " << name << " -image <image1> -points <pointset1> [-dof <dof1>]... [options]" << endl;
#endif
}

/// Print brief program usage information
void PrintUsage(const char* name)
{
  cout << endl;
  PrintSynopsis(name);
  cout << endl;
  cout << "Required options:" << endl;
  cout << "  -dofout <file>          Write transformation to specified file." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -model  <name>          Transformation model(s). (default: Rigid+Affine+FFD)" << endl;
  cout << "  -image  <file>...       Input image(s) to be registered." << endl;
#if MIRTK_Registration_WITH_PointSet
  cout << "  -pset <file>...         Input points, curve(s), surface(s), and/or other simplicial complex(es)." << endl;
#endif
  cout << "  -dof    <file>          Affine transformation to be applied to the preceeding image and/or polydata." << endl;
  cout << "  -dof_i  <file>          Apply inverse of pre-computed affine transformation instead (cf. -dof)." << endl;
  cout << "  -mask   <file>          Reference mask which defines the domain within which to evaluate the" << endl;
  cout << "                          energy function (i.e. data fidelity terms). (default: none)" << endl;
  cout << "  -dofin  <file>          Initial transformation estimate. (default: align centroids)" << endl;
  cout << "  -par <name> <value>    Specify parameter value directly as command argument." << endl;
  cout << "  -parin  <file>          Read parameters from configuration file. If \"stdin\" or \"cin\"," << endl;
  cout << "                          the parameters are read from standard input instead. (default: none)" << endl;
  cout << "  -parout <file>          Write parameters to the named configuration file. (default: none)" << endl;
  cout << "  -v -verbose [n]         Increase/Set verbosity of output messages. (default: " << verbose << ")" << endl;
  cout << "  -h -[-]help             Print complete help and exit." << endl;
  cout << endl;
}

/// Print full program usage information
void PrintHelp(const char* name)
{
  cout << endl;
  PrintSynopsis(name);
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Registers a set of images, polygonal surface meshes, and/or point clouds (e.g. fiducial markers)." << endl;
  cout << "  The set of input images can be comprised of multiple channels (e.g., acquired with different imaging" << endl;
  cout << "  modalities) at different time points. For longitudinal data, the temporal origin in the NIfTI header" << endl;
  cout << "  identifies the time point that each input image belongs to. How all input images and polydata sets are" << endl;
  cout << "  registered with one another is determined by an energy function. This energy function is formulated in" << endl;
  cout << "  a simplified math expression using MATLAB-style indexing for the individual input files, i.e.," << endl;
  cout << endl;
  cout << "  Energy function =      SIM[Image dissimilarity](I(1), I(2:end) o T)" << endl;
  cout << "                  +      PDM[Point set distance](T o P(1), P(2:end))" << endl;
  cout << "                  + 1e-3 BE [Bending energy](T)" << endl;
  cout << "                  +    0 VP [Volume preservation](T)" << endl;
  cout << "                  +    0 JAC[Jacobian penalty](T)" << endl;
  cout << "                  +    0 Sparsity(T)" << endl;
  cout << endl;
  cout << "  where only energy terms with non-zero weights are active during the registration." << endl;
  cout << "  The image dissimilarity term is only added if at least two input images are given." << endl;
  cout << "  Similarly, only with at least two input polydata sets, the PDM term is added." << endl;
  cout << "  These energy terms are referenced in the configuration file using their respective" << endl;
  cout << "  identifier in square brackets. For example, to change the weight of the bending energy" << endl;
  cout << "  smoothness term, add \"Bending energy weight = 0.01\" to the configuration file or" << endl;
  cout << "  use the -par option on the command-line. To enable volume preservation, set the" << endl;
  cout << "  parameter \"Volume preservation weight\" to a positive value." << endl;
  cout << endl;
  cout << "Required arguments:" << endl;
  cout << "  -dofout <file>               Write transformation to specified file." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -model <m1>[+<m2>...]        Transformation model(s). (default: Rigid+Affine+FFD)" << endl;
  cout << "                               Alternatively, use \"-par 'Transformation model=<name>'\" (see :option:`-par`)." << endl;
  cout << "  -image <file>...             Input image to be registered. The order in which the images are" << endl;
  cout << "                               given on the command-line corresponds to the image indices used" << endl;
  cout << "                               in the registration energy formulation. The first input image has" << endl;
  cout << "                               index 1 and is denoted by the symbol I(1) in the image similarity term." << endl;
#if MIRTK_Registration_WITH_PointSet
  cout << "  -pset <file>...              Input points, curve(s), surface(s), or other simplicial complex(es)." << endl;
  cout << "                               The order in which the point set files are given on the command-line" << endl;
  cout << "                               corresponds to the point set indices used in the registration energy" << endl;
  cout << "                               formulation. The first input point set has index 1 and is denoted by" << endl;
  cout << "                               the symbol P(1), C(1), or S(1), respectively, in the point set distance" << endl;
  cout << "                               measure (PDM) term." << endl;
#endif
  cout << "  -dof <file>                  Specifies a pre-computed affine transformation to be applied to" << endl;
  cout << "                               the preceeding image and/or polydata. Note the difference to :option:`-dofin`," << endl;
  cout << "                               which specifies an initial estimate of the :option:`-dofout` transformation." << endl;
  cout << "                               (default: none)" << endl;
  cout << "  -dof_i <file>                Specifies a pre-computed affine transformation whose inverse is to be" << endl;
  cout << "                               applied to the preceeding image and/or polydata (see :option:`-dof`). (default: none)" << endl;
  cout << "  -images <file>               Text file with N lines containing the file names of each" << endl;
  cout << "                               2D/3D input image, optionally followed by a space character and the" << endl;
  cout << "                               associated (relative) acquisition time (in ms). If no acquisition" << endl;
  cout << "                               time is given, the temporal offset of the image header is used." << endl;
  cout << "                               Images specified as positional arguments are read before the images" << endl;
  cout << "                               named in this (then optional) text file and the index of the images" << endl;
  cout << "                               given in the text file are offset by the number of images named" << endl;
  cout << "                               as positional arguments or using :option:`-image`, respectively." << endl;
  cout << "  -psets <file>                Equivalent to :option:`-images`, but for point sets (:option:`-pset`)." << endl;
  cout << "  -dofs <file>                 Read names of affine transformations of input images from" << endl;
  cout << "                               specified text file. The affine transformations are set as" << endl;
  cout << "                               if the images had been transformed by these before. Hence," << endl;
  cout << "                               these transformations are not part of the output transformation." << endl;
  cout << "                               The text file must contain on each line the basename of the input image" << endl;
  cout << "                               followed by the corresponding affine transformation file. Input images" << endl;
  cout << "                               not included are assumed to not be further pre-transformed. (default: none)" << endl;
  cout << "  -mask <file>                 Reference mask which defines the domain within which to evaluate the" << endl;
  cout << "                               energy function (i.e. image similarity). The registered images will" << endl;
  cout << "                               thus be resampled within the corresponding domain of the world system." << endl;
  cout << "                               By default, the foreground of the target image defines this domain." << endl;
  cout << "  -dofins <file>               Read pairwise transformations from the files specified" << endl;
  cout << "                               in the given text file which lists for each pair-wise" << endl;
  cout << "                               transformation the corresponding target image, followed" << endl;
  cout << "                               by the source image and the name of the transformation file." << endl;
  cout << "                               Note that the target and source image names must be listed" << endl;
  cout << "                               in the :option:`-images` list file. (default: none)" << endl;
  cout << "  -dofin  <file>               Read initial transformation from file if :option:`-dofins` not specified." << endl;
  cout << "                               Otherwise, writes the initial transformation obtained by approximating" << endl;
  cout << "                               the pairwise transformations to the named file." << endl;
  cout << "                               If the given transformation cannot be used directly as starting" << endl;
  cout << "                               point of the registration, it will be approximated by an instance" << endl;
  cout << "                               of the chosen transformation model at the initial resolution level." << endl;
  cout << "                               The input transformation may thus be of different type than the" << endl;
  cout << "                               output transformation of the registration. (default: none)" << endl;
  cout << "  -par '<name>=<value>'        Specify parameter value directly as command argument." << endl;
  cout << "  -parin  <file>               Read parameters from configuration file. If \"stdin\" or \"cin\"," << endl;
  cout << "                               the parameters are read from standard input instead. (default: none)" << endl;
  cout << "  -parout <file>               Write parameters to the named configuration file. Note that after" << endl;
  cout << "                               initializaton of the registration, an interim configuration file is" << endl;
  cout << "                               written. This file is overwritten once the registration finished with" << endl;
  cout << "                               final configuration used during the course of the registration." << endl;
  cout << "                               (default: none)" << endl;
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
      unique_ptr<BaseImage> &image = (*_Image)[n];
      if (image.get() == NULL) {
        image.reset(BaseImage::New(_ImageName[n].c_str()));
        if (!IsIdentity(_DoFName[n])) {
          unique_ptr<Transformation> dof(Transformation::New(_DoFName[n].c_str()));
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
                  Array<unique_ptr<BaseImage> > &image,
                  Error                         *error)
  {
    ConcurrentImageReader body;
    body._ImageName = fname;
    body._DoFName   = tname;
    body._DoFInvert = tinv;
    body._Image     = &image;
    body._Error     = error;
    blocked_range<size_t> idx(0, fname.size());
//    parallel_for(idx, body);
    body(idx);
  }

private:

  Array<string>                  _ImageName;
  Array<string>                  _DoFName;
  Array<bool>                    _DoFInvert;
  Array<unique_ptr<BaseImage> > *_Image;
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
            unique_ptr<Transformation> t(Transformation::New(_DoFName[n].c_str()));
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
  InitializeImageIOLibrary();
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
    int len = strlen(POSARG(ARGIDX));
    if (len > 4 && (strcmp(POSARG(ARGIDX) + len - 4, ".vtk") == 0 ||
                    strcmp(POSARG(ARGIDX) + len - 4, ".VTK") == 0)) {
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
  const char *image_list_name    = NULL;
  const char *dofin_list_name    = NULL;
  const char *pset_list_name     = NULL;
  const char *dofin_name         = NULL;
  const char *dofout_name        = NULL;
  const char *parin_name         = NULL;
  const char *parout_name        = NULL;
  const char *mask_name          = NULL;
  stringstream params;

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
    else if (OPTION("-images")) image_list_name = ARGUMENT;
    else if (OPTION("-pointsets") || OPTION("-psets") || OPTION("-meshes") || OPTION("-polys")) {
      pset_list_name = ARGUMENT;
    }
    else if (OPTION("-dofins")) dofin_list_name = ARGUMENT;
    else if (OPTION("-dofin" )) dofin_name      = ARGUMENT;
    else if (OPTION("-dofout")) dofout_name     = ARGUMENT;
    else if (OPTION("-mask"))   mask_name       = ARGUMENT;
    else if (OPTION("-nodebug-level-prefix")) debug_output_level_prefix = false;
    // Parameter
    else if (OPTION("-par"))    params << ARGUMENT << " = " << ARGUMENT << endl;
    else if (OPTION("-parin" )) parin_name      = ARGUMENT;
    else if (OPTION("-parout")) parout_name     = ARGUMENT;
    // Shortcuts for often used -par "<parameter> = <value>"
    else if (OPTION("-model"))  params << "Transformation model = " << ARGUMENT << endl;
    else if (OPTION("-sim"))    params << "Image (dis-)similarity measure = " << ARGUMENT << endl;
    else if (OPTION("-ds"))     params << "Control point spacing [mm] = " << ARGUMENT << endl;
    else if (OPTION("-be"))     params << "Bending energy weight = " << ARGUMENT << endl;
    else if (OPTION("-vp"))     params << "Volume preservation weight = " << ARGUMENT << endl;
    else if (OPTION("-jl") ||
             OPTION("-jac"))    params << "Jacobian penalty weight = " << ARGUMENT << endl;
    // Unknown option
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (!dofout_name) {
    cerr << "Error: Missing -dofout argument" << endl;
    exit(1);
  }

  // TODO: Read initial transformations from list file (cf. obsolete tdreg)
  if (dofin_list_name) {
    cerr << "Option -dofins currently not implemented" << endl;
    exit(1);
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
  if (!registration.Read(params, verbose > 2)) {
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
  Array<double>                 image_times;
  Array<unique_ptr<BaseImage> > images;

  if (image_names.size() == 1) {
    unique_ptr<BaseImage> sequence(BaseImage::New(image_names[0].c_str()));
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
  // Read initial transformation
  unique_ptr<Transformation> dofin;
  if (dofin_name) {
    if (IsIdentity(dofin_name)) dofin.reset(new RigidTransformation());
    else                        dofin.reset(Transformation::New(dofin_name));
    registration.InitialGuess(dofin.get());
  }

  // ---------------------------------------------------------------------------
  // Read mask which defines domain on which similarity is evaluated
  unique_ptr<BinaryImage> mask;
  if (mask_name) {
    mask.reset(new BinaryImage(mask_name));
    registration.Domain(mask.get());
  }

  MIRTK_DEBUG_TIMING(1, "reading input data");

  // ---------------------------------------------------------------------------
  // Run registration
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

  Transformation *dofout = NULL;
  registration.Output(&dofout);

  // Write initial parameters to file
  // Note: Will be overwritten once the registration finished successfully
  //       with the actual parameters used below.
  if (parout_name) {
    registration.GuessParameter();
    registration.Write(parout_name);
  }

#ifdef HAVE_TBB
  tbb::tick_count start_time = tbb::tick_count::now();
#else
  clock_t start_time = clock();
#endif

  registration.Run();

  if (verbose) {
    double elapsed_time;
#ifdef HAVE_TBB
    elapsed_time = (tbb::tick_count::now() - start_time).seconds();
#else
    elapsed_time = static_cast<double>(clock() - start_time)
                 / static_cast<double>(CLOCKS_PER_SEC);
#endif
    int m = floor(elapsed_time / 60.0);
    int s = round(elapsed_time - m * 60);
    if (s == 60) m += 1, s = 0;
    cout << "\nFinished in " << m << " min " << s << " sec" << endl;
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

  // Clean up
  delete dofout;
  registration.DeleteObserver(logger);
  registration.DeleteObserver(debugger);

  return 0;
}

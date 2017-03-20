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
#include "mirtk/GenericImage.h"
#include "mirtk/Transformation.h"
#include "mirtk/AffineTransformation.h"
#include "mirtk/FreeFormTransformation.h"


using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " [<dof>...] [options]" << endl;
  cout << "\n";
  cout << "Description:" << endl;
  cout << "  Calculates registration transformation quality measures such as the voxel-wise\n";
  cout << "  cumulative or mean inverse-consistency error (CICE/MICE) for pairs of forward\n";
  cout << "  transformation from target to source and backward transformation from source\n";
  cout << "  to target image. Another voxel-wise measure that can be computed using this\n";
  cout << "  program are the cumulative or mean transitivity error (CTE/MTE) given three\n";
  cout << "  transformations, from target (A) to B, from B to C, and from C to A again.\n";
  cout << "\n";
  cout << "  This tool can further be used to evaluate the deviation of N affine transformations\n";
  cout << "  which map a template to each one of N images from the barycenter corresponding\n";
  cout << "  to the identity transformation. A larger deviation indicates a stronger bias of\n";
  cout << "  the template towards a subset of the images and/or unbalanced registration errors.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  <dof>   File path of input transformation. Depending on the evaluation measure,\n";
  cout << "          multiple consecutive transformations are processed in tuples of fixed\n";
  cout << "          size and an error measure is computed for each such tuple of transformations.\n";
  cout << "          In the default case of 1-tuples, each transformation is evaluated separately.\n";
  cout << "          For other evaluation measures which are based on more than just one transformation\n";
  cout << "          see the list of optional arguments below.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -dofs <file>\n";
  cout << "      Text file with (additional) input transformation paths. (default: none)\n";
  cout << "  -output <file>\n";
  cout << "      Name of voxel-wise output error map. (default: none)\n";
  cout << "  -target <file>\n";
  cout << "      Target image. (default: FFD lattice)\n";
  cout << "  -mask <file>\n";
  cout << "      Target region-of-interest mask. (default: none)\n";
  cout << "  -padding <value>\n";
  cout << "      Target background value, unused if -mask specified. (default: none)\n";
  cout << "  -cumulative\n";
  cout << "      Whether to compute the cumulative error. (default: mean)\n";
  cout << "  -inverse-consistency\n";
  cout << "      Request evaluation of inverse-consistency error for each\n";
  cout << "      pair of input transformations when applied in the given order.\n";
  cout << "  -transitivity\n";
  cout << "      Request evaluation of transitivity error for each triple\n";
  cout << "      of input transformations when applied in the given order.\n";
  cout << "  -barycentricity\n";
  cout << "      Evaluate norm of the average log map of all input transformations.\n";
  cout << "      When the transformations map each image to the common barycenter,\n";
  cout << "      the reported value should be (close to) zero.\n";
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Types
// =============================================================================

// -----------------------------------------------------------------------------
/// Enumeration of implemented evaluation metrics
enum Metric
{
  None,
  InverseConsistency,
  Transitivity,
  Barycentricity
};

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
/// Read input file names from text file
int read_list_file(const char *list_name, Array<string>& paths)
{
  string   base_dir = Directory(list_name); // Base directory containing files
  ifstream iff(list_name);                  // List input file stream
  string   line;                            // Input line
  int      l = 0;                           // Line number
  size_t   pos;

  // Read base directory for relative file paths
  if (!getline(iff, line)) {
    cerr << "Error: Cannot parse list file " << list_name << endl;
    return 0;
  }

  if (!base_dir.empty() && line[0] != PATHSEP) {
    if (line != ".") base_dir += PATHSEP + line;
  } else {
    base_dir = line;
  }
  l++;

  while (getline(iff, line)) {
    l++;
    // Trim line string
    pos = line.find_first_not_of(" \t");
    if (pos != string::npos) line = line.substr(pos);
    pos = line.find_last_not_of(" \t");
    if (pos != string::npos) line = line.substr(0, pos + 1);
    // Ignore blank lines and comment lines starting with # character
    if (line.empty() || line[0] == '#') continue;
    // Make relative file paths absolute
    if (!base_dir.empty() && line[0] != PATHSEP) line = base_dir + PATHSEP + line;
    paths.push_back(line);
  }

  return static_cast<int>(paths.size());
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  // Positional arguments
  REQUIRES_POSARGS(0);
  Array<string> dofin_name;
  for (ALL_POSARGS) dofin_name.push_back(ARGUMENT);

  // Optional arguments
  const char *list_name     = NULL;
  const char *output_name   = NULL;
  const char *target_name   = NULL;
  const char *mask_name     = NULL;
  double      padding_value = -numeric_limits<double>::infinity();
  bool        cumulative    = false;
  bool        squared       = false;
  Metric      metric        = None;

  for (ALL_OPTIONS) {
    if      (OPTION("-output")) output_name = ARGUMENT;
    else if (OPTION("-dofs"))   list_name   = ARGUMENT;
    else if (OPTION("-target")) target_name = ARGUMENT;
    else if (OPTION("-mask"))   mask_name   = ARGUMENT;
    else if (OPTION("-padding")) {
      const char *arg = ARGUMENT;
      if (!FromString(arg, padding_value)) {
        cerr << "Invalid -padding value" << endl;
        exit(1);
      }
    }
    else if (OPTION("-cumulative")) cumulative = true;
    else if (OPTION("-squared"))    squared    = true;
    else if (OPTION("-inverse-consistency") || OPTION("-ice") || OPTION("-ic")) {
      metric = InverseConsistency;
    }
    else if (OPTION("-transitivity") || OPTION("-te")) metric = Transitivity;
    else if (OPTION("-mice")) { metric = InverseConsistency; cumulative = false; }
    else if (OPTION("-cice")) { metric = InverseConsistency; cumulative = true;  }
    else if (OPTION("-mte"))  { metric = Transitivity;       cumulative = false; }
    else if (OPTION("-cte"))  { metric = Transitivity;       cumulative = true;  }
    else if (OPTION("-barycentric-equation")) metric = Barycentricity;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (list_name) {
    read_list_file(list_name, dofin_name);
  }
  if (dofin_name.empty()) {
    PrintHelp(EXECNAME);
    exit(1);
  }
  if (metric == None) {
    metric = Barycentricity;
  }

  if (metric == Barycentricity) {

    Matrix residual(4, 4);
    AffineTransformation dof;
    for (auto fname : dofin_name) {
      dof.Read(fname.c_str());
      residual += dof.GetMatrix().Log();
    }
    residual /= dofin_name.size();
    cout << residual.Norm() << endl;

  } else {

    // Shared pointers for input tuples
    Array<SharedPtr<Transformation> > dofs;
    if      (metric == InverseConsistency) dofs.resize(2);
    else if (metric == Transitivity)       dofs.resize(3);
    else                                   dofs.resize(1);

    if (dofin_name.size() % dofs.size() != 0) {
      cerr << "Invalid number of input transformations. Require " << dofs.size() << "-tuples for the evaluation!" << endl;
      exit(1);
    }

    // Read first transformation
    SharedPtr<Transformation> dof(Transformation::New(dofin_name[0].c_str()));

    // Read target, foreground mask, and/or determine attributes for voxel-wise metrics
    GreyImage        target;
    BinaryImage      mask;
    ImageAttributes  attr;
    WorldCoordsImage i2w;

    if (mask_name) {
      InitializeIOLibrary();
      mask.Read(mask_name);
      attr = mask.Attributes();
    }
    if (target_name) {
      InitializeIOLibrary();
      target.Read(target_name);
      if (!IsInf(padding_value)) target.PutBackgroundValueAsDouble(padding_value, true);
      attr = target.Attributes();
    }
    if (!attr) {
      FreeFormTransformation *ffd = dynamic_cast<FreeFormTransformation *>(dof.get());
      if (!ffd) {
        cerr << "First input transformation is no 3D FFD. Either -target or -mask image is required!" << endl;
        exit(1);
      }
      attr = ffd->Attributes();
    }
    if (!mask.IsEmpty() && !mask.Attributes().EqualInSpace(attr)) {
      cerr << "Mask and target image must have identical spatial attributs!" << endl;
      exit(1);
    }

    // Pre-compute world coordinates of reference voxels
    if (verbose) cout << "Pre-computing world coordinates of input points" << endl;
    i2w.Initialize(attr, 3);
    double *x = i2w.Data(0, 0, 0, 0);
    double *y = i2w.Data(0, 0, 0, 1);
    double *z = i2w.Data(0, 0, 0, 2);
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i) {
      *x = i, *y = j, *z = k;
      i2w.ImageToWorld(*x, *y, *z);
      ++x, ++y, ++z;
    }

    // Compute inverse-consistency error for each pair of transformations
    GenericImage<float> output(attr);
    const int nvox = output.NumberOfVoxels();
    int ntuple = 0;
    int ntotal = static_cast<int>(dofin_name.size() / dofs.size());

    size_t i, j;
    double x2, y2, z2, d;
    for (i = 0; i < dofin_name.size(); i += dofs.size()) {
      ++ntuple;
      if (verbose) {
        const char *suffix = "th";
        if      (ntuple == 1) suffix = "st";
        else if (ntuple == 2) suffix = "nd";
        else if (ntuple == 3) suffix = "rd";
        cout << "Performing evaluation of " << ntuple << suffix << " " << dofs.size() << "-tuple" << " out of " << ntotal << endl;
      }
      for (j = 0; j < dofs.size(); ++j) {
        if (i == 0 && j == 0) dofs[j] = dof;
        else dofs[j].reset(Transformation::New(dofin_name[i + j].c_str()));
      }
      const double *x1 = i2w.Data(0, 0, 0, 0);
      const double *y1 = i2w.Data(0, 0, 0, 1);
      const double *z1 = i2w.Data(0, 0, 0, 2);
      float *err = output.Data();
      for (int n = 0; n < nvox; ++n) {
        if ((target.IsEmpty() || target.IsForeground(n)) && (mask.IsEmpty() || mask.Get(n))) {
          x2 = x1[n], y2 = y1[n], z2 = z1[n];
          for (j = 0; j < dofs.size(); ++j) dofs[j]->Transform(x2, y2, z2);
          x2 -= x1[n], y2 -= y1[n], z2 -= z1[n];
          d = x2 * x2 + y2 * y2 + z2 * z2;
          err[n] += static_cast<float>(squared ? d : sqrt(d));
        }
      }
    }

    // Calculate mean inverse-consistency error
    if (!cumulative && ntuple > 0) output /= ntuple;

    // Write voxel-wise error map
    if (output_name) {
      InitializeIOLibrary();
      output.Write(output_name);
    }

  }

  return 0;
}

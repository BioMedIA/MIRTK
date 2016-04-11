/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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
#include "mirtk/Indent.h"
#include "mirtk/GenericImage.h"
#include "mirtk/Transformation.h"
#include "mirtk/RegisteredImage.h"
#include "mirtk/SimilarityMeasure.h"
#include "mirtk/ImageSimilarity.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <target> <source>... [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Computes the (dis-)similarity of two intensity images.\n";
  cout << "  If more than one source image is given, the (dis-)similarity between\n";
  cout << "  each of these and the target is evaluated. By default, common image\n";
  cout << "  image (dis-)similarity metrics are reported. One or more other metrics\n";
  cout << "  can be chosen using :option:`-metric` which can be given multiple times.\n";
  cout << "  The input source image can either be pre-aligned with the target image\n";
  cout << "  or the output of a previous registration can be specified using\n";
  cout << "  :option:`-dofin`.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  target   Target image.\n";
  cout << "  source   (Transformed) source image.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -dofin <file>      Source image transformation. (default: Id)\n";
  cout << "  -interp <mode>     Interpolation mode. (default: linear with padding)\n";
  cout << "  -metric <sim>...   Image (dis-)similarity measure(s). (default: AvgSI)\n";
  cout << "  -Tp <value>        Target image background value/threshold. (default: NaN)\n";
  cout << "  -Sp <value>        Source image background value/threshold. (default: NaN)\n";
  cout << "  -Rx1 <int>         Leftmost  target voxel index along x axis.\n";
  cout << "  -Rx2 <int>         Rightmost target voxel index along x axis.\n";
  cout << "  -Ry1 <int>         Leftmost  target voxel index along y axis.\n";
  cout << "  -Ry2 <int>         Rightmost target voxel index along y axis.\n";
  cout << "  -Rz1 <int>         Leftmost  target voxel index along z axis.\n";
  cout << "  -Rz2 <int>         Rightmost target voxel index along z axis.\n";
  cout << "\n";
  cout << "Output format options:\n";
  cout << "  -precision <int>   Number of significant digits. (default: 5)\n";
  cout << "  -delim <char>      Delimiter for output of multiple metric values. (default: ,)\n";
  cout << "  -table             Output in non-verbose tabular format, i.e., \"-v 0\".\n";
  cout << "  -csv               Output as comma separated values, i.e., \"-v 0 -delim ','\"\n";
  cout << "  -tsv               Output as tab   separated values, i.e., \"-v 0 -delim '\\t'\"\n";
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  // Check command line
  REQUIRES_POSARGS(2);

  // Parse positional arguments
  const char *target_name = POSARG(1);
  Array<const char *> source_name;
  source_name.reserve(NUM_POSARGS-1);
  source_name.push_back(POSARG(2));
  for (OPTIONAL_POSARGS) {
    source_name.push_back(POSARG(ARGIDX));
  }

  // Parse optional arguments
  Array<SimilarityMeasure> metric;
  InterpolationMode interp   = Interpolation_LinearWithPadding;
  const char *dofin_name     = nullptr;
  double      target_padding = numeric_limits<double>::quiet_NaN();
  double      source_padding = numeric_limits<double>::quiet_NaN();
  int         digits         = 5;
  string      delim;

  verbose = 1; // by default, include textual description of output value
               // can be disabled by caller using "-v 0" option in order to
               // only append the overlap value (excl. newline) to a file

  int i1 =  0, j1 =  0, k1 =  0;
  int i2 = -1, j2 = -1, k2 = -1;

  for (ALL_OPTIONS) {
    if      (OPTION("-dofin")) dofin_name = ARGUMENT;
    else if (OPTION("-interp")) PARSE_ARGUMENT(interp);
    else if (OPTION("-Tp")) PARSE_ARGUMENT(target_padding);
    else if (OPTION("-Sp")) PARSE_ARGUMENT(source_padding);
    else if (OPTION("-metric")) {
      do {
        SimilarityMeasure m;
        PARSE_ARGUMENT(m);
        metric.push_back(m);
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-precision")) PARSE_ARGUMENT(digits);
    else if (OPTION("-delim")) {
      delim = ARGUMENT;
      size_t pos;
      while ((pos = delim.find("\\n")) != string::npos) {
        delim.replace(pos, 2, "\n");
      }
      while ((pos = delim.find("\\t")) != string::npos) {
        delim.replace(pos, 2, "\t");
      }
    }
    else if (OPTION("-table")) verbose = 0;
    else if (OPTION("-csv")) verbose = 0, delim = ',';
    else if (OPTION("-tsv")) verbose = 0, delim = '\t';
    else if (OPTION("-Rx1")) PARSE_ARGUMENT(i1);
    else if (OPTION("-Rx2")) PARSE_ARGUMENT(i2);
    else if (OPTION("-Ry1")) PARSE_ARGUMENT(j1);
    else if (OPTION("-Ry2")) PARSE_ARGUMENT(j2);
    else if (OPTION("-Rz1")) PARSE_ARGUMENT(k1);
    else if (OPTION("-Rz2")) PARSE_ARGUMENT(k2);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (metric.empty()) {
    FatalError("Default similarity metric AvgSI not implemented yet, use -metric option");
  }

  if (delim.empty()) {
    delim = ',';
    if (verbose > 0) delim += ' ';
  }

  // Initialize I/O module
  InitializeIOLibrary();

  // Input and (transformed) image instances
  GenericImage<double> target_image, source_image;
  RegisteredImage      target,       source;

  target.InputImage(&target_image);
  source.InputImage(&source_image);
  source.InterpolationMode(interp);

  target.SelfUpdate(false);
  source.SelfUpdate(false);

  // Read target image
  target_image.Read(target_name);
  if (target_image.T() > 1) {
    FatalError("Target image must be 2D or 3D image!");
  }
  target_image.PutBackgroundValueAsDouble(target_padding, true);

  // Read input transformation
  unique_ptr<Transformation> dofin;
  if (dofin_name) {
    dofin.reset(Transformation::New(dofin_name));
    source.Transformation(dofin.get());
  }

  // Target region of interest
  if (i2 < 0) i2 = target_image.X();
  if (j2 < 0) j2 = target_image.Y();
  if (k2 < 0) k2 = target_image.Z();

  ImageAttributes target_region = target_image.Attributes();
  target_region._x       = i2 - i1;
  target_region._y       = j2 - j1;
  target_region._z       = k2 - k1;
  target_region._xorigin = .0;
  target_region._yorigin = .0;
  target_region._zorigin = .0;

  // Adjust origin
  double x1 = i1, y1 = j1, z1 = k1;
  target_image.ImageToWorld(x1, y1, z1);
  double x2 = .0, y2 = .0, z2 = .0;
  target_region.LatticeToWorld(x2, y2, z2);
  target_region._xorigin = x1 - x2;
  target_region._yorigin = y1 - y2;
  target_region._zorigin = z1 - z2;

  // Initialize (transformed) images
  target.Initialize(target_region);
  source.Initialize(target_region);

  target.Update();

  // Print header
  if (verbose < 1) {
    cout << "Target" << delim << "Source";
    for (size_t i = 0; i < metric.size(); ++i) {
      cout << delim << ToString(metric[i]);
    }
    cout << endl;
  }

  // Evaluate similarity of (transformed) source image(s)
  const string target_id = FileName(target_name);
  for (size_t n = 0; n < source_name.size(); ++n) {

    if (verbose > 0) cout << "Target = ";
    cout << target_id << delim;
    if (verbose > 0) cout << "Source = ";
    cout << FileName(source_name[n]);

    // Read source image
    source_image.Read(source_name[n]);
    source_image.PutBackgroundValueAsDouble(source_padding, true);
    source.Update();

    // Compute similarity measure(s)
    for (size_t i = 0; i < metric.size(); ++i) {
      cout << delim;
      if (verbose > 0) cout << ToString(metric[i]) << " = ";
      unique_ptr<ImageSimilarity> sim(ImageSimilarity::New(metric[i]));
      sim->Target(&target);
      sim->Source(&source);
      sim->Domain(target_region);
      sim->DivideByInitialValue(false);
      sim->SkipTargetInitialization(true);
      sim->SkipSourceInitialization(true);
      sim->Initialize();
      sim->Update(false);
      cout << setprecision(digits) << sim->RawValue();
      sim->ReleaseTarget();
      sim->ReleaseSource();
    }

    cout << endl;
  }

  return 0;
}

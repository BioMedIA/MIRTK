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

#include <mirtkImageIOConfig.h>
#include <mirtkGenericImage.h>
#include <mirtkImageFunction.h>
#include <mirtkHistogram1D.h>

using namespace mirtk;

// TODO: Add -labels option to evaluate overlap of all labels at once

// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <target> <source>... [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Computes the overlap of either two intensity images (average SI)" << endl;
  cout << "  or two segmentations (see :option:`-label`). If more than one source image" << endl;
  cout << "  is given, the overlap between each of these and the target is evaluated." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  target   Target image/segmentation." << endl;
  cout << "  source   Transformed source image/segmentation." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -Rx1 <int>         Region of interest" << endl;
  cout << "  -Ry1 <int>         Region of interest" << endl;
  cout << "  -Rz1 <int>         Region of interest" << endl;
  cout << "  -Rx2 <int>         Region of interest" << endl;
  cout << "  -Ry2 <int>         Region of interest" << endl;
  cout << "  -Rz2 <int>         Region of interest" << endl;
  cout << "  -Tp <value>        Padding value in target." << endl;
  cout << "  -label <int>       Segmentation label. (default: none)" << endl;
  cout << "  -precision <int>   Number of significant digits. (default: 2)" << endl;
  cout << "  -metric <Sensitivity|Specificity|Dice|Jaccard>   Overlap metric. (default: Dice)" << endl;
  cout << "  -delim <char>      Delimiter for output of multiple overlap values. (default: ,)" << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// Implemented (segmentation) overlap metrics
enum OverlapMetric { Unknown, Sensitivity, Specificity, Dice, Jaccard };

// -----------------------------------------------------------------------------
istream &operator >>(istream &is, OverlapMetric &metric)
{
  string str;
  is >> str;
  if      (str == "Sensitivity") metric = Sensitivity;
  else if (str == "Specificity") metric = Specificity;
  else if (str == "Dice coefficient"   || str == "Dice")    metric = Dice;
  else if (str == "Jaccard similarity" || str == "Jaccard") metric = Jaccard;
  else {
    metric = Unknown;
    is.setstate(ios::failbit);
  }
  return is;
}

// -----------------------------------------------------------------------------
ostream &operator <<(ostream &os, const OverlapMetric &metric)
{
  switch (metric) {
    case Sensitivity: os << "Sensitivity";        break;
    case Specificity: os << "Specificity";        break;
    case Dice:        os << "Dice coefficient";   break;
    case Jaccard:     os << "Jaccard similarity"; break;
    default:          os << "Unknown metric";     break;
  }
  return os;
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
  source_name.push_back(POSARG(2));

  int nposarg;
  for (nposarg = 3; nposarg < argc; ++nposarg) {
    if (argv[nposarg][0] == '-') break;
    source_name.push_back(argv[nposarg]);
  }
  --nposarg;

  // Parse optional arguments
  OverlapMetric metric  = Unknown;
  GreyPixel     padding = MIN_GREY;
  GreyPixel     label   = -1;
  const char   *delim   = ",";
  int           digits  = 2;

  verbose = 1; // by default, include textual description of output value
               // can be disabled by caller using "-v 0" option in order to
               // only append the overlap value (excl. newline) to a file

  int i1 =  0, j1 =  0, k1 =  0;
  int i2 = -1, j2 = -1, k2 = -1;

  for (ARGUMENTS_AFTER(nposarg)) {
    if      (OPTION("-Tp"))  padding = atof(ARGUMENT);
    else if (OPTION("-label")) label = atof(ARGUMENT);
    else if (OPTION("-metric")) {
      const char *arg = ARGUMENT;
      if (!FromString(arg, metric)) {
        FatalError("Invalid -metric argument: " << arg);
      }
    }
    else if (OPTION("-precision")) digits = atoi(ARGUMENT);
    else if (OPTION("-delim"))     delim  = ARGUMENT;
    else if (OPTION("-Rx1")) i1 = atoi(ARGUMENT);
    else if (OPTION("-Rx2")) i2 = atoi(ARGUMENT);
    else if (OPTION("-Ry1")) j1 = atoi(ARGUMENT);
    else if (OPTION("-Ry2")) j2 = atoi(ARGUMENT);
    else if (OPTION("-Rz1")) k1 = atoi(ARGUMENT);
    else if (OPTION("-Rz2")) k2 = atoi(ARGUMENT);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read target image
  InitializeImageIOLibrary();
  GreyImage target(target_name);
  const ImageAttributes target_attributes = target.Attributes();

  // Extract region of interest
  if (i2 < 0) i2 = target.X();
  if (j2 < 0) j2 = target.Y();
  if (k2 < 0) k2 = target.Z();
  if ((i1 != 0) || (i2 != target.GetX()) ||
      (j1 != 0) || (j2 != target.GetY()) ||
      (k1 != 0) || (k2 != target.GetZ())) {
    target = target.GetRegion(i1, j1, k1, i2, j2, k2);
  }

  for (size_t n = 0; n < source_name.size(); ++n) {

    if (verbose) {
      if (source_name.size() > 1) cout << source_name[n] << ": ";
    } else if (n > 0) cout << delim;

    // Read source image
    GreyImage source(source_name[n]);
    if (source.Attributes() != target_attributes) {
      FatalError("Input images must be sampled on same lattice!");
    }

    // Extract region of interest
    if ((i1 != 0) || (i2 != source.GetX()) ||
        (j1 != 0) || (j2 != source.GetY()) ||
        (k1 != 0) || (k2 != source.GetZ())) {
      source = source.GetRegion(i1, j1, k1, i2, j2, k2);
    }

    // ---------------------------------------------------------------------------
    // Segmentation overlap
    if (label > -1) {

      // Default metric
      if (metric == Unknown) metric = Dice;

      // Determine TP, FP, TN, FN
      int tp = 0, fp = 0, tn = 0, fn = 0;
      for (int k = 0; k < target.Z(); ++k) {
        for (int j = 0; j < target.Y(); ++j) {
          for (int i = 0; i < target.X(); ++i) {
            if (target(i, j, k) == label) {
              if (source(i, j, k) == label) ++tp;
              else                          ++fn;
            }
            else if (source(i, j, k) == label) ++fp;
            else                               ++tn;
          }
        }
      }

      // Compute overlap metric
      double overlap = .0;
      switch (metric) {
        case Sensitivity:
          overlap = double(tp) / double(tp + fn);
          break;
        case Specificity:
          overlap = double(tn) / double(tn + fp);
          break;
        case Dice:
          overlap = 2.0 * double(tp) / double((fp + tp) + (tp + fn));
          break;
        case Jaccard:
          overlap = double(tp) / double(fp + tp + fn);
          break;
        default:
          FatalError(metric << " not implemented");
          exit(1);
      }

      if (verbose) {
        const int d = (fequal(overlap, 1.0, 1e-6) ? 3 : digits);
        cout << metric << " = " << setprecision(d) << (100.0 * overlap) << "%" << endl;
      } else {
        cout << setprecision(digits) << overlap;
      }

    // ---------------------------------------------------------------------------
    // Intensity image overlap
    } else {

      // Default metric
      if (metric != Unknown) {
        FatalError(metric << " not suitable for intensity images");
        exit(1);
      }

      // Determine maximum intensity
      int max = 0;
      for (int k = 0; k < target.Z(); ++k) {
        for (int j = 0; j < target.Y(); ++j) {
          for (int i = 0; i < target.X(); ++i) {
            if (target(i, j, k) > padding) {
              if (target(i, j, k) > max) max = target(i, j, k);
              if (source(i, j, k) > max) max = source(i, j, k);
            }
          }
        }
      }

      // Compute image overlap
      Histogram1D<int> histogramA (max);
      Histogram1D<int> histogramB (max);
      Histogram1D<int> histogramAB(max);

      for (int k = 0; k < target.Z(); ++k)
      for (int j = 0; j < target.Y(); ++j)
      for (int i = 0; i < target.X(); ++i) {
        if (target(i, j, k) > padding) {
          histogramA.Add(target(i, j, k));
          histogramB.Add(source(i, j, k));
          if (target(i, j, k) == source(i, j, k)) histogramAB.Add(target(i, j, k));
        }
      }
      if (histogramA.NumberOfSamples() == 0) {
        FatalError("No samples in target histogram");
      }

      // Calculate average SI
      double w, si = .0;
      for (int i = 1; i < histogramA.NumberOfBins(); ++i) {
        w = histogramA(i) / double(histogramA.NumberOfSamples() - histogramA(0));
        if (histogramA(i) + histogramB(i) != 0) {
          si += w * 2.0 * histogramAB(i) / double(histogramA(i) + histogramB(i));
        }
      }

      if (verbose) cout << "Average SI = ";
      cout << si;
      if (verbose) cout << endl;

    }
  }

  return 0;
}

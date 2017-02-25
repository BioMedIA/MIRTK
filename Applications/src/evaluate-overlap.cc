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

#include "mirtk/IOConfig.h"
#include "mirtk/GenericImage.h"
#include "mirtk/ImageFunction.h"
#include "mirtk/Histogram1D.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <target> [<source>...] [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Computes the overlap of either two intensity images (average SI)\n";
  cout << "  or two segmentations (see :option:`-label`). If more than one source image\n";
  cout << "  is given, the overlap between each of these and the target is evaluated.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  target   Target image/segmentation.\n";
  cout << "  source   Transformed source image/segmentation.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -images <file>\n";
  cout << "      Text file with source image file paths.\n";
  cout << "  -images-delim, -images-delimiter <c>\n";
  cout << "      Delimiter used in :option:`-images` file.\n";
  cout << "      (default: ',' for .csv, '\\t' for .tsv, and ' ' otherwise)\n";
  cout << "  -Rx1 <int>\n";
  cout << "      Region of interest lower index threshold in x dimension.\n";
  cout << "  -Ry1 <int>\n";
  cout << "      Region of interest lower index threshold in y dimension.\n";
  cout << "  -Rz1 <int>\n";
  cout << "      Region of interest lower index threshold in z dimension.\n";
  cout << "  -Rx2 <int>\n";
  cout << "      Region of interest upper index threshold in x dimension.\n";
  cout << "  -Ry2 <int>\n";
  cout << "      Region of interest upper index threshold in y dimension.\n";
  cout << "  -Rz2 <int>\n";
  cout << "      Region of interest upper index threshold in z dimension.\n";
  cout << "  -Tp <value>\n";
  cout << "      Padding value in target.\n";
  cout << "  -label <int>...\n";
  cout << "      Segmentation labels. When more than one specified, the segments are merged.\n";
  cout << "      This option can be given multiple times to evaluate the overlap of more than\n";
  cout << "      one (merged) segment at once. (default: none)\n";
  cout << "  -labels\n";
  cout << "      All (positive) segmentation labels individually.\n";
  cout << "  -precision <int>\n";
  cout << "      Number of significant digits. (default: 2)\n";
  cout << "  -metric <Sensitivity|Specificity|Dice|Jaccard>\n";
  cout << "      Overlap metric. (default: Dice)\n";
  cout << "  -delimiter, -delim <char>\n";
  cout << "      Delimiter for output of multiple overlap values. (default: ,)\n";
  cout << "  -one-row-per-source [yes|no]\n";
  cout << "      Print table with overlap measures for all segments in a row,\n";
  cout << "      one row per input source image which is compared to the target\n";
  cout << "      segmentation. By default, the output table with verbosity 0 is\n";
  cout << "      the transpose table and printed without table header. (default: no)\n";
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// Implemented (segmentation) overlap metrics
enum OverlapMetric { UnknownMetric, Sensitivity, Specificity, Dice, Jaccard };

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
    metric = UnknownMetric;
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

// -----------------------------------------------------------------------------
/// Read image file names from input list file
int ReadImageListFile(const char *image_list_name, Array<string> &names, const char *delim = ",")
{
  if (delim == nullptr) {
    const string ext = Extension(image_list_name);
    if      (ext == ".tsv") delim = "\t";
    else if (ext == ".csv") delim = ",";
    else delim = " ";
  }

  string   base_dir = Directory(image_list_name); // Base directory for relative paths
  ifstream iff(image_list_name);                  // List input file stream
  string   line;                                  // Input line
  int      l = 0;                                 // Line number

  names.clear();

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
  }

  return static_cast<int>(names.size());
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  // Check command line
  REQUIRES_POSARGS(1);

  // Parse positional arguments
  const char *target_name = POSARG(1);
  Array<string> source_name;

  int nposarg;
  for (nposarg = 2; nposarg < argc; ++nposarg) {
    if (argv[nposarg][0] == '-') break;
    source_name.push_back(argv[nposarg]);
  }
  --nposarg;

  // Parse optional arguments
  const char    *image_list  = nullptr;
  const char    *image_delim = nullptr;
  OverlapMetric  metric      = UnknownMetric;
  GreyPixel      padding     = MIN_GREY;
  const char    *delim       = ",";
  int            digits      = 2;
  bool           all_labels  = false;
  bool           row_per_seg = false;
  Array<OrderedSet<GreyPixel> > segments;

  verbose = 1; // by default, include textual description of output value
               // can be disabled by caller using "-v 0" option in order to
               // only append the overlap value (excl. newline) to a file

  int i1 =  0, j1 =  0, k1 =  0;
  int i2 = -1, j2 = -1, k2 = -1;

  for (ARGUMENTS_AFTER(nposarg)) {
    if (OPTION("-Tp")) PARSE_ARGUMENT(padding);
    else if (OPTION("-label") || OPTION("-segment")) {
      GreyPixel label;
      OrderedSet<GreyPixel> segment;
      do {
        PARSE_ARGUMENT(label);
        segment.insert(label);
      } while (HAS_ARGUMENT);
      segments.push_back(segment);
    }
    else if (OPTION("-images")) image_list = ARGUMENT;
    else if (OPTION("-images-delimiter") || OPTION("-images-delim")) image_delim = ARGUMENT;
    else if (OPTION("-labels")) all_labels = true;
    else if (OPTION("-metric")) PARSE_ARGUMENT(metric);
    else if (OPTION("-precision")) PARSE_ARGUMENT(digits);
    else if (OPTION("-delim")) delim = ARGUMENT;
    else if (OPTION("-Rx1")) PARSE_ARGUMENT(i1);
    else if (OPTION("-Rx2")) PARSE_ARGUMENT(i2);
    else if (OPTION("-Ry1")) PARSE_ARGUMENT(j1);
    else if (OPTION("-Ry2")) PARSE_ARGUMENT(j2);
    else if (OPTION("-Rz1")) PARSE_ARGUMENT(k1);
    else if (OPTION("-Rz2")) PARSE_ARGUMENT(k2);
    else HANDLE_BOOLEAN_OPTION("one-row-per-source", row_per_seg);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read source image file paths from text file
  if (image_list) {
    if (ReadImageListFile(image_list, source_name, image_delim) == 0) {
      Warning("No source image file paths read from file: " << image_list);
    }
  }
  if (source_name.empty()) {
    FatalError("No source images specified!");
  }

  // Read target image
  InitializeIOLibrary();
  GreyImage target(target_name);
  const ImageAttributes target_attributes = target.Attributes();

  if (all_labels) {
    OrderedSet<int> labels;
    bool cont;
    do {
      cont = false;
      for (auto it = segments.begin(); it != segments.end(); ++it) {
        if (it->size() == 1) {
          labels.insert(*(it->begin()));
          segments.erase(it);
          cont = true;
          break;
        }
      }
    } while (cont);
    for (int i = 0; i < target.NumberOfVoxels(); ++i) {
      GreyPixel label = target(i);
      if (label > 0) {
        labels.insert(label);
      }
    }
    OrderedSet<GreyPixel> segment;
    for (auto label : labels) {
      segment.clear();
      segment.insert(label);
      segments.push_back(segment);
    }
  }

  // Extract region of interest
  if (i2 < 0) i2 = target.X();
  if (j2 < 0) j2 = target.Y();
  if (k2 < 0) k2 = target.Z();
  if ((i1 != 0) || (i2 != target.GetX()) ||
      (j1 != 0) || (j2 != target.GetY()) ||
      (k1 != 0) || (k2 != target.GetZ())) {
    target = target.GetRegion(i1, j1, k1, i2, j2, k2);
  }

  // -------------------------------------------------------------------------------
  // Intensity image overlap
  if (segments.empty()) {

    for (size_t n = 0; n < source_name.size(); ++n) {

      if (verbose) {
        if (source_name.size() > 1) cout << source_name[n] << ": ";
      } else if (n > 0) cout << delim;

      // Read source image
      GreyImage source(source_name[n].c_str());
      if (source.Attributes() != target_attributes) {
        FatalError("Input images must be sampled on same lattice!");
      }

      // Extract region of interest
      if ((i1 != 0) || (i2 != source.GetX()) ||
          (j1 != 0) || (j2 != source.GetY()) ||
          (k1 != 0) || (k2 != source.GetZ())) {
        source = source.GetRegion(i1, j1, k1, i2, j2, k2);
      }

      // Default metric
      if (metric != UnknownMetric) {
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

  // -------------------------------------------------------------------------------
  // Segmentation overlap
  } else {

    // Default metric
    if (metric == UnknownMetric) metric = Dice;

    if (row_per_seg) {

      for (size_t roi = 0; roi < segments.size(); ++roi) {
        if (roi > 0) cout << delim;
        const auto &labels = segments[roi];
        for (auto it = labels.begin(); it != labels.end(); ++it) {
          if (it != labels.begin()) cout << "+";
          cout << *it;
        }
      }
      cout << "\n";

      for (size_t n = 0; n < source_name.size(); ++n) {

        // Read source image and extract region of interest
        GreyImage source(source_name[n].c_str());
        if (source.Attributes() != target_attributes) {
          FatalError("Input images must be sampled on same lattice!");
        }
        if ((i1 != 0) || (i2 != source.X()) ||
            (j1 != 0) || (j2 != source.Y()) ||
            (k1 != 0) || (k2 != source.Z())) {
          source = source.GetRegion(i1, j1, k1, i2, j2, k2);
        }

        // Iterative over segments
        for (size_t roi = 0; roi < segments.size(); ++roi) {
          const auto &labels = segments[roi];

          // Determine TP, FP, TN, FN
          int tp = 0, fp = 0, tn = 0, fn = 0;
          for (int k = 0; k < target.Z(); ++k)
          for (int j = 0; j < target.Y(); ++j)
          for (int i = 0; i < target.X(); ++i) {
            if (labels.find(target(i, j, k)) != labels.end()) {
              if (labels.find(source(i, j, k)) != labels.end()) ++tp;
              else                                              ++fn;
            }
            else if (labels.find(source(i, j, k)) != labels.end()) ++fp;
            else                                                   ++tn;
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

          if (roi > 0) cout << delim;
          cout << setprecision(digits) << overlap;
        }
      }

    } else {

      // Pre-load first n source images and extract region of interest
      GreyImage source(target_attributes);
      Array<GreyImage> sources(min(size_t(10u), source_name.size()));
      for (size_t n = 0; n < sources.size(); ++n) {
        sources[n].Read(source_name[n].c_str());
        if (sources[n].Attributes() != target_attributes) {
          FatalError("Input images must all be sampled on same lattice!");
        }
        if ((i1 != 0) || (i2 != sources[n].X()) ||
            (j1 != 0) || (j2 != sources[n].Y()) ||
            (k1 != 0) || (k2 != sources[n].Z())) {
          sources[n] = sources[n].GetRegion(i1, j1, k1, i2, j2, k2);
        }
      }

      // Evaluate overlap for each segment
      for (size_t roi = 0; roi < segments.size(); ++roi) {
        const auto &labels = segments[roi];
        if (segments.size() > 1) {
          if (verbose) {
            if (roi > 0 && source_name.size() > 1) {
              cout << "\n";
            }
            cout << "Label ";
            for (auto it = labels.begin(); it != labels.end(); ++it) {
              if (it != labels.begin()) cout << "+";
              cout << *it;
            }
            if (source_name.size() > 1) {
              cout << ":\n";
            } else {
              cout << ", ";
            }
          } else {
            for (auto it = labels.begin(); it != labels.end(); ++it) {
              if (it != labels.begin()) cout << "+";
              cout << *it;
            }
            cout << delim;
          }
        }
        for (size_t n = 0; n < source_name.size(); ++n) {

          if (verbose) {
            if (source_name.size() > 1) {
              if (segments.size() > 1) {
                cout << "  ";
              }
              cout << source_name[n] << ": ";
            }
          } else if (n > 0) cout << delim;

          // Read source image and extract region of interest
          if (n < sources.size()) {
            source = sources[n];
          } else {
            source.Read(source_name[n].c_str());
            if (source.Attributes() != target_attributes) {
              FatalError("Input images must be sampled on same lattice!");
            }
            if ((i1 != 0) || (i2 != source.GetX()) ||
                (j1 != 0) || (j2 != source.GetY()) ||
                (k1 != 0) || (k2 != source.GetZ())) {
              source = source.GetRegion(i1, j1, k1, i2, j2, k2);
            }
          }

          // Determine TP, FP, TN, FN
          int tp = 0, fp = 0, tn = 0, fn = 0;
          for (int k = 0; k < target.Z(); ++k) {
            for (int j = 0; j < target.Y(); ++j) {
              for (int i = 0; i < target.X(); ++i) {
                if (labels.find(target(i, j, k)) != labels.end()) {
                  if (labels.find(source(i, j, k)) != labels.end()) ++tp;
                  else                                              ++fn;
                }
                else if (labels.find(source(i, j, k)) != labels.end()) ++fp;
                else                                                   ++tn;
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
            cout << metric << " = " << setprecision(d) << (100.0 * overlap) << "%\n";
          } else {
            cout << setprecision(digits) << overlap;
          }
        }

        if (verbose == 0) cout << "\n";
      }

    }
  }

  return 0;
}

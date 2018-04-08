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
  cout << "  -label, -segment <value|lower..upper>...\n";
  cout << "      Segmentation labels. When more than one specified, the segments are merged.\n";
  cout << "      This option can be given multiple times to evaluate the overlap of more than\n";
  cout << "      one (merged) segment at once. Use '-label 0' to evaluate overlap of all labels\n";
  cout << "      merged into one segment (i.e., union of all labels). (default: none)\n";
  cout << "  -labels\n";
  cout << "      All (positive) segmentation labels individually. Can be combined with :option:`-label`.\n";
  cout << "  -probs [<threshold>]\n";
  cout << "      Evaluate overlap of class probability maps (fuzzy membership functions).\n";
  cout << "      The values of the input images are rescaled to [0, 1]. When no <threshold>\n";
  cout << "      is specified, a value of 0.5 is used to binarize the fuzzy segmentation masks.\n";
  cout << "      When no <threshold> is specified and the metric is either Dice (DSC) or Jaccard (JSC),\n";
  cout << "      these metrics are evaluated using their respective definition based on scalar products.\n";
  cout << "  -precision <int>\n";
  cout << "      Number of significant digits. (default: 2)\n";
  cout << "  -metric <name>\n";
  cout << "      Segmentation overlap metric (see https://en.wikipedia.org/wiki/Confusion_matrix).\n";
  cout << "      (default: Dice / F1-score)\n";
  cout << "  -delimiter, -delim <char>\n";
  cout << "      Delimiter for output of multiple overlap values. (default: ,)\n";
  cout << "  -table [<file>|stdout|cout]\n";
  cout << "      Write table with overlap measures for all segments in a row,\n";
  cout << "      one row per input source image which is compared to the target\n";
  cout << "      segmentation. By default, the output to STDIN with verbosity 0\n";
  cout << "      is the transpose table and excludes a table header. (default: off)\n";
  cout << "  -[no]header [on|off]\n";
  cout << "      Whether to output a header row when :option:`-table` is used. (default: on)\n";
  cout << "  -id [<name>]\n";
  cout << "      Print ID column with given <name> when :option:`-table` is used and multiple source\n";
  cout << "      images are given. When only one image is given, the table is normally transposed with\n";
  cout << "      the label value in the first column. With this option, the table is always such that\n";
  cout << "      each row corresponds to one input source image, also when only one is given. (default: off)\n";
  cout << "  -noid\n";
  cout << "      Do not print ID column when :option:`-table` is used.\n";
  cout << "  -[no]id-path [on|off]\n";
  cout << "      Use full input source image file path as ID column entry. (default: off)\n";
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// Implemented segmentation overlap metrics
// (https://en.wikipedia.org/wiki/Confusion_matrix)
enum OverlapMetric {
  UnknownMetric,
  TruePositives,
  TrueNegatives,
  FalsePositives,
  FalseNegatives,
  Sensitivity,
  Specificity,
  PositivePredictiveValue,
  NegativePredictiveValue,
  FalsePositiveRate,
  FalseDiscoveryRate,
  FalseNegativeRate,
  Accuracy,
  F1Score,  // same as Dice
  MatthewsCorrelation,
  Informedness,
  Markedness,
  Dice,
  Jaccard,
  // Aliases
  Recall = Sensitivity,
  HitRate = Sensitivity,
  TruePositiveRate = Sensitivity,
  TrueNegativeRate = Specificity,
  Precision = PositivePredictiveValue,
  FallOut = FalsePositiveRate,
  MissRate = FalseNegativeRate
};

// -----------------------------------------------------------------------------
istream &operator >>(istream &is, OverlapMetric &metric)
{
  string str;
  is >> str;
  str = ToLower(str);
  if (str == "true positives" || str == "truepositives" || str == "tp") {
    metric = TruePositives;
  } else if (str == "true negatives" || str == "truenegatives" || str == "tn") {
    metric = TrueNegatives;
  } else if (str == "false negatives" || str == "falsenegatives" || str == "fn") {
    metric = FalseNegatives;
  } else if (str == "false positives" || str == "falsepositives" || str == "fp") {
    metric = FalsePositives;
  } else if (str == "sensitivity") {
    metric = Sensitivity;
  } else if (str == "specificity") {
    metric = Specificity;
  } else if (str == "ppv" || str == "positivepredictivevalue" || str == "positive predictive value" || str == "precision") {
    metric = PositivePredictiveValue;
  } else if (str == "npv" || str == "negativepredictivevalue" || str == "negative predictive value") {
    metric = NegativePredictiveValue;
  } else if (str == "fpr" || str == "falsepositiverate" || str == "false positive rate" || str == "fallout" || str == "fall-out" || str == "fall out") {
    metric = FalsePositiveRate;
  } else if (str == "fdr" || str == "falsediscoveryrate" || str == "false discovery rate") {
    metric = FalseDiscoveryRate;
  } else if (str == "fnr" || str == "falsenegativerate" || str == "false negative rate" || str == "missrate" || str == "miss rate") {
    metric = FalseNegativeRate;
  } else if (str == "accuracy") {
    metric = Accuracy;
  } else if (str == "f1score" || str == "fscore" || str == "f-score" || str == "fmeasure" || str == "f-measure") {
    metric = F1Score;
  } else if (str == "matthewscorrelationcoefficient" || str == "matthews correlation coefficient" ||
             str == "matthewscorrelation" || str == "matthews correlation" || str == "mcc") {
    metric = MatthewsCorrelation;
  } else if (str == "bm" || str == "informedness" || str == "bookmakerinformedness" || str == "bookmaker informedness") {
    metric = Informedness;
  } else if (str == "mk" || str == "markedness") {
    metric = Markedness;
  } else if (str == "dsc" || str == "dice coefficient" || str == "dicecoefficient" || str == "dice") {
    metric = Dice;
  } else if (str == "jsc" || str == "jaccard similarity" || str == "jacccardsimilarity" || str == "jaccard") {
    metric = Jaccard;
  } else {
    metric = UnknownMetric;
    is.setstate(ios::failbit);
  }
  return is;
}

// -----------------------------------------------------------------------------
ostream &operator <<(ostream &os, const OverlapMetric &metric)
{
  switch (metric) {
    case TruePositives:           os << "True positives"; break;
    case TrueNegatives:           os << "True negatives"; break;
    case FalsePositives:          os << "False positives"; break;
    case FalseNegatives:          os << "False negatives"; break;
    case Sensitivity:             os << "Sensitivity"; break;
    case Specificity:             os << "Specificity"; break;
    case PositivePredictiveValue: os << "Positive predictive value"; break;
    case NegativePredictiveValue: os << "Negative predictive value"; break;
    case FalsePositiveRate:       os << "False positive rate"; break;
    case FalseDiscoveryRate:      os << "False discovery rate"; break;
    case FalseNegativeRate:       os << "False negative rate"; break;
    case Accuracy:                os << "Accuracy"; break;
    case MatthewsCorrelation:     os << "Matthews correlation coefficient"; break;
    case Informedness:            os << "Bookmaker informedness"; break;
    case Markedness:              os << "Markedness"; break;
    case F1Score:                 os << "F1-score"; break;
    case Dice:                    os << "Dice similarity coefficient"; break;
    case Jaccard:                 os << "Jaccard similarity coefficient"; break;
    default:                      os << "Unknown metric"; break;
  }
  return os;
}

// -----------------------------------------------------------------------------
string Abbreviation(const OverlapMetric &metric)
{
  switch (metric) {
    case TruePositives:           return "TP";
    case TrueNegatives:           return "TN";
    case FalsePositives:          return "FP";
    case FalseNegatives:          return "FN";
    case Sensitivity:             return "TPR";
    case Specificity:             return "TNR";
    case PositivePredictiveValue: return "PPV";
    case NegativePredictiveValue: return "NPV";
    case FalsePositiveRate:       return "FPR";
    case FalseDiscoveryRate:      return "FDR";
    case FalseNegativeRate:       return "FNR";
    case Accuracy:                return "ACC";
    case MatthewsCorrelation:     return "MCC";
    case Informedness:            return "BM";
    case Markedness:              return "MK";
    case Dice:                    return "DSC";
    case Jaccard:                 return "JSC";
    default:                      return "Unknown";
  }
}

// -----------------------------------------------------------------------------
double Evaluate(OverlapMetric metric, int tp, int fn, int fp, int tn)
{
  switch (metric) {
    case TruePositives:
      return double(tp);
    case TrueNegatives:
      return double(tn);
    case FalsePositives:
      return double(fp);
    case FalseNegatives:
      return double(fn);
    case Sensitivity:
      return double(tp) / double(tp + fn);
    case Specificity:
      return double(tn) / double(fp + tn);
    case PositivePredictiveValue:
      return double(tp) / double(tp + fp);
    case NegativePredictiveValue:
      return double(tn) / double(tn + fn);
    case FalsePositiveRate:
      return double(fp) / double(fp + tn);
    case FalseDiscoveryRate:
      return double(fp) / double(fp + tp);
    case FalseNegativeRate:
      return double(fn) / double(tp + fn);
    case Accuracy:
      return double(tp + tn) / double(tp + fn + fp + tn);
    case Informedness:
      return double(tp) / double(tp + fn) + double(tn) / double(fp + tn) - 1.;
    case MatthewsCorrelation: {
      double denom = sqrt(double(tp + fp) * double(tp + fn) * double(tn + fp) * double(tn + fn));
      return (double(tp) * double(tn) - double(fp) * double(fn)) / denom;
    } break;
    case Markedness:
      return double(tp) / double(tp + fp) + double(tn) / double(tn + fn) - 1.;
    case F1Score: case Dice:
      return double(2 * tp) / double((fp + tp) + (tp + fn));
    case Jaccard:
      return double(tp) / double(fp + tp + fn);
    default:
      FatalError("Unknown overlap metric: " << metric);
  }
}

// -----------------------------------------------------------------------------
double Evaluate(OverlapMetric metric, const BaseImage *target, const BaseImage *source, double min_value = NaN)
{
  const int num = target->NumberOfVoxels();
  if (source->NumberOfVoxels() != num) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Both images must have the same number of voxels");
  }
  if (IsNaN(min_value) && (metric == Dice || metric == Jaccard)) {
    double t2 = 0., s2 = 0., ts = 0.;
    for (int idx = 0; idx < num; ++idx) {
      t2 += target->GetAsDouble(idx) * target->GetAsDouble(idx);
      s2 += source->GetAsDouble(idx) * source->GetAsDouble(idx);
      ts += target->GetAsDouble(idx) * source->GetAsDouble(idx);
    }
    if (metric == Jaccard) {
      return ts / (t2 + s2 - ts);
    }
    return 2.0 * ts / (t2 + s2);
  } else {
    if (IsNaN(min_value)) min_value = 0.5;
    int tp = 0, tn = 0, fp = 0, fn = 0;
    for (int idx = 0; idx < num; ++idx) {
      if (target->GetAsDouble(idx) >= min_value) {
        if (source->GetAsDouble(idx) >= min_value) ++tp;
        else                                       ++fn;
      } else {
        if (source->GetAsDouble(idx) >= min_value) ++fp;
        else                                       ++tn;
      }
    }
    return Evaluate(metric, tp, fn, fp, tn);
  }
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
  auto columns = Split(Trim(line), delim, 0, discard_empty, handle_quotes);
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
  const char *table_name  = nullptr;
  const char *image_list  = nullptr;
  const char *image_delim = nullptr;
  const char *idcol_name  = nullptr;
  bool        idcol_path  = false;
  GreyPixel   padding     = MIN_GREY;
  bool        header      = true;
  const char *delim       = ",";
  int         digits      = 2;
  bool        all_labels  = false;
  bool        pbmaps      = false;
  double      pbmap_min   = NaN;
  Array<OrderedSet<GreyPixel> > segments;
  Array<OverlapMetric>          metrics;

  verbose = 1; // by default, include textual description of output value
               // can be disabled by caller using "-v 0" option in order to
               // only append the overlap value (excl. newline) to a file
               // or even better, use the -table output option

  int i1 =  0, j1 =  0, k1 =  0;
  int i2 = -1, j2 = -1, k2 = -1;

  for (ARGUMENTS_AFTER(nposarg)) {
    if (OPTION("-Tp")) PARSE_ARGUMENT(padding);
    else if (OPTION("-label") || OPTION("-segment")) {
      GreyPixel a, b;
      OrderedSet<GreyPixel> segment;
      do {
        bool invalid_a_or_b = false;
        const char *arg = ARGUMENT;
        const Array<string> parts = Split(arg, "..");
        if (parts.size() == 1) {
          if (FromString(parts[0], a)) {
            b = a;
          } else {
            invalid_a_or_b = true;
          }
        } else if (parts.size() == 2) {
          if (!FromString(parts[0], a) || !FromString(parts[1], b)) {
            invalid_a_or_b = true;
          }
        } else {
          invalid_a_or_b = true;
        }
        if (invalid_a_or_b) {
          FatalError("Invalid -label, -segment argument: " << arg);
        }
        if (a > b) swap(a, b);
        for (GreyPixel label = a; label <= b; ++label) {
          segment.insert(label);
        }
      } while (HAS_ARGUMENT);
      if (segment.size() == 0) {
        FatalError("Failed to parse -label, -segment option argument");
      }
      segments.push_back(segment);
    }
    else if (OPTION("-table")) {
      if (HAS_ARGUMENT) table_name = ARGUMENT;
      else table_name = "cout";
    }
    else if (OPTION("-images")) image_list = ARGUMENT;
    else if (OPTION("-images-delimiter") || OPTION("-images-delim")) image_delim = ARGUMENT;
    else if (OPTION("-labels")) all_labels = true;
    else if (OPTION("-metric")) {
      OverlapMetric metric;
      do {
        const char *arg = ARGUMENT;
        if (!FromString(arg, metric) || metric == UnknownMetric) {
          FatalError("Unknown overlap metric: " << arg);
        }
        metrics.push_back(metric);
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-precision")) PARSE_ARGUMENT(digits);
    else if (OPTION("-delim")) delim = ARGUMENT;
    else if (OPTION("-Rx1")) PARSE_ARGUMENT(i1);
    else if (OPTION("-Rx2")) PARSE_ARGUMENT(i2);
    else if (OPTION("-Ry1")) PARSE_ARGUMENT(j1);
    else if (OPTION("-Ry2")) PARSE_ARGUMENT(j2);
    else if (OPTION("-Rz1")) PARSE_ARGUMENT(k1);
    else if (OPTION("-Rz2")) PARSE_ARGUMENT(k2);
    else if (OPTION("-probs")) {
      pbmaps = true;
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(pbmap_min);
      } else {
        pbmap_min = NaN;
      }
    }
    else if (OPTION("-id")) {
      if (HAS_ARGUMENT) {
        idcol_name = ARGUMENT;
        bool flag;
        if (FromString(idcol_name, flag)) {
          idcol_name = (flag ? "ID" : nullptr);
        }
      } else {
        idcol_name = "ID";
      }
    }
    else if (OPTION("-noid")) {
      idcol_name = nullptr;
    }
    else HANDLE_BOOLEAN_OPTION("id-path", idcol_path);
    else HANDLE_BOOL_OPTION(header);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }
  if (pbmaps && (all_labels || !segments.empty())) {
    FatalError("Options -label[s] and -pbmaps are mutually exclusive");
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

  // Initialize image reader/writer factories
  InitializeIOLibrary();

  // Read target image and extract region of interest
  UniquePtr<BaseImage> target(BaseImage::New(target_name));
  const ImageAttributes target_attributes = target->Attributes();

  // Determine set of hard segmentation labels
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
    for (int i = 0; i < target->NumberOfVoxels(); ++i) {
      GreyPixel label = voxel_cast<GreyPixel>(target->GetAsDouble(i));
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
  if (i2 < 0) i2 = target->X();
  if (j2 < 0) j2 = target->Y();
  if (k2 < 0) k2 = target->Z();
  if ((i1 != 0) || (i2 != target->X()) ||
      (j1 != 0) || (j2 != target->Y()) ||
      (k1 != 0) || (k2 != target->Z())) {
    BaseImage *region = nullptr;
    target->GetRegion(region, i1, j1, k1, i2, j2, k2);
    target.reset(region);
  }

  // Number of voxels
  const int num = target->NumberOfVoxels();

  // -------------------------------------------------------------------------------
  // Intensity image overlap
  if (segments.empty() && !pbmaps) {

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
      if ((i1 != 0) || (i2 != source.X()) ||
          (j1 != 0) || (j2 != source.Y()) ||
          (k1 != 0) || (k2 != source.Z())) {
        source = source.GetRegion(i1, j1, k1, i2, j2, k2);
      }

      // Default image overlap metric
      if (!metrics.empty()) {
        FatalError("Segmentation overlap metrics not suitable for intensity images, use default SI!");
      }

      // Determine maximum intensity
      int max = 0;
      for (int idx = 0; idx < num; ++idx) {
        auto tgt = voxel_cast<GreyPixel>(target->GetAsDouble(idx));
        if (tgt > padding) {
          auto src = source(idx);
          if (tgt > max) max = tgt;
          if (src > max) max = src;
        }
      }

      // Compute image overlap
      Histogram1D<int> histogramA (max);
      Histogram1D<int> histogramB (max);
      Histogram1D<int> histogramAB(max);

      for (int idx = 0; idx < num; ++idx) {
        auto tgt = voxel_cast<GreyPixel>(target->GetAsDouble(idx));
        if (tgt > padding) {
          auto src = source(idx);
          histogramA.Add(tgt);
          histogramB.Add(src);
          if (tgt == src) histogramAB.Add(tgt);
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

    // Default segmentation overlap metric
    if (metrics.empty()) {
      if (source_name.size() == 1) {
        metrics.push_back(TruePositives);
        metrics.push_back(FalseNegatives);
        metrics.push_back(FalsePositives);
        metrics.push_back(TrueNegatives);
        metrics.push_back(Sensitivity);
        metrics.push_back(Specificity);
        metrics.push_back(Dice);
        metrics.push_back(Jaccard);
      } else {
        metrics.push_back(Dice);
      }
    }

    if (table_name || pbmaps) {

      ostream *os = &cout;
      ofstream ofs;

      if (table_name && ToLower(table_name) != "stdout" && ToLower(table_name) != "cout") {
        ofs.open(table_name);
        if (!ofs) FatalError("Failed to open output file: " << table_name);
        os = &ofs;
      }

      bool transposed = (!idcol_name && source_name.size() == 1);

      if (table_name && header) {
        if (idcol_name) {
          *os << idcol_name << delim;
        }
        if (pbmaps) {
          for (size_t m = 0; m < metrics.size(); ++m) {
            if (m > 0) *os << delim;
            *os << Abbreviation(metrics[m]);
          }
        } else if (transposed) {
          *os << "Label";
          for (size_t m = 0; m < metrics.size(); ++m) {
            *os << delim << Abbreviation(metrics[m]);
          }
        } else {
          for (size_t roi = 0; roi < segments.size(); ++roi)
          for (size_t m = 0; m < metrics.size(); ++m) {
            if (roi > 0 || m > 0) {
              *os << delim;
            }
            if (metrics.size() > 1) {
              *os << Abbreviation(metrics[m]);
              if (strcmp(delim, " ") != 0) {
                *os << ' ';
              }
            }
            const auto &labels = segments[roi];
            for (auto it = labels.begin(); it != labels.end(); ++it) {
              if (it != labels.begin()) *os << "+";
              *os << *it;
            }
          }
        }
        *os << "\n";
      }

      for (size_t n = 0; n < source_name.size(); ++n) {

        // Print source image ID
        if (idcol_name) {
          if (idcol_path) {
            *os << source_name[n] << delim;
          } else {
            *os << FileName(source_name[n]) << delim;
          }
        }

        // Read source image and extract region of interest
        UniquePtr<BaseImage> source(BaseImage::New(source_name[n].c_str()));
        if (source->Attributes() != target_attributes) {
          FatalError("Input images must be sampled on same lattice!");
        }
        if ((i1 != 0) || (i2 != source->X()) ||
            (j1 != 0) || (j2 != source->Y()) ||
            (k1 != 0) || (k2 != source->Z())) {
          BaseImage *region = nullptr;
          source->GetRegion(region, i1, j1, k1, i2, j2, k2);
          source.reset(region);
        }

        if (pbmaps) {
          for (size_t m = 0; m < metrics.size(); ++m) {
            if (m > 0) {
              *os << delim;
            }
            *os << setprecision(digits) << Evaluate(metrics[m], target.get(), source.get(), pbmap_min);
          }
        } else {

          // Iterate over segments
          for (size_t roi = 0; roi < segments.size(); ++roi) {
            // Print segment ID
            const auto &labels = segments[roi];
            if (transposed) {
              for (auto it = labels.begin(); it != labels.end(); ++it) {
                if (it != labels.begin()) *os << "+";
                *os << *it;
              }
            }
            // Determine TP, FP, TN, FN
            int tp = 0, fp = 0, tn = 0, fn = 0;
            if (labels.size() == 0 && (*labels.begin()) == 0) {
              for (int idx = 0; idx < num; ++idx) {
                auto tgt = voxel_cast<GreyPixel>(target->GetAsDouble(idx));
                auto src = voxel_cast<GreyPixel>(source->GetAsDouble(idx));
                tgt = (tgt > 0 ? 1 : 0);
                src = (src > 0 ? 1 : 0);
                if      (tgt == 1 && src == 1) ++tp;
                else if (tgt == 1 && src == 0) ++fn;
                else if (tgt == 0 && src == 1) ++fp;
                else                           ++tn;
              }
            } else {
              for (int idx = 0; idx < num; ++idx) {
                auto tgt = voxel_cast<GreyPixel>(target->GetAsDouble(idx));
                auto src = voxel_cast<GreyPixel>(source->GetAsDouble(idx));
                if (labels.find(tgt) != labels.end()) {
                  if (labels.find(src) != labels.end()) ++tp;
                  else                                  ++fn;
                }
                else if (labels.find(src) != labels.end()) ++fp;
                else                                       ++tn;
              }
            }
            // Compute overlap metrics
            for (size_t m = 0; m < metrics.size(); ++m) {
              if (transposed || roi > 0 || m > 0) {
                *os << delim;
              }
              switch (metrics[m]) {
                case TruePositives:
                  *os << tp;
                  break;
                case TrueNegatives:
                  *os << tn;
                  break;
                case FalsePositives:
                  *os << fp;
                  break;
                case FalseNegatives:
                  *os << fn;
                  break;
                default:
                  *os << setprecision(digits) << Evaluate(metrics[m], tp, fn, fp, tn);
              }
            }

            if (transposed) *os << "\n";
          }
        }
        if (!transposed) *os << "\n";
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
          } else {
            if (n > 0) cout << delim;
          }

          // Read source image and extract region of interest
          if (n < sources.size()) {
            source = sources[n];
          } else {
            source.Read(source_name[n].c_str());
            if (source.Attributes() != target_attributes) {
              FatalError("Input images must be sampled on same lattice!");
            }
            if ((i1 != 0) || (i2 != source.X()) ||
                (j1 != 0) || (j2 != source.Y()) ||
                (k1 != 0) || (k2 != source.Z())) {
              source = source.GetRegion(i1, j1, k1, i2, j2, k2);
            }
          }

          // Determine TP, FP, TN, FN
          int tp = 0, fp = 0, tn = 0, fn = 0;
          if (labels.size() == 0 && (*labels.begin()) == 0) {
            for (int idx = 0; idx < num; ++idx) {
              auto tgt = voxel_cast<GreyPixel>(target->GetAsDouble(idx));
              auto src = source(idx);
              tgt = (tgt > 0 ? 1 : 0);
              src = (src > 0 ? 1 : 0);
              if      (tgt == 1 && src == 1) ++tp;
              else if (tgt == 1 && src == 0) ++fn;
              else if (tgt == 0 && src == 1) ++fp;
              else                           ++tn;
            }
          } else {
            for (int idx = 0; idx < num; ++idx) {
              auto tgt = voxel_cast<GreyPixel>(target->GetAsDouble(idx));
              if (labels.find(tgt) != labels.end()) {
                if (labels.find(source(idx)) != labels.end()) ++tp;
                else                                          ++fn;
              }
              else if (labels.find(source(idx)) != labels.end()) ++fp;
              else                                               ++tn;
            }
          }

          // Compute overlap metrics
          for (size_t m = 0; m < metrics.size(); ++m) {
            double overlap = Evaluate(metrics[m], tp, fn, fp, tn);
            if (verbose) {
              if (m > 0) {
                cout << ", ";
              }
              if (verbose == 1) {
                cout << Abbreviation(metrics[m]);
              } else {
                cout << metrics[m];
              }
              cout << " = ";
              if (metrics[m] == TruePositives || metrics[m] == FalsePositives ||
                  metrics[m] == TrueNegatives || metrics[m] == FalseNegatives) {
                cout << int(overlap);
              } else {
                const int d = (fequal(overlap, 1., 1e-6) ? 3 : digits);
                cout << setprecision(d) << (100. * overlap) << "%";
              }
            } else {
              if (m > 0) {
                cout << delim;
              }
              if (metrics[m] == TruePositives || metrics[m] == FalsePositives ||
                  metrics[m] == TrueNegatives || metrics[m] == FalseNegatives) {
                cout << int(overlap);
              } else {
                cout << setprecision(digits) << overlap;
              }
            }
          }
          if (verbose) cout << "\n";
        }
        if (verbose == 0) cout << "\n";
      }

    }
  }

  return 0;
}

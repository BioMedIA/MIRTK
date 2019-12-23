/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2019 Andreas Schuh
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
  cout << "Usage: " << name << " [<image>...] [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Computes the overlap of one or more segmentation labels. If more than\n";
  cout << "  one image is given, either the overlap between each of these and a given\n";
  cout << "  :option:`-target` image or the overlap between all unordered pairs of\n";
  cout << "  images is evaluated.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  image   Hard segmentation with class labels or\n";
  cout << "          soft segmentation with class probabilities.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -target <file>\n";
  cout << "       Target image with respect to which overlap is evaluated.\n";
  cout << "       When no target image is given, the overlap between all unique\n";
  cout << "       unordered pairs of images is evaluated.\n";
  cout << "  -images <file>\n";
  cout << "      Text file with input image file paths.\n";
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
  cout << "  -probs, -pbmaps [<threshold>]\n";
  cout << "      Evaluate overlap of class probability maps (fuzzy membership functions).\n";
  cout << "      The values of the input images are rescaled to [0, 1]. When no <threshold>\n";
  cout << "      is specified, a value of 0.5 is used to binarize the fuzzy segmentation masks.\n";
  cout << "      When no <threshold> is specified and the metric is either Dice (DSC) or Jaccard (JSC),\n";
  cout << "      these metrics are evaluated using their respective definition based on scalar products.\n";
  cout << "  -precision <int>\n";
  cout << "      Number of significant digits. (default: 2)\n";
  cout << "  -metric <name>...\n";
  cout << "      Segmentation overlap metrics (see https://en.wikipedia.org/wiki/Confusion_matrix).\n";
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
  cout << "      Do not print ID column.\n";
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
  } else if (str == "tpr" || str == "sensitivity") {
    metric = Sensitivity;
  } else if (str == "tnr" || str == "specificity") {
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
  } else if (str == "acc" || str == "accuracy") {
    metric = Accuracy;
  } else if (str == "f1score" || str == "fscore" || str == "f-score" || str == "fmeasure" || str == "f-measure") {
    metric = F1Score;
  } else if (str == "mcc" || str == "matthewscorrelationcoefficient" || str == "matthews correlation coefficient" ||
             str == "matthewscorrelation" || str == "matthews correlation") {
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
struct Arguments
{
  ImageAttributes domain;       ///< Common image attributes.
  GreyPixel       padding;      ///< Target padding value.
  int             i1, i2;
  int             j1, j2;
  int             k1, k2;
  int             digits;       ///< Number of decimal digits.
  const char     *delim;        ///< Column delimiter.
  bool            table;        ///< Whether to print results in table format.
  bool            header;       ///< Whether to include table header.
  bool            idcol_flag;   ///< Whether to include index column(s).
  bool            idcol_path;   ///< Whether to use full file path as image ID.
  bool            all_labels;   ///< Whether to evaluate overlap for all labels.
  bool            pbmaps;       ///< Whether input is class probabilities.
  double          pbmap_min;    ///< Class probability threshold.

  Array<OrderedSet<GreyPixel> > segments;   ///< Merged label segments.
  Array<OverlapMetric>          metrics;    ///< Overlap metrics.

  Arguments()
  :
    padding(MIN_GREY),
    i1(0), i2(-1),
    j1(0), j2(-1),
    k1(0), k2(-1),
    digits(2),
    delim(","),
    header(true),
    idcol_flag(true),
    idcol_path(false),
    all_labels(false),
    pbmaps(false),
    pbmap_min(NaN)
  {}
};

// -----------------------------------------------------------------------------
void AppendSingleLabelSegments(const BaseImage *image, Array<OrderedSet<GreyPixel> > &segments)
{
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
  for (int i = 0; i < image->NumberOfVoxels(); ++i) {
    GreyPixel label = voxel_cast<GreyPixel>(image->GetAsDouble(i));
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

// -----------------------------------------------------------------------------
UniquePtr<BaseImage> ReadImage(const char *name, Arguments &args)
{
  // Read image from file
  UniquePtr<BaseImage> image(BaseImage::New(name));
  // Store target / initial image attributes
  if (!args.domain) {
    args.domain = image->Attributes();
    if (args.all_labels) {
      AppendSingleLabelSegments(image.get(), args.segments);
    }
  } else {
    if (image->Attributes() != args.domain) {
      FatalError("Input images must be sampled on the same lattice!");
    }
  }
  if (args.i2 < 0) args.i2 = image->X();
  if (args.j2 < 0) args.j2 = image->Y();
  if (args.k2 < 0) args.k2 = image->Z();
  // Extract region of interest
  if ((args.i1 != 0) || (args.i2 != image->X()) ||
      (args.j1 != 0) || (args.j2 != image->Y()) ||
      (args.k1 != 0) || (args.k2 != image->Z())) {
    BaseImage *region = nullptr;
    image->GetRegion(region, args.i1, args.j1, args.k1, args.i2, args.j2, args.k2);
    image.reset(region);
  }
  return image;
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
    // Trim leading and trailing whitespace characters
    line = Trim(line);
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

// -----------------------------------------------------------------------------
// Print header row
void PrintHeader(ostream *os, const char *target_name, const Arguments &args)
{
  if (args.table && args.header) {
    if (args.idcol_flag) {
      if (!target_name) {
        *os << "target" << args.delim;
      }
      *os << "source" << args.delim;
    }
    if (!args.pbmaps) {
      *os << "label";
    }
    for (size_t m = 0; m < args.metrics.size(); ++m) {
      *os << args.delim << ToLower(Abbreviation(args.metrics[m]));
    }
    *os << "\n";
  }
}

// -----------------------------------------------------------------------------
// Print index column(s)
void PrintIndex(ostream *os,
                const char *target_name,
                const char *source_name,
                const Arguments &args)
{
  string target_id, source_id;
  if (args.idcol_flag) {
    if (args.idcol_path) {
      if (target_name) {
        target_id = target_name;
      }
      source_id = source_name;
    } else {
      if (target_name) {
        target_id = FileName(target_name);
      }
      source_id = FileName(source_name);
    }
  }
  if (args.table) {
    if (!target_id.empty()) *os << target_id << args.delim;
    if (!source_id.empty()) *os << source_id << args.delim;
  } else {
    if (!target_id.empty()) *os << "Target = " << target_id << args.delim;
    if (!source_id.empty()) *os << "Source = " << source_id << args.delim;
  }
}

// -----------------------------------------------------------------------------
// Evaluate overlap between pair of images
void AppendOverlap(ostream *os,
                   const BaseImage *target,
                   const BaseImage *source,
                   const char *target_name,
                   const char *source_name,
                   const Arguments &args)
{
  if (args.pbmaps) {

    PrintIndex(os, target_name, source_name, args);
    for (size_t m = 0; m < args.metrics.size(); ++m) {
      if (m > 0) {
        *os << args.delim;
      }
      if (!args.table) *os << args.metrics[m] << " = ";
      *os << setprecision(args.digits) << Evaluate(args.metrics[m], target, source, args.pbmap_min);
    }
    *os << "\n";

  } else {

    const int num = target->NumberOfVoxels();

    // Iterate over segments
    for (size_t roi = 0; roi < args.segments.size(); ++roi) {
      // Print index column(s)
      PrintIndex(os, target_name, source_name, args);
      // Print segment ID
      if (!args.table) *os << "Label = ";
      const auto &labels = args.segments[roi];
      for (auto it = labels.begin(); it != labels.end(); ++it) {
        if (it != labels.begin()) *os << "+";
        *os << *it;
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
      for (size_t m = 0; m < args.metrics.size(); ++m) {
        *os << args.delim;
        if (!args.table) *os << args.metrics[m] << " = ";
        switch (args.metrics[m]) {
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
            *os << setprecision(args.digits) << Evaluate(args.metrics[m], tp, fn, fp, tn);
        }
      }
      *os << "\n";
    }
  }
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  // Check command line
  REQUIRES_POSARGS(0);

  // Parse positional arguments
  int nposarg;
  Array<string> image_names;
  for (nposarg = 1; nposarg < argc; ++nposarg) {
    if (argv[nposarg][0] == '-') break;
    image_names.push_back(argv[nposarg]);
  }
  --nposarg;

  // Parse optional arguments
  Arguments args;
  const char *target_name  = nullptr;
  const char *table_name   = nullptr;
  const char *image_list   = nullptr;
  const char *image_delim  = nullptr;

  verbose = 1; // by default, include textual description of output value
               // can be disabled by caller using "-v 0" option in order to
               // only append the overlap value (excl. newline) to a file
               // or even better, use the -table output option

  for (ARGUMENTS_AFTER(nposarg)) {
    if (OPTION("-Tp")) PARSE_ARGUMENT(args.padding);
    else if (OPTION("-labels")) args.all_labels = true;
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
      args.segments.push_back(segment);
    }
    else if (OPTION("-table")) {
      if (HAS_ARGUMENT) {
        table_name = ARGUMENT;
        auto larg = ToLower(table_name);
        if (larg == "cout" || larg == "stdout") {
          table_name = "cout";
        }
      } else {
        table_name = "cout";
      }
      args.table = true;
    }
    else if (OPTION("-target")) {
      target_name = ARGUMENT;
    }
    else if (OPTION("-images")) {
      image_list = ARGUMENT;
    }
    else if (OPTION("-images-delimiter") || OPTION("-images-delim")) {
      image_delim = ARGUMENT;
    }
    else if (OPTION("-metric")) {
      OverlapMetric metric;
      do {
        const char *arg = ARGUMENT;
        if (!FromString(arg, metric) || metric == UnknownMetric) {
          FatalError("Unknown overlap metric: " << arg);
        }
        args.metrics.push_back(metric);
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-precision")) PARSE_ARGUMENT(args.digits);
    else if (OPTION("-delim")) args.delim = ARGUMENT;
    else if (OPTION("-Rx1")) PARSE_ARGUMENT(args.i1);
    else if (OPTION("-Rx2")) PARSE_ARGUMENT(args.i2);
    else if (OPTION("-Ry1")) PARSE_ARGUMENT(args.j1);
    else if (OPTION("-Ry2")) PARSE_ARGUMENT(args.j2);
    else if (OPTION("-Rz1")) PARSE_ARGUMENT(args.k1);
    else if (OPTION("-Rz2")) PARSE_ARGUMENT(args.k2);
    else if (OPTION("-probs") || OPTION("-pbmaps")) {
      args.pbmaps = true;
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(args.pbmap_min);
      } else {
        args.pbmap_min = NaN;
      }
    }
    else if (OPTION("-index") || OPTION("-id")) {
      args.idcol_flag = true;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(args.idcol_flag);
    }
    else if (OPTION("-noindex") || OPTION("-noid")) {
      args.idcol_flag = false;
    }
    else HANDLE_BOOLEAN_OPTION("id-path", args.idcol_path);
    else HANDLE_BOOLEAN_OPTION("index-by-path", args.idcol_path);
    else HANDLE_BOOL_OPTION(args.header);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }
  if (args.pbmaps && (args.all_labels || !args.segments.empty())) {
    FatalError("Options -label[s] and -probs are mutually exclusive");
  }

  // By default, assume hard segmentations and evalute overlap of all non-zero labels
  if (args.segments.empty() && !args.pbmaps) {
    args.all_labels = true;
  }

  // Read source image file paths from text file
  if (image_list) {
    if (ReadImageListFile(image_list, image_names, image_delim) == 0) {
      Warning("No source image file paths read from file: " << image_list);
    }
  }
  const int required_images = (target_name ? 1 : 2);
  if (image_names.size() < required_images) {
    FatalError("At least " << required_images << " image"
               << (required_images > 1 ? "s" : "") << " must be given!");
  }

  // Initialize image reader/writer factories
  InitializeIOLibrary();

  // Default segmentation overlap metric
  if (args.metrics.empty()) {
    args.metrics.push_back(TruePositives);
    args.metrics.push_back(FalseNegatives);
    args.metrics.push_back(FalsePositives);
    args.metrics.push_back(TrueNegatives);
    args.metrics.push_back(FallOut);
    args.metrics.push_back(Sensitivity);
    args.metrics.push_back(Specificity);
    args.metrics.push_back(Dice);
    args.metrics.push_back(Jaccard);
  }

  // Open output stream
  ostream *os = &cout;

  ofstream ofs;
  if (args.table && strcmp(table_name, "cout") != 0) {
    ofs.open(table_name);
    if (!ofs) FatalError("Failed to open output file: " << table_name);
    os = &ofs;
  }

  // Print evaluated overlap
  PrintHeader(os, target_name, args);

  if (target_name) {

    auto target = ReadImage(target_name, args);
    for (size_t n = 0; n < image_names.size(); ++n) {
      const char *source_name = image_names[n].c_str();
      if (os != &cout) {
        cout << "Evaluating target overlap with " << FileName(source_name) << endl;
      }
      auto source = ReadImage(source_name, args);
      AppendOverlap(os, target.get(), source.get(), target_name, source_name, args);
    }

  } else {

    const char *source_name;
    for (size_t i = 0; i < image_names.size(); ++i) {
      target_name = image_names[i].c_str();
      auto target = ReadImage(target_name, args);
      for (size_t j = i + 1; j < image_names.size(); ++j) {
        source_name = image_names[j].c_str();
        if (os != &cout) {
          cout << "Evaluating overlap between " << FileName(target_name) << " and " << FileName(source_name) << endl;
        }
        auto source = ReadImage(source_name, args);
        AppendOverlap(os, target.get(), source.get(), target_name, source_name, args);
      }
    }

  }

  return 0;
}

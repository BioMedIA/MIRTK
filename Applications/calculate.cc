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

#include <mirtkImageConfig.h>
#include <mirtkIOConfig.h>

#include <mirtkDataOp.h>
#include <mirtkDataStatistics.h>
#include <mirtkDataFunctions.h>

#if MIRTK_Image_WITH_VTK
#  include <vtkDataSet.h>
#  include <vtkSmartPointer.h>
#endif

using namespace mirtk;
using namespace mirtk::data;
using namespace mirtk::data::op;
using namespace mirtk::data::statistic;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <input> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  This tool can be used for basic calculations from a sequence of data values read" << endl;
  cout << "  either from an image or a VTK pointset. It can be used, for example, to add two" << endl;
  cout << "  data sequences and to divide the result by a constant. The current sequence can" << endl;
  cout << "  be written to an output file again using the -out option. Additionally, statistics" << endl;
  cout << "  of the current data sequence can be computed such as the mean or variance." << endl;
  cout << "  The order of the data transformations and calculation of statistics is determined" << endl;
  cout << "  by the order of the command-line arguments."  << endl;
  cout << endl;
  cout << "  By default, data statistics are printed to STDOUT in a human readable format." << endl;
  cout << "  This output can be appended to a text file using the -append option instead." << endl;
  cout << "  For a more machine readable output, e.g., as comma separated values (CSV)," << endl;
  cout << "  specify a delimiting string using the -d option. In this case, a header line" << endl;
  cout << "  is also printed when the -header option is given with optional user specified" << endl;
  cout << "  column names for the individual output values." << endl;
  cout << endl;
  cout << "Input options:" << endl;
  cout << "  -scalars <name>   Name of input point data array. (default: active SCALARS array)" << endl;
  cout << endl;
  cout << "Data manipulation options:" << endl;
  cout << "  -mask <float|file>                Exclude values equal a given threshold or with zero input mask value." << endl;
  cout << "                                    Note that this does not modify the data values, but only marks them to be" << endl;
  cout << "                                    ignored from now on. Use -outside/-pad following this operation to set an outside value." << endl;
  cout << "  -mask-inside  <lower> <upper>     Exclude values which are inside a given closed interval." << endl;
  cout << "  -mask-outside <lower> <upper>     Exclude values which are outside a given open interval." << endl;
  cout << "  -mask-below <value>               Exclude values less or equal a given threshold." << endl;
  cout << "  -mask-above <value>               Exclude values greator or equal a given threshold." << endl;
  cout << "  -even                             Exclude values which are not an even number when cast to an integer." << endl;
  cout << "  -odd                              Exclude values which are not an odd number when cast to an integer." << endl;
  cout << "  -threshold <lower> [<upper>]      Clamp values which are less or equal and optionally greater or equal a threshold." << endl;
  cout << "  -threshold-below <value>          Clamp values less or equal a given threshold." << endl;
  cout << "  -threshold-above <value>          Clamp values greator or equal a given threshold." << endl;
  cout << "  -rescale <min> <max>              Linearly rescale values to the interval [min, max]." << endl;
  cout << "  -inside <float>                   Set new value for all currently unmasked data values." << endl;
  cout << "  -outside, -pad <float>            Set new value for all currently masked data values." << endl;
  cout << "  -out <file> [<type>]              Write current data sequence to file in the format of the input file." << endl;
  cout << "                                    Output data type can be: uchar, short, ushort, int, uint, float, double." << endl;
  cout << "  -add <float|file>                 Add constant value or data sequence read from specified file." << endl;
  cout << "  -sub <float|file>                 Subtract constant value or data sequence read from specified file." << endl;
  cout << "  -mul <float|file>                 Multiply by constant value or data sequence read from specified file." << endl;
  cout << "  -div <float|file>                 Divide by constant value or data sequence read from specified file." << endl;
  cout << "                                    When dividing by zero values in the input file, the result is NaN." << endl;
  cout << "                                    Use :option:`-mask` with argument NaN and :option:`-outside` to replace" << endl;
  cout << "                                    these undefined values by a constant such as zero." << endl;
  cout << "  -div-with-zero <file>             Same as :option:`-div` but set result to zero instead of NaN in case of" << endl;
  cout << "                                    a division by zero due to a zero value in the specified file." << endl;
  cout << "  -abs                              Take absolute values." << endl;
  cout << "  -pow <exponent>                   Raise values to the power of the given exponent." << endl;
  cout << "  -sq -square                       Raise values to the power of 2 (i.e, -pow 2)." << endl;
  cout << "  -sqrt                             Calculate square root of each value (i.e, -pow .5)." << endl;
  cout << "  -exp                              Calculate exponential of data sequence." << endl;
  cout << "  -log [<threshold>] [<base>]       Compute logarithm after applying an optional threshold." << endl;
  cout << "                                    (default threshold: min double, default base: e)" << endl;
  cout << "  -lb [<threshold>]                 Compute binary logarithm, alias for -log <threshold> 2." << endl;
  cout << "  -log2 [<threshold>]               Compute binary logarithm, alias for -log <threshold> 2." << endl;
  cout << "  -ln [<threshold>]                 Compute natural logarithm, alias for -log [<threshold>]." << endl;
  cout << "  -lg [<threshold>]                 Compute logarithm to base 10, alias for -log <threshold> 10." << endl;
  cout << "  -log10 [<threshold>]              Compute logarithm to base 10, alias for -log <threshold> 10." << endl;
  cout << endl;
  cout << "Data statistics options:" << endl;
  cout << "  -append <file>                    Append output to a file. (default: STDOUT)" << endl;
  cout << "  -delimiter | -delim | -d | -sep   Delimiting character(s). (default: '')" << endl;
  cout << "  -header [<name>...]               Request output of header line if delimiter was" << endl;
  cout << "                                    specified as well. If the output is appended to" << endl;
  cout << "                                    a text file, the header is only printed if it does not exist." << endl;
  cout << "                                    If no or fewer custom column names are given, the default" << endl;
  cout << "                                    names for each statistic are printed. (default: none)" << endl;
  cout << "  -prefix <str>...                  One or more prefix strings to print. If no delimiter" << endl;
  cout << "                                    is specified, the concatenated strings are printed before" << endl;
  cout << "                                    each line of the output. Otherwise, each prefix string" << endl;
  cout << "                                    is printed as entry for the first columns in the delimited" << endl;
  cout << "                                    output row, separated by the specified delimiter. (default: none)" << endl;
  cout << "  -precision <int>                  Number of significant digits. (default: 5)" << endl;
  cout << "  -mean | -avg                      Print mean intensity. (default: on)" << endl;
  cout << "  -variance | -var                  Print variance of intensity values. (default: off)" << endl;
  cout << "  -sigma | -std                     Print standard deviation of intensity values. (default: on)" << endl;
  cout << "  -minimum | -min                   Print minimum intensity value. (default: off)" << endl;
  cout << "  -maximum | -max                   Print maximum intensity value. (default: off)" << endl;
  cout << "  -extrema | -minmax                Print minimum and maximum intensity value. (default: on)" << endl;
  cout << "  -range                            Print range of intensity values (i.e., max - min). (default: off)" << endl;
  cout << "  -percentile <n> | -pct <n>        Print n-th percentile. (default: none)" << endl;
  PrintCommonOptions(cout);
  cout << endl;
  cout << "Examples:" << endl;
  cout << endl;
  cout << "  " << name << " mni305.nii.gz" << endl;
  cout << "      Mean = 26.9753" << endl;
  cout << "      Standard deviation = 50.3525" << endl;
  cout << "      Extrema = [0, 254]" << endl;
  cout << "      Range = 254" << endl;
  cout << endl;
  cout << "  " << name << " mni305.nii.gz -pct 77" << endl;
  cout << "      77th percentile = 25" << endl;
  cout << endl;
  cout << "  " << name << " mni305.nii.gz -padding 25 -range -percentile 25 50 75 -prefix MNI305 '[>25]'" << endl;
  cout << "      MNI305 [>25] range = 254" << endl;
  cout << "      MNI305 [>25] 25th percentile = 69" << endl;
  cout << "      MNI305 [>25] 50th percentile = 113" << endl;
  cout << "      MNI305 [>25] 75th percentile = 150" << endl;
  cout << endl;
  cout << "  " << name << " mni305.nii.gz -d , -prefix MNI305" << endl;
  cout << "      MNI305,26.9753,50.3525,0,254,254 [no newline at end of line]" << endl;
  cout << endl;
  cout << "  " << name << " mni305.nii.gz -d , -prefix MNI305 -header" << endl;
  cout << "      ,Mean,Sigma,Min,Max,Range" << endl;
  cout << "      MNI305,26.9753,50.3525,0,254,254" << endl;
  cout << endl;
  cout << "  " << name << " mni305.nii.gz -d , -prefix MNI305 -header ID Mean SD" << endl;
  cout << "      ID,Mean,SD,Min,Max,Range" << endl;
  cout << "      MNI305,26.9753,50.3525,0,254,254" << endl;
  cout << endl;
  cout << "  " << name << " a.nii.gz + b.nii.gz = c.nii.gz" << endl;
  cout << endl;
  cout << "  " << name << " a.vtk + b.nii.gz - 10 / c.nii = d.vtk" << endl;
  cout << "      Adds data values at identical sequential memory indices in a and b," << endl;
  cout << "      subtracts the constant 10, and then divides by the values in image c." << endl;
  cout << endl;
  cout << "      Note: Operations are always executed from left to right," << endl;
  cout << "            i.e., no mathematical operator precedence is considered!" << endl;
  cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
// Some special options do not start with a '-' as otherwise required
#undef  HAS_ARGUMENT
#define HAS_ARGUMENT                                                           \
  _IsArgument(ARGIDX, argc, argv) &&                                           \
  strcmp(argv[ARGIDX+1], "+") != 0 &&                                          \
  strcmp(argv[ARGIDX+1], "/") != 0 &&                                          \
  strcmp(argv[ARGIDX+1], "=") != 0

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  InitializeIOLibrary();

  // Initial data values
  REQUIRES_POSARGS(1);

  double *data = NULL;
  int datatype = MIRTK_VOXEL_DOUBLE;
  ImageAttributes attr;

#if MIRTK_Image_WITH_VTK
  const char *scalars_name = NULL;
  for (ARGUMENTS_AFTER(1)) {
    if (OPTION("-scalars")) scalars_name = ARGUMENT;
  }
  vtkSmartPointer<vtkDataSet> dataset;
  int n = Read(POSARG(1), data, &datatype, &attr, &dataset, scalars_name);
#else // MIRTK_Image_WITH_VTK
  int n = Read(POSARG(1), data, &datatype, &attr);
#endif // MIRTK_Image_WITH_VTK

  // Optional arguments
  const char *append_name   = NULL;
  const char *delimiter     = NULL;
  bool        print_header  = false;
  int         digits        = 5;

  Array<string> header;
  Array<string> prefix;

  Array<unique_ptr<Op> > ops;

  for (ARGUMENTS_AFTER(1)) {
    if (OPTION("-append")) append_name = ARGUMENT;
#if MIRTK_Image_WITH_VTK
    else if (OPTION("-scalars")) scalars_name = ARGUMENT;
#endif // MIRTK_Image_WITH_VTK
    else if (OPTION("-prefix")) {
      do {
        prefix.push_back(ARGUMENT);
      } while (HAS_ARGUMENT);
    } else if (OPTION("-header")) {
      print_header = true;
      while (HAS_ARGUMENT) header.push_back(ARGUMENT);
    } else if (OPTION("-mask")) {
      double c;
      do {
        const char *arg = ARGUMENT;
        if (FromString(arg, c)) ops.push_back(unique_ptr<Op>(new Mask(c)));
        else                    ops.push_back(unique_ptr<Op>(new Mask(arg)));
      } while (HAS_ARGUMENT);
    } else if (OPTION("-mask-outside")) {
      const char *arg = ARGUMENT;
      double lower, upper;
      if (!FromString(arg, lower)) {
        cerr << "Invalid -mask-outside lower argument, must be a number!" << endl;
      }
      arg = ARGUMENT;
      if (!FromString(arg, upper)) {
        cerr << "Invalid -mask-outside upper argument, must be a number!" << endl;
      }
      ops.push_back(unique_ptr<Op>(new MaskOutsideOpenInterval(lower, upper)));
    } else if (OPTION("-mask-inside")) {
      const char *arg = ARGUMENT;
      double lower, upper;
      if (!FromString(arg, lower)) {
        cerr << "Invalid -mask-inside lower argument, must be a number!" << endl;
      }
      arg = ARGUMENT;
      if (!FromString(arg, upper)) {
        cerr << "Invalid -mask-inside upper argument, must be a number!" << endl;
      }
      ops.push_back(unique_ptr<Op>(new MaskInsideInterval(lower, upper)));
    } else if (OPTION("-mask-below")) {
      const char *arg = ARGUMENT;
      double lower;
      if (!FromString(arg, lower)) {
        cerr << "Invalid -mask-below argument, must be a number!" << endl;
      }
      ops.push_back(unique_ptr<Op>(new MaskOutsideOpenInterval(lower, numeric_limits<double>::infinity())));
    } else if (OPTION("-mask-above")) {
      const char *arg = ARGUMENT;
      double upper;
      if (!FromString(arg, upper)) {
        cerr << "Invalid -mask-above argument, must be a number!" << endl;
      }
      ops.push_back(unique_ptr<Op>(new MaskOutsideOpenInterval(-numeric_limits<double>::infinity(), upper)));
    } else if (OPTION("-even")) {
      ops.push_back(unique_ptr<Op>(new MaskOddValues()));
    } else if (OPTION("-odd")) {
      ops.push_back(unique_ptr<Op>(new MaskEvenValues()));
    } else if (OPTION("-threshold")) {
      double lower, upper;
      if (!FromString(ARGUMENT, lower)) {
        cerr << "Invalid -threshold argument, must be a number!" << endl;
        exit(1);
      }
      if (HAS_ARGUMENT) {
        if (!FromString(ARGUMENT, upper)) {
          cerr << "Invalid -threshold upper argument, must be a number!" << endl;
          exit(1);
        }
        ops.push_back(unique_ptr<Op>(new Clamp(lower, upper)));
      } else {
        ops.push_back(unique_ptr<Op>(new LowerThreshold(lower)));
      }
    } else if (OPTION("-threshold-below") || OPTION("-lower-threshold")) {
      double threshold;
      if (!FromString(ARGUMENT, threshold)) {
        cerr << "Invalid -threshold-below, must be a number!" << endl;
        exit(1);
      }
      ops.push_back(unique_ptr<Op>(new LowerThreshold(threshold)));
    } else if (OPTION("-threshold-above") || OPTION("-upper-threshold")) {
      double threshold;
      if (!FromString(ARGUMENT, threshold)) {
        cerr << "Invalid -threshold-above, must be a number!" << endl;
        exit(1);
      }
      ops.push_back(unique_ptr<Op>(new UpperThreshold(threshold)));
    } else if (OPTION("-rescale")) {
      double min, max;
      if (!FromString(ARGUMENT, min)) {
        cerr << "Invalid -rescale minimum, must be a number!" << endl;
        exit(1);
      }
      if (!FromString(ARGUMENT, max)) {
        cerr << "Invalid -rescale maximum, must be a number!" << endl;
        exit(1);
      }
      ops.push_back(unique_ptr<Op>(new Rescale(min, max)));
    } else if (OPTION("-inside")) {
      double inside_value;
      if (!FromString(ARGUMENT, inside_value)) {
        cerr << "Invalid -inside value, must be a number!" << endl;
        exit(1);
      }
      ops.push_back(unique_ptr<Op>(new SetInsideValue(inside_value)));
    } else if (OPTION("-outside") || OPTION("-pad")) {
      double outside_value;
      if (!FromString(ARGUMENT, outside_value)) {
        cerr << "Invalid -outside value, must be a number!" << endl;
        exit(1);
      }
      ops.push_back(unique_ptr<Op>(new SetOutsideValue(outside_value)));
    // Data transformations
    } else if (OPTION("-add") || OPTION("-plus") || OPTION("+")) {
      const char *arg = ARGUMENT;
      double c;
      if (FromString(arg, c)) {
        ops.push_back(unique_ptr<Op>(new Add(c)));
      } else {
        ops.push_back(unique_ptr<Op>(new Add(arg)));
      }
    } else if (OPTION("-sub") || OPTION("-subtract") || OPTION("-minus") || OPTION("-")) {
      const char *arg = ARGUMENT;
      double c;
      if (FromString(arg, c)) {
        ops.push_back(unique_ptr<Op>(new Sub(c)));
      } else {
        ops.push_back(unique_ptr<Op>(new Sub(arg)));
      }
    } else if (OPTION("-mul") || OPTION("-multiply-by") || OPTION("-times") || OPTION("*")) {
      const char *arg = ARGUMENT;
      double c;
      if (FromString(arg, c)) {
        ops.push_back(unique_ptr<Op>(new Mul(c)));
      } else {
        ops.push_back(unique_ptr<Op>(new Mul(arg)));
      }
    } else if (OPTION("-div") || OPTION("-divide-by") || OPTION("-over") || OPTION("/")) {
      const char *arg = ARGUMENT;
      double c;
      if (FromString(arg, c)) {
        if (fequal(c, .0)) {
          cerr << "Invalid -div argument, value must not be zero!" << endl;
          exit(1);
        }
        ops.push_back(unique_ptr<Op>(new Div(c)));
      } else {
        ops.push_back(unique_ptr<Op>(new Div(arg)));
      }
    } else if (OPTION("-div-with-zero")) {
      ops.push_back(unique_ptr<Op>(new DivWithZero(ARGUMENT)));
    } else if (OPTION("-abs")) {
      ops.push_back(unique_ptr<Op>(new Abs()));
    } else if (OPTION("-pow") || OPTION("-power")) {
      const char *arg = ARGUMENT;
      double exponent;
      if (!FromString(arg, exponent)) {
        cerr << "Invalid -power value, must be a number!" << endl;
        exit(1);
      }
      ops.push_back(unique_ptr<Op>(new Pow(exponent)));
    } else if (OPTION("-sqrt")) {
      ops.push_back(unique_ptr<Op>(new Pow(.5)));
    } else if (OPTION("-square") || OPTION("-sq")) {
      ops.push_back(unique_ptr<Op>(new Pow(2.0)));
    } else if (OPTION("-exp")) {
      ops.push_back(unique_ptr<Op>(new Exp()));
    } else if (OPTION("-log") || OPTION("-log2") || OPTION("-log10") || OPTION("-lb") || OPTION("-ln") || OPTION("-lg")) {
      double threshold = numeric_limits<double>::min();
      if (HAS_ARGUMENT) {
        if (!FromString(ARGUMENT, threshold) || threshold <= .0) {
          cerr << "Invalid -log threshold argument, must be a positive number" << endl;
          exit(1);
        }
      }
      Op *op = NULL;
      if (strcmp(OPTNAME, "-log") == 0) {
        if (HAS_ARGUMENT) {
          double base;
          if (!FromString(ARGUMENT, base)) {
            cerr << "Invalid -log base argument, must be a positive number" << endl;
            exit(1);
          }
          op = new Log(base, threshold);
        } else {
          op = new Ln(threshold);
        }
      } else if (strcmp(OPTNAME, "-log2") == 0 || strcmp(OPTNAME, "-lb") == 0) {
        op = new Lb(threshold);
      } else if (strcmp(OPTNAME, "-log10") == 0 || strcmp(OPTNAME, "-lg") == 0) {
        op = new Lg(threshold);
      } else if (strcmp(OPTNAME, "-ln") == 0) {
        op = new Ln(threshold);
      }
      ops.push_back(unique_ptr<Op>(op));
    } else if (OPTION("-o") || OPTION("-out") || OPTION("-output") || OPTION("=")) {
      const char *fname = ARGUMENT;
      int dtype = datatype;
      if (HAS_ARGUMENT) {
        const char *arg = ARGUMENT;
        dtype = ToDataType(arg);
        if (dtype == MIRTK_VOXEL_UNKNOWN) {
          cerr << "Invalid -out data type " << arg << endl;
          exit(1);
        }
      }
#if MIRTK_Image_WITH_VTK
      ops.push_back(unique_ptr<Op>(new Write(fname, dtype, attr, dataset, scalars_name)));
#else
      ops.push_back(unique_ptr<Op>(new Write(fname, dtype, attr)));
#endif
    // Data statistics
    } else if (OPTION("-mean") || OPTION("-average") || OPTION("-avg")) {
      ops.push_back(unique_ptr<Op>(new Mean()));
    } else if (OPTION("-sigma") || OPTION("-stddev") || OPTION("-stdev") || OPTION("-std") || OPTION("-sd")) {
      ops.push_back(unique_ptr<Op>(new StDev()));
    } else if (OPTION("-normal-distribution") ||
               OPTION("-mean+sigma") || OPTION("-mean+stddev") || OPTION("-mean+stdev") || OPTION("-mean+std") || OPTION("-mean+sd") ||
               OPTION("-avg+sigma")  || OPTION("-avg+stddev")  || OPTION("-avg+stdev")  || OPTION("-avg+std")  || OPTION("-avg+sd")) {
      ops.push_back(unique_ptr<Op>(new NormalDistribution()));
    } else if (OPTION("-variance") || OPTION("-var")) {
      ops.push_back(unique_ptr<Op>(new Var()));
    } else if (OPTION("-minimum")  || OPTION("-min")) {
      ops.push_back(unique_ptr<Op>(new Min()));
    } else if (OPTION("-maximum")  || OPTION("-max")) {
      ops.push_back(unique_ptr<Op>(new Max()));
    } else if (OPTION("-extrema") || OPTION("-minmax")) {
      ops.push_back(unique_ptr<Op>(new Extrema()));
    } else if (OPTION("-range")) {
      ops.push_back(unique_ptr<Op>(new Range()));
    } else if (OPTION("-percentile") || OPTION("-pct") || OPTION("-p")) {
      do {
        int p;
        if (FromString(ARGUMENT, p) && 0 <= p && p <= 100) {
          ops.push_back(unique_ptr<Op>(new Percentile(p)));
        } else {
          cerr << "Invalid -percentile value, must be integer in the range [0, 100]!" << endl;
          exit(1);
        }
      } while (HAS_ARGUMENT);
    } else if (OPTION("-lower-percentile-mean") || OPTION("-lpctavg")) {
      do {
        int p;
        if (FromString(ARGUMENT, p) && 0 <= p && p <= 100) {
          ops.push_back(unique_ptr<Op>(new LowerPercentileMean(p)));
        } else {
          cerr << "Invalid -lower-percentile-mean value, must be integer in the range [0, 100]!" << endl;
          exit(1);
        }
      } while (HAS_ARGUMENT);
    } else if (OPTION("-upper-percentile-mean") || OPTION("-upctavg")) {
      do {
        int p;
        if (FromString(ARGUMENT, p) && 0 <= p && p <= 100) {
          ops.push_back(unique_ptr<Op>(new UpperPercentileMean(p)));
        } else {
          cerr << "Invalid -upper-percentile-mean value, must be integer in the range [0, 100]!" << endl;
          exit(1);
        }
      } while (HAS_ARGUMENT);
    } else if (OPTION("-delimiter") || OPTION("-delim") || OPTION("-d") || OPTION("-sep")) {
      delimiter = ARGUMENT;
    } else if (OPTION("-precision") || OPTION("-digits")) {
      if (!FromString(ARGUMENT, digits) || digits < 0) {
        cerr << "Invalid -precision argument, value must be non-negative integer!" << endl;
        exit(1);
      }
    } else {
      HANDLE_COMMON_OR_UNKNOWN_OPTION();
    }
  }

  // If delimiter explicitly set to empty string, use none
  if (delimiter && delimiter[0] == '\0') delimiter = NULL;

  // Default statistics to compute
  if (ops.empty()) {
    ops.push_back(unique_ptr<Statistic>(new Mean()));
    ops.push_back(unique_ptr<Statistic>(new StDev()));
    ops.push_back(unique_ptr<Statistic>(new Extrema()));
    ops.push_back(unique_ptr<Statistic>(new Range()));
  }

  // Initial data mask
  bool *mask = new bool[n];
  for (int i = 0; i < n; ++i) {
    if (IsNaN(data[i])) {
      mask[i] = false;
    } else {
      mask[i] = true;
    }
  }

  // Open output file to append to or use STDOUT if none specified
  ofstream ofs;
  if (append_name) {
    if (print_header) {
      ifstream ifs(append_name);
      if (ifs.is_open()) {
        print_header = false;
        ifs.close();
      }
    }
    ofs.open(append_name, ios_base::app);
    if (!ofs.is_open()) {
      cerr << "Cannot append to file " << append_name << endl;
      delete[] mask;
      exit(1);
    }
  }
  ostream &out = (ofs.is_open() ? ofs : cout);

  // Print column names if requested
  if (delimiter && print_header) {
    size_t c = 0;
    for (size_t i = 0; i < prefix.size(); ++i, ++c) {
      if (c > 0) out << delimiter;
      if (c < header.size()) out << header[c];
    }
    for (size_t i = 0; i < ops.size(); ++i) {
      Statistic *stat = dynamic_cast<Statistic *>(ops[i].get());
      if (stat) {
        for (size_t j = 0; j < stat->Names().size(); ++j, ++c) {
          if (c > 0) out << delimiter;
          if (c < header.size()) out << header[c];
          else out << stat->Names()[j];
        }
      }
    }
    out << endl;
  }

  // Process input data, either transform it or compute statistics from it
  for (size_t i = 0; i < ops.size(); ++i) ops[i]->Process(n, data, mask);
  delete[] mask;

  // Print image statistics
  if (delimiter) {
    for (size_t i = 0; i < prefix.size(); ++i) {
      if (i > 0) out << delimiter;
      out << prefix[i];
    }
    bool first = prefix.empty();
    for (size_t i = 0; i < ops.size(); ++i) {
      Statistic *stat = dynamic_cast<Statistic *>(ops[i].get());
      if (stat) {
        if (!first) out << delimiter;
        else first = false;
        stat->PrintValues(out, digits, delimiter);
      }
    }
    // No newline at end of row if printing results to STDOUT which in this
    // case is usually assigned to a string in a calling script
    if (print_header || ofs.is_open()) out << endl;
  } else {
    string prefix_string;
    for (size_t i = 0; i < prefix.size(); ++i) {
      if (i > 0) prefix_string += ' ';
      prefix_string += prefix[i];
    }
    for (size_t i = 0; i < ops.size(); ++i) {
      Statistic *stat = dynamic_cast<Statistic *>(ops[i].get());
      if (stat) stat->Print(out, digits, prefix_string.c_str());
    }
  }

  ofs.close();
  return 0;
}

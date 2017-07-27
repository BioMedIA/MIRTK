/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2015-2017 Imperial College London
 * Copyright 2015-2017 Andreas Schuh
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

#include "mirtk/ImageConfig.h"
#include "mirtk/IOConfig.h"

#include "mirtk/DataOp.h"
#include "mirtk/DataStatistics.h"
#include "mirtk/DataFunctions.h"

#if MIRTK_Image_WITH_VTK
  #include "vtkDataSet.h"
  #include "vtkSmartPointer.h"
  #include "vtkPointData.h"
  #include "vtkCellData.h"
  #include "vtkDataArray.h"
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
  cout << "\n";
  cout << "Usage: " << name << " <input> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  This tool can be used for basic calculations from a sequence of data values read\n";
  cout << "  either from an image or a VTK pointset. It can be used, for example, to add two\n";
  cout << "  data sequences and to divide the result by a constant. The current sequence can\n";
  cout << "  be written to an output file again using :option:`-out`. Additionally, statistics\n";
  cout << "  of the current data sequence can be computed such as the mean or variance.\n";
  cout << "  The order of the data transformations and calculation of statistics is determined\n";
  cout << "  by the order of the command-line arguments.\n";
  cout << "\n";
  cout << "  The data mask is used to include/exclude values from subsequent operations.\n";
  cout << "  Initially, all NaN values in the input data sequence are excluded.\n";
  cout << "  Further values can be excluded using one or more of the masking operations.\n";
  cout << "  Using the mask, operations can be performed on only a subset of the data,\n";
  cout << "  and the mask then reset using :option:`-reset-mask`.\n";
  cout << "\n";
  cout << "  By default, data statistics are printed to STDOUT in a human readable format.\n";
  cout << "  This output can be appended to a text file using :option:`-append` instead.\n";
  cout << "  For a more machine readable output, e.g., as comma separated values (CSV),\n";
  cout << "  specify a delimiting string using :option:`-delimiter`. In this case, a header\n";
  cout << "  line is also printed when :option:`-header` is given with optional user\n";
  cout << "  specified column names for the individual output values.\n";
  cout << "\n";
  cout << "Input options:\n";
  cout << "  -pd, -point-data, -scalars <name>   Name of input point data array. (default: active SCALARS array)\n";
  cout << "  -cd, -cell-data <name>              Name of input cell  data array. Overrides :option:`-pd`.\n";
  cout << "\n";
  cout << "Data masking options:\n";
  cout << "  -even\n";
  cout << "      Exclude values which are not an even number when cast to an integer.\n";
  cout << "  -odd\n";
  cout << "      Exclude values which are not an odd number when cast to an integer.\n";
  cout << "  -label <value|lower..upper>...\n";
  cout << "      Include data points with a value equal to either one of the given values.\n";
  cout << "      Closed intervals of values can be specified as \"lower..upper\".\n";
  cout << "      For example, \"-label 1 3 5..6 10 20..50\". This option is a shorthand for\n";
  cout << "        :option:`-mask-all` :option:`-threshold-inside` <lower> <upper> :option:`-invert-mask`\n";
  cout << "      where one :option:`-threshold-inside` operation is performed for each argument.\n";
  cout << "  -mask <value>... | <file> [<scalars>] [<value>]\n";
  cout << "      Exclude values equal a given threshold or with specified input mask <value>.\n";
  cout << "      The default mask value of values to be excluded is zero. When the input file\n";
  cout << "      is a point set file (e.g., .vtk, .vtp), the optional <scalars> argument can be\n";
  cout << "      used to specify the name of the point/cell data array to use as mask.\n";
  cout << "      Note that this operation does not modify the data values, but only marks them\n";
  cout << "      to be ignored from now on. Use :option:`-pad` following this operation to\n";
  cout << "      replace these values by a constant background value.\n";
  cout << "  -mask-all\n";
  cout << "      Exclude all values.\n";
  cout << "  -reset-mask\n";
  cout << "      Reset mask to include all values again.\n";
  cout << "  -invert-mask\n";
  cout << "      Invert mask to include all values that where excluded before and\n";
  cout << "      exclude all values that were included before.\n";
  cout << "  -set, -inside <value>\n";
  cout << "      Set new value for all currently included data values.\n";
  cout << "  -pad, -outside <value>\n";
  cout << "      Set new value for all currently excluded data values.\n";
  cout << "\n";
  cout << "Data thresholding options:\n";
  cout << "  -threshold <lower> [<upper>]\n";
  cout << "      This masking operation is equivalent to :option:`-threshold-outside`.\n";
  cout << "      When no upper threshold is specified, it defaults to +inf. Therefore,\n";
  cout << "      \"-threshold 0\" will exclude all negative values.\n";
  cout << "  -percentile-threshold, -pct-threshold <lower>\n";
  cout << "      This masking operation is equivalent to :option:`-threshold-outside-percentiles`.\n";
  cout << "      with an upper threshold of +inf. Therefore, \"-threshold 0\" excludes all negative values.\n";
  cout << "  -threshold-percentiles, -threshold-pcts <lower> <upper>\n";
  cout << "      This masking operation is equivalent to :option:`-threshold-outside-percentiles`.\n";
  cout << "  -threshold-inside, -mask-inside <lower> <upper>\n";
  cout << "      Exclude values which are inside a given closed interval.\n";
  cout << "      When the lower threshold is greater than the upper threshold,\n";
  cout << "      values less than or equal to the upper threshold and values greater\n";
  cout << "      than or equal to the lower threshold are excluded.\n";
  cout << "  -threshold-inside-percentiles, -threshold-inside-pcts, -mask-inside-percentiles, -mask-inside-pct <lower> <upper>\n";
  cout << "      Exclude values which are inside a given closed interval of percentiles.\n";
  cout << "      When the lower percentile is greater than the upper percentile,\n";
  cout << "      values less than or equal to the upper percentile and values greater\n";
  cout << "      than or equal to the lower percentile are excluded.\n";
  cout << "  -threshold-outside, -mask-outside <lower> <upper>\n";
  cout << "      Exclude values which are outside a given open interval.\n";
  cout << "      When the lower threshold is greater than the upper threshold,\n";
  cout << "      values inside the closed interval <upper>..<lower> are excluded.\n";
  cout << "  -threshold-outside-percentiles, -threshold-outside-pcts, -mask-outside-percentiles, -mask-outside-pcts <lower> <upper>\n";
  cout << "      Exclude values which are outside a given open interval of percentiles.\n";
  cout << "      When the lower percentile is greater than the upper percentile,\n";
  cout << "      values inside the closed interval <upper>..<lower> are excluded.\n";
  cout << "  -threshold-lt, -lower-threshold, -mask-lt <value>\n";
  cout << "      Exclude values less than a given threshold.\n";
  cout << "  -threshold-lt-percentile, -threshold-lt-pct, -lower-percentile-threshold, -lower-pct-threshold, -mask-lt-percentile, -mask-lt-pct <value>\n";
  cout << "      Exclude values less than a given precentile.\n";
  cout << "  -threshold-le, -mask-le, -mask-below <value>\n";
  cout << "      Exclude values less than or equal to a given threshold.\n";
  cout << "  -threshold-le-percentile, -threshold-le-pct, -mask-le-percentile, -mask-le-pct, -mask-below-percentile, -mask-below-pct <value>\n";
  cout << "      Exclude values less than or equal to a given percentile.\n";
  cout << "  -threshold-ge, -mask-ge, -mask-above <value>\n";
  cout << "      Exclude values greater than or equal to a given threshold.\n";
  cout << "  -threshold-ge-percentile, -threshold-ge-pct, -mask-ge-percentile, -mask-ge-pct, -mask-above-percentile, -mask-above-pct <value>\n";
  cout << "      Exclude values greater than or equal to a given percentile.\n";
  cout << "  -threshold-gt, -upper-threshold, -mask-gt <value>\n";
  cout << "      Exclude values greater than a given threshold.\n";
  cout << "  -threshold-gt-percentile, -threshold-gt-pct, -upper-percentile-threshold, -upper-pct-threshold, -mask-gt-percentile, -mask-gt-pct <value>\n";
  cout << "      Exclude values greater than a given percentile.\n";
  cout << "\n";
  cout << "Data rescaling options:\n";
  cout << "  -binarize <lower> [<upper>]\n";
  cout << "      Set values inside the closed interval <lower>..<upper> to one,\n";
  cout << "      and all other values to zero. The default upper threshold is +inf.\n";
  cout << "      When the lower threshold is greater than the upper threshold,\n";
  cout << "      values inside the closed interval <upper>..<lower> are set to zero\n";
  cout << "      and all other values to one instead. This operation is short for:\n";
  cout << "      :option:`-threshold-inside` <lower> <upper> :option:`-set` 1 :option:`-pad` 0\n";
  cout << "  -clamp <lower> <upper>\n";
  cout << "      Clamp values which are less than a lower or greater than an upper threshold.\n";
  cout << "  -clamp-percentiles, -clamp-pcts <lower> <upper>\n";
  cout << "      Clamp values which are less than a lower percentile or greater than an upper percentile.\n";
  cout << "  -clamp-below, -clamp-lt <value>\n";
  cout << "      Clamp values less than a given threshold.\n";
  cout << "  -clamp-below-percentile, -clamp-below-pct, -clamp-lt-percentile, -clamp-lt-pct <value>\n";
  cout << "      Clamp values less than a given percentile.\n";
  cout << "  -clamp-above, -clamp-gt <value>\n";
  cout << "      Clamp values greater than a given threshold.\n";
  cout << "  -clamp-above-percentile, -clamp-above-pct, -clamp-gt-percentile, -clamp-gt-pct <value>\n";
  cout << "      Clamp values greater than a given percentile.\n";
  cout << "  -rescale <min> <max>\n";
  cout << "      Linearly rescale values to the interval [min, max].\n";
  cout << "  -map <from> <to>...\n";
  cout << "      Replaces values equal to <from> by the specified <to> value. Multiple pairs of <from>\n";
  cout << "      and <to> value replacements can be specified in order to perform the substitutions in\n";
  cout << "      one step. For example, to swap the two values 1 and 2, use ``-map 1 2 2 1``.\n";
  cout << "\n";
  cout << "Arithmetic operation options:\n";
  cout << "  -add, -plus <value> | <file> [<scalars>]\n";
  cout << "      Add constant value or data sequence read from specified file.\n";
  cout << "      Another name for this option is the '+' sign, see Examples.\n";
  cout << "  -sub, -subtract, -minus <value> | <file> [<scalars>]\n";
  cout << "      Subtract constant value or data sequence read from specified file.\n";
  cout << "      Another name for this option is the '-' sign, see Examples.\n";
  cout << "  -mul, -multiply-with, -times <value> | <file> [<scalars>]\n";
  cout << "      Multiply by constant value or data sequence read from specified file.\n";
  cout << "      Another name for this option is the '*' sign, see Examples.\n";
  cout << "  -div, -divide-by, -over <value> | sum | <file> [<scalars>]\n";
  cout << "      Divide by constant value or data sequence read from specified file.\n";
  cout << "      When the argument is \"sum\", the divisor is the sum of the values.\n";
  cout << "      When dividing by zero values in the input file, the result is NaN.\n";
  cout << "      Use :option:`-mask` with argument NaN and :option:`-pad` to replace\n";
  cout << "      these undefined values by a constant such as zero.\n";
  cout << "      Another name for this option is the '/' sign, see Examples.\n";
  cout << "  -div-with-zero <value> | sum | <file> [<scalars>]\n";
  cout << "      Same as :option:`-div`, but set result to zero in case of division by zero.\n";
  cout << "  -abs\n";
  cout << "      Replace values by their respective absolute value.\n";
  cout << "  -pow, -power <exponent>\n";
  cout << "      Raise values to the power of the given exponent.\n";
  cout << "  -sq, -square\n";
  cout << "      Raise values to the power of 2 (i.e, -pow 2).\n";
  cout << "  -sqrt\n";
  cout << "      Calculate square root of each value (i.e, -pow .5).\n";
  cout << "  -exp\n";
  cout << "      Calculate exponential of data sequence.\n";
  cout << "  -log [<threshold>] [<base>]\n";
  cout << "      Compute logarithm after applying an optional threshold.\n";
  cout << "      (default threshold: min double, default base: e)\n";
  cout << "  -lb, -log2 [<threshold>]\n";
  cout << "      Compute binary logarithm, alias for :option:`-log` with base 2.\n";
  cout << "  -ln, -loge [<threshold>]\n";
  cout << "      Compute natural logarithm, alias for :option:`-log` with base e.\n";
  cout << "  -lg, -log10 [<threshold>]\n";
  cout << "      Compute logarithm to base 10, alias for :option:`-log` with base 10.\n";
  cout << "\n";
  cout << "Data output options:\n";
  cout << "  -out, -o, -output <file> [<type>] [<name>]\n";
  cout << "      Write current data sequence to file in the format of the input file.\n";
  cout << "      Output data type can be: uchar, short, ushort, int, uint, float, double.\n";
  cout << "      The optional <name> argument can be used to save the modified data\n";
  cout << "      of an input point set data array with a different name along with the\n";
  cout << "      input data. Otherwise, the input data values are replaced by the modified\n";
  cout << "      values and stored with point data array name is unchanged.\n";
  cout << "      Another name for this option is the '=' sign, but the optional arguments are\n";
  cout << "      are not supported by this alternative notation. See Examples for usage.\n";
  cout << "\n";
  cout << "Data statistics options:\n";
  cout << "  -append <file>\n";
  cout << "      Append output to a file. (default: STDOUT)\n";
  cout << "  -delimiter, -delim, -d, -sep\n";
  cout << "      Delimiting character(s). (default: '')\n";
  cout << "  -header [<name>...]\n";
  cout << "      Request output of header line if delimiter was specified as well.\n";
  cout << "      If the output is appended to a text file, the header is only printed\n";
  cout << "      if it does not exist. If no or fewer custom column names are given,\n";
  cout << "      the default names for each statistic are printed. (default: none)\n";
  cout << "  -prefix <str>...\n";
  cout << "      One or more prefix strings to print. If no delimiter is specified,\n";
  cout << "      the concatenated strings are printed before each line of the output.\n";
  cout << "      Otherwise, each prefix string is printed as entry for the first columns\n";
  cout << "      in the delimited output row, separated by the specified delimiter. (default: none)\n";
  cout << "  -precision, -digits <int>\n";
  cout << "      Number of significant digits. (default: 5)\n";
  cout << "  -median\n";
  cout << "      Print median value, i.e., 50th percentile. (default: off)\n";
  cout << "  -mean, -avg, -average\n";
  cout << "      Print mean value. (default: on)\n";
  cout << "  -variance, -var\n";
  cout << "      Print variance of values. (default: off)\n";
  cout << "  -sigma, -std, -stddev, -stdev, -sd\n";
  cout << "      Print standard deviation of values. (default: on)\n";
  cout << "  -normal-distribution\n";
  cout << "      Print mean and standard deviation of values.\n";
  cout << "      Other option names: -mean+sigma, -mean+sd, -avg+std,... (default: off)\n";
  cout << "  -mad, -mean-absolute-difference, -mean-absolute-deviation\n";
  cout << "      Print mean absolute difference/deviation around the mean. (default: off)\n";
  cout << "  -mad-median, -median-absolute-difference, -median-absolute-deviation\n";
  cout << "      Print mean absolute difference/deviation around the median. (default: off)\n";
  cout << "  -minimum, -min\n";
  cout << "      Print minimum value. (default: off)\n";
  cout << "  -maximum, -max\n";
  cout << "      Print maximum value. (default: off)\n";
  cout << "  -extrema, -minmax\n";
  cout << "      Print minimum and maximum value. (default: on)\n";
  cout << "  -range\n";
  cout << "      Print range of values (i.e., max - min). (default: off)\n";
  cout << "  -percentile, -pct, -p <n>...\n";
  cout << "      Print n-th percentile. (default: none)\n";
  cout << "  -lower-percentile-mean, -lpctavg <n>\n";
  cout << "      Print mean intensity of values less than or equal to the n-th percentile. (default: off)\n";
  cout << "  -upper-percentile-mean, -upctavg <n>\n";
  cout << "      Print mean intensity of values greater than or equal to the n-th percentile. (default: off)\n";
  cout << "  -sum\n";
  cout << "      Print sum of values. Can be used to count values within a certain range using a thresholding\n";
  cout << "      followed by :option:`-set` 1 before summing these values. (default: off)\n";
  cout << "  -count\n";
  cout << "      Print number of values inside the mask, i.e., values not currently excluded. (default: off)\n";
  PrintCommonOptions(cout);
  cout << "\n";
  cout << "Examples:\n";
  cout << "\n";
  cout << "  " << name << " mni305.nii.gz\n";
  cout << "      Mean = 26.9753\n";
  cout << "      Standard deviation = 50.3525\n";
  cout << "      Extrema = [0, 254]\n";
  cout << "      Range = 254\n";
  cout << "\n";
  cout << "  " << name << " mni305.nii.gz -pct 77\n";
  cout << "      77th percentile = 25\n";
  cout << "\n";
  cout << "  " << name << " mni305.nii.gz -padding 25 -range -percentile 25 50 75 -prefix MNI305 '[>25]'\n";
  cout << "      MNI305 [>25] range = 254\n";
  cout << "      MNI305 [>25] 25th percentile = 69\n";
  cout << "      MNI305 [>25] 50th percentile = 113\n";
  cout << "      MNI305 [>25] 75th percentile = 150\n";
  cout << "\n";
  cout << "  " << name << " mni305.nii.gz -d , -prefix MNI305\n";
  cout << "      MNI305,26.9753,50.3525,0,254,254 [no newline at end of line]\n";
  cout << "\n";
  cout << "  " << name << " mni305.nii.gz -d , -prefix MNI305 -header\n";
  cout << "      ,Mean,Sigma,Min,Max,Range\n";
  cout << "      MNI305,26.9753,50.3525,0,254,254\n";
  cout << "\n";
  cout << "  " << name << " mni305.nii.gz -d , -prefix MNI305 -header ID Mean SD\n";
  cout << "      ID,Mean,SD,Min,Max,Range\n";
  cout << "      MNI305,26.9753,50.3525,0,254,254\n";
  cout << "\n";
  cout << "  " << name << " a.nii.gz + b.nii.gz = c.nii.gz\n";
  cout << "\n";
  cout << "  " << name << " a.vtk + b.nii.gz - 10 / c.nii = d.vtk\n";
  cout << "      Adds data values at identical sequential memory indices in a and b,\n";
  cout << "      subtracts the constant 10, and then divides by the values in image c.\n";
  cout << "\n";
  cout << "      Note: Operations are always executed from left to right,\n";
  cout << "            i.e., no mathematical operator precedence is considered!\n";
  cout << "\n";
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

  const char *input_name = POSARG(1);

  UniquePtr<double[]> data;
  int datatype = MIRTK_VOXEL_DOUBLE;
  ImageAttributes attr;

#if MIRTK_Image_WITH_VTK
  const char *scalars_name = nullptr;
  bool        cell_data    = false;
  for (ARGUMENTS_AFTER(1)) {
    if (OPTION("-point-data") || OPTION("-pointdata") || OPTION("-pd") || OPTION("-scalars")) {
      scalars_name = ARGUMENT;
      cell_data    = false;
    }
    else if (OPTION("-cell-data") || OPTION("-celldata") || OPTION("-cd")) {
      scalars_name = ARGUMENT;
      cell_data    = true;
    }
  }
  vtkSmartPointer<vtkDataSet> dataset;
  vtkSmartPointer<vtkDataSetAttributes> arrays;
  int n = Read(input_name, data, &datatype, &attr, &dataset, scalars_name, cell_data);
  if (dataset) {
    if (cell_data) {
      arrays = dataset->GetCellData();
    } else {
      arrays = dataset->GetPointData();
    }
  }
#else // MIRTK_Image_WITH_VTK
  int n = Read(input_name, data, &datatype, &attr);
#endif // MIRTK_Image_WITH_VTK

  // Optional arguments
  const double inf = numeric_limits<double>::infinity();
  const double nan = numeric_limits<double>::quiet_NaN();

  double      a, b;
  int         p;
  const char *append_name   = NULL;
  const char *delimiter     = NULL;
  bool        print_header  = false;
  int         digits        = 5;

  Array<string> header;
  Array<string> prefix;

  Array<UniquePtr<Op> > ops;

  for (ARGUMENTS_AFTER(1)) {
    if (OPTION("-append")) {
      append_name = ARGUMENT;
    } else if (OPTION("-point-data") || OPTION("-pointdata") || OPTION("-pd") || OPTION("-scalars")) {
      #if MIRTK_Image_WITH_VTK
        // Parsed before Read above
        scalars_name = ARGUMENT;
        cell_data    = false;
      #else
        FatalError("Cannot process -point-data of VTK file because MIRTK Image library was built without VTK!");
      #endif // MIRTK_Image_WITH_VTK
    } else if (OPTION("-cell-data") || OPTION("-celldata") || OPTION("-cd")) {
      #if MIRTK_Image_WITH_VTK
        // Parsed before Read above
        scalars_name = ARGUMENT;
        cell_data    = true;
      #else
        FatalError("Cannot process -cell-data of VTK file because MIRTK Image library was built without VTK!");
      #endif // MIRTK_Image_WITH_VTK
    } else if (OPTION("-prefix")) {
      do {
        prefix.push_back(ARGUMENT);
      } while (HAS_ARGUMENT);
    } else if (OPTION("-header")) {
      print_header = true;
      while (HAS_ARGUMENT) header.push_back(ARGUMENT);
    // Masking
    } else if (OPTION("-label")) {
      ops.push_back(UniquePtr<Op>(new ResetMask(true)));
      do {
        const char *arg = ARGUMENT;
        const Array<string> parts = Split(arg, "..");
        if (parts.size() == 1) {
          if (!FromString(parts[0], a)) a = nan;
          b = a;
        } else if (parts.size() == 2) {
          if (!FromString(parts[0], a) || !FromString(parts[1], b)) {
            a = b = nan;
          }
        } else {
          a = b = nan;
        }
        if (IsNaN(a) || IsNaN(b)) {
          FatalError("Invalid -label argument: " << arg);
        }
        ops.push_back(UniquePtr<Op>(new MaskInsideInterval(a, b)));
      } while (HAS_ARGUMENT);
      ops.push_back(UniquePtr<Op>(new InvertMask()));
    } else if (OPTION("-mask-all")) {
      ops.push_back(UniquePtr<Op>(new ResetMask(false)));
    } else if (OPTION("-reset-mask")) {
      ops.push_back(UniquePtr<Op>(new ResetMask(true)));
    } else if (OPTION("-invert-mask")) {
      ops.push_back(UniquePtr<Op>(new InvertMask()));
    } else if (OPTION("-mask")) {
      double c;
      do {
        const char *arg = ARGUMENT;
        if (FromString(arg, c)) {
          ops.push_back(UniquePtr<Op>(new Mask(c)));
        } else {
          const char *fname = arg;
          const char *aname = nullptr;
          if (HAS_ARGUMENT) {
            arg = ARGUMENT;
            if (HAS_ARGUMENT) {
              aname = arg;
              PARSE_ARGUMENT(c);
            } else if (!FromString(arg, c)) {
              aname = arg, c = 0.;
            }
          } else {
            c = 0.;
            #if MIRTK_Image_WITH_VTK
              if (dataset && arrays->HasArray(fname)) {
                aname = fname;
                fname = input_name;
              }
            #endif
          }
          UniquePtr<Mask> op(new Mask(fname, c));
          if (aname) {
            #if MIRTK_Image_WITH_VTK
              op->ArrayName(aname);
              op->IsCellData(cell_data);
            #else
              FatalError("Cannot read point set files when build without VTK or wrong usage!");
            #endif
          }
          ops.push_back(UniquePtr<Op>(op.release()));
          break;
        }
      } while (HAS_ARGUMENT);
    } else if (OPTION("-threshold-outside") || OPTION("-mask-outside")) {
      PARSE_ARGUMENT(a);
      PARSE_ARGUMENT(b);
      ops.push_back(UniquePtr<Op>(new MaskOutsideOpenInterval(a, b)));
    } else if (OPTION("-threshold-outside-percentiles") || OPTION("-threshold-outside-pcts") ||
               OPTION("-mask-outside-percentiles")      || OPTION("-mask-outside-pcts")) {
      PARSE_ARGUMENT(p);
      Statistic *a = new Percentile(p);
      a->Hidden(verbose < 1);
      ops.push_back(UniquePtr<Op>(a));
      PARSE_ARGUMENT(p);
      Statistic *b = new Percentile(p);
      b->Hidden(verbose < 1);
      ops.push_back(UniquePtr<Op>(b));
      Op *op = new MaskOutsideOpenInterval(&a->Value(), &b->Value());
      ops.push_back(UniquePtr<Op>(op));
    } else if (OPTION("-threshold")) {
      PARSE_ARGUMENT(a);
      if (HAS_ARGUMENT) PARSE_ARGUMENT(b);
      else b = inf;
      ops.push_back(UniquePtr<Op>(new MaskOutsideInterval(a, b)));
    } else if (OPTION("-percentile-threshold") || OPTION("-pct-threshold")) {
      PARSE_ARGUMENT(p);
      Statistic *a = new Percentile(p);
      a->Hidden(verbose < 1);
      ops.push_back(UniquePtr<Op>(a));
      Op *op = new MaskOutsideInterval(&a->Value(), inf);
      ops.push_back(UniquePtr<Op>(op));
    } else if (OPTION("-threshold-percentiles") || OPTION("-threshold-pcts")) {
      PARSE_ARGUMENT(p);
      Statistic *a = new Percentile(p);
      a->Hidden(verbose < 1);
      ops.push_back(UniquePtr<Op>(a));
      PARSE_ARGUMENT(p);
      Statistic *b = new Percentile(p);
      b->Hidden(verbose < 1);
      ops.push_back(UniquePtr<Op>(b));
      Op *op = new MaskOutsideInterval(&a->Value(), &b->Value());
      ops.push_back(UniquePtr<Op>(op));
    } else if (OPTION("-threshold-inside") || OPTION("-mask-inside")) {
      PARSE_ARGUMENT(a);
      PARSE_ARGUMENT(b);
      ops.push_back(UniquePtr<Op>(new MaskInsideInterval(a, b)));
    } else if (OPTION("-threshold-inside-percentiles") || OPTION("-threshold-inside-pcts") ||
               OPTION("-mask-inside-percentiles")      || OPTION("-mask-inside-pcts")) {
      PARSE_ARGUMENT(p);
      Statistic *a = new Percentile(p);
      a->Hidden(verbose < 1);
      ops.push_back(UniquePtr<Op>(a));
      PARSE_ARGUMENT(p);
      Statistic *b = new Percentile(p);
      b->Hidden(verbose < 1);
      ops.push_back(UniquePtr<Op>(b));
      Op *op = new MaskInsideInterval(&a->Value(), &b->Value());
      ops.push_back(UniquePtr<Op>(op));
    } else if (OPTION("-threshold-lt") || OPTION("-lower-threshold") || OPTION("-mask-lt")) {
      PARSE_ARGUMENT(a);
      ops.push_back(UniquePtr<Op>(new MaskOutsideInterval(a, inf)));
    } else if (OPTION("-threshold-lt-percentile")    || OPTION("-threshold-lt-pct") ||
               OPTION("-lower-percentile-threshold") || OPTION("-lower-pct-threshold") ||
               OPTION("-mask-lt-percentile")         || OPTION("-mask-lt-pct")) {
      PARSE_ARGUMENT(p);
      Statistic *a = new Percentile(p);
      a->Hidden(verbose < 1);
      ops.push_back(UniquePtr<Op>(a));
      ops.push_back(UniquePtr<Op>(new MaskOutsideInterval(&a->Value(), inf)));
    } else if (OPTION("-threshold-le") || OPTION("-mask-below") || OPTION("-mask-le")) {
      PARSE_ARGUMENT(a);
      ops.push_back(UniquePtr<Op>(new MaskOutsideOpenInterval(a, inf)));
    } else if (OPTION("-threshold-le-percentile") || OPTION("-threshold-le-pct") ||
               OPTION("-mask-below-percentile")   || OPTION("-mask-below-pct") ||
               OPTION("-mask-le-percentile")      || OPTION("-mask-le-pct")) {
      PARSE_ARGUMENT(p);
      Statistic *a = new Percentile(p);
      a->Hidden(verbose < 1);
      ops.push_back(UniquePtr<Op>(a));
      ops.push_back(UniquePtr<Op>(new MaskOutsideOpenInterval(&a->Value(), inf)));
    } else if (OPTION("-threshold-ge") || OPTION("-mask-above") || OPTION("-mask-ge")) {
      PARSE_ARGUMENT(b);
      ops.push_back(UniquePtr<Op>(new MaskOutsideOpenInterval(-inf, b)));
    } else if (OPTION("-threshold-ge-percentile") || OPTION("-threshold-ge-pct") ||
               OPTION("-mask-above-percentile")   || OPTION("-mask-above-pct") ||
               OPTION("-mask-ge-percentile")      || OPTION("-mask-ge-pct")) {
      PARSE_ARGUMENT(p);
      Statistic *b = new Percentile(p);
      b->Hidden(verbose < 1);
      ops.push_back(UniquePtr<Op>(b));
      ops.push_back(UniquePtr<Op>(new MaskOutsideOpenInterval(-inf, &b->Value())));
    } else if (OPTION("-threshold-gt") || OPTION("-upper-threshold") || OPTION("-mask-gt")) {
      PARSE_ARGUMENT(b);
      ops.push_back(UniquePtr<Op>(new MaskOutsideInterval(-inf, b)));
    } else if (OPTION("-threshold-gt-percentile")    || OPTION("-threshold-gt-pct") ||
               OPTION("-upper-percentile-threshold") || OPTION("-upper-pct-threshold") ||
               OPTION("-mask-gt-percentile")         || OPTION("-mask-gt-pct")) {
      PARSE_ARGUMENT(p);
      Statistic *b = new Percentile(p);
      b->Hidden(verbose < 1);
      ops.push_back(UniquePtr<Op>(b));
      ops.push_back(UniquePtr<Op>(new MaskOutsideInterval(-inf, &b->Value())));
    } else if (OPTION("-even")) {
      ops.push_back(UniquePtr<Op>(new MaskOddValues()));
    } else if (OPTION("-odd")) {
      ops.push_back(UniquePtr<Op>(new MaskEvenValues()));
    // Clamping
    } else if (OPTION("-clamp")) {
      PARSE_ARGUMENT(a);
      PARSE_ARGUMENT(b);
      ops.push_back(UniquePtr<Op>(new Clamp(a, b)));
    } else if (OPTION("-clamp-percentiles") || OPTION("-clamp-pcts")) {
      PARSE_ARGUMENT(p);
      Statistic *a = new Percentile(p);
      a->Hidden(verbose < 1);
      ops.push_back(UniquePtr<Op>(a));
      PARSE_ARGUMENT(p);
      Statistic *b = new Percentile(p);
      b->Hidden(verbose < 1);
      ops.push_back(UniquePtr<Op>(b));
      ops.push_back(UniquePtr<Op>(new Clamp(&a->Value(), &b->Value())));
    } else if (OPTION("-clamp-lt") || OPTION("-clamp-below")) {
      PARSE_ARGUMENT(a);
      ops.push_back(UniquePtr<Op>(new LowerThreshold(a)));
    } else if (OPTION("-clamp-lt-percentile")    || OPTION("-clamp-lt-pct") ||
               OPTION("-clamp-below-percentile") || OPTION("-clamp-below-pct")) {
      PARSE_ARGUMENT(p);
      Statistic *a = new Percentile(p);
      a->Hidden(verbose < 1);
      ops.push_back(UniquePtr<Op>(a));
      ops.push_back(UniquePtr<Op>(new LowerThreshold(&a->Value())));
    } else if (OPTION("-clamp-gt") || OPTION("-clamp-above")) {
      PARSE_ARGUMENT(b);
      ops.push_back(UniquePtr<Op>(new UpperThreshold(b)));
    } else if (OPTION("-clamp-gt-percentile")    || OPTION("-clamp-gt-pct") ||
               OPTION("-clamp-above-percentile") || OPTION("-clamp-above-pct")) {
      PARSE_ARGUMENT(p);
      Statistic *b = new Percentile(p);
      b->Hidden(verbose < 1);
      ops.push_back(UniquePtr<Op>(b));
      ops.push_back(UniquePtr<Op>(new UpperThreshold(&b->Value())));
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
      ops.push_back(UniquePtr<Op>(new Rescale(min, max)));
    } else if (OPTION("-set") || OPTION("-inside")) {
      double inside_value;
      if (!FromString(ARGUMENT, inside_value)) {
        cerr << "Invalid -inside value, must be a number!" << endl;
        exit(1);
      }
      ops.push_back(UniquePtr<Op>(new SetInsideValue(inside_value)));
    } else if (OPTION("-pad") || OPTION("-outside")) {
      double outside_value;
      if (!FromString(ARGUMENT, outside_value)) {
        cerr << "Invalid -outside value, must be a number!" << endl;
        exit(1);
      }
      ops.push_back(UniquePtr<Op>(new SetOutsideValue(outside_value)));
    // Data transformations
    } else if (OPTION("-binarize")) {
      PARSE_ARGUMENT(a);
      if (HAS_ARGUMENT) PARSE_ARGUMENT(b);
      else b = inf;
      ops.push_back(UniquePtr<Op>(new Binarize(a, b)));
    } else if (OPTION("-map")) {
      UniquePtr<Map> map(new Map());
      do {
        const char * const arg1 = ARGUMENT;
        const char * const arg2 = ARGUMENT;
        if (!FromString(arg1, a) || !FromString(arg2, b)) {
          FatalError("Arguments of -map option must be pairs of two numbers (i.e., number of arguments must be even)!");
        }
        map->Insert(a, b);
      } while (HAS_ARGUMENT);
      ops.push_back(UniquePtr<Op>(map.release()));
    } else if (OPTION("-add") || OPTION("-plus") || OPTION("+")) {
      const char *arg = ARGUMENT;
      double c;
      if (FromString(arg, c)) {
        ops.push_back(UniquePtr<Op>(new Add(c)));
      } else {
        const char *fname = arg;
        const char *aname = nullptr;
        if (HAS_ARGUMENT) {
          aname = ARGUMENT;
        } else {
          #if MIRTK_Image_WITH_VTK
            if (dataset && arrays->HasArray(fname)) {
              aname = fname;
              fname = input_name;
            }
          #endif
        }
        UniquePtr<Add> op(new Add(fname));
        if (aname) {
          #if MIRTK_Image_WITH_VTK
            op->ArrayName(aname);
            op->IsCellData(cell_data);
          #else
            FatalError("Cannot read scalars from point set file when build without VTK or wrong usage!");
          #endif
        }
        ops.push_back(UniquePtr<Op>(op.release()));
      }
    } else if (OPTION("-sub") || OPTION("-subtract") || OPTION("-minus") || OPTION("-")) {
      const char *arg = ARGUMENT;
      double c;
      if (FromString(arg, c)) {
        ops.push_back(UniquePtr<Op>(new Sub(c)));
      } else {
        const char *fname = arg;
        const char *aname = nullptr;
        if (HAS_ARGUMENT) {
          aname = ARGUMENT;
        } else {
          #if MIRTK_Image_WITH_VTK
            if (dataset && arrays->HasArray(fname)) {
              aname = fname;
              fname = input_name;
            }
          #endif
        }
        UniquePtr<Sub> op(new Sub(fname));
        if (aname) {
          #if MIRTK_Image_WITH_VTK
            op->ArrayName(aname);
            op->IsCellData(cell_data);
          #else
            FatalError("Cannot read point set files when build without VTK or wrong usage!");
          #endif
        }
        ops.push_back(UniquePtr<Op>(op.release()));
      }
    } else if (OPTION("-mul") || OPTION("-multiply-by") || OPTION("-times") || OPTION("*")) {
      const char *arg = ARGUMENT;
      double c;
      if (FromString(arg, c)) {
        ops.push_back(UniquePtr<Op>(new Mul(c)));
      } else {
        const char *fname = arg;
        const char *aname = nullptr;
        if (HAS_ARGUMENT) {
          aname = ARGUMENT;
        } else {
          #if MIRTK_Image_WITH_VTK
            if (dataset && arrays->HasArray(fname)) {
              aname = fname;
              fname = input_name;
            }
          #endif
        }
        UniquePtr<Mul> op(new Mul(fname));
        if (aname) {
          #if MIRTK_Image_WITH_VTK
            op->ArrayName(aname);
            op->IsCellData(cell_data);
          #else
            FatalError("Cannot read point set files when build without VTK or wrong usage!");
          #endif
        }
        ops.push_back(UniquePtr<Op>(op.release()));
      }
    } else if (OPTION("-div") || OPTION("-divide-by") || OPTION("-over") || OPTION("/")) {
      const char *arg = ARGUMENT;
      double c;
      if (ToLower(arg) == "sum") {
        Statistic *a = new Sum();
        a->Hidden(verbose < 1);
        ops.push_back(UniquePtr<Op>(a));
        ops.push_back(UniquePtr<Op>(new Div(&a->Value())));
      } else if (FromString(arg, c)) {
        if (fequal(c, .0)) {
          cerr << "Invalid -div argument, value must not be zero!" << endl;
          exit(1);
        }
        ops.push_back(UniquePtr<Op>(new Div(c)));
      } else {
        const char *fname = arg;
        const char *aname = nullptr;
        if (HAS_ARGUMENT) {
          aname = ARGUMENT;
        } else {
          #if MIRTK_Image_WITH_VTK
            if (dataset && arrays->HasArray(fname)) {
              aname = fname;
              fname = input_name;
            }
          #endif
        }
        UniquePtr<Div> op(new Div(fname));
        if (aname) {
          #if MIRTK_Image_WITH_VTK
            op->ArrayName(aname);
            op->IsCellData(cell_data);
          #else
            FatalError("Cannot read point set files when build without VTK or wrong usage!");
          #endif
        }
        ops.push_back(UniquePtr<Op>(op.release()));
      }
    } else if (OPTION("-div-with-zero")) {
      const char *arg = ARGUMENT;
      double c;
      if (ToLower(arg) == "sum") {
        Statistic *a = new Sum();
        a->Hidden(verbose < 1);
        ops.push_back(UniquePtr<Op>(a));
        ops.push_back(UniquePtr<Op>(new DivWithZero(&a->Value())));
      } else if (FromString(arg, c)) {
        ops.push_back(UniquePtr<Op>(new DivWithZero(c)));
      } else {
        const char *fname = arg;
        const char *aname = nullptr;
        if (HAS_ARGUMENT) {
          aname = ARGUMENT;
        } else {
          #if MIRTK_Image_WITH_VTK
            if (dataset && arrays->HasArray(fname)) {
              aname = fname;
              fname = input_name;
            }
          #endif
        }
        UniquePtr<DivWithZero> op(new DivWithZero(fname));
        if (aname) {
          #if MIRTK_Image_WITH_VTK
            op->ArrayName(aname);
            op->IsCellData(cell_data);
          #else
            FatalError("Cannot read point set files when build without VTK or wrong usage!");
          #endif
        }
        ops.push_back(UniquePtr<Op>(op.release()));
      }
    } else if (OPTION("-abs")) {
      ops.push_back(UniquePtr<Op>(new Abs()));
    } else if (OPTION("-pow") || OPTION("-power")) {
      const char *arg = ARGUMENT;
      double exponent;
      if (!FromString(arg, exponent)) {
        cerr << "Invalid -power value, must be a number!" << endl;
        exit(1);
      }
      ops.push_back(UniquePtr<Op>(new Pow(exponent)));
    } else if (OPTION("-sqrt")) {
      ops.push_back(UniquePtr<Op>(new Pow(.5)));
    } else if (OPTION("-square") || OPTION("-sq")) {
      ops.push_back(UniquePtr<Op>(new Pow(2.0)));
    } else if (OPTION("-exp")) {
      ops.push_back(UniquePtr<Op>(new Exp()));
    } else if (OPTION("-log") || OPTION("-log2") || OPTION("-loge") || OPTION("-log10") || OPTION("-lb") || OPTION("-ln") || OPTION("-lg")) {
      a = numeric_limits<double>::min();
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(a);
        if (a <= .0) {
          cerr << "Invalid -log threshold argument, must be a positive number" << endl;
          exit(1);
        }
      }
      Op *op = nullptr;
      if (strcmp(OPTNAME, "-log") == 0) {
        if (HAS_ARGUMENT) {
          double base;
          if (!FromString(ARGUMENT, base)) {
            char c;
            if (!FromString(ARGUMENT, c) || c != 'e') {
              cerr << "Invalid -log base argument, must be a positive number or character e" << endl;
              exit(1);
            }
            op = new Ln(a);
          } else {
            op = new Log(base, a);
          }
        } else {
          op = new Ln(a);
        }
      } else if (strcmp(OPTNAME, "-log2")  == 0 || strcmp(OPTNAME, "-lb") == 0) {
        op = new Lb(a);
      } else if (strcmp(OPTNAME, "-log10") == 0 || strcmp(OPTNAME, "-lg") == 0) {
        op = new Lg(a);
      } else if (strcmp(OPTNAME, "-loge")  == 0 || strcmp(OPTNAME, "-ln") == 0) {
        op = new Ln(a);
      }
      ops.push_back(UniquePtr<Op>(op));
    } else if (OPTION("=")) {
      const char *fname = ARGUMENT;
      #if MIRTK_Image_WITH_VTK
        ops.push_back(UniquePtr<Op>(new Write(fname, datatype, attr, dataset, scalars_name, scalars_name)));
      #else
        ops.push_back(UniquePtr<Op>(new Write(fname, datatype, attr)));
      #endif
    } else if (OPTION("-o") || OPTION("-out") || OPTION("-output")) {
      const char *fname = ARGUMENT;
      int         dtype = datatype;
      #if MIRTK_Image_WITH_VTK
        const char *output_scalars_name = scalars_name;
      #endif
      if (HAS_ARGUMENT) {
        const char *arg = ARGUMENT;
        dtype = ToDataType(arg);
        if (dtype == MIRTK_VOXEL_UNKNOWN) {
          cerr << "Invalid -out data type " << arg << endl;
          exit(1);
        }
        if (HAS_ARGUMENT) {
          #if MIRTK_Image_WITH_VTK
            output_scalars_name = ARGUMENT;
          #else
            Warning("Output scalars array name argument of -output option ignored");
          #endif
        }
      }
      #if MIRTK_Image_WITH_VTK
        ops.push_back(UniquePtr<Op>(new Write(fname, dtype, attr, dataset, scalars_name, output_scalars_name, cell_data)));
      #else
        ops.push_back(UniquePtr<Op>(new Write(fname, dtype, attr)));
      #endif
    // Data statistics
    } else if (OPTION("-median")) {
      ops.push_back(UniquePtr<Op>(new Median()));
    } else if (OPTION("-mean") || OPTION("-average") || OPTION("-avg")) {
      ops.push_back(UniquePtr<Op>(new Mean()));
    } else if (OPTION("-sigma") || OPTION("-stddev") || OPTION("-stdev") || OPTION("-std") || OPTION("-sd")) {
      ops.push_back(UniquePtr<Op>(new StDev()));
    } else if (OPTION("-normal-distribution") ||
               OPTION("-mean+sigma") || OPTION("-mean+stddev") || OPTION("-mean+stdev") || OPTION("-mean+std") || OPTION("-mean+sd") ||
               OPTION("-avg+sigma")  || OPTION("-avg+stddev")  || OPTION("-avg+stdev")  || OPTION("-avg+std")  || OPTION("-avg+sd")) {
      ops.push_back(UniquePtr<Op>(new NormalDistribution()));
    } else if (OPTION("-variance") || OPTION("-var")) {
      ops.push_back(UniquePtr<Op>(new Var()));
    } else if (OPTION("-mean-absolute-difference") || OPTION("-mean-absolute-deviation") || OPTION("-mad") || OPTION("-mad-mean")) {
      ops.push_back(UniquePtr<Op>(new MeanAbsoluteDifference()));
    } else if (OPTION("-median-absolute-difference") || OPTION("-median-absolute-deviation") || OPTION("-mad-median")) {
      ops.push_back(UniquePtr<Op>(new MedianAbsoluteDifference()));
    } else if (OPTION("-minimum")  || OPTION("-min")) {
      ops.push_back(UniquePtr<Op>(new Min()));
    } else if (OPTION("-maximum")  || OPTION("-max")) {
      ops.push_back(UniquePtr<Op>(new Max()));
    } else if (OPTION("-extrema") || OPTION("-minmax")) {
      ops.push_back(UniquePtr<Op>(new Extrema()));
    } else if (OPTION("-range")) {
      ops.push_back(UniquePtr<Op>(new Range()));
    } else if (OPTION("-percentile") || OPTION("-pct") || OPTION("-p")) {
      do {
        int p;
        if (FromString(ARGUMENT, p) && 0 <= p && p <= 100) {
          ops.push_back(UniquePtr<Op>(new Percentile(p)));
        } else {
          cerr << "Invalid -percentile value, must be integer in the range [0, 100]!" << endl;
          exit(1);
        }
      } while (HAS_ARGUMENT);
    } else if (OPTION("-lower-percentile-mean") || OPTION("-lpctavg")) {
      do {
        int p;
        if (FromString(ARGUMENT, p) && 0 <= p && p <= 100) {
          ops.push_back(UniquePtr<Op>(new LowerPercentileMean(p)));
        } else {
          cerr << "Invalid -lower-percentile-mean value, must be integer in the range [0, 100]!" << endl;
          exit(1);
        }
      } while (HAS_ARGUMENT);
    } else if (OPTION("-upper-percentile-mean") || OPTION("-upctavg")) {
      do {
        int p;
        if (FromString(ARGUMENT, p) && 0 <= p && p <= 100) {
          ops.push_back(UniquePtr<Op>(new UpperPercentileMean(p)));
        } else {
          cerr << "Invalid -upper-percentile-mean value, must be integer in the range [0, 100]!" << endl;
          exit(1);
        }
      } while (HAS_ARGUMENT);
    } else if (OPTION("-sum")) {
      ops.push_back(UniquePtr<Op>(new Sum()));
    } else if (OPTION("-count")) {
      ops.push_back(UniquePtr<Op>(new Count()));
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
    ops.push_back(UniquePtr<Statistic>(new Mean()));
    ops.push_back(UniquePtr<Statistic>(new StDev()));
    ops.push_back(UniquePtr<Statistic>(new Extrema()));
    ops.push_back(UniquePtr<Statistic>(new Range()));
  }

  // Initial data mask
  UniquePtr<bool[]> mask(new bool[n]);
  for (int i = 0; i < n; ++i) {
    if (IsNaN(data[i])) {
      mask[i] = false;
    } else {
      mask[i] = true;
    }
  }

  // Process input data, either transform it or compute statistics from it
  for (size_t i = 0; i < ops.size(); ++i) {
    ops[i]->Process(n, data.get(), mask.get());
  }
  mask.reset();

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
      FatalError("Cannot append to file " << append_name);
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
      if (stat != nullptr && !stat->Hidden()) {
        for (size_t j = 0; j < stat->Names().size(); ++j, ++c) {
          if (c > 0) out << delimiter;
          if (c < header.size()) out << header[c];
          else out << stat->Names()[j];
        }
      }
    }
    out << endl;
  }

  // Print image statistics
  if (delimiter) {
    for (size_t i = 0; i < prefix.size(); ++i) {
      if (i > 0) out << delimiter;
      out << prefix[i];
    }
    bool first = prefix.empty();
    for (size_t i = 0; i < ops.size(); ++i) {
      Statistic *stat = dynamic_cast<Statistic *>(ops[i].get());
      if (stat != nullptr && !stat->Hidden() && !stat->Names().empty()) {
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
      if (stat != nullptr && !stat->Hidden()) {
        stat->Print(out, digits, prefix_string.c_str());
      }
    }
  }

  ofs.close();
  return 0;
}

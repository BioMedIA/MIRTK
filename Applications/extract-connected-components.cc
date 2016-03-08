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
#include <mirtkConnectedComponents.h>

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <input> <output> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Extracts connected components from input segmentation label image.\n";
  cout << "  By default, the largest connected component is extracted.\n";
  cout << "  In case of :option:`-output-component-labels`, the labels\n";
  cout << "  of the components are written instead.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Input segmentation image.\n";
  cout << "  output   Extracted components or components label image.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -m <m>\n";
  cout << "      Zero-based index of first connected component to extract. (default: 0)\n";
  cout << "\n";
  cout << "  -n <n>\n";
  cout << "      Extract (at most) n components. (default: 1)\n";
  cout << "\n";
  cout << "  -all\n";
  cout << "      Extract all components.\n";
  cout << "\n";
  cout << "  -connectivity <num>\n";
  cout << "      Type of voxel connectivity (4, 6, 18, or 26). (default: 26)\n";
  cout << "\n";
  cout << "  -ordering (none|largest|smallest)\n";
  cout << "      Ordering of connected component labels based on component size:\n";
  cout << "      - ``none``:     Label components in no specific order.\n";
  cout << "      - ``largest``:  Label largest component first. (default)\n";
  cout << "      - ``smallest``: Label smallest component first.\n";
  cout << "\n";
  cout << "  -output-component-labels\n";
  cout << "      Write component labels instead of input labels.\n";
  cout << "\n";
  cout << "  -output-component-mask\n";
  cout << "      Write binary mask with non-zero value at extracted components.\n";
  PrintStandardOptions(cout);
  cout << "\n";
  cout.flush();
}

// =============================================================================
// Main
// =============================================================================

enum OutputType
{
  ComponentLabels,
  ComponentMask,
  InputLabels
};

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  EXPECTS_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  int                         m = 0, n = 1;
  ConnectedComponentsOrdering ordering     = CC_LargestFirst;
  ConnectivityType            connectivity = CONNECTIVITY_26;
  OutputType                  output_type  = InputLabels;

  for (ALL_OPTIONS) {
    if      (OPTION("-m"))            PARSE_ARGUMENT(m);
    else if (OPTION("-n"))            PARSE_ARGUMENT(n);
    else if (OPTION("-all"))          n = 0;
    else if (OPTION("-connectivity")) PARSE_ARGUMENT(connectivity);
    else if (OPTION("-ordering"))     PARSE_ARGUMENT(ordering);
    else if (OPTION("-output-component-labels")) output_type = ComponentLabels;
    else if (OPTION("-output-component-mask"))   output_type = ComponentMask;
    else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
  }

  InitializeImageIOLibrary();
 
  GreyImage labels(input_name);
  GreyImage components;

  if (verbose) cout << "Labeling connected components...", cout.flush();
  ConnectedComponents<GreyPixel> cc(ordering, connectivity);
  cc.Input(&labels);
  cc.Output(&components);
  cc.Run();
  if (verbose) cout << " done" << endl;

  const int nregions = cc.NumberOfComponents();

  if (verbose) {
    int nvoxels = 0;
    for (int idx = 0; idx < labels.NumberOfVoxels(); ++idx) {
      if (labels(idx) != 0) ++nvoxels;
    }
    cout << "No. of labeled voxels:       " << nvoxels << "\n";
    cout << "No. of connected components: " << nregions << "\n";
    if (verbose > 1) {
      for (int i = 0; i < nregions; ++i) {
        cout << "  Size of component " << ToString(i+1, 3) << ": "
             << ToString(cc.ComponentSize()[i], 6) << "\n";
      }
    }
    cout.flush();
  }

  if (m >= nregions) {
    FatalError("Start index -m is greater or equal than the total no. of components!");
  }
  if (n <= 0) n = nregions;
  if (m + n > nregions) n = nregions - m;

  const GreyPixel min_label = m + 1; // component labels are 1-based
  const GreyPixel max_label = min_label + n - 1;

  GreyImage output(labels.Attributes());
  GreyPixel component_label;

  if (verbose) {
    switch (output_type) {
      case ComponentLabels: cout << "Relabeling"; break;
      case ComponentMask:   cout << "Masking";    break;
      case InputLabels:     cout << "Extracting"; break;
    }
    cout << " " << n << " components starting with component " << m << "...";
    cout.flush();
  }
  for (int idx = 0; idx < output.NumberOfVoxels(); ++idx) {
    component_label = cc.Output()->Get(idx);
    if (min_label <= component_label && component_label <= max_label) {
      switch (output_type) {
        case ComponentLabels: output(idx) = component_label - min_label + 1; break;
        case ComponentMask:   output(idx) = 1; break;
        case InputLabels:     output(idx) = labels(idx); break;
      }
    }
  }
  if (verbose) cout << " done" << endl;

  if (verbose) cout << "Writing output to file " << output_name << "...", cout.flush();
  output.Write(output_name);
  if (verbose) cout << " done" << endl;

  return 0;
}


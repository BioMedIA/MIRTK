/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#include "mirtk/PointSetIO.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/EdgeConnectivity.h"
#include "mirtk/ClosePointData.h"


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
  cout << "  Closes scalar point data of an input point set by perfoming a dilation\n";
  cout << "  followed by the same number of erosions. When the input data array has more\n";
  cout << "  than one component, each component is processed separately.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Input  point set.\n";
  cout << "  output   Output point set with closed scalars.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -a, -array, -scalars <name>\n";
  cout << "      Name of point data array. (default: active scalars)\n";
  cout << "  -c, -connectivity, -neighborhood <int>\n";
  cout << "      Maximum node edge-connectivity of neighbors. (default: 1)\n";
  cout << "  -r, -radius <float>\n";
  cout << "      Maximum point distance of neighboring points. Used instead of :option:`-c` when positive. (default: 0)\n";
  cout << "  -n, -iterations, -iter <n>\n";
  cout << "      Number of iterations. (default: 1)\n";
  cout << "\n";
  cout << "Output format options:\n";
  cout << "  -ascii, -nobinary" << endl;
  cout << "      Write legacy VTK in ASCII format. (default: input)" << endl;
  cout << "  -binary, -noascii" << endl;
  cout << "      Write legacy VTK in binary format. (default: input)" << endl;
  cout << "  -[no]compress" << endl;
  cout << "      Write XML VTK file with or without compression. (default: on)" << endl;
  PrintStandardOptions(cout);
  cout << "\n";
  cout.flush();
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  EXPECTS_POSARGS(2);

  const char *input_name        = POSARG(1);
  const char *output_name       = POSARG(2);
  const char *array_name        = nullptr;
  const char *output_array_name = nullptr;

  int iterations    = 1;
  int connectivity  = 1;
  double     radius = 0.;
  FileOption fopt   = FO_Default;

  for (ALL_OPTIONS) {
    if (OPTION("-a") || OPTION("-array") || OPTION("-scalars") || OPTION("-name")) {
      array_name = ARGUMENT;
    }
    else if (OPTION("-o") || OPTION("-output-array") || OPTION("-output-scalars") || OPTION("-output-name")) {
      output_array_name = ARGUMENT;
    }
    else if (OPTION("-n") || OPTION("-iterations") || OPTION("-iter")) {
      PARSE_ARGUMENT(iterations);
    }
    else if (OPTION("-c") || OPTION("-connectivity") || OPTION("-neighborhood") || OPTION("-neighbourhood")) {
      PARSE_ARGUMENT(connectivity);
    }
    else if (OPTION("-r") || OPTION("-radius")) {
      PARSE_ARGUMENT(radius);
    }
    else HANDLE_POINTSETIO_OPTION(fopt);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (verbose) cout << "Reading " << input_name << "...", cout.flush();
  vtkSmartPointer<vtkPolyData> input = ReadPolyData(input_name);
  if (verbose) cout << " done" << endl;

  SharedPtr<EdgeTable> edges(new EdgeTable(input));
  SharedPtr<EdgeConnectivity> neighbors;
  if (radius > 0.) {
    neighbors = NewShared<EdgeConnectivity>(input, radius, edges.get());
  } else {
    neighbors = NewShared<EdgeConnectivity>(input, connectivity, edges.get());
  }
  if (verbose) cout << "Closing with c=" << neighbors->Maximum() << "...", cout.flush();
  ClosePointData filter;
  filter.Input(input);
  filter.DataName(array_name);
  filter.EdgeTable(edges);
  filter.Neighbors(neighbors);
  filter.Connectivity(connectivity);
  filter.Radius(radius);
  filter.Iterations(iterations);
  filter.Run();
  if (output_array_name && (!array_name || strcmp(array_name, output_array_name) != 0)) {
    filter.OutputData()->SetName(output_array_name);
    if (filter.InputData()->GetName()) {
      filter.Output()->GetPointData()->AddArray(filter.InputData());
    }
  }
  if (verbose) cout << " done" << endl;

  if (!WritePolyData(output_name, filter.Output(), fopt)) {
    FatalError("Failed to write result to " << output_name << "!");
  }

  return 0;
}

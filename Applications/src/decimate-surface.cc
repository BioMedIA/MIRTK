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

#include "mirtk/PointSetIO.h"
#include "mirtk/Vtk.h"

#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkTriangleFilter.h"
#include "vtkDecimatePro.h"
#include "vtkSmartPointer.h"

#include "vtkQuadricDecimation.h"

using namespace mirtk;


// Global filter objects to enable output of default settings in help
vtkSmartPointer<vtkQuadricDecimation> quadric   = vtkSmartPointer<vtkQuadricDecimation>::New();
vtkSmartPointer<vtkDecimatePro>       decimater = vtkSmartPointer<vtkDecimatePro>      ::New();

// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <input> <output> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Decimates a (triangular) mesh using VTK's ``vtkQuadricDecimation`` filter." << endl;
  cout << "  In case of :option:`-pro`, the ``vtkDecimatePro`` filter is used instead." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -reduceby <percentage>       Target reduction value in percentage. (default: " << (100.0 * quadric->GetTargetReduction()) << ")" << endl;
  cout << "  -pro                         Use vtkDecimatePro filter. (default: off)" << endl;
  cout << "  -[no]compress                Whether to compress output .vtp file. (default: on)" << endl;
  cout << "  -binary                      Write binary data when output file name extension is .vtk. (default: on)" << endl;
  cout << "  -ascii                       Write ASCII  data when output file name extension is .vtk. (default: off)" << endl;
  cout << endl;
  cout << "Options for vtkDecimatePro filter (see :option:`-pro`):" << endl;
  cout << "  -preservetopology [on|off]   Whether to preserve topology. (default: " << (decimater->GetPreserveTopology() ? "on" : "off") << ")" << endl;
  cout << "  -splitangle <angle>          Specify the mesh splitting angle. A negative value turns splitting off. (default: "
                                          << (decimater->GetSplitting() ? decimater->GetSplitAngle() : -1.0) << ")" << endl;
  cout << "  -featureangle <angle>        Specify the mesh feature angle. (default: " << decimater->GetFeatureAngle() << ")" << endl;
  cout << "  -maxerror <float>            Specify the largest decimation error that is allowed. (default: ";
  if      (decimater->GetErrorIsAbsolute()) cout << "-abserror";
  else if (decimater->GetAccumulateError()) cout << "-accerror";
  else                                      cout << decimater->GetMaximumError();
  cout << ")" << endl;
  cout << "  -abserror <float>            Specify absolute decimation error. (default: ";
  if      (decimater->GetErrorIsAbsolute()) cout << decimater->GetAbsoluteError();
  else if (decimater->GetAccumulateError()) cout << ":option:`-accerror`";
  else                                      cout << ":option:`-maxerror`";
  cout << ")" << endl;
  cout << "  -accerror [on|off]           Accumulate decimation error. (default: " << (decimater->GetAccumulateError() ? "on" : "off") << ")" << endl;
  cout << "  -boundarydeletion [on|off]   Whether to allow deletion of vertices on the boundary of the mesh."
                                          << " (default: " << (decimater->GetBoundaryVertexDeletion() ? "on" : "off") << ")" << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
// Auxiliary function to convert string to boolean value
bool atob(const char *value)
{
  bool b = false;
  FromString(value, b);
  return b;
}

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  // Set/Check number of required positional arguments
  EXPECTS_POSARGS(2);

  quadric->AttributeErrorMetricOn(); // Otherwise point data is not copied

  // Change decimater settings to user specified values
  FileOption            fopt   = FO_Default;
  vtkPolyDataAlgorithm *filter = quadric;

  for (ALL_OPTIONS) {
    if (OPTION("-reduceby") || OPTION("-targetreduction")) {
      double ratio = atof(ARGUMENT) / 100.0;
      decimater->SetTargetReduction(ratio);
      quadric  ->SetTargetReduction(ratio);
    }
    else if (OPTION("-maxerror") || OPTION("-maximumerror")) {
      decimater->SetMaximumError(atof(ARGUMENT));
      decimater->SetErrorIsAbsolute(false);
    }
    else if (OPTION("-abserror") || OPTION("-absoluteerror")) {
      decimater->SetAbsoluteError(atof(ARGUMENT));
      decimater->SetErrorIsAbsolute(true);
    }
    else if (OPTION("-accerror") || OPTION("-accumulateerror")) {
      if (HAS_ARGUMENT) decimater->SetAccumulateError(atob(ARGUMENT));
      else              decimater->SetAccumulateError(true);
    }
    else if (OPTION("-preservetopology")) {
      if (HAS_ARGUMENT) decimater->SetPreserveTopology(atob(ARGUMENT));
      else              decimater->SetPreserveTopology(true);
    }
    else if (OPTION("-featureangle")) decimater->SetFeatureAngle(atof(ARGUMENT));
    else if (OPTION("-splitangle")) {
      double arg = atof(ARGUMENT);
      if (arg < .0) decimater->SplittingOff();
      else          decimater->SplittingOn(), decimater->SetSplitAngle(arg);
    }
    else if (OPTION("-boundarydeletion") || OPTION("-boundaryvertexdeletion")) {
      if (HAS_ARGUMENT) decimater->SetBoundaryVertexDeletion(atob(ARGUMENT));
      else              decimater->SetBoundaryVertexDeletion(true);
    }
    else if (OPTION("-pro")) {
      filter = decimater;
    }
    else HANDLE_POINTSETIO_OPTION(fopt);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read input mesh
  vtkSmartPointer<vtkPolyData> input = ReadPolyData(POSARG(1), fopt);

  // Execute pipeline
  vtkSmartPointer<vtkTriangleFilter> triangulator;
  triangulator = vtkSmartPointer<vtkTriangleFilter>::New();

  SetVTKInput(triangulator, input);
  SetVTKConnection(filter, triangulator);
  filter->Update();

  vtkSmartPointer<vtkPolyData> output = filter->GetOutput();

  // Write output mesh
  WritePolyData(POSARG(2), output, fopt);

  // Some verbose reporting
  if (verbose) {
    cout << "Number of output vertices: " << output->GetNumberOfPoints();
    cout << " (" << setw(3) << round(100.0 * output->GetNumberOfPoints()
                                     / input->GetNumberOfPoints()) << "%)" << endl;
    cout << "Number of output cells:    " << output->GetNumberOfCells();
    cout << " (" << setw(3) << round(100.0 * output->GetNumberOfCells()
                                     / input->GetNumberOfCells()) << "%)" << endl;
  }

  exit(0);
}

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

#include "mirtk/EdgeConnectivity.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"

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
  cout << "  Dilates scalar point data of an input point set by replacing a point's\n";
  cout << "  value by the maximum of the values of its neighboring points. When the\n";
  cout << "  input data array has more than one component, each component is processed\n";
  cout << "  separately.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Input  point set.\n";
  cout << "  output   Output point set with dilated scalars.\n";
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
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
struct DilateScalars
{
  vtkDataArray           *_Input;
  vtkDataArray           *_Output;
  const EdgeConnectivity *_Neighbors;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    int        nbrPts;
    const int *nbrIds;
    double     value;

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      for (int j = 0; j < _Input->GetNumberOfComponents(); ++j) {
        value = _Input->GetComponent(ptId, j);
        _Neighbors->GetConnectedPoints(ptId, nbrPts, nbrIds);
        for (int i = 0; i < nbrPts; ++i) {
          value = max(value, _Input->GetComponent(nbrIds[i], j));
        }
        _Output->SetComponent(ptId, j, value);
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Dilate point data
///
/// \todo Implement this as MeshFilter as part of the PointSet module.
vtkSmartPointer<vtkPointSet> Dilate(vtkSmartPointer<vtkPointSet> input,
                                    const EdgeConnectivity &neighbors,
                                    vtkSmartPointer<vtkDataArray> arr, int niter = 1)
{
  const int npoints = static_cast<int>(input->GetNumberOfPoints());

  vtkSmartPointer<vtkPointSet> output;
  output.TakeReference(input->NewInstance());
  output->ShallowCopy(input);

  int attr = -1;
  vtkPointData * const pd = output->GetPointData();
  for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
    if (pd->GetArray(i) == arr) {
      attr = pd->IsArrayAnAttribute(i);
      pd->RemoveArray(i);
      break;
    }
  }

  // Attention: vtkDataArray::DeepCopy (of older VTK versions?) does not
  //            copy array name (and possibly not array component names).
  vtkSmartPointer<vtkDataArray> res;
  res.TakeReference(arr->NewInstance());
  res->SetNumberOfComponents(arr->GetNumberOfComponents());
  res->SetNumberOfTuples(arr->GetNumberOfTuples());
  res->SetName(arr->GetName());
  for (int j = 0; j < arr->GetNumberOfComponents(); ++j) {
    res->SetComponentName(j, arr->GetComponentName(j));
  }

  DilateScalars body;
  body._Input     = arr;
  body._Output    = res;
  body._Neighbors = &neighbors;
  for (int iter = 0; iter < niter; ++iter) {
    if (iter == 1) {
      arr.TakeReference(res->NewInstance());
      arr->SetNumberOfComponents(res->GetNumberOfComponents());
      arr->SetNumberOfTuples(res->GetNumberOfTuples());
      arr->SetName(res->GetName());
      for (int j = 0; j < res->GetNumberOfComponents(); ++j) {
        arr->SetComponentName(j, res->GetComponentName(j));
        arr->CopyComponent(j, res, j);
      }
      body._Input = arr;
    } else if (iter > 1) {
      swap(body._Input, body._Output);
    }
    parallel_for(blocked_range<int>(0, npoints), body);
  }

  const int idx = pd->AddArray(body._Output);
  if (attr >= 0) pd->SetActiveAttribute(idx, attr);

  return output;
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
  vtkSmartPointer<vtkPointSet> pset = ReadPointSet(input_name);
  if (verbose) cout << " done" << endl;

  vtkSmartPointer<vtkDataArray> arr;
  if (array_name) {
    arr = GetArrayByCaseInsensitiveName(pset->GetPointData(), array_name);
    if (arr == nullptr) {
      FatalError("Input point set has not point data array named " << array_name << "!");
    }
  } else {
    arr = pset->GetPointData()->GetScalars();
    if (arr == nullptr) {
      FatalError("Input point set has no active scalars, use -array, -scalars option!");
    }
  }
  if (output_array_name) {
    if (arr->GetName()) {
      if (strcmp(output_array_name, arr->GetName()) != 0) {
        vtkSmartPointer<vtkDataArray> in = arr;
        arr.TakeReference(in->NewInstance());
        arr->DeepCopy(in);
        arr->SetName(output_array_name);
        for (int j = 0; j < in->GetNumberOfComponents(); ++j) {
          arr->SetComponentName(j, in->GetComponentName(j));
        }
        pset->GetPointData()->AddArray(arr);
      }
    } else {
      arr->SetName(output_array_name);
    }
  }

  EdgeConnectivity neighbors;
  if (radius > 0.) {
    neighbors.Initialize(pset, radius);
  } else {
    neighbors.Initialize(pset, connectivity);
  }

  if (verbose) cout << "Dilating with c=" << neighbors.Maximum() << "...", cout.flush();
  pset = Dilate(pset, neighbors, arr, iterations);
  arr  = nullptr;
  if (verbose) cout << " done" << endl;

  if (!WritePointSet(output_name, pset, fopt)) {
    FatalError("Failed to write result to " << output_name << "!");
  }

  return 0;
}

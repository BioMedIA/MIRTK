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

#include "mirtk/DataSelection.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/MeshSmoothing.h"

#include "mirtk/Vtk.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkPoints.h"
#include "vtkPolyDataNormals.h"

using namespace mirtk;
using namespace mirtk::data;
using namespace mirtk::data::select;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <input1> <input2> <output> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Blends two surface meshes with index based one-to-one point correspondences\n";
  cout << "  by taking the point coordinates for selected points from the first surface\n";
  cout << "  mesh and all others from the second surface mesh. Points at the boundary of\n";
  cout << "  the point selection mask are blended between the two mesh positions of the\n";
  cout << "  corresponding points. For other points, i.e., those further away from the\n";
  cout << "  boundary, the blending factor is either 0 or 1.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input1   First  input surface.\n";
  cout << "  input2   Second input surface.\n";
  cout << "  output   Output surface.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
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
UnorderedSet<int> RegionBoundaryIds(vtkDataArray *labels, const EdgeTable &edgeTable)
{
  int        adjPts;
  const int *adjIds;
  double     label;
  UnorderedSet<int> ptIds;
  for (int ptId = 0; ptId < edgeTable.NumberOfPoints(); ++ptId) {
    label = labels->GetComponent(ptId, 0);
    edgeTable.GetAdjacentPoints(ptId, adjPts, adjIds);
    for (int i = 0; i < adjPts; ++i) {
      if (labels->GetComponent(adjIds[i], 0) != label) {
        ptIds.insert(ptId);
        break;
      }
    }
  }
  return ptIds;
}

// -----------------------------------------------------------------------------
void DilateLabel(vtkDataArray *labels, const EdgeTable &edgeTable, double label)
{
  int        adjPts;
  const int *adjIds;
  vtkSmartPointer<vtkDataArray> output;
  output.TakeReference(labels->NewInstance());
  output->DeepCopy(labels);
  if (IsNaN(label)) {
    for (int ptId = 0; ptId < edgeTable.NumberOfPoints(); ++ptId) {
      edgeTable.GetAdjacentPoints(ptId, adjPts, adjIds);
      for (int i = 0; i < adjPts; ++i) {
        if (IsNaN(labels->GetComponent(adjIds[i], 0))) {
          output->SetComponent(ptId, 0, label);
          break;
        }
      }
    }
  } else {
    for (int ptId = 0; ptId < edgeTable.NumberOfPoints(); ++ptId) {
      edgeTable.GetAdjacentPoints(ptId, adjPts, adjIds);
      for (int i = 0; i < adjPts; ++i) {
        if (labels->GetComponent(adjIds[i], 0) == label) {
          output->SetComponent(ptId, 0, label);
          break;
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
/// Average adjacent non-invalid values at invalid data points
double ApplyLaplaceOperator(vtkDataArray *input, vtkDataArray *output,
                            const EdgeTable &edgeTable, vtkDataArray *mask,
                            double invalid, int *num_invalid = nullptr)
{
  double     val, label, delta = 0.;
  int        num, adjPts, n_invalid = 0;
  const int *adjIds;

  if (IsNaN(invalid)) {
    for (int ptId = 0; ptId < edgeTable.NumberOfPoints(); ++ptId) {
      if (mask->GetComponent(ptId, 0) == 0.) {
        output->SetComponent(ptId, 0, input->GetComponent(ptId, 0));
      } else {
        val = 0., num = 0;
        edgeTable.GetAdjacentPoints(ptId, adjPts, adjIds);
        for (int i = 0; i < adjPts; ++i) {
          label = input->GetComponent(adjIds[i], 0);
          if (!IsNaN(label)) {
            val += label, ++num;
          }
        }
        if (num == 0) {
          output->SetComponent(ptId, 0, invalid);
          ++n_invalid;
        } else {
          val /= num;
          label = input->GetComponent(ptId, 0);
          if (IsNaN(label)) label = 0.;
          delta = max(delta, abs(label - val));
          output->SetComponent(ptId, 0, val);
        }
      }
    }
  } else {
    for (int ptId = 0; ptId < edgeTable.NumberOfPoints(); ++ptId) {
      if (mask->GetComponent(ptId, 0) == 0.) {
        output->SetComponent(ptId, 0, input->GetComponent(ptId, 0));
      } else {
        val = 0., num = 0;
        edgeTable.GetAdjacentPoints(ptId, adjPts, adjIds);
        for (int i = 0; i < adjPts; ++i) {
          label = input->GetComponent(adjIds[i], 0);
          if (label != invalid) {
            val += label, ++num;
          }
        }
        if (num == 0) {
          output->SetComponent(ptId, 0, invalid);
          ++n_invalid;
        } else {
          val /= num;
          label = input->GetComponent(ptId, 0);
          delta = max(delta, abs(label - val));
          output->SetComponent(ptId, 0, val);
        }
      }
    }
  }

  if (num_invalid) *num_invalid = n_invalid;
  return delta;
}

// -----------------------------------------------------------------------------
void LaplacianInterpolation(vtkDataArray *input, vtkDataArray *mask,
                            const EdgeTable &edgeTable, double invalid,
                            int max_iter, double max_delta)
{
  int n_invalid = 0;
  double delta  = inf;
  vtkSmartPointer<vtkDataArray> buffer1 = input, buffer2;
  buffer2.TakeReference(buffer1->NewInstance());
  buffer2->SetNumberOfComponents(buffer1->GetNumberOfComponents());
  buffer2->SetNumberOfTuples(buffer1->GetNumberOfTuples());
  for (int iter = 0; n_invalid > 0 || (iter < max_iter && delta > max_delta); ++iter) {
    delta = ApplyLaplaceOperator(buffer1, buffer2, edgeTable, mask, invalid, &n_invalid);
    swap(buffer1, buffer2);
  }
  if (buffer1 != input) input->CopyComponent(0, buffer1, 0);
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  EXPECTS_POSARGS(3);

  const char *input1_name = POSARG(1);
  const char *input2_name = POSARG(2);
  const char *output_name = POSARG(3);

  double      value;
  const char *alpha_name      = nullptr;
  const char *scalars_name    = nullptr;
  const char *source_name     = nullptr;
  int         comp_index      = 0;
  int         radius          = 0;
  int         max_iter        = 0;
  double      max_delta       = 1e-6;
  int         update_pnormals = -1;
  int         update_cnormals = -1;
  FileOption  fopt            = FO_Default;

  SharedPtr<LogicalOp> selector(new LogicalAnd());

  for (ALL_OPTIONS) {
    if (OPTION("-a") || OPTION("-array") || OPTION("-scalars") || OPTION("-where") || OPTION("-select")) {
      scalars_name = ARGUMENT;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(comp_index);
    }
    else if (OPTION("-source-array")) {
      source_name = ARGUMENT;
    }
    else if (OPTION("-alpha-array")) {
      alpha_name = ARGUMENT;
    }
    else if (OPTION("-transition-extent")) {
      PARSE_ARGUMENT(radius);
    }
    else if (OPTION("-max-iter") || OPTION("-smooth-iter") || OPTION("-smooth-iterations")) {
      PARSE_ARGUMENT(max_iter);
    }
    else if (OPTION("-max-delta")) {
      PARSE_ARGUMENT(max_delta);
    }
    else if (OPTION("-and")) {
      if (dynamic_cast<LogicalAnd *>(selector.get()) == nullptr) {
        SharedPtr<LogicalAnd> op(new LogicalAnd());
        if (selector->NumberOfCriteria() > 1) {
          op->Push(selector);
        } else if (selector->NumberOfCriteria() == 1) {
          op->Push(selector->Criterium(0));
        }
        selector = op;
      }
    }
    else if (OPTION("-or")) {
      if (dynamic_cast<LogicalOr *>(selector.get()) == nullptr) {
        SharedPtr<LogicalOr> op(new LogicalOr());
        if (selector->NumberOfCriteria() > 1) {
          op->Push(selector);
        } else if (selector->NumberOfCriteria() == 1) {
          op->Push(selector->Criterium(0));
        }
        selector = op;
      }
    }
    else if (OPTION("-eq")) {
      PARSE_ARGUMENT(value);
      selector->Push(NewShared<Equal>(value));
    }
    else if (OPTION("-ne")) {
      PARSE_ARGUMENT(value);
      selector->Push(NewShared<NotEqual>(value));
    }
    else if (OPTION("-lt")) {
      PARSE_ARGUMENT(value);
      selector->Push(NewShared<LessThan>(value));
    }
    else if (OPTION("-le")) {
      PARSE_ARGUMENT(value);
      selector->Push(NewShared<LessOrEqual>(value));
    }
    else if (OPTION("-gt")) {
      PARSE_ARGUMENT(value);
      selector->Push(NewShared<GreaterThan>(value));
    }
    else if (OPTION("-ge")) {
      PARSE_ARGUMENT(value);
      selector->Push(NewShared<GreaterOrEqual>(value));
    }
    else HANDLE_POINTSETIO_OPTION(fopt);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read input surfaces
  if (verbose) cout << "Reading " << input1_name << "...", cout.flush();
  vtkSmartPointer<vtkPolyData> input1 = ReadPolyData(input1_name);
  if (verbose) cout << " done" << endl;

  if (verbose) cout << "Reading " << input2_name << "...", cout.flush();
  vtkSmartPointer<vtkPolyData> input2 = ReadPolyData(input2_name);
  if (verbose) cout << " done" << endl;

  const vtkIdType npoints = input1->GetNumberOfPoints();
  if (input2->GetNumberOfPoints() != npoints) {
    FatalError("Input surface meshes must have equal number of points!");
  }

  // Select which points of second surface to blend to
  vtkPointData * const pd = input2->GetPointData();
  vtkDataArray *scalars = pd->GetScalars();
  if (scalars_name) {
    scalars = GetArrayByCaseInsensitiveName(pd, scalars_name);
    if (scalars == nullptr) {
      FatalError("Second input point set has no point data array named: " << scalars_name);
    }
  } else if (!scalars && !alpha_name) {
    FatalError("Second input point set has no SCALARS point data array! Use -where option.");
  }

  Selection selection;
  vtkSmartPointer<vtkPolyData> output;
  vtkSmartPointer<vtkDataArray> mask;

  if (!alpha_name) {
    if (comp_index < 0 || comp_index >= scalars->GetNumberOfComponents()) {
      FatalError("Point data array has only " << scalars->GetNumberOfComponents()
                 << " component" << (scalars->GetNumberOfComponents() == 1 ? "" : "s") << "!");
    }
    Array<double> values(npoints);
    for (vtkIdType ptId = 0; ptId < npoints; ++ptId) {
      values[ptId] = scalars->GetComponent(ptId, comp_index);
    }
    selection = selector->Evaluate(values);
    values.clear();
    if (selection.empty()) {
      output = input1;
      if (source_name) {
        mask = NewVtkDataArray(VTK_UNSIGNED_CHAR, npoints, 1, source_name);
        mask->FillComponent(0, 0.);
        output->GetPointData()->AddArray(mask);
      }
    } else if (selection.size() == static_cast<size_t>(npoints)) {
      output = input2;
      if (source_name) {
        mask = NewVtkDataArray(VTK_UNSIGNED_CHAR, npoints, 1, source_name);
        mask->FillComponent(0, 1.);
        output->GetPointData()->AddArray(mask);
      }
    }
  }

  if (!output) {
    // Initialize output surface
    vtkSmartPointer<vtkPoints> points;
    points = vtkSmartPointer<vtkPoints>::New();
    points->SetDataType(input1->GetPoints()->GetDataType());
    points->DeepCopy(input1->GetPoints());

    output.TakeReference(input1->NewInstance());
    output->ShallowCopy(input1);
    output->SetPoints(points);

    // Calculate blending weights
    vtkSmartPointer<vtkDataArray> alpha;
    if (alpha_name) {
      alpha = GetArrayByCaseInsensitiveName(input2->GetPointData(), alpha_name);
      if (alpha == nullptr) {
        FatalError("Second input surface has no point data array named: " << alpha_name);
      }
    } else {
      alpha = NewVtkDataArray(VTK_FLOAT, npoints, 1, "BlendAlpha");
      for (vtkIdType ptId = 0; ptId < npoints; ++ptId) {
        alpha->SetComponent(ptId, 0, selection.find(ptId) == selection.end() ? 0. : 1.);
      }
      output->GetPointData()->RemoveArray("BlendAlpha");
      output->GetPointData()->AddArray(alpha);

      output->BuildLinks();
      SharedPtr<EdgeTable> edgeTable(new EdgeTable(output));
      if (radius > 0) {
        mask = NewVtkDataArray(VTK_UNSIGNED_CHAR, npoints, 1, "BlendMask");
        mask->FillComponent(0, 0.);
        const double invalid = NaN;
        for (auto ptId : RegionBoundaryIds(alpha, *edgeTable)) {
          mask->SetComponent(ptId, 0, 1.);
        }
        for (int iter = 1; iter < radius; ++iter) {
          DilateLabel(mask, *edgeTable, 1.);
        }
        if (max_iter <= 0) max_iter = 100;
        LaplacianInterpolation(alpha, mask, *edgeTable, invalid, max_iter, max_delta);
      } else {
        if (max_iter <= 0) max_iter = 10;
        MeshSmoothing smoother;
        smoother.Input(output);
        smoother.Mask(alpha);
        smoother.EdgeTable(edgeTable);
        smoother.SmoothArray(alpha->GetName());
        smoother.AdjacentValuesOnlyOn();
        smoother.Weighting(MeshSmoothing::Combinatorial);
        smoother.Lambda(1.);
        smoother.NumberOfIterations(max_iter);
        smoother.Run();
        alpha->CopyComponent(0, smoother.Output()->GetPointData()->GetArray(alpha->GetName()), 0);
      }

      if (debug <= 0) {
        output->GetPointData()->RemoveArray("BlendAlpha");
      }
    }

    // Blend between point positions
    Point  p1, p2;
    double w1, w2;
    for (vtkIdType ptId = 0; ptId < npoints; ++ptId) {
      w2 = clamp(alpha->GetComponent(ptId, 0), 0., 1.);
      w1 = 1. - w2;
      input1->GetPoint(ptId, p1);
      input2->GetPoint(ptId, p2);
      points->SetPoint(ptId, (p1 *= w1) += (p2 *= w2));
    }
    if (source_name) {
      if (!mask) {
        mask = NewVtkDataArray(VTK_UNSIGNED_CHAR, npoints, 1, source_name);
      } else {
        mask->SetName(source_name);
      }
      for (vtkIdType ptId = 0; ptId < npoints; ++ptId) {
        mask->SetComponent(ptId, 0, alpha->GetComponent(ptId, 0) <= .5 ? 0. : 1.);
      }
      output->GetPointData()->AddArray(mask);
    }
  }

  // Recompute normals
  if (update_pnormals == -1) {
    update_pnormals = (output->GetPointData()->GetNormals() == nullptr ? 0 : 1);
  }
  if (update_cnormals == -1) {
    update_cnormals = (output->GetCellData ()->GetNormals() == nullptr ? 0 : 1);
  }
  if (update_pnormals || update_cnormals) {
    if (verbose) cout << "\nRecalculating normals...", cout.flush();
    vtkNew<vtkPolyDataNormals> filter;
    filter->SplittingOff();
    filter->NonManifoldTraversalOff();
    filter->ConsistencyOff();
    filter->AutoOrientNormalsOff();
    filter->SetComputePointNormals(update_pnormals);
    filter->SetComputeCellNormals(update_cnormals);
    SetVTKInput(filter, output);
    filter->Update();
    output = filter->GetOutput();
    if (verbose) cout << " done" << endl;
  }

  // Write output surface
  if (!WritePolyData(output_name, output, fopt)) {
    FatalError("Failed to write output surface to file " << output_name);
  }

  return 0;
}

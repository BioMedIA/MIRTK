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
#include "mirtk/PointSetUtils.h"
#include "mirtk/MeshSmoothing.h"

#include "mirtk/Vtk.h"
#include "mirtk/VtkMath.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkPolyDataNormals.h"
#include "vtkImplicitModeller.h"
#include "vtkContourFilter.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkCellLocator.h"
#include "vtkQuadricDecimation.h"
#include "vtkCellDataToPointData.h"

using namespace mirtk;


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
  cout << "  Displaces surface mesh points by a given amount along the surface normal" << endl;
  cout << "  to create an offset surface mesh. To prevent self-intersections, a dense" << endl;
  cout << "  offset surface can optionally first be sampled from an implicit offset" << endl;
  cout << "  surface model onto which each point of the input surface is projected." << endl;
  cout << endl;
  cout << "  Another use of this tool is to enforce a minimum distance between two" << endl;
  cout << "  surface meshes such as in particular the inner (WM/cGM) and outer (cGM/CSF)" << endl;
  cout << "  cortical surfaces. Points of the <input> surface are displaced if they" << endl;
  cout << "  are closer than the allowed offset distance to the second surface." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input    Input surface mesh file." << endl;
  cout << "  output   Output surface mesh file." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -offset <distance>              Offset distance. (default: 0)" << endl;
  cout << "  -relative                       Multiply given distance value by length of bounding box diagonal." << endl;
  cout << "  -implicit                       Use offset surface reconstructed from implicit surface model." << endl;
  cout << "  -size <nx> [<ny> <nz>]          Number of voxels to use for implicit surface model." << endl;
  cout << "                                  (default: bounding box divided by :option:`-voxelsize`)" << endl;
  cout << "  -voxelsize <dx> [<dy> <dz>]     Size of voxels to use for implicit surface model." << endl;
  cout << "                                  This option is ignored when image :option:`-size` is specified." << endl;
  cout << "                                  (default: length of bounding box diagonal divided by 512)" << endl;
  cout << "  -surface <file>                 Reference surface to which a minimum offset distance" << endl;
  cout << "                                  to the input surface is required. (default: <input>)" << endl;
  cout << "  -scalars <name> <min> [<max>]   Only move points whose given scalars point data" << endl;
  cout << "                                  is within the given range." << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  // Parse arguments
  REQUIRES_POSARGS(2);

  const char *input_name     = POSARG(1);
  const char *output_name    = POSARG(2);
  const char *surface_name   = nullptr;
  const char *scalars_name   = nullptr;
  double      scalars_min    = -inf;
  double      scalars_max    = +inf;

  int    smooth_iter      = 0;
  double target_reduction = .0;
  bool   relative         = false;
  bool   implicit         = false;
  double offset           = 0.;
  bool   along_normal     = false;

  FileOption fopt = FO_Default;

  if (NUM_POSARGS == 3) {
    if (!FromString(POSARG(3), offset)) {
      FatalError("Invalid <offset> argument, must be floating point number!");
    }
  } else if (NUM_POSARGS > 3) {
    FatalError("Too many positional arguments!");
  }

  // resolution of implicit surface model
  int    nx =  0, ny =  0, nz =  0;
  double dx = .0, dy = .0, dz = .0;

  for (ALL_OPTIONS) {
    if (OPTION("-voxelsize") || OPTION("-voxel-size")) {
      PARSE_ARGUMENT(dx);
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(dy);
        PARSE_ARGUMENT(dz);
      } else {
        dy = dz = dx;
      }
    }
    else if (OPTION("-size")) {
      PARSE_ARGUMENT(nx);
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(ny);
        PARSE_ARGUMENT(nz);
      } else {
        ny = nz = nx;
      }
    }
    else if (OPTION("-offset")) PARSE_ARGUMENT(offset);
    else if (OPTION("-relative")) relative = true;
    else if (OPTION("-implicit")) implicit = true;
    else if (OPTION("-surface")) surface_name = ARGUMENT;
    else if (OPTION("-smooth")) PARSE_ARGUMENT(smooth_iter);
    else if (OPTION("-decimate")) PARSE_ARGUMENT(target_reduction);
    else if (OPTION("-scalars")) {
      scalars_name = ARGUMENT;
      PARSE_ARGUMENT(scalars_min);
      if (HAS_ARGUMENT) PARSE_ARGUMENT(scalars_max);
      else scalars_max = scalars_min;
    }
    else HANDLE_POINTSETIO_OPTION(fopt);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read input point set
  vtkSmartPointer<vtkPointSet> pointset = ReadPointSet(input_name, fopt);
  vtkSmartPointer<vtkPolyData> output   = DataSetSurface(pointset);
  vtkSmartPointer<vtkPolyData> surface  = output;
  vtkSmartPointer<vtkPolyData> offset_surface;

  if (surface_name) {
    surface = DataSetSurface(ReadPointSet(surface_name));
  }

  vtkSmartPointer<vtkDataArray> scalars;
  if (scalars_name) {
    scalars = output->GetPointData()->GetArray(scalars_name);
    if (!scalars) {
      vtkDataArray *cell_scalars = output->GetCellData()->GetArray(scalars_name);
      if (cell_scalars) {
        vtkSmartPointer<vtkPolyData> copy;
        copy.TakeReference(output->NewInstance());
        copy->ShallowCopy(output);
        copy->GetPointData()->Initialize();
        copy->GetCellData()->Initialize();
        copy->GetCellData()->AddArray(cell_scalars);
        vtkNew<vtkCellDataToPointData> converter;
        SetVTKInput(converter, copy);
        converter->PassCellDataOff();
        converter->Update();
        scalars = converter->GetOutput()->GetPointData()->GetArray(scalars_name);
        if (!scalars) {
          FatalError("Failed to convert cell -scalars to point data!");
        }
        output->GetPointData()->AddArray(scalars);
      }
    }
    if (!scalars) {
      FatalError("Input surface mesh has no data array named: " << scalars_name);
    }
  }

  double center[3], bounds[6];
  surface->GetCenter(center);
  surface->GetBounds(bounds);
  if (relative) offset *= surface->GetLength();

  if (implicit) {

    const bool inside = (offset < .0);
    offset = abs(offset);

    // Adjust implicit surface model bounds
    const double margin = 1.15 * offset;
    bounds[0] -= margin;
    bounds[1] += margin;
    bounds[2] -= margin;
    bounds[3] += margin;
    bounds[4] -= margin;
    bounds[5] += margin;

    // Sampling of implicit surface model
    if (nx <= 0 || ny <= 0 || nz <= 0) {
      if (dx <= .0 && nx > 0) dx = (bounds[1] - bounds[0]) / nx;
      if (dy <= .0 && ny > 0) dy = (bounds[3] - bounds[2]) / ny;
      if (dz <= .0 && nz > 0) dz = (bounds[5] - bounds[4]) / nz;
      double ds = min(min(dx, dy), dz);
      if (ds <= .0) ds = surface->GetLength() / 512;
      if (dx <= .0) dx = ds;
      if (dy <= .0) dy = ds;
      if (dz <= .0) dz = ds;
      if (nx <=  0) nx = iceil((bounds[1] - bounds[0]) / dx);
      if (ny <=  0) ny = iceil((bounds[3] - bounds[2]) / dy);
      if (nz <=  0) nz = iceil((bounds[5] - bounds[4]) / dz);
    }

    // Create implicit surface model from input surface
    vtkNew<vtkImplicitModeller> model;
    model->SetSampleDimensions(nx, ny, nz);
    model->SetMaximumDistance(1.1 * offset);
    model->SetModelBounds(bounds);
    model->AdjustBoundsOff();
    model->ScaleToMaximumDistanceOff();
    model->SetOutputScalarTypeToFloat();
    model->CappingOn();
    model->SetCapValue(margin);
    model->SetProcessModeToPerVoxel();
    SetVTKInput(model, surface);

    // Extract inside/outside offset surfaces
    // Note: Unfortunately, vtkImplicitModeller creates an unsigned distance field.
    vtkNew<vtkContourFilter> contours;
    contours->UseScalarTreeOn();
    contours->SetNumberOfContours(1);
    contours->SetValue(0, offset);
    SetVTKConnection(contours,  model);
    contours->Update();
    if (contours->GetOutput()->GetNumberOfPoints() == 0) {
      FatalError("Failed to contour offset surfaces");
    }

    // Only keep offset surface closest to bounding box corner (i.e., outside)
    vtkNew<vtkPolyDataConnectivityFilter> extract;
    if (inside) extract->SetClosestPoint(center[0], center[1], center[2]);
    else        extract->SetClosestPoint(bounds[0], bounds[2], bounds[4]);
    extract->SetExtractionModeToClosestPointRegion();
    SetVTKConnection(extract, contours);
    extract->Update();
    if (extract->GetOutput()->GetNumberOfPoints() == 0) {
      FatalError("Failed to extract offset surface");
    }

    extract->GetOutput()->GetPointData()->RemoveArray("ImageScalars");
    offset_surface = extract->GetOutput();

    // Smooth offset surface to reduce sampling artifacts
    if (smooth_iter > 0) {
      MeshSmoothing smoother;
      smoother.Input(offset_surface);
      smoother.NumberOfIterations(2 * smooth_iter);
      smoother.Weighting(MeshSmoothing::Combinatorial);
      smoother.SmoothPointsOn();
      smoother.AdjacentValuesOnlyOn();
      smoother.Lambda(.33);
      smoother.Mu(-.34);
      smoother.Run();
      offset_surface = smoother.Output();
    }

    // Decimate offset surface mesh
    if (target_reduction > .0) {
      vtkNew<vtkQuadricDecimation> decimate;
      SetVTKInput(decimate, offset_surface);
      decimate->SetTargetReduction(target_reduction);
      decimate->Update();
      offset_surface = decimate->GetOutput();
    }

  } else {

    // Make shallow copy to not modify/add normals of/to input surface
    vtkSmartPointer<vtkPolyData> copy;
    copy = vtkSmartPointer<vtkPolyData>::NewInstance(surface);
    copy->ShallowCopy(surface);

    // Compute outwards point normals
    vtkNew<vtkPolyDataNormals> calc_normals;
    SetVTKInput(calc_normals, copy);
    calc_normals->SplittingOff();
    calc_normals->ConsistencyOff();
    calc_normals->AutoOrientNormalsOff();
    calc_normals->ComputeCellNormalsOff();
    calc_normals->ComputePointNormalsOn();
    calc_normals->NonManifoldTraversalOff();
    calc_normals->Update();
    vtkSmartPointer<vtkDataArray> normals;
    normals = calc_normals->GetOutput()->GetPointData()->GetNormals();

    // Move points of input surface mesh
    double p[3], n[3], value;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(surface->GetNumberOfPoints());

    if (surface == output && scalars != nullptr) {
      for (vtkIdType ptId = 0; ptId < points->GetNumberOfPoints(); ++ptId) {
        surface->GetPoint(ptId, p);
        value = scalars->GetComponent(ptId, 0);
        if (scalars_min <= value && value <= scalars_max) {
          normals->GetTuple(ptId, n);
          p[0] += offset * n[0];
          p[1] += offset * n[1];
          p[2] += offset * n[2];
        }
        points->SetPoint(ptId, p);
      }
    } else {
      for (vtkIdType ptId = 0; ptId < points->GetNumberOfPoints(); ++ptId) {
        surface->GetPoint(ptId, p);
        normals->GetTuple(ptId, n);
        p[0] += offset * n[0];
        p[1] += offset * n[1];
        p[2] += offset * n[2];
        points->SetPoint(ptId, p);
      }
    }

    // Offset surface
    offset_surface = vtkSmartPointer<vtkPolyData>::NewInstance(surface);
    offset_surface->ShallowCopy(surface);
    offset_surface->SetPoints(points);

  }

  // Displace input surface points whose scalar values are within the specified
  // range and that are optionally too close to the reference surface.
  // If no reference surface is given, all points whose scalar values are
  // in the given range are displaced. When neither a reference surface nor
  // a scalar range is given, the output is the offset surface itself.
  if (surface == output && (!implicit || scalars == nullptr)) {

    output = offset_surface;

  } else {

    const double mindist2 = offset * offset;
    const double tol      = 1e-6;
    double       value, p[3], n[3], q[3], x[3], t, pcoords[3], dist2;
    vtkIdType    cellId;
    int          subId, nmoved = 0;

    // Make shallow copy to not modify/add normals of/to input surface
    vtkSmartPointer<vtkPolyData> copy;
    copy = vtkSmartPointer<vtkPolyData>::NewInstance(output);
    copy->ShallowCopy(output);

    // Compute outwards point normals
    vtkNew<vtkPolyDataNormals> calc_normals;
    SetVTKInput(calc_normals, copy);
    calc_normals->SplittingOff();
    calc_normals->ConsistencyOff();
    calc_normals->AutoOrientNormalsOff();
    calc_normals->ComputeCellNormalsOff();
    calc_normals->ComputePointNormalsOn();
    calc_normals->NonManifoldTraversalOff();
    calc_normals->Update();
    vtkSmartPointer<vtkDataArray> normals;
    normals = calc_normals->GetOutput()->GetPointData()->GetNormals();

    vtkSmartPointer<vtkCellLocator> locator;
    if (surface != output) {
      locator = vtkSmartPointer<vtkCellLocator>::New();
      locator->SetDataSet(surface);
      locator->BuildLocator();
    }

    vtkSmartPointer<vtkCellLocator> offset_locator;
    offset_locator = vtkSmartPointer<vtkCellLocator>::New();
    offset_locator->SetDataSet(offset_surface);
    offset_locator->BuildLocator();

    vtkPoints * const points = output->GetPoints();
    for (vtkIdType ptId = 0; ptId < output->GetNumberOfPoints(); ++ptId) {
      if (scalars) {
        value = scalars->GetComponent(ptId, 0);
        if (value < scalars_min || value > scalars_max) continue;
      }
      points->GetPoint(ptId, p);
      if (locator) {
        locator->FindClosestPoint(p, x, cellId, subId, dist2);
        if (dist2 >= mindist2) continue;
      }
      if (along_normal) {
        normals->GetTuple(ptId, n);
        q[0] = p[0] + 1.1 * offset * n[0];
        q[1] = p[1] + 1.1 * offset * n[1];
        q[2] = p[2] + 1.1 * offset * n[2];
        if (offset_locator->IntersectWithLine(p, q, tol, t, x, pcoords, subId)) {
          points->SetPoint(ptId, x);
          ++nmoved;
        }
      } else {
        offset_locator->FindClosestPoint(p, x, cellId, subId, dist2);
        points->SetPoint(ptId, x);
        ++nmoved;
      }
    }

    if (verbose) {
      cout << "Modified " << nmoved << " point" << (nmoved == 1 ? "" : "s") << endl;
    }
  }

  // Write output surface mesh
  if (!WritePolyData(output_name, output, fopt)) {
    FatalError("Failed to write offset surface to " << output_name);
  }
  return 0;
}

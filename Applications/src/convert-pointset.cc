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

#include "mirtk/Point.h"
#include "mirtk/PointSet.h"
#include "mirtk/GenericImage.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/EuclideanDistanceTransform.h"

#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkAppendFilter.h"
#include "vtkAppendPolyData.h"
#include "vtkCleanPolyData.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <input>... <output> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Convert point set from one (file) format to another.\n";
  cout << "\n";
  cout << "  If more than one input point sets is given, all input point sets\n";
  cout << "  are appended into a single point set before conversion. When all\n";
  cout << "  input point sets are of type vtkPolyData, the output point set\n";
  cout << "  is of type vtkPolyData. If more than one input point set is given\n";
  cout << "  and not all are of type vtkPolyData, the output point set is of type\n";
  cout << "  vtkUnstructuredGrid.\n";
  cout << "\n";
  cout << "  The current implementation can only convert between different\n";
  cout << "  point set file formats based on the file name extension.\n";
  cout << "  Besides the common formats supported by VTK, it can also read/write\n";
  cout << "  BrainSuite .dfs files and write a Piecewise Linear Complex (PLC)\n";
  cout << "  B-Rep description in the TetGen formats .poly and .smesh. It also\n";
  cout << "  supports the Object File Format (.off) used by the CGAL library.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Input  point set file (.vtk, .vtp, .vtu, .stl, .ply, .off, .dfs, .obj).\n";
  cout << "  output   Output point set file (.vtk, .vtp, .vtu, .stl, .ply, .off, .dfs,\n";
  cout << "                                  .node, .poly, .smesh, .gii, .csv, .tsv, .txt).\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -merge [<tol>]   Merge points of non-polygonal input point sets.\n";
  cout << "                   When only vtkPolyData are merged, the maximum distance between\n";
  cout << "                   points to be merged can be optionally given as argument. (default: off)\n";
  cout << "  -holes           Add a point inside the input surface meshes, except off the first input mesh,\n";
  cout << "                   to the hole list of the output .poly or .smesh file. (default: no holes)\n";
  cout << "  -nocelldata      Do not write cell  data to output file. (default: off)\n";
  cout << "  -nopointdata     Do not write point data to output file. (default: off)\n";
  cout << "  -nopoints        Do not write point coordinates to output .csv, .tsv, or .txt file.\n";
  cout << "                   This option is implicit when the output file name ends with\n";
  cout << "                   .attr.csv, .attr.tsv, or .attr.txt (default: off)\n";
  cout << "  -ascii           Write legacy VTK files encoded in ASCII. (default: off)\n";
  cout << "  -binary          Write legacy VTK files in binary form. (default: on)\n";
  cout << "  -compress        Compress XML VTK files. (default: on)\n";
  cout << "  -nocompress      Do not compress XML VTK files. (default: off)\n";
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Get image attributes for implicit surface representation
ImageAttributes ImplicitSurfaceAttributes(vtkSmartPointer<vtkPolyData> surface)
{
  ImageAttributes attr;
  double bounds[6];
  surface->GetBounds(bounds);
  double ds = sqrt(pow(bounds[1] - bounds[0], 2) +
                   pow(bounds[3] - bounds[2], 2) +
                   pow(bounds[5] - bounds[4], 2)) / 128;
  attr._xorigin = bounds[0] + .5 * (bounds[1] - bounds[0]);
  attr._yorigin = bounds[2] + .5 * (bounds[3] - bounds[2]);
  attr._zorigin = bounds[4] + .5 * (bounds[5] - bounds[4]);
  attr._x       = int(ceil((bounds[1] - bounds[0]) / ds));
  attr._y       = int(ceil((bounds[3] - bounds[2]) / ds));
  attr._z       = int(ceil((bounds[5] - bounds[4]) / ds));
  attr._dx      = ds;
  attr._dy      = ds;
  attr._dz      = ds;
  return attr;
}

// -----------------------------------------------------------------------------
/// Get binary mask of region enclosed by surface
void GetSurfaceMask(vtkSmartPointer<vtkPolyData> surface, BinaryImage &mask)
{
  surface = vtkPolyData::SafeDownCast(WorldToImage(surface, &mask));
  vtkSmartPointer<vtkImageData>        vtkmask = NewVtkMask(mask.X(), mask.Y(), mask.Z());
  vtkSmartPointer<vtkImageStencilData> stencil = ImageStencil(vtkmask, surface);
  ImageStencilToMask(stencil, vtkmask);
  mask.CopyFrom(reinterpret_cast<BinaryPixel *>(vtkmask->GetScalarPointer()));
}

// -----------------------------------------------------------------------------
/// Compute signed distance map
void ComputeDistanceMap(vtkSmartPointer<vtkPolyData> surface, const BinaryImage &mask, bool invert, RealImage &dmap)
{
  const int nvox = mask.NumberOfSpatialVoxels();

  // Initialize distance map
  dmap.Initialize(mask.Attributes(), 1);
  if (invert) {
    for (int vox = 0; vox < nvox; ++vox) {
      if (mask(vox)) dmap(vox) = 0.0;
      else           dmap(vox) = 1.0;
    }
  } else {
    for (int vox = 0; vox < nvox; ++vox) {
      if (mask(vox)) dmap(vox) = 1.0;
      else           dmap(vox) = 0.0;
    }
  }

  // Calculate Euclidean distance transforms
  typedef EuclideanDistanceTransform<RealPixel> DistanceTransform;
  DistanceTransform edt(DistanceTransform::DT_3D);
  edt.Input (&dmap);
  edt.Output(&dmap);
  edt.Run();

  for (int vox = 0; vox < nvox; ++vox) {
    dmap(vox) = sqrt(dmap(vox));
  }
  if (invert) dmap *= -1.0;
}

// -----------------------------------------------------------------------------
/// Get a point inside with maximum distance from the implicit surface
bool GetPointInside(const RealImage &dmap, Point &p)
{
  const int nvox = dmap.NumberOfSpatialVoxels();
  int    idx  = -1, i, j, k;
  double dist = .0;
  for (int vox = 0; vox < nvox; ++vox) {
    if (dmap(vox) < dist) {
      idx  = vox;
      dist = dmap(vox);
    }
  }
  if (idx == -1) return false;
  dmap.IndexToVoxel(idx, i, j, k);
  p._x = i, p._y = j, p._z = k;
  dmap.ImageToWorld(p);
  return true;
}

// -----------------------------------------------------------------------------
/// Get a point inside with maximum distance from the surface mesh
bool GetPointInside(vtkSmartPointer<vtkPointSet> &pointset, Point &p)
{
  vtkSmartPointer<vtkPolyData> surface = DataSetSurface(pointset);
  BinaryImage mask(ImplicitSurfaceAttributes(surface));
  RealImage   dmap;
  GetSurfaceMask(surface, mask);
  ComputeDistanceMap(surface, mask, true, dmap);
  return GetPointInside(dmap, p);
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  REQUIRES_POSARGS(2);

  // Positional arguments
  Array<const char *> input_names(NUM_POSARGS - 1, NULL);
  for (int i = 1; i < NUM_POSARGS; ++i) input_names[i-1] = POSARG(i);
  const char *output_name = POSARG(NUM_POSARGS);

  // Optional arguments
  bool pointdata   = true;
  bool celldata    = true;
  bool pointcoords = true;
  bool pointids    = true;
  bool merge       = false;
  bool with_holes  = false;
  double merge_tol = 1e-6;

  FileOption fopt = FO_Default;

  for (ALL_OPTIONS) {
    if      (OPTION("-nopointdata")) pointdata   = false;
    else if (OPTION("-nocelldata"))  celldata    = false;
    else if (OPTION("-nopoints"))    pointcoords = false;
    else if (OPTION("-nopointids"))  pointids    = false;
    else if (OPTION("-merge")) {
      merge = true;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(merge_tol);
    }
    else if (OPTION("-holes")) with_holes = true;
    else HANDLE_POINTSETIO_OPTION(fopt);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Output file name extension
  const string ext = Extension(output_name);

  // Read input point sets
  Array<vtkSmartPointer<vtkPointSet> > pointsets(input_names.size());
  for (size_t i = 0; i < input_names.size(); ++i) {
    FileOption opt;
    pointsets[i] = ReadPointSet(input_names[i], opt);
    if (fopt == FO_Default) fopt = opt;
  }

  // Combine point sets into single point set
  vtkSmartPointer<vtkPointSet> output;

  if (pointsets.size() == 1) {
    output = pointsets[0];
  } else {
    bool output_polydata = true;
    for (size_t i = 0; i < pointsets.size(); ++i) {
      if (vtkPolyData::SafeDownCast(pointsets[i]) == NULL) {
        output_polydata = false;
        break;
      }
    }
    if (output_polydata) {
      vtkNew<vtkAppendPolyData> append;
      for (size_t i = 0; i < pointsets.size(); ++i) {
        AddVTKInput(append, vtkPolyData::SafeDownCast(pointsets[i]));
      }
      append->Update();
      output = append->GetOutput();
    } else {
      vtkNew<vtkAppendFilter> append;
      for (size_t i = 0; i < pointsets.size(); ++i) {
        AddVTKInput(append, pointsets[i]);
      }
      append->SetMergePoints(merge);
      append->Update();
      output = append->GetOutput();
    }
  }

  // Merge points/cells
  vtkSmartPointer<vtkPolyData> surface = vtkPolyData::SafeDownCast(output);
  if (merge && surface) {
    // Merge points
    vtkNew<vtkCleanPolyData> merger;
    merger->ConvertLinesToPointsOff();
    merger->ConvertPolysToLinesOff();
    merger->ConvertStripsToPolysOn();
    merger->PointMergingOn();
    merger->ToleranceIsAbsoluteOn();
    merger->SetAbsoluteTolerance(merge_tol);
    SetVTKInput(merger, output);
    merger->Update();
    surface = merger->GetOutput();
    // Remove duplicate cells
    surface->BuildLinks();
    vtkNew<vtkIdList> ptIds1, ptIds2, cellIds;
    for (vtkIdType cellId = 0; cellId < output->GetNumberOfCells(); ++cellId) {
      if (surface->GetCellType(cellId) != VTK_EMPTY_CELL) {
        surface->GetCellPoints(cellId, ptIds1.GetPointer());
        for (vtkIdType i = 0; i < ptIds1->GetNumberOfIds(); ++i) {
          surface->GetPointCells(ptIds1->GetId(i), cellIds.GetPointer());
          for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
            if (cellIds->GetId(j) > cellId && surface->GetCellType(cellIds->GetId(j)) != VTK_EMPTY_CELL) {
              surface->GetCellPoints(cellIds->GetId(j), ptIds2.GetPointer());
              if (ptIds1->GetNumberOfIds() == ptIds2->GetNumberOfIds()) {
                ptIds2->IntersectWith(ptIds1.GetPointer());
                if (ptIds1->GetNumberOfIds() == ptIds2->GetNumberOfIds()) {
                  surface->DeleteCell(cellIds->GetId(j));
                }
              }
            }
          }
        }
      }
    }
    surface->RemoveDeletedCells();
    output = surface;
  }

  // Reset point/cell data
  if (!pointdata) output->GetPointData()->Initialize();
  if (!celldata ) output->GetCellData ()->Initialize();

  // Write .csv, .tsv, and .txt with or without point coordinates
  if (ext == ".csv" || ext == ".tsv" || ext == ".txt") {
    char sep = ',';
    if (ext == ".tsv") sep = '\t';
    string type = output_name;
    type = type.substr(0, type.length() - ext.length());
    type = Extension(type, EXT_Last);
    if (type == ".attr") pointcoords = false;
    if (!WritePointSetTable(output_name, output, sep, pointids, pointcoords)) {
      FatalError("Failed to write output point set to " << output_name);
    }
    return 0;
  }

  // Write TetGen .poly/.smesh with holes list
  if (with_holes && input_names.size() > 1) {
    if (ext == ".poly" || ext == ".smesh") {
      vtkSmartPointer<vtkPolyData> surface = vtkPolyData::SafeDownCast(output);
      if (surface == NULL) {
        FatalError("Can only save surface meshes (vtkPolyData) to " << ext << " file");
      }
      int ok = 0;
      Point    hole;
      PointSet holes;
      holes.Reserve(static_cast<int>(pointsets.size()) - 1);
      for (size_t i = 1; i < pointsets.size(); ++i) {
        if (GetPointInside(pointsets[i], hole)) {
          if (verbose) cout << "Add hole at [" << hole << "]" << endl;
          holes.Add(hole);
        }
      }
      if (ext == ".poly") ok = WriteTetGenPoly (output_name, surface, &holes);
      else                ok = WriteTetGenSMesh(output_name, surface, &holes);
      if (!ok) {
        FatalError("Failed to write output point set to " << output_name);
      }
      return 0;
    }
  }

  // Write output point set
  if (!WritePointSet(output_name, output, fopt)) {
    FatalError("Failed to write output point set to " << output_name);
  }
  return 0;
}

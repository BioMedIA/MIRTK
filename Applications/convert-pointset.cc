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

#include <mirtkPoint.h>
#include <mirtkPointSet.h>
#include <mirtkGenericImage.h>
#include <mirtkPointSetUtils.h>
#include <mirtkEuclideanDistanceTransform.h>

#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkAppendFilter.h>
#include <vtkAppendPolyData.h>

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <input>... <output> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Convert point set from one (file) format to another." << endl;
  cout << endl;
  cout << "  If more than one input point sets is given, all input point sets" << endl;
  cout << "  are appended into a single point set before conversion. When all" << endl;
  cout << "  input point sets are of type vtkPolyData, the output point set" << endl;
  cout << "  is of type vtkPolyData. If more than one input point set is given" << endl;
  cout << "  and not all are of type vtkPolyData, the output point set is of type" << endl;
  cout << "  vtkUnstructuredGrid." << endl;
  cout << endl;
  cout << "  The current implementation can only convert between different" << endl;
  cout << "  point set file formats based on the file name extension." << endl;
  cout << "  Besides the common formats supported by VTK, it can also read/write" << endl;
  cout << "  BrainSuite .dfs files and write a Piecewise Linear Complex (PLC)" << endl;
  cout << "  B-Rep description in the TetGen formats .poly and .smesh." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input    Input  point set file (.vtk, .vtp, .vtu, .stl, .ply, .dfs, .obj)." << endl;
  cout << "  output   Output point set file (.vtk, .vtp, .vtu, .stl, .ply, .dfs, .node, .poly, .smesh)." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -merge         Merge points of non-polygonal input point sets. (default: off)" << endl;
  cout << "  -holes         Add a point inside the input surface meshes, except off the first input mesh," << endl;
  cout << "                 to the hole list of the output .poly or .smesh file. (default: no holes)" << endl;
  cout << "  -nocelldata    Do not write cell  data to output file. (default: off)" << endl;
  cout << "  -nopointdata   Do not write point data to output file. (default: off)" << endl;
  cout << "  -ascii         Write legacy VTK files encoded in ASCII. (default: off)" << endl;
  cout << "  -binary        Write legacy VTK files in binary form. (default: on)" << endl;
  cout << "  -compress      Compress XML VTK files. (default: on)" << endl;
  cout << "  -nocompress    Do not compress XML VTK files. (default: off)" << endl;
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
  bool pointdata  = true;
  bool celldata   = true;
  bool ascii      = false;
  bool compress   = true;
  bool merge      = false;
  bool with_holes = false;

  for (ALL_OPTIONS) {
    if      (OPTION("-nopointdata")) pointdata = false;
    else if (OPTION("-nocelldata"))  celldata  = false;
    else if (OPTION("-merge"))       merge = true;
    else if (OPTION("-holes"))       with_holes = true;
    else if (OPTION("-ascii"))       ascii = true;
    else if (OPTION("-binary"))      ascii = false;
    else if (OPTION("-compress"))    compress = true;
    else if (OPTION("-nocompress"))  compress = false;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read input point sets
  Array<vtkSmartPointer<vtkPointSet> > pointsets(input_names.size());
  for (size_t i = 0; i < input_names.size(); ++i) {
    pointsets[i] = ReadPointSet(input_names[i], NULL, true);
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

  // Reset point/cell data
  if (!pointdata) output->GetPointData()->Initialize();
  if (!celldata ) output->GetCellData ()->Initialize();

  // Write TetGen .poly/.smesh with holes list
  if (with_holes && input_names.size() > 1) {
    const std::string ext = Extension(output_name);
    if (ext == ".poly" || ext == ".smesh") {
      vtkSmartPointer<vtkPolyData> surface = vtkPolyData::SafeDownCast(output);
      if (surface == NULL) {
        FatalError("Can only save surface meshes (vtkPolyData) to " << ext << " file");
      }
      int ok = 0;
      Point    hole;
      PointSet holes;
      holes.Reserve(pointsets.size() - 1);
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
  if (!WritePointSet(output_name, output, compress, ascii)) {
    FatalError("Failed to write output point set to " << output_name);
    exit(1);
  }
  return 0;
}

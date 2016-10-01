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
#include "mirtk/GenericImage.h"
#include "mirtk/EdgeConnectivity.h"
#include "mirtk/MeshSmoothing.h"
#include "mirtk/SurfaceBoundary.h"
#include "mirtk/SurfaceCollisions.h"
#include "mirtk/EuclideanDistanceTransform.h"

#include "mirtk/InterpolateImageFunction.h"
#include "mirtk/LinearInterpolateImageFunction.hxx"

#include "mirtk/Vtk.h"

#include "vtkSmartPointer.h"
#include "vtkPolygon.h"
#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkPolyData.h"
#include "vtkImageData.h"
#include "vtkImageStencilData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSphereSource.h"
#include "vtkDelaunay2D.h"
#include "vtkDelaunay3D.h"
#include "vtkBooleanOperationPolyDataFilter.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkFillHolesFilter.h"
#include "vtkPointLocator.h"
#include "vtkPolyDataNormals.h"
#include "vtkCleanPolyData.h"
#include "vtkGenericCell.h"


using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " [options]" << endl;
  cout << "       " << name << " [<input>...] <output> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Extract surface of point set. If more than one input point set is given," << endl;
  cout << "  it computes the boundary of the union, intersection, or difference volume" << endl;
  cout << "  computed from the volumes defined by the individual input surfaces." << endl;
  cout << endl;
  cout << "Input options:" << endl;
  cout << "  -i -input <file>       Input point set." << endl;
  cout << "  -triangulate           Triangulate input point set surface." << endl;
  cout << "  -delaunay              Use surface of Delaunay tesselation of input point set." << endl;
  cout << "  -hull                  Use convex hull of input point set (implies :option:`-delaunay`)." << endl;
  cout << endl;
  cout << "Boolean operations:" << endl;
  cout << "  -union                 Compute boundary of union of input volumes. (default)" << endl;
  cout << "  -difference            Compute boundary of difference of input volumes." << endl;
  cout << "  -intersection          Compute boundary of intersection of input volumes." << endl;
  cout << endl;
  cout << "Output options:" << endl;
  cout << "  -source-array <name>    Add point/cell data array with the specified name with one-based" << endl;
  cout << "                          labels corresponding to the input point set from which an output" << endl;
  cout << "                          point/cell originates from. When the first input point set has" << endl;
  cout << "                          a scalar array with the specified name, the labels of this first" << endl;
  cout << "                          surface are preserved, while successive labels are offset by the" << endl;
  cout << "                          maximum integer value of the input data array. This is useful when" << endl;
  cout << "                          successively merging surface meshes instead of with a single execution" << endl;
  cout << "                          of this command. (default: none)" << endl;
  cout << "  -merge [<float>]        Merge points closer than the specified distance, default is 1e-6. (default: off)" << endl;
  cout << "  -tolerance <float>      Distance tolerance value to use for boolean operations. (default: 1e-6)" << endl;
  cout << "  -fill-holes [<size>]    Fill holes of given maximum size (i.e., circumsphere radius)." << endl;
  cout << "  -largest [<n>]          Output largest n (default 1 if not specified) components." << endl;
  cout << "  -insphere               Output sphere which is inscribed the surface mesh." << endl;
  cout << "  -bounding-sphere        Output minimum sphere which fully contains the surface mesh." << endl;
  cout << "  -o -surface <file>      Write output surface mesh to named file." << endl;
  cout << "  -mask <file>            Write binary inside/outside mask to named image file." << endl;
  cout << "  -implicit <file>        Write signed implicit surface distance to named image file." << endl;
  cout << "  -[no]outside [on|off]   Swap inside/outside of output :option:`-mask` and :option:`-implicit` surface. (default: off)" << endl;
  cout << "  -reference <image>      Reference image for :option:`-mask` or :option:`-implicit` output. (default: bounding box)" << endl;
  cout << "  -resolution <float>     Resolution of :option:`-insphere`, :option:`-mask`, or :option:`-implicit` output." << endl;
  cout << "                          (default: fraction of length of bounding box diagonal)" << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Shock removal (experimental)
// =============================================================================

// -----------------------------------------------------------------------------
// Typedefs
typedef RealImage                                            DistanceImage;
typedef GenericLinearInterpolateImageFunction<DistanceImage> DistanceFunction;

// -----------------------------------------------------------------------------
/// Get distance value at given world position
inline double Evaluate(const DistanceFunction *dist, double p[3])
{
  double d, x = p[0], y = p[1], z = p[2];
  dist->WorldToImage(x, y, z);
  if (dist->IsInside(x, y, z)) {
    d = dist->GetInside(x, y, z);
  } else {
    double xmin, ymin, zmin, xmax, ymax, zmax;
    dist->Inside(xmin, ymin, zmin, xmax, ymax, zmax);
    x = clamp(x, xmin, xmax);
    y = clamp(y, ymin, ymax);
    z = clamp(z, zmin, zmax);
    d = dist->GetInside(x, y, z);
    dist->ImageToWorld(x, y, z);
    x -= p[0], y -= p[1], z -= p[2];
    d += sqrt(x*x + y*y + z*z);
  }
  return d;
}

// -----------------------------------------------------------------------------
/// Remove "shocks" from cortical surface mesh and fix topology afterwards
///
/// These shocks may occur when an initial surface mesh such as the convex hull
/// or a bounding sphere is deformed towards an implicit surface given by the
/// zero level isosurface of a discrete distance function. A hard
/// non-self-intersection constraint prevents the surface from intersecting
/// itself. Additionally, a repulsive node force may favor a certain minimum
/// distance between non-neighboring nodes. The surface deformation may depending
/// on the topology of the input isosurface/segmentation, or simply due to a
/// non-smooth propagation of the surface mesh and the aforementioned constraints,
/// end up with a number of "shocks", i.e., moving fronts which nearly touch
/// each other and thus their further deformation is prevented. This function
/// aims to fixup such output surface (cf. deformmesh tool).
void RemoveShocks(vtkSmartPointer<vtkPolyData> &surface, double radius,
                  const DistanceImage *dmap = NULL)
{

  const int ncells = surface->GetNumberOfCells();
  surface->BuildLinks();

  SurfaceCollisions collisions;
  collisions.Input(surface);
  collisions.MinBackfaceDistance(radius);
  collisions.MinFrontfaceDistance(radius);
  collisions.AdjacentIntersectionTestOff();
  collisions.NonAdjacentIntersectionTestOff();
  collisions.FrontfaceCollisionTestOn();
  collisions.BackfaceCollisionTestOn();
  collisions.Run();

  int    subId;
  double pcoords[3], p[3], *weights = new double[surface->GetMaxCellSize()];
  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();

  double mind = .0;
  DistanceFunction distance;
  if (dmap) {
    distance.DefaultValue(numeric_limits<double>::infinity());
    distance.Input(dmap);
    distance.Initialize();
    mind = .1 * min(dmap->GetXSize(), min(dmap->GetYSize(), dmap->GetZSize()));
  }

  int ncollisions = 0;
  for (int cellId = 0; cellId < ncells; ++cellId) {
    if (collisions.GetCollisionType(cellId) != SurfaceCollisions::NoCollision) {
      if (dmap) {
        surface->GetCell(cellId, cell);
        subId = cell->GetParametricCenter(pcoords);
        cell->EvaluateLocation(subId, pcoords, p, weights);
        if (Evaluate(&distance, p) <= mind) continue;
      }
      surface->DeleteCell(cellId);
      ++ncollisions;
    }
  }
  delete[] weights;
  surface->RemoveDeletedCells();

  if (verbose > 1) cout << "  No. of deleted cells  = " << ncollisions << endl;

  vtkNew<vtkCleanPolyData> cleaner;
  SetVTKInput(cleaner, surface);
  cleaner->ToleranceIsAbsoluteOn();
  cleaner->SetAbsoluteTolerance(.0);
  cleaner->Update();
  surface = cleaner->GetOutput();
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
enum OutputType
{
  OutputSurface,
  OutputInSphere,
  OutputBoundingSphere
};

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // Print help if called with no arguments
  if (argc == 1) {
    PrintHelp(EXECNAME);
    exit(1);
  }

  // ---------------------------------------------------------------------------
  // Parse arguments
  REQUIRES_POSARGS(0);

  Array<const char *> input_names;
  const char *output_name       = nullptr;
  const char *mask_name         = nullptr;
  const char *dmap_name         = nullptr;
  int         operation         = vtkBooleanOperationPolyDataFilter::VTK_UNION;
  int         hull_levels       = -1;
  bool        delaunay          = false;
  bool        triangulate       = (input_names.size() > 1);
  double      merge_points_tol  = -1.;
  double      max_hole_size     = 0.;
  double      max_shock_dist    = NaN;
  const char *shock_dmap_name   = NULL;
  int         nlargest          = 0;
  OutputType  type              = OutputSurface;
  const char *ref_name          = NULL;
  double      resolution        = NaN;
  double      tolerance         = 1e-6;
  const char *source_array_name = nullptr;
  bool        inside_out        = false;

  if (NUM_POSARGS > 0) {
    if (NUM_POSARGS > 1) {
      input_names.resize(NUM_POSARGS - 1, nullptr);
      for (int i = 1; i < NUM_POSARGS; ++i) input_names[i-1] = POSARG(i);
    }
    output_name = POSARG(NUM_POSARGS);
  }

  for (ALL_OPTIONS) {
    if (OPTION("-i") || OPTION("-input")) {
      do {
        input_names.push_back(ARGUMENT);
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-o") || OPTION("-output") || OPTION("-surface")) {
      if (output_name) input_names.push_back(output_name);
      output_name = ARGUMENT;
    }
    else if (OPTION("-mask")) {
      mask_name = ARGUMENT;
    }
    else if (OPTION("-implicit")) {
      dmap_name = ARGUMENT;
    }
    else if (OPTION("-hull")){
      if (HAS_ARGUMENT) PARSE_ARGUMENT(hull_levels);
      else hull_levels = 3;
    }
    else if (OPTION("-merge")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(merge_points_tol);
      else merge_points_tol = 1e-6;
    }
    else if (OPTION("-nomerge")) {
      merge_points_tol = -1.;
    }
    else if (OPTION("-remove-shocks") || OPTION("-removeshocks")) {
      shock_dmap_name = nullptr;
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(max_shock_dist);
        if (HAS_ARGUMENT) shock_dmap_name = ARGUMENT;
      }
      else max_shock_dist = -1.0;
    }
    else if (OPTION("-fill-holes") || OPTION("-fillholes")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(max_hole_size);
      else max_hole_size = inf;
    }
    else if (OPTION("-largest")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(nlargest);
      else nlargest = 1;
    }
    else if (OPTION("-delaunay"))    delaunay    = true;
    else if (OPTION("-triangulate")) triangulate = true;
    else if (OPTION("-insphere")) type = OutputInSphere;
    else if (OPTION("-bounding-sphere")) type = OutputBoundingSphere;
    else if (OPTION("-resolution") || OPTION("-res")) PARSE_ARGUMENT(resolution);
    else if (OPTION("-reference")  || OPTION("-ref")) ref_name = ARGUMENT;
    else if (OPTION("-union"))        operation = vtkBooleanOperationPolyDataFilter::VTK_UNION;
    else if (OPTION("-difference"))   operation = vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE;
    else if (OPTION("-intersection")) operation = vtkBooleanOperationPolyDataFilter::VTK_INTERSECTION;
    else if (OPTION("-tolerance") || OPTION("-tol")) PARSE_ARGUMENT(tolerance);
    else if (OPTION("-source-array")) source_array_name = ARGUMENT;
    else HANDLE_BOOLEAN_OPTION("outside", inside_out);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (input_names.empty()) {
    FatalError("No input point sets given!");
  }
  if (!output_name && !mask_name && !dmap_name) {
    FatalError("Missing output option (-o/-surface, -mask, and/or -implicit)!");
  }

  // ---------------------------------------------------------------------------
  // Read input point sets and extract their surfaces
  Array<vtkSmartPointer<vtkPolyData> > surfaces(input_names.size());
  for (size_t i = 0; i < input_names.size(); ++i) {
    if (verbose) cout << "Reading surface from " << input_names[i] << endl;
    vtkSmartPointer<vtkPointSet> pointset = ReadPointSet(input_names[i]);
    if (hull_levels >= 0) {
      surfaces[i] = ConvexHull(pointset, hull_levels);
    } else {
      if (delaunay) {
        vtkNew<vtkDelaunay3D> tesselation;
        SetVTKInput(tesselation, pointset);
        tesselation->Update();
        pointset = tesselation->GetOutput();
      }
      surfaces[i] = DataSetSurface(pointset);
    }
  }

  // ---------------------------------------------------------------------------
  // Compute union, difference, or intersection of input surfaces
  if (verbose) cout << "Combining surface meshes...", cout.flush();

  vtkSmartPointer<vtkDataArray> source_array;
  const char *temp_source_array_name = source_array_name;
  if (source_array_name != nullptr) {
    if (strcmp(source_array_name, "PointSource") == 0) {
      temp_source_array_name = "_PointSource";
    }
    size_t offset = 1;
    source_array = surfaces[0]->GetPointData()->GetArray(source_array_name);
    if (source_array != nullptr && source_array->GetNumberOfComponents() == 1) {
      offset = static_cast<size_t>(ceil(source_array->GetRange(0)[1]));
      if (offset > 1) source_array->SetName(temp_source_array_name);
    }
    for (size_t i = 0; i < surfaces.size(); ++i) {
      const auto &surface = surfaces[i];
      surface->GetCellData()->RemoveArray(source_array_name);
      if (i > 0 || offset <= 1) {
        surface->GetPointData()->RemoveArray(source_array_name);
        source_array = NewVTKDataArray(VTK_UNSIGNED_SHORT);
        source_array->SetName(temp_source_array_name);
        source_array->SetNumberOfComponents(1);
        source_array->SetNumberOfTuples(surface->GetNumberOfPoints());
        source_array->FillComponent(0, static_cast<double>(offset + i));
        surface->GetPointData()->AddArray(source_array);
      }
    }
    source_array = nullptr;
  }

  vtkSmartPointer<vtkPolyData> surface = surfaces[0];
  for (size_t i = 1; i < surfaces.size(); ++i) {
    vtkNew<vtkBooleanOperationPolyDataFilter> joiner;
    SetNthVTKInput(joiner, 0, surface);
    SetNthVTKInput(joiner, 1, surfaces[i]);
    joiner->SetOperation(operation);
    joiner->SetTolerance(tolerance);
    joiner->ReorientDifferenceCellsOn();
    joiner->Update();

    surface = joiner->GetOutput();
    surface->GetPointData()->RemoveArray("PointSource");
    surface->GetCellData ()->RemoveArray("CellSource");
    surface->GetPointData()->RemoveArray("Distance");
    surface->GetCellData ()->RemoveArray("Distance");

    vtkNew<vtkCleanPolyData> cleaner;
    cleaner->ConvertPolysToLinesOff();
    cleaner->ConvertLinesToPointsOff();
    cleaner->ConvertStripsToPolysOff();
    if (merge_points_tol >= 0.) {
      cleaner->PointMergingOn();
      cleaner->ToleranceIsAbsoluteOn();
      cleaner->SetAbsoluteTolerance(merge_points_tol);
    } else {
      cleaner->PointMergingOff();
    }
    SetVTKInput(cleaner, surface);
    cleaner->Update();
    surface = cleaner->GetOutput();
  }

  if (temp_source_array_name != source_array_name) {
    surface->GetPointData()->RemoveArray(source_array_name);
    source_array = surface->GetPointData()->GetArray(temp_source_array_name);
    source_array->SetName(source_array_name);
    temp_source_array_name = nullptr;
  }

  surfaces.clear();
  if (verbose) cout << endl;

  // ---------------------------------------------------------------------------
  // Remove "shocks" from surface mesh
  if (!IsNaN(max_shock_dist)) {
    if (verbose) {
      cout << "Removing shocks...";
      if (verbose > 1) cout << "\n";
      cout.flush();
    }
    if (max_shock_dist < .0) {
      max_shock_dist = abs(max_shock_dist) * AverageEdgeLength(surface);
    }
    if (verbose > 1) cout << "  Maximum shock distance = " << max_shock_dist << endl;
    UniquePtr<DistanceImage> dmap;
    if (shock_dmap_name) dmap.reset(new DistanceImage(shock_dmap_name));
    RemoveShocks(surface, max_shock_dist, dmap.get());
    if (verbose) {
      if (verbose > 1) cout << "Removing shocks...";
      cout << " done" << endl;
    }
  }

  // ---------------------------------------------------------------------------
  // Fill holes
  if (max_hole_size > .0) {
    if (verbose) cout << "Filling holes...", cout.flush();
    vtkNew<vtkFillHolesFilter> fill;
    SetVTKInput(fill, surface);
    fill->SetHoleSize(max_hole_size);
    fill->Update();
    surface = fill->GetOutput();
    if (verbose) cout << " done" << endl;
  }

  // ---------------------------------------------------------------------------
  // Extract largest n components
  if (nlargest > 0) {
    if (verbose) {
      cout << "Extracting largest ";
      if (nlargest > 1) cout << nlargest << " ";
      cout << "component";
      if (nlargest > 1) cout << "s";
      cout << "...";
      cout.flush();
    }
    vtkNew<vtkPolyDataConnectivityFilter> conn;
    SetVTKInput(conn, surface);
    conn->SetExtractionModeToAllRegions();
    conn->Update();
    conn->SetExtractionModeToSpecifiedRegions();
    for (int i = 0; i < conn->GetNumberOfExtractedRegions(); ++i) {
      conn->DeleteSpecifiedRegion(i);
    }
    nlargest = min(nlargest, conn->GetNumberOfExtractedRegions());
    for (int i = 0; i < nlargest; ++i) {
      conn->AddSpecifiedRegion(i);
    }
    conn->Update();
    surface = conn->GetOutput();
    vtkNew<vtkCleanPolyData> cleaner;
    cleaner->ConvertPolysToLinesOff();
    cleaner->ConvertLinesToPointsOff();
    cleaner->ConvertStripsToPolysOff();
    cleaner->PointMergingOff();
    SetVTKInput(cleaner, surface);
    cleaner->Update();
    surface = cleaner->GetOutput();
    if (verbose) cout << " done" << endl;
  }

  // ---------------------------------------------------------------------------
  // Generate desired output from given surface mesh
  switch (type) {
    case OutputInSphere:
    case OutputBoundingSphere: {
      double c[3], p[3], d;
      double rmin = +numeric_limits<double>::infinity();
      double rmax = -numeric_limits<double>::infinity();
      surface->GetCenter(c);
      for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId) {
        surface->GetPoint(ptId, p);
        d = sqrt(pow(p[0] - c[0], 2) + pow(p[1] - c[1], 2) + pow(p[2] - c[2], 2));
        if (d < rmin) rmin = d;
        if (d > rmax) rmax = d;
      }
      double r = (type == OutputBoundingSphere ? rmax : rmin);
      const double l = two_pi * r;
      if (IsNaN(resolution)) resolution = r / 8;
      vtkNew<vtkSphereSource> sphere;
      sphere->SetCenter(c);
      sphere->SetRadius(r);
      sphere->SetThetaResolution(round(l / resolution));
      sphere->SetPhiResolution(sphere->GetThetaResolution());
      sphere->Update();
      surface = sphere->GetOutput();
    } break;

    default: {
      if (triangulate) {
        if (verbose) cout << "Triangulating surface...";
        surface = Triangulate(surface);
        if (verbose) cout << " done" << endl;
      }
    } break;
  }

  // ---------------------------------------------------------------------------
  // Write polygonal surface mesh
  if (output_name && !WritePointSet(output_name, surface)) {
    FatalError("Failed to write surface to " << output_name);
  }

  // ---------------------------------------------------------------------------
  // Write binary mask and/or implicit surface representation
  //
  // Note: vtkImplicitModeller would be suitable to generate implicit
  //       surface distance function, but unfortunately it only produces
  //       an unsigned distance map...
  if (mask_name || dmap_name) {

    ImageAttributes attr;
    if (ref_name) {
      BinaryImage ref(ref_name);
      attr = ref.Attributes();
      if (!IsNaN(resolution)) {
        attr._x  = int(ceil(attr._x * attr._dx / resolution));
        attr._y  = int(ceil(attr._y * attr._dy / resolution));
        attr._z  = int(ceil(attr._z * attr._dz / resolution));
        attr._dx = resolution;
        attr._dy = resolution;
        attr._dz = resolution;
      }
    } else {
      double bounds[6], ds = resolution;
      surface->GetBounds(bounds);
      if (IsNaN(ds)) {
        ds = sqrt(pow(bounds[1] - bounds[0], 2) +
                  pow(bounds[3] - bounds[2], 2) +
                  pow(bounds[5] - bounds[4], 2)) / 256;
      }
      attr._xorigin = bounds[0] + .5 * (bounds[1] - bounds[0]);
      attr._yorigin = bounds[2] + .5 * (bounds[3] - bounds[2]);
      attr._zorigin = bounds[4] + .5 * (bounds[5] - bounds[4]);
      attr._x       = int(ceil((bounds[1] - bounds[0]) / ds));
      attr._y       = int(ceil((bounds[3] - bounds[2]) / ds));
      attr._z       = int(ceil((bounds[5] - bounds[4]) / ds));
      attr._dx      = ds;
      attr._dy      = ds;
      attr._dz      = ds;
    }
    const int nvox = attr.NumberOfLatticePoints();

    // Generate binary inside/outside mask
    BinaryImage mask(attr);
    surface = vtkPolyData::SafeDownCast(WorldToImage(surface, &mask));
    vtkSmartPointer<vtkImageData>        vtkmask = NewVtkMask(attr._x, attr._y, attr._z);
    vtkSmartPointer<vtkImageStencilData> stencil = ImageStencil(vtkmask, surface);
    ImageStencilToMask(stencil, vtkmask);
    mask.CopyFrom(reinterpret_cast<BinaryPixel *>(vtkmask->GetScalarPointer()));
    vtkmask = NULL, stencil = NULL;

    // Invert inside/outside mask
    if (inside_out) {
      for (int vox = 0; vox < nvox; ++vox) {
        mask(vox) = BinaryPixel(mask(vox) == 0 ? 1 : 0);
      }
    }

    // Write binary inside/outside mask
    if (mask_name) mask.Write(mask_name);

    // Compute signed distance map
    if (dmap_name) {
    	// Create separate inside/outside masks
      RealImage inside_mask (attr);
      RealImage outside_mask(attr);
      for (int vox = 0; vox < nvox; ++vox) {
        if (mask(vox)) {
          inside_mask (vox) = 0.0;
          outside_mask(vox) = 1.0;
        } else {
          inside_mask (vox) = 1.0;
          outside_mask(vox) = 0.0;
        }
      }

      // Calculate Euclidean distance transforms
      typedef EuclideanDistanceTransform<RealPixel> DistanceTransform;
      DistanceTransform edt(DistanceTransform::DT_3D);
      RealImage inside_dmap, outside_dmap;

      edt.Input (&inside_mask);
      edt.Output(&inside_dmap);
      edt.Run();

      edt.Input (&outside_mask);
      edt.Output(&outside_dmap);
      edt.Run();

      for (int vox = 0; vox < nvox; ++vox) {
        inside_dmap(vox) = sqrt(outside_dmap(vox)) - sqrt(inside_dmap(vox));
      }

      // Write signed distance map
      inside_dmap.Write(dmap_name);
    }
  }

  return 0;
}

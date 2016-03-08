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

#include <mirtkGenericImage.h>
#include <mirtkPointSetUtils.h>
#include <mirtkEdgeConnectivity.h>
#include <mirtkSurfaceCollisions.h>
#include <mirtkEuclideanDistanceTransform.h>

#include <mirtkInterpolateImageFunction.h>
#include <mirtkLinearInterpolateImageFunction.hxx>

#include <vtkSmartPointer.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkImageStencilData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSphereSource.h>
#include <vtkDelaunay2D.h>
#include <vtkDelaunay3D.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkFillHolesFilter.h>
#include <vtkPointLocator.h>
#include <vtkPolyDataNormals.h>
#include <vtkCleanPolyData.h>
#include <vtkGenericCell.h>

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
  cout << "  -fillholes [<size>]    Fill holes of given maximum size (i.e., circumsphere radius)." << endl;
  cout << "  -largest [<n>]         Output largest n (default 1 if not specified) components." << endl;
  cout << "  -insphere              Output sphere which is inscribed the surface mesh." << endl;
  cout << "  -boundingsphere        Output minimum sphere which fully contains the surface mesh." << endl;
  cout << "  -o -surface <file>     Write output surface mesh to named file." << endl;
  cout << "  -mask <file>           Write binary inside/outside mask to named image file." << endl;
  cout << "  -implicit <file>       Write signed implicit surface distance to named image file." << endl;
  cout << "  -reference <image>     Reference image for :option:`-mask` or :option:`-implicit` output. (default: bounding box)" << endl;
  cout << "  -resolution <float>    Resolution of :option:`-insphere`, :option:`-mask`, or :option:`-implicit` output." << endl;
  cout << "                         (default: fraction of length of bounding box diagonal)" << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
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
  const char *output_name     = NULL;
  const char *mask_name       = NULL;
  const char *dmap_name       = NULL;
  int         operation       = vtkBooleanOperationPolyDataFilter::VTK_UNION;
  int         hull_levels     = -1;
  bool        delaunay        = false;
  bool        triangulate     = (input_names.size() > 1);
  double      max_hole_size   = .0;
  double      max_shock_dist  = numeric_limits<double>::quiet_NaN();
  const char *shock_dmap_name = NULL;
  int         nlargest        = 0;
  OutputType  type            = OutputSurface;
  const char *ref_name        = NULL;
  double      resolution      = numeric_limits<double>::quiet_NaN();
  double      tolerance       = 1e-6;

  if (NUM_POSARGS > 0) {
    if (NUM_POSARGS > 1) {
      input_names.resize(NUM_POSARGS - 1, NULL);
      for (int i = 1; i < NUM_POSARGS; ++i) input_names[i-1] = POSARG(i);
    }
    output_name = POSARG(NUM_POSARGS);
  }

  for (ALL_OPTIONS) {
    if (OPTION("-i") || OPTION("-input")) input_names.push_back(ARGUMENT);
    else if (OPTION("-o") || OPTION("-output") || OPTION("-surface")) {
      output_name = ARGUMENT;
    }
    else if (OPTION("-mask"))     mask_name = ARGUMENT;
    else if (OPTION("-implicit")) dmap_name = ARGUMENT;
    else if (OPTION("-hull")){
      if (HAS_ARGUMENT) hull_levels = atoi(ARGUMENT);
      else              hull_levels = 3;
    }
    else if (OPTION("-removeshocks")) {
      shock_dmap_name = NULL;
      if (HAS_ARGUMENT) {
        max_shock_dist = atof(ARGUMENT);
        if (HAS_ARGUMENT) shock_dmap_name = ARGUMENT;
      }
      else max_shock_dist = -1.0;
    }
    else if (OPTION("-fillholes")) {
      if (HAS_ARGUMENT) max_hole_size = atof(ARGUMENT);
      else              max_hole_size = numeric_limits<double>::infinity();
    }
    else if (OPTION("-largest")) {
      if (HAS_ARGUMENT) nlargest = atoi(ARGUMENT);
      else              nlargest = 1;
    }
    else if (OPTION("-delaunay"))    delaunay    = true;
    else if (OPTION("-triangulate")) triangulate = true;
    else if (OPTION("-insphere")) type = OutputInSphere;
    else if (OPTION("-boundingsphere")) type = OutputBoundingSphere;
    else if (OPTION("-resolution") || OPTION("-res")) resolution = atof(ARGUMENT);
    else if (OPTION("-reference")  || OPTION("-ref")) ref_name   = ARGUMENT;
    else if (OPTION("-union"))        operation = vtkBooleanOperationPolyDataFilter::VTK_UNION;
    else if (OPTION("-difference"))   operation = vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE;
    else if (OPTION("-intersection")) operation = vtkBooleanOperationPolyDataFilter::VTK_INTERSECTION;
    else if (OPTION("-tolerance") || OPTION("-tol")) tolerance = atof(ARGUMENT);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (input_names.empty()) {
    FatalError("No input point sets given!");
  }
  if (!output_name && !mask_name && !dmap_name) {
    FatalError("Missing output option (-o/-surface, -mask, and/or -implicit)!");
  }

  // ---------------------------------------------------------------------------
  // Read input point sets and extract (triangulated) surface
  Array<vtkSmartPointer<vtkPolyData> > surfaces(input_names.size());
  for (size_t i = 0; i < input_names.size(); ++i) {
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
      if (triangulate) surfaces[i] = Triangulate(surfaces[i]);
    }
  }

  // ---------------------------------------------------------------------------
  // Compute union, difference, or intersection of input surfaces
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
  }
  surfaces.clear();

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
    unique_ptr<DistanceImage> dmap;
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
    if (verbose) cout << "Extracting largest ...", cout.flush();
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

    default: break;
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

    // Generate binary inside/outside mask
    BinaryImage mask(attr);
    surface = vtkPolyData::SafeDownCast(WorldToImage(surface, &mask));
    vtkSmartPointer<vtkImageData>        vtkmask = NewVtkMask(attr._x, attr._y, attr._z);
    vtkSmartPointer<vtkImageStencilData> stencil = ImageStencil(vtkmask, surface);
    ImageStencilToMask(stencil, vtkmask);
    mask.CopyFrom(reinterpret_cast<BinaryPixel *>(vtkmask->GetScalarPointer()));
    vtkmask = NULL, stencil = NULL;

    // Write binary inside/outside mask
    if (mask_name) mask.Write(mask_name);

    // Compute signed distance map
    if (dmap_name) {
      const int nvox = attr.NumberOfLatticePoints();

    	// Invert inside mask and convert to real type
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

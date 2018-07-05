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
#include "mirtk/PointSetUtils.h"
#include "mirtk/GenericImage.h"
#include "mirtk/Stripper.h"
#include "mirtk/SurfaceBoundary.h"
#include "mirtk/ConnectedComponents.h"
#include "mirtk/SurfacePatches.h"
#include "mirtk/SurfaceRemeshing.h"
#include "mirtk/MeshSmoothing.h"
#include "mirtk/Matrix3x3.h"
#include "mirtk/Vector3.h"
#include "mirtk/Triangle.h"

#include "mirtk/Vtk.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkAppendPolyData.h"
#include "vtkImageData.h"
#include "vtkImageStencilData.h"
#include "vtkIdList.h"
#include "vtkCellArray.h"
#include "vtkPolygon.h"
#include "vtkPolyDataNormals.h"
#include "vtkCellDataToPointData.h"
#include "vtkPointLocator.h"
#include "vtkCellLocator.h"
#include "vtkPlaneSource.h"
#include "vtkIntersectionPolyDataFilter.h"
#include "vtkGenericCell.h"
#include "vtkMergePoints.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkLine.h"
#include "vtkDelaunay2D.h"
#include "vtkMatrix4x4.h"
#include "vtkMatrixToLinearTransform.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " [options]\n";
  cout << "       " << name << " [<input>...] <output> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Merge surface meshes at segmentation label boundaries. The input\n";
  cout << "  surfaces must follow closely the boundary of a segment in the given\n";
  cout << "  image segmentation. Two surfaces are then merged at the longest\n";
  cout << "  common intersection boundary. When the two surfaces share no such\n";
  cout << "  segmentation boundary, the surfaces are not connected.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    File name of input surface mesh.\n";
  cout << "  output   File name of output surface mesh.\n";
  cout << "\n";
  cout << "Options:\n";
  cout << "  -i, -input <file>...\n";
  cout << "      Input surfaces. These are added to the list of input surfaces\n";
  cout << "      to be merged after the <input> arguments.\n";
  cout << "  -o, -output <file>\n";
  cout << "      Output surface file. When this option is given, the <output>\n";
  cout << "      argument is appended to the <input> surfaces before any surface\n";
  cout << "      specified using :option:`-i`. Either this option or an <output>\n";
  cout << "      argument is required.\n";
  cout << "  -labels <file>\n";
  cout << "      Segmentation labels image. Currently required.\n";
  cout << "  -labels-percentage <value>\n";
  cout << "      Percentage of voxels within an input surface that must have the same label\n";
  cout << "      to be considered for determining boundaries between segments.\n";
  cout << "  -source-array <name>\n";
  cout << "      Add point/cell data array with the specified name with one-based\n";
  cout << "      labels corresponding to the input point set from which an output\n";
  cout << "      point/cell originates from. When the first input point set has\n";
  cout << "      a scalar array with the specified name, the labels of this first\n";
  cout << "      surface are preserved, while successive labels are offset by the\n";
  cout << "      maximum integer value of the input data array. This is useful when\n";
  cout << "      successively merging surface meshes instead of with a single execution\n";
  cout << "      of this command. (default: none)\n";
  cout << "  -point-source-array <name>\n";
  cout << "      Name of output point data :option:`-source-array`.\n";
  cout << "  -cell-source-array <name>\n";
  cout << "      Name of output cell data :option:`-source-array`.\n";
  cout << "  -join [on|off], -nojoin\n";
  cout << "      Whether to join surfaces at intersection boundaries.\n";
  cout << "      When off, the output data set has unconnected components.\n";
  cout << "      (default: on)\n";
  cout << "  -largest [on|off]\n";
  cout << "      Retain only largest surface component after :option:`-join` of boundaries. (default: off)\n";
  cout << "  -tolerance, -tol <float>\n";
  cout << "      Maximum distance of an input cell from the segmentation\n";
  cout << "      boundary to be removed. The resulting intersection boundary\n";
  cout << "      edges are joined again after the removal of these cells.\n";
  cout << "      (default: max voxel size)\n";
  cout << "  -smooth-boundaries, -boundary-smoothing [<n>]\n";
  cout << "      Number of intersection boundary edge smoothing iterations.\n";
  cout << "      (default: 3)\n";
  cout << "  -nosmooth-boundaries, -noboundary-smoothing\n";
  cout << "      Disable smoothing of boundaries, same as :option:`-smooth-boundaries` 0.\n";
  cout << "  -remeshing, -remesh [<n>]\n";
  cout << "      Enable remeshing of newly inserted triangles at intersection\n";
  cout << "      boundaries. The default number of remeshing iterations <n> is 3.\n";
  cout << "      A different value can be specified using this option or\n";
  cout << "      :option:`-remeshing-iterations`. (default: 3)\n";
  cout << "  -remeshing-iterations, -remesh-iterations\n";
  cout << "      Number of intersection triangle remeshing iterations. (default: 3)\n";
  cout << "  -noremeshing, -noremesh\n";
  cout << "      Disable remeshing of newly inserted intersection triangles,\n";
  cout << "      same as :option:`-remeshing` or :option:`-remeshing-iterations` 0.\n";
  cout << "  -edge-length <min> [<max>]\n";
  cout << "      Edge length range for remeshing newly inserted triangles at\n";
  cout << "      intersection boundaries. When min and/or max edge length not\n";
  cout << "      specified, the range is chosen based on the mean local edge\n";
  cout << "      length of the input meshes at the intersection boundaries\n";
  cout << "      minus/plus :option:`-edge-length-sigma` times the local\n";
  cout << "      standard deviation of input mesh edge lengths.\n";
  cout << "  -min-edge-length <min>\n";
  cout << "      Minimum edge length for remeshing of newly inserted triangles.\n";
  cout << "  -max-edge-length <min>\n";
  cout << "      Maximum edge length for remeshing of newly inserted triangles.\n";
  cout << "  -edge-length-sigma <scale>\n";
  cout << "      Standard deviation scaling factor of local edge length standard\n";
  cout << "      deviation used to set min/max edge length range for local remeshing\n";
  cout << "      of newly inserted triangles at intersection boundaries. (default: 2)\n";
  cout << "  -smoothing, -smooth [<n>]\n";
  cout << "      Enable smoothing of merged surface points nearby intersection boundaries.\n";
  cout << "      When :option:`-smoothing-mu` is zero or NaN, Gaussian weighting function\n";
  cout << "      is used for the mesh smoothing. Otherwise, a uniform weighting is used. (default: 0)\n";
  cout << "  -smoothing-iterations, -smooth-iterations\n";
  cout << "      Number of intersection points smoothing iterations. (default: 0)\n";
  cout << "  -nosmoothing, -nosmooth\n";
  cout << "      Disable smoothing of nearby intersection points, same as\n";
  cout << "      :option:`-smoothing` or :option:`-smoothing-iterations` 0.\n";
  cout << "  -smoothing-lambda, -smooth-lambda <float>\n";
  cout << "      Lambda parameter used for odd iterations of Laplacian point :option:`-smoothing`.\n";
  cout << "      The default value is set based on the :option:`-smooth-mu` parameter. (default: abs(mu) - .01)\n";
  cout << "      When only :option:`-smoothing-lambda` is given, :option:`-smoothing-mu` is set to\n";
  cout << "      to the same value and a Gaussian weighting function is used.\n";
  cout << "  -smoothing-mu, -smooth-mu <float>\n";
  cout << "      Lambda parameter used for even iterations of Laplacian point :option:`-smoothing`.\n";
  cout << "      The default smoothing parameters avoid shrinking of the surface mesh\n";
  cout << "      and serve mainly to relax any small self-intersections that may have\n";
  cout << "      been introduced while joining the split input surfaces. (default: -.75)\n";
  cout << "  -neighborhood, -neighbourhood, -radius <int>\n";
  cout << "      Edge connectivity radius around intersection border cells/points used\n";
  cout << "      to define the local neighborhood for which to apply the remeshing and\n";
  cout << "      smoothing operations. See also :option:`-remeshing-neighborhood`.\n";
  cout << "  -remeshing-neighborhood, -remeshing-radius <int>\n";
  cout << "      Edge connectivity radius around intersection border cells/points used\n";
  cout << "      to define neighborhood for which to apply the :option:`-remeshing`. (default: 3)\n";
  cout << "  -smoothing-neighborhood, -smoothing-radius <int>\n";
  cout << "      Edge connectivity radius around intersection border cells/points used\n";
  cout << "      to define neighborhood for which to apply the :option:`-smoothing`. (default: 2)\n";
  cout << "  -smooth-source, source-smoothing [<n>]\n";
  cout << "      Number of iterations for which to smooth the :option:`-source-array`\n";
  cout << "      labels to form smoother boundaries between intersection borders.\n";
  cout << "  -dividers [on|off], -nodividers\n";
  cout << "      Enable/disable insertion of segmentation cutting planes at mesh intersections.\n";
  cout << "      The :option:`-source-array` has unique negative labels for each such inserted\n";
  cout << "      divider surface patch which enable the extraction or removal of these\n";
  cout << "      dividers again when a genus-0 surface is required again. (default: off)\n";
  cout << "  -normals [on|off], -nonormals\n";
  cout << "      Enable/disable output of surface point and cell normals. (default: on)\n";
  cout << "  -point-normals [on|off], -nopoint-normals\n";
  cout << "      Enable/disable output of surface point normals. (default: on)\n";
  cout << "  -cell-normals [on|off], -nocell-normals\n";
  cout << "      Enable/disable output of surface cell normals. (default: on)\n";
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Constants
// =============================================================================

const char * const SOURCE_ARRAY_NAME          = "_SourceId";
const char * const INTERSECTION_ARRAY_NAME    = "_IntersectionMask";
const char * const MIN_EDGE_LENGTH_ARRAY_NAME = "_MinEdgeLength";
const char * const MAX_EDGE_LENGTH_ARRAY_NAME = "_MaxEdgeLength";

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Get binary mask of region enclosed by surface
BinaryImage InsideMask(const ImageAttributes &attr, vtkSmartPointer<vtkPolyData> surface)
{
  BinaryImage mask(attr, 1);
  surface = vtkPolyData::SafeDownCast(WorldToImage(surface, &mask));
  vtkSmartPointer<vtkImageData>        vtkmask = NewVtkMask(mask.X(), mask.Y(), mask.Z());
  vtkSmartPointer<vtkImageStencilData> stencil = ImageStencil(vtkmask, surface);
  ImageStencilToMask(stencil, vtkmask);
  mask.CopyFrom(reinterpret_cast<BinaryPixel *>(vtkmask->GetScalarPointer()));
  return mask;
}

// -----------------------------------------------------------------------------
/// Determine dominant surface label
UnorderedSet<int> InsideLabels(const GreyImage &labels, double percentage, vtkPolyData *surface)
{
  OrderedMap<int, int> hist;
  BinaryImage mask = InsideMask(labels.Attributes(), surface);
  const int nvox = labels.NumberOfSpatialVoxels();
  int size = 0;
  for (int vox = 0; vox < nvox; ++vox) {
    if (mask(vox) != 0) {
      ++hist[labels(vox)];
      ++size;
    }
  }
  UnorderedSet<int> label_set;
  if (size == 0) return label_set;
  const int min_count = ifloor(percentage * static_cast<double>(size) / 100.);
  for (const auto &bin : hist) {
    if (bin.second > min_count) {
      label_set.insert(bin.first);
    }
  }
  return label_set;
}

// -----------------------------------------------------------------------------
/// Add set of segmentation boundary points
vtkSmartPointer<vtkPolyData>
LabelBoundary(const GreyImage &labels, const UnorderedSet<int> &label_set1, const UnorderedSet<int> &label_set2)
{
  const NeighborhoodOffsets offsets(&labels, labels.Z() > 1 ? CONNECTIVITY_6 : CONNECTIVITY_4);
  const int nvox = labels.NumberOfSpatialVoxels();

  const GreyPixel *l1, *l2;
  const GreyPixel * const begin = labels.Data();
  const GreyPixel * const end   = labels.Data() + nvox;

  l1 = begin;

  bool is_boundary_voxel;
  ByteImage boundary_voxels(labels.Attributes());
  for (int vox = 0; vox < nvox; ++vox, ++l1) {
    if (label_set1.find(*l1) != label_set1.end()) {
      is_boundary_voxel = false;
      for (int n = 0; n < offsets.Size(); ++n) {
        l2 = l1 + offsets(n);
        if (begin <= l2 && l2 < end) {
          if (*l2 != *l1 && label_set2.find(*l2) != label_set2.end()) {
            is_boundary_voxel = true;
            break;
          }
        }
      }
      if (is_boundary_voxel) {
        boundary_voxels(vox) = 1;
      }
    }
  }

  ConnectedComponents<BytePixel> cc;
  cc.Input (&boundary_voxels);
  cc.Output(&boundary_voxels);
  cc.Run();

  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkCellArray> polys;
  points = vtkSmartPointer<vtkPoints>::New();
  polys  = vtkSmartPointer<vtkCellArray>::New();

  l1 = begin;
  double p[3], c[3];
  vtkIdType ptIds[4];
  int i1, j1, k1, i2, j2, k2;
  for (int vox = 0; vox < nvox; ++vox, ++l1) {
    if (boundary_voxels(vox) == 1) {
      for (int n = 0; n < offsets.Size(); ++n) {
        l2 = l1 + offsets(n);
        if (begin <= l2 && l2 < end) {
          if (*l2 != *l1 && label_set2.find(*l2) != label_set2.end()) {
            labels.IndexToVoxel(static_cast<int>(l1 - begin), i1, j1, k1);
            labels.IndexToVoxel(static_cast<int>(l2 - begin), i2, j2, k2);
            c[0] = i1, c[1] = j1, c[2] = k1;
            if (i1 != i2) {
              p[0] = .5 * (i1 + i2);
              p[1] = c[1] - .5;
              p[2] = c[2] - .5;
              labels.ImageToWorld(p[0], p[1], p[2]);
              ptIds[0] = points->InsertNextPoint(p);
              p[0] = .5 * (i1 + i2);
              p[1] = c[1] - .5;
              p[2] = c[2] + .5;
              labels.ImageToWorld(p[0], p[1], p[2]);
              ptIds[1] = points->InsertNextPoint(p);
              p[0] = .5 * (i1 + i2);
              p[1] = c[1] + .5;
              p[2] = c[2] + .5;
              labels.ImageToWorld(p[0], p[1], p[2]);
              ptIds[2] = points->InsertNextPoint(p);
              p[0] = .5 * (i1 + i2);
              p[1] = c[1] + .5;
              p[2] = c[2] - .5;
              labels.ImageToWorld(p[0], p[1], p[2]);
              ptIds[3] = points->InsertNextPoint(p);
            } else if (j1 != j2) {
              p[0] = c[0] - .5;
              p[1] = .5 * (j1 + j2);
              p[2] = c[2] - .5;
              labels.ImageToWorld(p[0], p[1], p[2]);
              ptIds[0] = points->InsertNextPoint(p);
              p[0] = c[0] - .5;
              p[1] = .5 * (j1 + j2);
              p[2] = c[2] + .5;
              labels.ImageToWorld(p[0], p[1], p[2]);
              ptIds[1] = points->InsertNextPoint(p);
              p[0] = c[0] + .5;
              p[1] = .5 * (j1 + j2);
              p[2] = c[2] + .5;
              labels.ImageToWorld(p[0], p[1], p[2]);
              ptIds[2] = points->InsertNextPoint(p);
              p[0] = c[0] + .5;
              p[1] = .5 * (j1 + j2);
              p[2] = c[2] - .5;
              labels.ImageToWorld(p[0], p[1], p[2]);
              ptIds[3] = points->InsertNextPoint(p);
            } else {
              p[0] = c[0] - .5;
              p[1] = c[1] - .5;
              p[2] = .5 * (k1 + k2);
              labels.ImageToWorld(p[0], p[1], p[2]);
              ptIds[0] = points->InsertNextPoint(p);
              p[0] = c[0] - .5;
              p[1] = c[1] + .5;
              p[2] = .5 * (k1 + k2);
              labels.ImageToWorld(p[0], p[1], p[2]);
              ptIds[1] = points->InsertNextPoint(p);
              p[0] = c[0] + .5;
              p[1] = c[1] + .5;
              p[2] = .5 * (k1 + k2);
              labels.ImageToWorld(p[0], p[1], p[2]);
              ptIds[2] = points->InsertNextPoint(p);
              p[0] = c[0] + .5;
              p[1] = c[1] - .5;
              p[2] = .5 * (k1 + k2);
              labels.ImageToWorld(p[0], p[1], p[2]);
              ptIds[3] = points->InsertNextPoint(p);
            }
            polys->InsertNextCell(4, ptIds);
          }
        }
      }
    }
  }

  vtkSmartPointer<vtkPolyData> boundary;
  boundary = vtkSmartPointer<vtkPolyData>::New();
  boundary->SetPoints(points);
  boundary->SetPolys(polys);

  vtkNew<vtkCleanPolyData> cleaner;
  SetVTKInput(cleaner, boundary);
  cleaner->ConvertPolysToLinesOff();
  cleaner->ConvertLinesToPointsOff();
  cleaner->ConvertStripsToPolysOff();
  cleaner->PointMergingOn();
  cleaner->ToleranceIsAbsoluteOn();
  cleaner->SetAbsoluteTolerance(1e-9);
  cleaner->Update();
  boundary = cleaner->GetOutput();

  if (debug) {
    static int ncall = 0; ++ncall;
    char fname[64];
    snprintf(fname, 64, "debug_boundary_voxels_%d.nii.gz", ncall);
    boundary_voxels.Write(fname);
    snprintf(fname, 64, "debug_boundary_%d.vtp", ncall);
    WritePolyData(fname, boundary);
  }

  return boundary;
}

// -----------------------------------------------------------------------------
/// Minimum squared distance of point to segmentation boundary surface
inline double MinSquaredDistance(vtkAbstractCellLocator *cut, const Point &p)
{
  double    x[3], dist2;
  vtkIdType cellId;
  int       subId;
  cut->FindClosestPoint(const_cast<Point &>(p), x, cellId, subId, dist2);
  return dist2;
}

// -----------------------------------------------------------------------------
/// Minimum squared distance of point to surface boundary
inline double MinSquaredDistance(const BoundarySegment &seg, const Point &p)
{
  const auto i = seg.FindClosestPoint(p);
  return p.SquaredDistance(seg.Point(i));
}

// -----------------------------------------------------------------------------
/// Minimum squared distance of surface boundary segment to segmentation boundary surface
double MinSquaredDistance(const BoundarySegment &seg, vtkAbstractCellLocator *cut)
{
  double    x[3], dist2, min_dist2 = inf;
  vtkIdType cellId;
  int       subId;
  for (int i = 0; i < seg.NumberOfPoints(); ++i) {
    cut->FindClosestPoint(seg.Point(i), x, cellId, subId, dist2);
    if (dist2 < min_dist2) min_dist2 = dist2;
  }
  return min_dist2;
}

// -----------------------------------------------------------------------------
/// Hausdorff distance of two boundary segments
double HausdorffDistance(const BoundarySegment &a, const BoundarySegment &b)
{
  double d_ab = 0., d_ba = 0.;
  for (int i = 0; i < a.NumberOfPoints(); ++i) {
    d_ab = max(d_ab, MinSquaredDistance(b, a.Point(i)));
  }
  for (int i = 0; i < b.NumberOfPoints(); ++i) {
    d_ba = max(d_ba, MinSquaredDistance(a, b.Point(i)));
  }
  return sqrt(max(d_ab, d_ba));
}

// -----------------------------------------------------------------------------
/// Delete boundary points with only one cell
void DeleteSingleBoundaryPointCells(SurfaceBoundary &boundary, vtkAbstractCellLocator *cut, double max_dist2)
{
  vtkPolyData * const surface = boundary.Surface();

  unsigned short ncells;
  vtkIdType      *cells;

  while (true) {
    int ndel = 0;
    surface->BuildLinks();
    for (int j = 0; j < boundary.NumberOfSegments(); ++j) {
      auto &seg = boundary.Segment(j);
      if (MinSquaredDistance(seg, cut) < max_dist2) {
        for (int i = 0; i < seg.NumberOfPoints(); ++i) {
          surface->GetPointCells(seg.PointId(i), ncells, cells);
          if (ncells == 1) {
            surface->RemoveCellReference(cells[0]);
            surface->DeleteCell(cells[0]);
            ++ndel;
          }
        }
      }
    }
    if (ndel == 0) break;
    surface->RemoveDeletedCells();
    boundary = SurfaceBoundary(surface);
  }
}

// -----------------------------------------------------------------------------
/// Smooth intersection boundaries
void SmoothBoundaries(SurfaceBoundary &boundary, vtkAbstractCellLocator *cut, double max_dist2, int niter = 1)
{
  if (niter <= 0) return;
  PointSet new_points;
  for (int j = 0; j < boundary.NumberOfSegments(); ++j) {
    auto &seg = boundary.Segment(j);
    if (MinSquaredDistance(seg, cut) < max_dist2) {
      new_points.Resize(seg.NumberOfPoints());
      for (int iter = 0; iter < niter; ++iter) {
        for (int i = 0; i < seg.NumberOfPoints(); ++i) {
          new_points(i) = (seg.Point(i-1) + seg.Point(i) + seg.Point(i+1)) / 3.;
        }
        for (int i = 0; i < seg.NumberOfPoints(); ++i) {
          seg.Point(i, new_points(i));
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
void InterpolateCellData(vtkCellData *cd, vtkPointData *pd, const Array<int> &type,
                         vtkIdType cellId, vtkIdType ptId1, vtkIdType ptId2, vtkIdType ptId3)
{
  double v1, v2, v3;
  Array<double> tuple;
  for (int i = 0; i < cd->GetNumberOfArrays(); ++i) {
    vtkDataArray * const dst = cd->GetArray(i);
    vtkDataArray * const src = pd->GetArray(i);
    tuple.resize(src->GetNumberOfComponents());
    if (type[i] == 2) {
      for (int j = 0; j < src->GetNumberOfComponents(); ++j) {
        tuple[j] = 0.;
      }
    } else if (type[i] == 1) {
      for (int j = 0; j < src->GetNumberOfComponents(); ++j) {
        v1 = src->GetComponent(ptId1, j);
        v2 = src->GetComponent(ptId2, j);
        v3 = src->GetComponent(ptId3, j);
        tuple[j] = (v2 == v3 ? v2 : v1);
      }
    } else {
      for (int j = 0; j < src->GetNumberOfComponents(); ++j) {
        v1 = src->GetComponent(ptId1, j);
        v2 = src->GetComponent(ptId2, j);
        v3 = src->GetComponent(ptId3, j);
        tuple[j] = (v1 + v2 + v3) / 3.;
      }
    }
    dst->InsertTuple(cellId, tuple.data());
  }
}

// -----------------------------------------------------------------------------
/// Convert cell data to point data
vtkSmartPointer<vtkPointData> CellToPointData(vtkDataSet *dataset, const Array<int> &type)
{
  vtkCellData * const cd = dataset->GetCellData();
  // Use VTK filter to convert cell data to point data using averaging
  vtkNew<vtkCellDataToPointData> cd_to_pd;
  SetVTKInput(cd_to_pd, dataset);
  cd_to_pd->PassCellDataOff();
  cd_to_pd->Update();
  vtkSmartPointer<vtkPointData> pd = cd_to_pd->GetOutput()->GetPointData();
  // Fix label data
  double v;
  UnorderedMap<double, int> bins;
  UnorderedMap<double, int>::iterator bin;
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
    if (type[i] == 1) {
      vtkDataArray * const src = cd->GetArray(i);
      vtkDataArray * const dst = pd->GetArray(i);
      for (vtkIdType ptId = 0; ptId < dataset->GetNumberOfPoints(); ++ptId) {
        dataset->GetPointCells(ptId, cellIds);
        for (int j = 0; j < src->GetNumberOfComponents(); ++j) {
          bins.clear();
          for (vtkIdType k = 0; k < cellIds->GetNumberOfIds(); ++k) {
            v = src->GetComponent(cellIds->GetId(k), j);
            bin = bins.find(v);
            if (bin == bins.end()) bins[v] = 1;
            else bin->second += 1;
          }
          int    max_cnt = 0;
          double max_val = 0.;
          for (bin = bins.begin(); bin != bins.end(); ++bin) {
            if (bin->second > max_cnt) {
              max_val = bin->first;
              max_cnt = bin->second;
            }
          }
          dst->SetComponent(ptId, j, max_val);
        }
      }
    }
  }
  return pd;
}

// -----------------------------------------------------------------------------
/// Join intersected components at the segmentation boundary
void JoinBoundaries(SurfaceBoundary &boundary, vtkAbstractCellLocator *cut, double max_dist2, double max_hdist)
{
  vtkIdType cellId, ptIds[3];

  vtkPolyData  * const surface = boundary.Surface();
  vtkCellArray * const polys   = surface->GetPolys();
  vtkCellData  * const cd      = surface->GetCellData();

  Array<int> cd_type(cd->GetNumberOfArrays(), 0);
  for (int i = 0; i < cd->GetNumberOfArrays(); ++i) {
    vtkDataArray * const arr = cd->GetArray(i);
    if (arr->GetName()) {
      if (strcmp(arr->GetName(), SOURCE_ARRAY_NAME) == 0) {
        cd_type[i] = 2;
      } else if (IsCategoricalArrayName(arr->GetName())) {
        cd_type[i] = 1;
      }
    }
  }
  vtkSmartPointer<vtkPointData> cd_as_pd = CellToPointData(surface, cd_type);

  while (boundary.NumberOfSegments() > 1) {

    Array<double> dist2(boundary.NumberOfSegments());
    for (int j = 0; j < boundary.NumberOfSegments(); ++j) {
      dist2[j] = MinSquaredDistance(boundary.Segment(j), cut);
    }

    int idx1, idx2;
    Array<Array<double>> hdist(boundary.NumberOfSegments());
    for (idx1 = 0; idx1 < boundary.NumberOfSegments(); ++idx1) {
      hdist[idx1].resize(boundary.NumberOfSegments(), inf);
      if (dist2[idx1] <= max_dist2) {
        const auto &seg1 = boundary.Segment(idx1);
        if (seg1.NumberOfPoints() > 1) {
          for (idx2 = idx1 + 1; idx2 < boundary.NumberOfSegments(); ++idx2) {
            if (dist2[idx2] <= max_dist2) {
              const auto &seg2 = boundary.Segment(idx2);
              if (seg2.NumberOfPoints() > 1) {
                hdist[idx1][idx2] = HausdorffDistance(seg1, seg2);
              }
            }
          }
        }
      }
    }
    Array<double> min_hdist(boundary.NumberOfSegments());
    Array<int>    min_idx2 (boundary.NumberOfSegments());
    for (idx1 = 0; idx1 < boundary.NumberOfSegments(); ++idx1) {
      min_idx2 [idx1] = IncreasingOrder(hdist[idx1]).front();
      min_hdist[idx1] = hdist[idx1][min_idx2[idx1]];
    }
    idx1 = IncreasingOrder(min_hdist).front();
    if (min_hdist[idx1] > max_hdist) break;
    idx2 = min_idx2[idx1];

    if (verbose > 1) {
      cout << "  Joining boundary segments " << idx1+1 << " and " << idx2+1
           << " out of " << boundary.NumberOfSegments() << " remaining surface boundaries"
           << " (Hausdorff distance = " << hdist[idx1][idx2] << ")" << endl;
    }

    // Get boundary segments to be joined
    auto &seg1 = boundary.Segment(idx1);
    auto &seg2 = boundary.Segment(idx2);

    const int i0 = 0;
    const int j0 = seg2.FindClosestPoint(seg1.Point(i0));

    // TODO: Determine in which direction to traverse each boundary segment
    // such that orientation of newly added triangles is consistent with the
    // orientation of the boundary triangles
    const int di = 1;

    // Construct list of new triangles based on either dj=-1 or dj=+1 traversal
    // direction and select those triangles with the least cost/error
    const int max_iter = seg1.NumberOfPoints() + seg2.NumberOfPoints();
    double l1, l2, cost1 = 0., cost2 = 0.;
    Array<Vector3D<vtkIdType>> tris1, tris2;
    for (int dj = -1; dj <= 1; dj += 2) {
      int i = i0, j = j0;
      seg1.ClearSelection();
      seg2.ClearSelection();
      auto &tris = (dj == -1 ? tris1 : tris2);
      auto &cost = (dj == -1 ? cost1 : cost2);
      tris.reserve(max_iter);
      for (int iter = 0; iter < max_iter; ++iter) {
        l1 = (seg1.IsSelected(i + di) ? inf : seg2.Point(j).SquaredDistance(seg1.Point(i + di)));
        l2 = (seg2.IsSelected(j + dj) ? inf : seg1.Point(i).SquaredDistance(seg2.Point(j + dj)));
        if (IsInf(l1) && IsInf(l2)) break;
        if (l1 <= l2) {
          cost += l1;
          ptIds[0] = seg1.PointId(i);
          ptIds[1] = seg1.PointId(i + di);
          ptIds[2] = seg2.PointId(j);
          i = seg1.IndexModuloNumberOfPoints(i + di);
          seg1.SelectPoint(i);
        } else {
          cost += l2;
          ptIds[0] = seg1.PointId(i);
          ptIds[1] = seg2.PointId(j + dj);
          ptIds[2] = seg2.PointId(j);
          j = seg2.IndexModuloNumberOfPoints(j + dj);
          seg2.SelectPoint(j);
        }
        tris.push_back(Vector3D<vtkIdType>(ptIds));
      }
    }

    // Add triangles joining the two boundary segments
    for (auto &tri : (cost1 <= cost2 ? tris1 : tris2)) {
      cellId = polys->InsertNextCell(3);
      polys->InsertCellPoint(tri._x);
      polys->InsertCellPoint(tri._y);
      polys->InsertCellPoint(tri._z);
      InterpolateCellData(cd, cd_as_pd, cd_type, cellId, tri._x, tri._y, tri._z);
    }

    surface->DeleteLinks();
    surface->DeleteCells();

    #if 0
      boundary = SurfaceBoundary(surface);
    #else
      break;
    #endif
  }

  // Free excessive memory
  polys->Squeeze();
  cd->Squeeze();
}

// -----------------------------------------------------------------------------
/// Merge two surface meshes at given segmentation boundary
vtkSmartPointer<vtkPolyData>
Merge(vtkPolyData *s1, vtkPolyData *s2, vtkPolyData *label_boundary, double tol, int smooth, bool join)
{
  const double tol2      = tol * tol;
  const double max_dist2 =  4. * tol2;
  const double max_hdist = 10. * tol;

  double         p[3], x[3], dist2;
  vtkIdType      cellId;
  int            subId;
  unsigned short ncells;
  vtkIdType      *cells;

  // Build links
  s1->BuildLinks();
  s2->BuildLinks();

  // Locator to find distance to segmentation boundary
  vtkSmartPointer<vtkAbstractCellLocator> cut;
  cut = vtkSmartPointer<vtkCellLocator>::New();
  cut->SetDataSet(label_boundary);
  cut->BuildLocator();

  // Mark cells close to segmentation boundary as deleted
  for (vtkIdType ptId = 0; ptId < s1->GetNumberOfPoints(); ++ptId) {
    s1->GetPoint(ptId, p);
    cut->FindClosestPoint(p, x, cellId, subId, dist2);
    if (dist2 < tol2) {
      s1->GetPointCells(ptId, ncells, cells);
      for (unsigned short i = 0; i < ncells; ++i) {
        s1->DeleteCell(cells[i]);
      }
      s1->DeletePoint(ptId);
    }
  }

  for (vtkIdType ptId = 0; ptId < s2->GetNumberOfPoints(); ++ptId) {
    s2->GetPoint(ptId, p);
    cut->FindClosestPoint(p, x, cellId, subId, dist2);
    if (dist2 < tol2) {
      s2->GetPointCells(ptId, ncells, cells);
      for (unsigned short i = 0; i < ncells; ++i) {
        s2->DeleteCell(cells[i]);
      }
      s2->DeletePoint(ptId);
    }
  }

  // Remove deleted cells
  s1->RemoveDeletedCells();
  s2->RemoveDeletedCells();

  // Combine surfaces into one polygonal data set
  vtkNew<vtkAppendPolyData> appender;
  AddVTKInput(appender, s1);
  AddVTKInput(appender, s2);
  appender->Update();

  // Get combined output data set
  vtkSmartPointer<vtkPolyData> output;
  output = appender->GetOutput();

  // Clean intersection boundaries
  //
  // This is in particular needed to avoid issues with these cells when
  // smoothing and/or joining the intersection boundaries...
  SurfaceBoundary boundary(output);
  DeleteSingleBoundaryPointCells(boundary, cut, max_dist2);

  // Smooth intersection boundaries
  SmoothBoundaries(boundary, cut, max_dist2, smooth);

  // Join intersection boundaries
  if (join) JoinBoundaries(boundary, cut, max_dist2, max_hdist);

  return output;
}

// -----------------------------------------------------------------------------
/// Add array with given point source label
int AddPointSourceArray(vtkPolyData * const surface, size_t idx)
{
  vtkSmartPointer<vtkDataArray> arr;
  arr = NewVtkDataArray(VTK_SHORT, surface->GetNumberOfPoints(), 1, SOURCE_ARRAY_NAME);
  arr->FillComponent(0, static_cast<double>(idx));
  surface->GetPointData()->RemoveArray(SOURCE_ARRAY_NAME);
  return surface->GetPointData()->AddArray(arr);
}

// -----------------------------------------------------------------------------
/// Add array with given cell source label
int AddCellSourceArray(vtkPolyData * const surface, size_t idx)
{
  vtkSmartPointer<vtkDataArray> arr;
  arr = NewVtkDataArray(VTK_SHORT, surface->GetNumberOfCells(), 1, SOURCE_ARRAY_NAME);
  arr->FillComponent(0, static_cast<double>(idx));
  surface->GetCellData()->RemoveArray(SOURCE_ARRAY_NAME);
  return surface->GetCellData()->AddArray(arr);
}

// -----------------------------------------------------------------------------
/// Calculate point normals of surface mesh while fixing vertex order if necessary
vtkSmartPointer<vtkPolyData>
CalculateNormals(vtkPolyData *surface, bool point_normals, bool cell_normals)
{
  vtkNew<vtkPolyDataNormals> normals;
  SetVTKInput(normals, surface);
  normals->AutoOrientNormalsOff();
  normals->SplittingOff();
  normals->NonManifoldTraversalOff();
  normals->ConsistencyOn();
  normals->SetComputePointNormals(point_normals);
  normals->SetComputeCellNormals(cell_normals);
  normals->Update();
  return normals->GetOutput();
}

// -----------------------------------------------------------------------------
/// Fix normals of border edge points
void FixBorderPointNormals(vtkPolyData *surface)
{
  vtkDataArray * const cnormals = surface->GetCellData ()->GetNormals();
  vtkDataArray * const pnormals = surface->GetPointData()->GetNormals();
  vtkDataArray * const csource  = surface->GetCellData ()->GetArray(SOURCE_ARRAY_NAME);
  vtkDataArray * const psource  = surface->GetPointData()->GetArray(SOURCE_ARRAY_NAME);

  unsigned short ncells, n;
  vtkIdType      *cells;
  Vector3        pn, cn;

  for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId) {
    if (psource->GetComponent(ptId, 0) == 0.) {
      surface->GetPointCells(ptId, ncells, cells);
      n  = 0;
      pn = 0.;
      for (unsigned short i = 0; i < ncells; ++i) {
        if (csource->GetComponent(cells[i], 0) >= 0.) {
          cnormals->GetTuple(cells[i], cn);
          pn += cn, ++n;
        }
      }
      if (n > 0) {
        pn /= n;
        pn.Normalize();
        pnormals->SetTuple(ptId, pn);
      }
    }
  }
}

// -----------------------------------------------------------------------------
/// Get cell data array with cells within vicinity of joint boundaries masked
vtkSmartPointer<vtkDataArray> IntersectionCellMask(vtkPolyData * const surface, int nconn = 3)
{
  const vtkIdType n = surface->GetNumberOfCells();

  vtkCellData  * const cd     = surface->GetCellData();
  vtkDataArray * const source = cd->GetArray(SOURCE_ARRAY_NAME);

  unsigned short ncells;
  vtkIdType npts, *pts, *cells;

  vtkSmartPointer<vtkDataArray> mask;
  mask = NewVtkDataArray(VTK_UNSIGNED_CHAR, n, 1, INTERSECTION_ARRAY_NAME);
  for (vtkIdType cellId = 0; cellId < n; ++cellId) {
    mask->SetComponent(cellId, 0, source->GetComponent(cellId, 0) == 0. ? 1. : 0.);
  }
  for (int iter = 0; iter < nconn; ++iter) {
    vtkSmartPointer<vtkDataArray> next;
    next.TakeReference(mask->NewInstance());
    next->DeepCopy(mask);
    for (vtkIdType cellId = 0; cellId < n; ++cellId) {
      if (mask->GetComponent(cellId, 0) != 0.) {
        surface->GetCellPoints(cellId, npts, pts);
        for (vtkIdType i = 0; i < npts; ++i) {
          surface->GetPointCells(pts[i], ncells, cells);
          for (unsigned short j = 0; j < ncells; ++j) {
            next->SetComponent(cells[j], 0, 1.);
          }
        }
      }
    }
    mask = next;
  }

  return mask;
}

// -----------------------------------------------------------------------------
/// Get point data array with points within vicinity of joint boundaries masked
vtkSmartPointer<vtkDataArray> IntersectionPointMask(vtkPolyData * const surface, int nconn = 3)
{
  const vtkIdType n = surface->GetNumberOfPoints();

  vtkCellData  * const cd     = surface->GetCellData();
  vtkDataArray * const source = cd->GetArray(SOURCE_ARRAY_NAME);

  unsigned short ncells;
  vtkIdType npts, *pts, *cells;

  vtkSmartPointer<vtkDataArray> mask;
  mask = NewVtkDataArray(VTK_UNSIGNED_CHAR, n, 1, INTERSECTION_ARRAY_NAME);
  mask->FillComponent(0, 0.);
  for (vtkIdType ptId = 0; ptId < n; ++ptId) {
    surface->GetPointCells(ptId, ncells, cells);
    for (unsigned short i = 0; i < ncells; ++i) {
      if (source->GetComponent(cells[i], 0) == 0.) {
        mask->SetComponent(ptId, 0, 1.);
        break;
      }
    }
  }
  for (int iter = 0; iter < nconn; ++iter) {
    vtkSmartPointer<vtkDataArray> next;
    next.TakeReference(mask->NewInstance());
    next->DeepCopy(mask);
    for (vtkIdType ptId = 0; ptId < n; ++ptId) {
      if (mask->GetComponent(ptId, 0) != 0.) {
        surface->GetPointCells(ptId, ncells, cells);
        for (unsigned short i = 0; i < ncells; ++i) {
          surface->GetCellPoints(cells[i], npts, pts);
          for (vtkIdType j = 0; j < npts; ++j) {
            next->SetComponent(pts[j], 0, 1.);
          }
        }
      }
    }
    mask = next;
  }

  return mask;
}

// -----------------------------------------------------------------------------
/// Determine edge length range within surface intersection
void GetEdgeLengthRange(vtkPolyData * const surface, vtkDataArray * const mask,
                        double &lmin, double &lmax, double sd = 2.5)
{
  vtkCellData  * const cd     = surface->GetCellData();
  vtkDataArray * const source = cd->GetArray(SOURCE_ARRAY_NAME);

  double min_edge_length = +inf;
  double max_edge_length = -inf;

  int    n = 0;
  double a[3], b[3], ab[3], l2, lsum = 0., l2sum = 0.;
  vtkIdType npts, *pts;
  for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
    if (mask->GetComponent(cellId, 0) != 0. && source->GetComponent(cellId, 0) != 0.) {
      surface->GetCellPoints(cellId, npts, pts);
      for (vtkIdType i = 0, j; i < npts; ++i) {
        j = (i + 1) % npts;
        surface->GetPoint(pts[i], a);
        surface->GetPoint(pts[j], b);
        vtkMath::Subtract(b, a, ab);
        l2 = vtkMath::Dot(ab, ab);
        lsum  += sqrt(l2);
        l2sum += l2;
        ++n;
      }
    }
  }
  if (n > 0) {
    const double mean  = lsum / n;
    const double sigma = sqrt(l2sum / n - mean * mean);
    min_edge_length = max(0., mean - sd * sigma);
    max_edge_length = mean + sd * sigma;
  }

  if (IsNaN(lmin)) {
    lmin = (min_edge_length > max_edge_length ? 0. : min_edge_length);
  }
  if (IsNaN(lmax)) {
    lmax = (min_edge_length > max_edge_length ? inf : max_edge_length);
  }
}

// -----------------------------------------------------------------------------
/// Principal component analysis of point set
bool PrincipalDirections(vtkPointSet * const points, Vector3 dir[3], bool twod = false)
{
  if (points->GetNumberOfPoints() == 0) return false;

  Point c, p;
  points->GetCenter(c);

  Matrix3x3 covar(0.);
  for (vtkIdType ptId = 0; ptId < points->GetNumberOfPoints(); ++ptId) {
    points->GetPoint(ptId, p);
    vtkMath::Subtract(p, c, p);
    for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c) {
      covar[r][c] += p[r] * p[c];
    }
  }

  Vector3       axis[3];
  Array<double> eigval(3);
  covar.EigenSolveSymmetric(eigval.data(), axis);

  eigval[0] = abs(eigval[0]);
  eigval[1] = abs(eigval[1]);
  eigval[2] = abs(eigval[2]);
  Array<int> order = DecreasingOrder(eigval);
  if (eigval[order[twod ? 1 : 2]] < 1e-6) return false;

  dir[0] = axis[order[0]], dir[0].Normalize();
  dir[1] = axis[order[1]], dir[1].Normalize();
  dir[2] = axis[order[2]], dir[2].Normalize();
  return true;
}

// -----------------------------------------------------------------------------
/// Structure with cutting plane parameters
struct PlaneAttributes
{
  Point   o; ///< Origin/Center point
  Vector3 n; ///< Normal vector
  Vector3 x; ///< In-plane x axis
  Vector3 y; ///< In-plane y axis
  Vector3 &z = n; ///< Out of plane z axis, i.e., plane normal
};

// -----------------------------------------------------------------------------
/// Get cutting plane from segmentation boundary
PlaneAttributes CuttingPlaneAttributes(vtkSmartPointer<vtkPointSet> points)
{
  MIRTK_START_TIMING();
  Vector3 dir[3], p;
  if (!PrincipalDirections(points, dir)) {
    FatalError("Failed to compute principal directions of point set");
  }
  double u, v, minu = +inf, maxu = -inf, minv = +inf, maxv = -inf;
  for (vtkIdType ptId = 0; ptId < points->GetNumberOfPoints(); ++ptId) {
    points->GetPoint(ptId, p);
    u = dir[0].Dot(p);
    v = dir[1].Dot(p);
    if (u < minu) minu = u;
    if (u > maxu) maxu = u;
    if (v < minv) minv = v;
    if (v > maxv) maxv = v;
  }
  dir[0] *= .5 * (maxu - minu);
  dir[1] *= .5 * (maxv - minv);
  dir[2].Normalize();
  PlaneAttributes attr;
  points->GetCenter(attr.o);
  attr.n = dir[2];
  attr.x = dir[0];
  attr.y = dir[1];
  MIRTK_DEBUG_TIMING(2, "CuttingPlaneAttributes");
  return attr;
}

// -----------------------------------------------------------------------------
/// Get cutting plane from segmentation boundary
vtkSmartPointer<vtkPolyData> CuttingPlane(const PlaneAttributes &attr, double tn, double rx, double ry, double margin)
{
  MIRTK_START_TIMING();

  // Calculate rotation matrices in local plane coordinates
  const double cx = cos(rx * rad_per_deg), sx = sin(rx * rad_per_deg);
  const double cy = cos(ry * rad_per_deg), sy = sin(ry * rad_per_deg);
  const Matrix3x3 Rx(1., 0., 0., 0., cx, -sx, 0., sx, cx);
  const Matrix3x3 Ry(cy, 0., sy, 0., 1., 0., -sy, 0., cy);

  // Rotate local plane axes vectors
  Vector3 dx(1., 0., 0.), dy(0., 1., 0.), n(0., 0., 1.);
  dx = Rx * Ry * dx;
  dy = Rx * Ry * dy;
  n  = Rx * Ry * n;

  // Map to world vectors of unit length
  Vector3 e1(attr.x), e2(attr.y), e3(attr.z);
  e1.Normalize(), e2.Normalize(), e3.Normalize();
  dx = dx._x * e1 + dx._y * e2 + dx._z * e3;
  dy = dy._x * e1 + dy._y * e2 + dy._z * e3;
  n  = n ._x * e1 + n ._y * e2 + n ._z * e3;

  // Add plane margin to initial plane extend
  dx *= (attr.x.Length() + 2. * margin);
  dy *= (attr.y.Length() + 2. * margin);

  // Translate plane along original normal vector
  Point o = attr.o + tn * attr.n;

  // vtkPlaneSource origin is in corner, not center
  o  -= dx;
  o  -= dy;
  dx *= 2.;
  dy *= 2.;

  // Tesselate finite plane
  vtkSmartPointer<vtkPolyData> plane;
  vtkNew<vtkPlaneSource> source;
  source->SetXResolution(1);
  source->SetYResolution(1);
  source->SetOrigin(o);
  source->SetPoint1(o + dx);
  source->SetPoint2(o + dy);
  source->Update();
  plane = Triangulate(source->GetOutput());

  MIRTK_DEBUG_TIMING(2, "CuttingPlane");
  return plane;
}

// -----------------------------------------------------------------------------
/// Intersect two polygonal data sets and return single intersection polygon
vtkSmartPointer<vtkPolyData> LargestClosedIntersection(vtkPolyData *s1, vtkPolyData *s2)
{
  MIRTK_START_TIMING();

  vtkSmartPointer<vtkPolyData> cut;
  {
    MIRTK_START_TIMING();
    vtkNew<vtkIntersectionPolyDataFilter> intersection;
    intersection->SplitFirstOutputOff();
    intersection->SplitSecondOutputOff();
    SetNthVTKInput(intersection, 0, s1);
    SetNthVTKInput(intersection, 1, s2);

    vtkNew<vtkCleanPolyData> merger;
    SetVTKConnection(merger, intersection);
    merger->ConvertStripsToPolysOff();
    merger->ConvertPolysToLinesOff();
    merger->ConvertLinesToPointsOff();
    merger->PointMergingOn();
    merger->ToleranceIsAbsoluteOn();
    merger->SetAbsoluteTolerance(1e-12);

    merger->Update();
    cut = merger->GetOutput();
    MIRTK_DEBUG_TIMING(3, "LargestClosedIntersection (intersect)");
  }

  if (debug > 2) {
    static int iter = 0; ++iter;
    char fname[64];
    snprintf(fname, 64, "debug_cutting_lines_%d.vtp", iter);
    WritePolyData(fname, cut);
  }

  {
    MIRTK_START_TIMING();
    cut->BuildLinks();

    unsigned short ncells;
    vtkIdType ptId, nbrId, cellId, *cells;
    vtkNew<vtkIdList> ptIds;
    ptIds->Allocate(2);

    Stack<vtkIdType> activePtIds;
    for (ptId = 0; ptId < cut->GetNumberOfPoints(); ++ptId) {
      cut->GetPointCells(ptId, ncells, cells);
      if (ncells == 1) activePtIds.push(ptId);
    }
    while (!activePtIds.empty()) {
      ptId = activePtIds.top(), activePtIds.pop();
      cut->GetPointCells(ptId, ncells, cells);
      if (ncells == 1) {
        cellId = cells[0];
        cut->GetCellPoints(cellId, ptIds.GetPointer());
        cut->RemoveCellReference(cellId);
        cut->DeleteCell(cellId);
        for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
          nbrId = ptIds->GetId(i);
          if (nbrId != ptId) {
            cut->GetPointCells(nbrId, ncells, cells);
            if (ncells == 1) activePtIds.push(nbrId);
          }
        }
      }
    }
    cut->RemoveDeletedCells();
    vtkNew<vtkCleanPolyData> merger;
    SetVTKInput(merger, cut);
    merger->ConvertStripsToPolysOff();
    merger->ConvertPolysToLinesOff();
    merger->ConvertLinesToPointsOn();
    merger->PointMergingOn();
    merger->ToleranceIsAbsoluteOn();
    merger->SetAbsoluteTolerance(1e-12);
    merger->Update();
    merger->GetOutput()->SetVerts(nullptr);
    cut = merger->GetOutput();
    MIRTK_DEBUG_TIMING(3, "LargestClosedIntersection (clean)");
  }

  {
    MIRTK_START_TIMING();
    vtkNew<vtkPolyDataConnectivityFilter> cc;
    SetVTKInput(cc, cut);
    cc->ScalarConnectivityOff();
    cc->SetExtractionModeToLargestRegion();
    cc->Update();
    cut = cc->GetOutput();
    MIRTK_DEBUG_TIMING(3, "LargestClosedIntersection (largest)");
  }

  MIRTK_DEBUG_TIMING(2, "LargestClosedIntersection");
  return cut;
}

// -----------------------------------------------------------------------------
/// Merge points and make line strips
vtkSmartPointer<vtkPolyData> LineStrips(vtkSmartPointer<vtkPolyData> cut)
{
  vtkNew<vtkCleanPolyData> cut_merger;
  SetVTKInput(cut_merger, cut);
  cut_merger->ConvertStripsToPolysOff();
  cut_merger->ConvertPolysToLinesOff();
  cut_merger->ConvertLinesToPointsOn();
  cut_merger->PointMergingOn();
  cut_merger->ToleranceIsAbsoluteOn();
  cut_merger->SetAbsoluteTolerance(1e-12);
  cut_merger->Update();
  cut_merger->GetOutput()->SetVerts(nullptr);
  cut = cut_merger->GetOutput();

  Stripper stripper;
  stripper.Input(cut);
  stripper.Run();
  return stripper.Output();
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkDataArray> NonBoundaryPointsMask(vtkSmartPointer<vtkPolyData> input)
{
  vtkSmartPointer<vtkDataArray> mask;
  mask = NewVtkDataArray(VTK_UNSIGNED_CHAR, input->GetNumberOfPoints(), 1);
  mask->FillComponent(0, 1.);
  UnorderedSet<int> ptIds = BoundaryPoints(input);
  for (auto ptId : ptIds) mask->SetComponent(ptId, 0, 0.);
  return mask;
}

// -----------------------------------------------------------------------------
/// Make divider surface from intersection curve
vtkSmartPointer<vtkPolyData> Divider(vtkSmartPointer<vtkPolyData> cut)
{
  vtkSmartPointer<vtkPolyData> divider = LineStrips(cut);
  if (debug) {
    static int callId = 0; ++callId;
    char fname[64];
    snprintf(fname, 64, "debug_split_surface_lines_%d.vtp", callId);
    WritePolyData(fname, divider);
  }
  vtkIdType npts, *pts, stripPts = 0, *stripIds = nullptr;
  for (vtkIdType cellId = 0; cellId < divider->GetLines()->GetNumberOfCells(); ++cellId) {
    divider->GetLines()->GetCell(cellId, npts, pts);
    if (npts > stripPts) {
      stripPts = npts;
      stripIds = pts;
    }
  }
  if (stripPts == 0) {
    Throw(ERR_LogicError, __FUNCTION__, "Expected at least one contiguous intersection line");
  }
  if (stripPts <= 2) {
    Throw(ERR_LogicError, __FUNCTION__, "Expected polygon with more than two points");
  }
  if (stripIds[0] != stripIds[stripPts-1]) {
    Throw(ERR_LogicError, __FUNCTION__, "Expected closed intersection polygon");
  }
  vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
  polys->Allocate(polys->EstimateSize(1, stripPts-1));
  polys->InsertNextCell(stripPts-1, stripIds);
  divider->SetLines(nullptr);
  divider->SetPolys(polys);
  divider->DeleteCells();
  return divider;
}

// -----------------------------------------------------------------------------
/// Compute Delaunay triangulation of interior of divider polygon after projection to 2D plane
vtkSmartPointer<vtkPolyData> TesselateDivider(vtkSmartPointer<vtkPolyData> divider, double ds)
{
  if (divider->GetNumberOfPoints() < 2) return divider;

  const double min_dist2 = pow(.2 * ds, 2);

  int       subId;
  double    dist2, reverse_dist2;
  vtkIdType otherId, npts, *pts;
  Vector3   c, p, q, x, dir[3];

  // Compute principle directions of divider polygon plane
  const bool twod = true;
  PrincipalDirections(divider, dir, twod);

  // Get center of divider polygon
  divider->GetCenter(c);
  double offset = 0.;
  for (vtkIdType ptId = 0; ptId < divider->GetNumberOfPoints(); ++ptId) {
    divider->GetPoint(ptId, p);
    offset += dir[2].Dot(p - c);
  }
  offset /= divider->GetNumberOfPoints();
  c += offset * dir[2];

  // Set up homogeneous transformation from 3D world to divider plane
  vtkSmartPointer<vtkMatrix4x4> matrix;
  matrix = vtkSmartPointer<vtkMatrix4x4>::New();
  matrix->Identity();
  for (int i = 0; i < 3; ++i) {
    matrix->SetElement(i, 0, dir[i]._x);
    matrix->SetElement(i, 1, dir[i]._y);
    matrix->SetElement(i, 2, dir[i]._z);
    matrix->SetElement(i, 3, -dir[i].Dot(c));
  }

  vtkSmartPointer<vtkMatrixToLinearTransform> transform;
  transform = vtkSmartPointer<vtkMatrixToLinearTransform>::New();
  transform->SetInput(matrix);
  transform->Update();

  // Determine extent of divider polygon within divider plane
  Point minq, maxq;
  divider->GetPoint(0, minq);
  divider->GetPoint(0, maxq);
  for (vtkIdType ptId = 1; ptId < divider->GetNumberOfPoints(); ++ptId) {
    divider->GetPoint(ptId, p);
    transform->TransformPoint(p, q);
    if (q._x < minq._x) minq._x = q._x;
    if (q._y < minq._y) minq._y = q._y;
    if (q._z < minq._z) minq._z = q._z;
    if (q._x > maxq._x) maxq._x = q._x;
    if (q._y > maxq._y) maxq._y = q._y;
    if (q._z > maxq._z) maxq._z = q._z;
  }

  const double margin = ds;
  minq._x -= margin, maxq._x += margin;
  minq._y -= margin, maxq._y += margin;

  int nx = iceil((maxq._x - minq._x) / ds) + 1;
  int ny = iceil((maxq._y - minq._y) / ds) + 1;

  if (nx % 2 == 0) ++nx;
  if (ny % 2 == 0) ++ny;

  minq._x = minq._x + .5 * (maxq._x - minq._x) - ((nx - 1) / 2) * ds;
  maxq._x = minq._x + (nx - 1) * ds + 1e-12;

  minq._y = minq._y + .5 * (maxq._y - minq._y) - ((ny - 1) / 2) * ds;
  maxq._y = minq._y + (ny - 1) * ds + 1e-12;

  // Add additional discrete divider polygon plane grid points
  vtkNew<vtkPointLocator> locator;
  locator->SetDataSet(divider);
  locator->BuildLocator();

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetDataTypeToDouble();
  points->Allocate(divider->GetNumberOfPoints() + nx * ny);

  for (vtkIdType ptId = 0; ptId < divider->GetNumberOfPoints(); ++ptId) {
    divider->GetPoint(ptId, p);
    points->InsertNextPoint(p);
  }

  q._z = 0.;
  for (q._y = minq._y; q._y <= maxq._y; q._y += ds)
  for (q._x = minq._x; q._x <= maxq._x; q._x += ds) {
    p = c + q._x * dir[0] + q._y * dir[1];
    otherId = locator->FindClosestPoint(p);
    if (otherId >= 0) {
      divider->GetPoint(otherId, x);
      if ((p - x).SquaredLength() < min_dist2) {
        p = x;
      }
    }
    points->InsertNextPoint(p);
  }

  // Compute constrained Delanauy triangulation of divider polygon
  vtkSmartPointer<vtkPolyData> pointset;
  pointset = vtkSmartPointer<vtkPolyData>::New();
  pointset->SetPoints(points);

  vtkSmartPointer<vtkPolyData> boundary;
  boundary = vtkSmartPointer<vtkPolyData>::New();
  boundary->SetPoints(points);
  boundary->SetPolys(divider->GetPolys());

  vtkNew<vtkDelaunay2D> delaunay;
  delaunay->SetInputData(pointset);
  delaunay->SetSourceData(boundary);
  delaunay->SetAlpha(0.);
  delaunay->SetOffset(1.);
  delaunay->SetTolerance(0.);
  delaunay->SetTransform(transform);
  delaunay->Update();

  vtkSmartPointer<vtkPolyData> output;
  output = vtkSmartPointer<vtkPolyData>::New();
  output->DeepCopy(delaunay->GetOutput());

  // TODO: Determine correct order of polygon points before vtkDelaunay2D::Update
  vtkNew<vtkCellLocator> cell_locator;
  cell_locator->SetDataSet(delaunay->GetOutput());
  cell_locator->BuildLocator();
  cell_locator->FindClosestPoint(c, x, otherId, subId, dist2);

  boundary->GetPolys()->GetCell(0, npts, pts);
  vtkNew<vtkIdList> reverse_pts;
  reverse_pts->Allocate(npts);
  for (vtkIdType i = npts-1; i >= 0; --i) {
    reverse_pts->InsertNextId(pts[i]);
  }
  vtkSmartPointer<vtkCellArray> reverse_polys;
  reverse_polys = vtkSmartPointer<vtkCellArray>::New();
  reverse_polys->Allocate(reverse_polys->EstimateSize(1, npts));
  reverse_polys->InsertNextCell(reverse_pts.GetPointer());

  vtkSmartPointer<vtkPolyData> reverse_boundary;
  reverse_boundary = vtkSmartPointer<vtkPolyData>::New();
  reverse_boundary->SetPoints(points);
  reverse_boundary->SetPolys(reverse_polys);

  delaunay->SetSourceData(reverse_boundary);
  delaunay->Update();

  cell_locator->SetDataSet(delaunay->GetOutput());
  cell_locator->BuildLocator();
  cell_locator->FindClosestPoint(c, x, otherId, subId, reverse_dist2);

  if (reverse_dist2 < dist2) {
    output = delaunay->GetOutput();
  }

  // Remove unused points
  vtkNew<vtkCleanPolyData> cleaner;
  SetVTKInput(cleaner, output);
  cleaner->ConvertStripsToPolysOff();
  cleaner->ConvertPolysToLinesOff();
  cleaner->ConvertLinesToPointsOff();
  cleaner->PointMergingOff();
  cleaner->Update();
  output = cleaner->GetOutput();

  // Smooth divider
  // (needed when boundary points were snapped to surface mesh)
  MeshSmoothing smoother;
  smoother.Input(output);
  smoother.Mask(NonBoundaryPointsMask(output));
  smoother.Weighting(MeshSmoothing::Combinatorial);
  smoother.AdjacentValuesOnlyOn();
  smoother.SmoothPointsOn();
  smoother.NumberOfIterations(3);
  smoother.Lambda(1.);
  smoother.Run();
  output = smoother.Output();

  return output;
}

// -----------------------------------------------------------------------------
/// Check if given intersection curve is acceptable
bool IsValidIntersection(vtkSmartPointer<vtkPolyData> cut, double tol)
{
  if (cut->GetNumberOfCells() == 0) return false;
  #if 1
    return true;
  #else
    vtkNew<vtkCleanPolyData> merger;
    SetVTKInput(merger, cut);
    merger->ConvertStripsToPolysOff();
    merger->ConvertPolysToLinesOff();
    merger->ConvertLinesToPointsOn();
    merger->PointMergingOn();
    merger->ToleranceIsAbsoluteOn();
    merger->SetAbsoluteTolerance(2. * tol);
    merger->Update();
    merger->GetOutput()->SetVerts(nullptr);
    cut = merger->GetOutput();

    cut->BuildLinks();

    unsigned short ncells;
    vtkIdType      *cells;
    for (vtkIdType ptId = 0; ptId < cut->GetNumberOfPoints(); ++ptId) {
      cut->GetPointCells(ptId, ncells, cells);
      if (ncells != 2 || cut->GetCellType(cells[0]) != VTK_LINE || cut->GetCellType(cells[1]) != VTK_LINE) {
        return false;
      }
    }

    return true;
  #endif
}

// -----------------------------------------------------------------------------
/// Smooth line strip
void SmoothLineStrip(vtkSmartPointer<vtkPolyData> cut, int niter = 1)
{
  int   i1, i2;
  Point p1, p2, p3;

  cut->BuildLinks();

  unsigned short ncells;
  vtkIdType      *cells, *pts1, *pts2, npts;

  vtkSmartPointer<vtkPoints> points;
  points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(cut->GetNumberOfPoints());

  for (int iter = 0; iter < niter; ++iter) {
    for (vtkIdType ptId = 0; ptId < cut->GetNumberOfPoints(); ++ptId) {
      cut->GetPointCells(ptId, ncells, cells);
      if (ncells == 2 && cut->GetCellType(cells[0]) == VTK_LINE && cut->GetCellType(cells[1]) == VTK_LINE) {
        cut->GetCellPoints(cells[0], npts, pts1);
        cut->GetCellPoints(cells[1], npts, pts2);
        i1 = (pts1[0] == ptId ? 1 : 0);
        i2 = (pts2[0] == ptId ? 1 : 0);
        cut->GetPoint(pts1[i1], p1);
        cut->GetPoint(ptId,     p2);
        cut->GetPoint(pts2[i2], p3);
        points->SetPoint(ptId, (p1 + p2 + p3) / 3.);
      }
    }
    vtkSmartPointer<vtkPoints> tmp = cut->GetPoints();
    cut->SetPoints(points);
    points = tmp;
  }
}

// -----------------------------------------------------------------------------
/// Intersection incircle radius
double LineStripIncircleRadius(vtkSmartPointer<vtkPolyData> cut)
{
  Point p, c;
  cut->GetCenter(c);
  double r = inf;
  for (vtkIdType ptId = 0; ptId < cut->GetNumberOfPoints(); ++ptId) {
    cut->GetPoint(ptId, p);
    r = min(r, p.Distance(c));
  }
  return r;
}

// -----------------------------------------------------------------------------
/// Intersection incircle radius
double AverageLineStripRadius(vtkSmartPointer<vtkPolyData> cut)
{
  Point p, c;
  cut->GetCenter(c);
  double r = 0.;
  for (vtkIdType ptId = 0; ptId < cut->GetNumberOfPoints(); ++ptId) {
    cut->GetPoint(ptId, p);
    r += p.Distance(c);
  }
  return r / cut->GetNumberOfPoints();
}

// -----------------------------------------------------------------------------
/// Length of line strip
double LineStripLength(vtkSmartPointer<vtkPolyData> cut)
{
  Point p1, p2;
  vtkIdType npts, *pts;
  double l = 0.;

  cut->BuildLinks();

  for (vtkIdType cellId = 0; cellId < cut->GetNumberOfCells(); ++cellId) {
    cut->GetCellPoints(cellId, npts, pts);
    if (cut->GetCellType(cellId) == VTK_LINE && npts == 2) {
      cut->GetPoint(pts[0], p1);
      cut->GetPoint(pts[1], p2);
      l += p1.Distance(p2);
    }
  }

  return l;
}

// -----------------------------------------------------------------------------
/// Average curvature of line strip
double LineStripCurvature(vtkSmartPointer<vtkPolyData> cut)
{
  int     i1, i2, n = 0;
  double  l1, l2, angle, curv = 0.;
  Vector3 e1, e2;
  Point   p1, p2, p3;

  cut->GetPoint(cut->GetNumberOfPoints()-1, p1);
  cut->GetPoint(0, p2);
  e1 = p2 - p1;
  l1 = e1.Length();

  cut->BuildLinks();

  unsigned short ncells;
  vtkIdType      *cells, *pts1, *pts2, npts;

  for (vtkIdType ptId = 0; ptId < cut->GetNumberOfPoints(); ++ptId) {
    cut->GetPointCells(ptId, ncells, cells);
    if (ncells == 2 && cut->GetCellType(cells[0]) == VTK_LINE && cut->GetCellType(cells[1]) == VTK_LINE) {
      cut->GetCellPoints(cells[0], npts, pts1);
      cut->GetCellPoints(cells[1], npts, pts2);
      i1 = (pts1[0] == ptId ? 1 : 0);
      i2 = (pts2[0] == ptId ? 1 : 0);
      cut->GetPoint(pts1[i1], p1);
      cut->GetPoint(ptId,     p2);
      cut->GetPoint(pts2[i2], p3);
      e1 = p2 - p1;
      e2 = p3 - p2;
      l1 = e1.Length();
      l2 = e2.Length();
      angle = acos(e1.Dot(e2) / (l1 * l2));
      curv += angle * angle / (l1 + l2);
      ++n;
    }
  }

  return curv / n;
}

// -----------------------------------------------------------------------------
/// Maximum angle made up by adjacent line segments
double MaxLineStripAngle(vtkSmartPointer<vtkPolyData> cut)
{
  double  angle = 0.;
  int     i1, i2;
  Vector3 e1, e2;
  Point   p1, p2, p3;

  cut->BuildLinks();

  unsigned short ncells;
  vtkIdType      *cells, *pts1, *pts2, npts;

  for (vtkIdType ptId = 0; ptId < cut->GetNumberOfPoints(); ++ptId) {
    cut->GetPointCells(ptId, ncells, cells);
    if (ncells == 2 && cut->GetCellType(cells[0]) == VTK_LINE && cut->GetCellType(cells[1]) == VTK_LINE) {
      cut->GetCellPoints(cells[0], npts, pts1);
      cut->GetCellPoints(cells[1], npts, pts2);
      i1 = (pts1[0] == ptId ? 1 : 0);
      i2 = (pts2[0] == ptId ? 1 : 0);
      cut->GetPoint(pts1[i1], p1);
      cut->GetPoint(ptId,     p2);
      cut->GetPoint(pts2[i2], p3);
      e1 = p2 - p1;
      e2 = p3 - p2;
      e1.Normalize();
      e2.Normalize();
      angle = max(angle, acos(e1.Dot(e2)));
    }
  }

  return angle;
}

// -----------------------------------------------------------------------------
double DividerArea(vtkSmartPointer<vtkPolyData> divider)
{
  vtkPolygon *polygon = vtkPolygon::SafeDownCast(divider->GetCell(0));
  return polygon->ComputeArea();
}

// -----------------------------------------------------------------------------
double DividerError(vtkSmartPointer<vtkPolyData> divider, vtkSmartPointer<vtkPoints> points)
{
  vtkIdType cellId;
  int subId;
  double sum = 0., dist2;
  Point p, x;
  vtkNew<vtkCellLocator> locator;
  locator->SetDataSet(divider);
  locator->BuildLocator();
  for (vtkIdType ptId = 0; ptId < points->GetNumberOfPoints(); ++ptId) {
    points->GetPoint(ptId, p);
    locator->FindClosestPoint(p, x, cellId, subId, dist2);
    sum += dist2;
  }
  return sqrt(sum / points->GetNumberOfPoints());
}

// -----------------------------------------------------------------------------
/// Find cutting plane which intersects surface with one closed intersection boundary
bool FindCuttingPlane(vtkSmartPointer<vtkPolyData> surface,
                      vtkSmartPointer<vtkPointSet> boundary,
                      vtkSmartPointer<vtkPolyData> &plane,
                      vtkSmartPointer<vtkPolyData> &cut,
                      double max_offset, double ds = 0., double tol = 0.)
{
  static int call = 0; ++call;

  MIRTK_START_TIMING();

  if (ds <= 0.) ds = AverageEdgeLength(surface);

  // Abort search after the specified number of suitable planes were found
  const int  max_suitable_planes  = 10;
  const bool prefer_default_plane = false;

  // Compute cutting plane from segmentation boundary points
  const PlaneAttributes attr = CuttingPlaneAttributes(boundary);

  const double max_angle   = 20.;
  const double delta_angle = 5.;

  const double delta_offset = max_offset / 4.;

  const double margin   =  5. * ds;
  const double bbmargin = 20. * ds;

  double    r, l;
  double    best_tn, best_rx, best_ry, best_l = inf;
  vtkIdType best_ct = 0;

  // Remove cells which are never cut to speed up vtkPolyDataIntersectionFilter
  Point center;
  double bounds[6], cell_bounds[6];
  plane = CuttingPlane(attr, 0., 0., 0., 0.);
  plane->GetBounds(bounds);
  plane->GetCenter(center);
  for (int i = 0; i < 6; i += 2) {
    bounds[i  ] -= bbmargin;
    bounds[i+1] += bbmargin;
  }
  plane = nullptr;

  vtkSmartPointer<vtkPolyData> cells;
  cells.TakeReference(surface->NewInstance());
  cells->ShallowCopy(surface);
  cells->DeleteCells();
  cells->BuildCells();
  for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
    surface->GetCell(cellId)->GetBounds(cell_bounds);
    if (cell_bounds[0] > bounds[1] || cell_bounds[1] < bounds[0] ||
        cell_bounds[2] > bounds[3] || cell_bounds[3] < bounds[2] ||
        cell_bounds[4] > bounds[5] || cell_bounds[5] < bounds[4]) {
      cells->DeleteCell(cellId);
    }
  }
  cells->RemoveDeletedCells();

  if (debug > 1) {
    char fname[64];
    snprintf(fname, 64, "debug_cutting_cells_%d.vtp", call);
    WritePolyData(fname, cells);
  }

  const auto prev_cout_precision = cout.precision();
  const auto prev_cout_flags     = cout.flags();

  const streamsize offset_precision = 3;
  const streamsize angle_precision  = 0;

  // Try different normal offsets and plane rotations in the order of preference
  int n = 0, m = 0;
  bool abort = false;
  for (double angle = 0., rx, ry; angle <= max_angle && !abort; angle += delta_angle) {
    for (int dim = (angle > 0. ? 0 : 1); dim < 2 && !abort; ++dim) {
      if (dim == 0) {
        rx = angle;
        ry = 0.;
      } else {
        rx = 0.;
        ry = angle;
      }
      for (double tn = 0.; tn <= max_offset && !abort; tn += delta_offset) {
        for (double sry = (ry > 0. ? -1. : 1.); sry <= 1. && !abort; sry += 2.)
        for (double srx = (rx > 0. ? -1. : 1.); srx <= 1. && !abort; srx += 2.)
        for (double stn = (tn > 0. ? -1. : 1.); stn <= 1. && !abort; stn += 2.) {
          if (verbose > 4) {
            cout << "    Cutting plane with: offset = " << fixed
                     << setprecision(offset_precision) << setw(offset_precision+3) << stn * tn
                     << ", alpha = " << setprecision(angle_precision) << setw(angle_precision+3) << srx * rx
                     << ", beta = "  << setprecision(angle_precision) << setw(angle_precision+3) << sry * ry
                     << " ...";
          }
          plane  = CuttingPlane(attr, stn * tn, srx * rx, sry * ry, margin);
          cut    = LargestClosedIntersection(cells, plane);
          if (IsValidIntersection(cut, tol)) {
            if (verbose == 4) {
              cout << "    Cutting plane with: offset = " << fixed
                   << setprecision(offset_precision) << setw(offset_precision+3) << stn * tn
                   << ", alpha = " << setprecision(angle_precision) << setw(angle_precision+3) << srx * rx
                   << ", beta = "  << setprecision(angle_precision) << setw(angle_precision+3) << sry * ry
                   << " ...";
            }
            bool accept = false;
            SmoothLineStrip(cut, 2);
            r = LineStripIncircleRadius(cut);
            if (center.Distance(cut->GetCenter()) < r) {
              ++m;
              if (verbose == 3) {
                cout << "    Cutting plane with: offset = " << fixed
                     << setprecision(offset_precision) << setw(offset_precision+3) << stn * tn
                     << ", alpha = " << setprecision(angle_precision) << setw(angle_precision+3) << srx * rx
                     << ", beta = "  << setprecision(angle_precision) << setw(angle_precision+3) << sry * ry
                     << " ...";
              }
              l = LineStripLength(cut);
              if (best_ct == 0 || l < .95 * best_l) {
                accept = true;
              }
              if (accept) {
                if (verbose > 2) cout << " accepted";
                best_tn = stn * tn;
                best_rx = srx * rx;
                best_ry = sry * ry;
                best_l  = l;
                best_ct = cut->GetNumberOfCells();
                if (prefer_default_plane && best_tn == 0. && best_rx == 0. && best_ry == 0.) {
                  abort = true;
                }
              } else {
                if (verbose > 2) cout << " rejected";
              }
              if (verbose > 2) {
                vtkSmartPointer<vtkPolyData> divider = Divider(cut);
                double area  = DividerArea(divider);
                double error = DividerError(divider, boundary->GetPoints());
                cout << ": perimeter = " << setprecision(5) << setw(8) << l;
                cout << ", area = "      << setprecision(5) << setw(8) << area;
                cout << ", RMS = "       << setprecision(5) << setw(8) << error;
                cout << endl;
              }
              if (m >= max_suitable_planes && best_ct > 0) {
                abort = true;
              }
            } else {
              if (verbose > 3) cout << " rejected" << endl;
            }
            if (debug > 1) {
              if (debug > 2 || center.Distance(cut->GetCenter()) < r) {
                ++n;
                char fname[64];
                const char *suffix = accept ? "accepted" : "rejected";
                snprintf(fname, 64, "debug_cutting_plane_%d_%d_%s.vtp", call, n, suffix);
                WritePolyData(fname, plane);
                snprintf(fname, 64, "debug_cutting_lines_%d_%d_%s.vtp", call, n, suffix);
                WritePolyData(fname, cut);
              }
            }
          } else {
            if (verbose > 4) cout << " invalid" << endl;
          }
        }
      }
    }
  }

  // Pick best plane
  if (best_ct > 0) {
    if (verbose > 1) {
      cout << "    Found cutting plane with: offset = " << fixed
           << setprecision(offset_precision) << setw(offset_precision+3) << best_tn
           << ", alpha = " << setprecision(angle_precision) << setw(angle_precision+3) << best_rx
           << ", beta = "  << setprecision(angle_precision) << setw(angle_precision+3) << best_ry << endl;
    }
    plane  = CuttingPlane(attr, best_tn, best_rx, best_ry, margin);
    cut    = LargestClosedIntersection(surface, plane);
  }

  cout.precision(prev_cout_precision);
  cout.flags(prev_cout_flags);

  MIRTK_DEBUG_TIMING(2, "FindCuttingPlane");
  return best_ct > 0;
}

// -----------------------------------------------------------------------------
struct EdgeIntersection
{
  vtkIdType i; // Edge intersection point index
  double    t; // Edge interpolation parameter

  EdgeIntersection() : i(-1), t(0.) {}
  EdgeIntersection(vtkIdType i, double t) : i(i), t(t) {}
};

// -----------------------------------------------------------------------------
struct TriangleIntersections
{
  EdgeIntersection edge[3];
};

typedef UnorderedMap<vtkIdType, TriangleIntersections> TriangleIntersectionsMap;

// -----------------------------------------------------------------------------
inline void InsertEdgeIntersection(TriangleIntersectionsMap &intersections,
                                   vtkPolyData *surface, vtkIdType cellId,
                                   int edge, Vector3 &x, double t)
{
  vtkIdType         cellPts, *cellPtIds, newPtId, nbrId, nbrPts, *nbrPtIds;
  int               nbrEdge;
  double            nbrParam;
  vtkNew<vtkIdList> cellIds;

  auto it = intersections.find(cellId);
  if (it == intersections.end()) {
    it = intersections.insert(MakePair(cellId, TriangleIntersections())).first;
  }
  auto &intersection = it->second.edge[edge];
  if (intersection.i < 0) {
    newPtId = surface->GetPoints()->InsertNextPoint(const_cast<Vector3 &>(x));
    intersection.i = newPtId;
    intersection.t = t;
  } else {
    newPtId = intersection.i;
    surface->GetPoint(newPtId, x);
  }

  surface->GetCellPoints(cellId, cellPts, cellPtIds);
  const vtkIdType ptId1 = cellPtIds[edge];
  const vtkIdType ptId2 = cellPtIds[(edge + 1) % cellPts];
  surface->GetCellEdgeNeighbors(cellId, ptId1, ptId2, cellIds.GetPointer());

  for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
    nbrId   = cellIds->GetId(j);
    nbrEdge = -1;
    surface->GetCellPoints(nbrId, nbrPts, nbrPtIds);
    for (vtkIdType k = 0; k < nbrPts; ++k) {
      if (nbrPtIds[k] == ptId1 && nbrPtIds[(k + 1) % nbrPts] == ptId2) {
        nbrEdge  = k;
        nbrParam = t;
        break;
      }
      if (nbrPtIds[k] == ptId2 && nbrPtIds[(k + 1) % nbrPts] == ptId1) {
        nbrEdge  = k;
        nbrParam = 1. - t;
        break;
      }
    }
    if (nbrEdge != -1) {
      auto it = intersections.find(nbrId);
      if (it == intersections.end()) {
        it = intersections.insert(MakePair(nbrId, TriangleIntersections())).first;
      }
      auto &intersection = it->second.edge[nbrEdge];
      if (intersection.i < 0) {
        intersection.i = newPtId;
        intersection.t = nbrParam;
      }
    }
  }
}

// -----------------------------------------------------------------------------
inline void InterpolateEdge(vtkPointData *outputPD, vtkPointData *inputPD,
                            vtkIdType ptId1, vtkIdType ptId2,
                            const EdgeIntersection &intersection)
{
  outputPD->InterpolateEdge(inputPD, intersection.i, ptId1, ptId2, intersection.t);
}

// -----------------------------------------------------------------------------
/// Bisect first edge of triangle
void Bisect(vtkPolyData *input, vtkIdType cellId,
            vtkIdType ptId1, vtkIdType ptId2, vtkIdType ptId3,
            const EdgeIntersection &intersection, vtkPolyData *output)
{
  vtkCellArray * const polys    = output->GetPolys();
  vtkPointData * const inputPD  = input ->GetPointData();
  vtkPointData * const outputPD = output->GetPointData();
  vtkCellData  * const inputCD  = input ->GetCellData();
  vtkCellData  * const outputCD = output->GetCellData();

  const vtkIdType &newPtId1 = ptId1;
  const vtkIdType &newPtId2 = intersection.i;
  const vtkIdType &newPtId3 = ptId2;
  const vtkIdType &newPtId4 = ptId3;
  vtkIdType       newCellId;

  InterpolateEdge(outputPD, inputPD, ptId1, ptId2, intersection);

  newCellId = polys->InsertNextCell(3);
  polys->InsertCellPoint(newPtId1);
  polys->InsertCellPoint(newPtId2);
  polys->InsertCellPoint(newPtId4);
  outputCD->CopyData(inputCD, cellId, newCellId);

  newCellId = polys->InsertNextCell(3);
  polys->InsertCellPoint(newPtId2);
  polys->InsertCellPoint(newPtId3);
  polys->InsertCellPoint(newPtId4);
  outputCD->CopyData(inputCD, cellId, newCellId);

  input->DeleteCell(cellId);
}

// -----------------------------------------------------------------------------
/// Bisect first and second edge of triangle
void Trisect(vtkPolyData *input, vtkIdType cellId,
             vtkIdType ptId1, vtkIdType ptId2, vtkIdType ptId3,
             const EdgeIntersection &intersection1,
             const EdgeIntersection &intersection2,
             vtkPolyData *output)
{
  vtkCellArray * const polys    = output->GetPolys();
  vtkPointData * const inputPD  = input ->GetPointData();
  vtkPointData * const outputPD = output->GetPointData();
  vtkCellData  * const inputCD  = input ->GetCellData();
  vtkCellData  * const outputCD = output->GetCellData();

  const vtkIdType &newPtId1 = ptId1;
  const vtkIdType &newPtId2 = intersection1.i;
  const vtkIdType &newPtId3 = ptId2;
  const vtkIdType &newPtId4 = intersection2.i;
  const vtkIdType &newPtId5 = ptId3;
  vtkIdType       newCellId;

  InterpolateEdge(outputPD, inputPD, ptId1, ptId2, intersection1);
  InterpolateEdge(outputPD, inputPD, ptId2, ptId3, intersection2);

  newCellId = polys->InsertNextCell(3);
  polys->InsertCellPoint(newPtId1);
  polys->InsertCellPoint(newPtId2);
  polys->InsertCellPoint(newPtId5);
  outputCD->CopyData(inputCD, cellId, newCellId);

  newCellId = polys->InsertNextCell(3);
  polys->InsertCellPoint(newPtId2);
  polys->InsertCellPoint(newPtId3);
  polys->InsertCellPoint(newPtId4);
  outputCD->CopyData(inputCD, cellId, newCellId);

  newCellId = polys->InsertNextCell(3);
  polys->InsertCellPoint(newPtId4);
  polys->InsertCellPoint(newPtId5);
  polys->InsertCellPoint(newPtId2);
  outputCD->CopyData(inputCD, cellId, newCellId);

  input->DeleteCell(cellId);
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData>
AddClosedIntersectionDivider(vtkPolyData *surface, vtkPolyData *cut, double tol = 0.)
{
  const double tol2 = tol * tol;

  int       bisect[3], subId, snapId;
  double    pcoords[3], weights[3], dist2, t1, t2;
  vtkIdType npts, *pts, lineId, cellId, ptId;
  Vector3   pt[3], x[2], p, q, closestPt;

  Array<double>                 dists;
  Array<int>                    order;
  vtkNew<vtkIdList>             ptIds, cellIds, edge1, edge2;
  vtkNew<vtkGenericCell>        cell;
  vtkSmartPointer<vtkDataArray> arr;

  ptIds->Allocate(10);
  edge1->Allocate(2);
  edge2->Allocate(2);

  // Determine edge length for tesselation of divider polygon
  const double mean_edge_length = AverageEdgeLength(surface);

  // Build needed auxiliary structures
  cut->BuildCells();
  surface->BuildLinks();

  // Polygonal data set of intersected surface cells
  vtkPoints * const points = surface->GetPoints();

  vtkSmartPointer<vtkPolyData>  split;
  vtkSmartPointer<vtkCellArray> polys;

  polys = vtkSmartPointer<vtkCellArray>::New();
  polys->Allocate(polys->EstimateSize(3 * cut->GetNumberOfCells(), 3));

  split = vtkSmartPointer<vtkPolyData>::New();
  split->SetPoints(points);
  split->SetPolys(polys);

  split->GetPointData()->InterpolateAllocate(surface->GetPointData());
  split->GetCellData()->CopyAllocate(surface->GetCellData());

  // Map (intersected) surface cell ID to intersection line ID(s)
  UnorderedMap<vtkIdType, List<vtkIdType>> cellIdToLineIdsMap;
  vtkDataArray * const inputIds = cut->GetCellData()->GetArray("Input0CellID");
  for (lineId = 0; lineId < cut->GetNumberOfCells(); ++lineId) {
    cellId = static_cast<vtkIdType>(inputIds->GetComponent(lineId, 0));
    cellIdToLineIdsMap[cellId].push_back(lineId);
  }
  cut->GetCellData()->RemoveArray("Input0CellID");
  cut->GetCellData()->RemoveArray("Input1CellID");

  // Determine and insert new edge intersection points
  TriangleIntersectionsMap intersections;
  for (auto &entry : cellIdToLineIdsMap) {

    // Get intersected surface cell
    cellId = entry.first;
    surface->GetCell(cellId, cell.GetPointer());

    // Merge intersection lines at non-edge cell points, these are resulting
    // from the triangular tesselation of the intersection plane needed by
    // the vtkIntersectionPolyDataFilter
    if (entry.second.size() > 1) {
      ptIds->Reset();
      for (const auto &lineId : entry.second) {
        cut->GetCellPoints(lineId, npts, pts);
        for (int i = 0; i < npts; ++i) {
          ptIds->InsertUniqueId(pts[i]);
        }
      }
      dists.resize(ptIds->GetNumberOfIds());
      for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
        ptId = ptIds->GetId(i);
        cut->GetPoint(ptId, p);
        cell->EvaluatePosition(p, nullptr, subId, pcoords, dist2, weights);
        cell->CellBoundary(subId, pcoords, edge1.GetPointer());
        surface->GetPoint(edge1->GetId(0), x[0]);
        surface->GetPoint(edge1->GetId(1), x[1]);
        dists[i] = vtkLine::DistanceToLine(p, x[0], x[1], t1, closestPt);
      }
      order = IncreasingOrder(dists);
      edge1->SetNumberOfIds(2);
      edge1->SetId(0, ptIds->GetId(order[0]));
      edge1->SetId(1, ptIds->GetId(order[1]));
      auto lineIdIt = entry.second.begin();
      cut->ReplaceCell(*lineIdIt, 2, edge1->GetPointer(0));
      ++lineIdIt;
      while (lineIdIt != entry.second.end()) {
        cut->DeleteCell(*lineIdIt);
        ++lineIdIt;
      }
      entry.second.resize(1);
    }

    // Determine which edges are being intersected
    surface->GetCellPoints(cellId, npts, pts);
    surface->GetPoint(pts[0], pt[0]);
    surface->GetPoint(pts[1], pt[1]);
    surface->GetPoint(pts[2], pt[2]);

    lineId = entry.second.front();
    cut->GetCellPoints(lineId, ptIds.GetPointer());
    cut->GetPoint(ptIds->GetId(0), p);
    cut->GetPoint(ptIds->GetId(1), q);

    // Snap first intersection line point to either edge point or node
    cell->EvaluatePosition(p, nullptr, subId, pcoords, dist2, weights);
    cell->CellBoundary(subId, pcoords, edge1.GetPointer());
    surface->GetPoint(edge1->GetId(0), x[0]);
    surface->GetPoint(edge1->GetId(1), x[1]);

    snapId = -1;
    if (tol2 > 0.) {
      dists.resize(2);
      dists[0] = (p - x[0]).SquaredLength();
      dists[1] = (p - x[1]).SquaredLength();
      snapId = (dists[0] <= dists[1] ? 0 : 1);
      if (dists[snapId] < tol2) {
        p = x[snapId], t1 = 0.;
        edge1->SetId(snapId == 0 ? 1 : 0, edge1->GetId(snapId));
      }
    }
    if (snapId == -1) {
      vtkLine::DistanceToLine(p, x[0], x[1], t1, closestPt);
      t1 = clamp(t1, 0., 1.);
      p  = closestPt;
    }

    // Snap second intersection line point to either edge point or node
    cell->EvaluatePosition(q, nullptr, subId, pcoords, dist2, weights);
    cell->CellBoundary(subId, pcoords, edge2.GetPointer());
    surface->GetPoint(edge2->GetId(0), x[0]);
    surface->GetPoint(edge2->GetId(1), x[1]);

    snapId = -1;
    if (tol2 > 0.) {
      dists.resize(2);
      dists[0] = (q - x[0]).SquaredLength();
      dists[1] = (q - x[1]).SquaredLength();
      snapId = (dists[0] <= dists[1] ? 0 : 1);
      if (dists[snapId] < tol2) {
        q = x[snapId], t2 = 0.;
        edge2->SetId(snapId == 0 ? 1 : 0, edge2->GetId(snapId));
      }
    }
    if (snapId == -1) {
      vtkLine::DistanceToLine(q, x[0], x[1], t2, closestPt);
      t2 = clamp(t2, 0., 1.);
      q  = closestPt;
    }

    // Insert (new) edge intersection points
    //
    // This also updates the edge intersections of neighboring cells,
    // including in particular those that are not themselves intersected
    // by any of the intersection lines. Edge intersection points are
    // inserted only once into the points list, thus ensuring that
    // intersection line segments share a common point.
    if (npts != 3) {
      Throw(ERR_LogicError, __FUNCTION__, "Surface must have triangular faces");
    }
    if (edge1->GetId(0) == pts[0] && edge1->GetId(1) == pts[1]) {
      InsertEdgeIntersection(intersections, surface, cellId, 0, p, t1);
    } else if (edge2->GetId(0) == pts[0] && edge2->GetId(1) == pts[1]) {
      InsertEdgeIntersection(intersections, surface, cellId, 0, q, t2);
    }
    if (edge1->GetId(0) == pts[1] && edge1->GetId(1) == pts[2]) {
      InsertEdgeIntersection(intersections, surface, cellId, 1, p, t1);
    } else if (edge2->GetId(0) == pts[1] && edge2->GetId(1) == pts[2]) {
      InsertEdgeIntersection(intersections, surface, cellId, 1, q, t2);
    }
    if (edge1->GetId(0) == pts[2] && edge1->GetId(1) == pts[0]) {
      InsertEdgeIntersection(intersections, surface, cellId, 2, p, t1);
    } else if (edge2->GetId(0) == pts[2] && edge2->GetId(1) == pts[0]) {
      InsertEdgeIntersection(intersections, surface, cellId, 2, q, t2);
    }
    if (debug && (t1 != 0. || t2 != 0.) && intersections.find(cellId) == intersections.end()) {
      Throw(ERR_LogicError, __FUNCTION__, "No edge intersection recorded for intersected cell ", cellId);
    }

    // Update line end points (**after** InsertEdgeIntersection)
    cut->GetPoints()->SetPoint(ptIds->GetId(0), p);
    cut->GetPoints()->SetPoint(ptIds->GetId(1), q);
  }
  cut->RemoveDeletedCells();

  // Split surface cells at edge intersection points
  for (const auto &intersection : intersections) {
    cellId           = intersection.first;
    const auto &edge = intersection.second.edge;
    bisect[0] = (edge[0].i >= 0 ? 1 : 0);
    bisect[1] = (edge[1].i >= 0 ? 1 : 0);
    bisect[2] = (edge[2].i >= 0 ? 1 : 0);
    surface->GetCellPoints(cellId, npts, pts);
    switch (bisect[0] + bisect[1] + bisect[2]) {
      case 0: {
        Throw(ERR_LogicError, __FUNCTION__, "Exptected at least one edge intersection");
      } break;
      case 1: {
        if      (bisect[0]) Bisect(surface, cellId, pts[0], pts[1], pts[2], edge[0], split);
        else if (bisect[1]) Bisect(surface, cellId, pts[1], pts[2], pts[0], edge[1], split);
        else                Bisect(surface, cellId, pts[2], pts[0], pts[1], edge[2], split);
      } break;
      case 2: {
        if      (bisect[0] && bisect[1]) Trisect(surface, cellId, pts[0], pts[1], pts[2], edge[0], edge[1], split);
        else if (bisect[1] && bisect[2]) Trisect(surface, cellId, pts[1], pts[2], pts[0], edge[1], edge[2], split);
        else                             Trisect(surface, cellId, pts[2], pts[0], pts[1], edge[2], edge[0], split);
      } break;
      case 3: {
        Throw(ERR_LogicError, __FUNCTION__, "Quadsection not possible when intersecting triangle with one line");
      } break;
    }
  }
  surface->RemoveDeletedCells();

  if (debug) {
    static int callId = 0; ++callId;
    char fname[64];
    snprintf(fname, 64, "debug_split_surface_other_%d.vtp", callId);
    WritePolyData(fname, surface);
    snprintf(fname, 64, "debug_split_surface_cells_%d.vtp", callId);
    WritePolyData(fname, split);
    snprintf(fname, 64, "debug_split_surface_cut_%d.vtp", callId);
    WritePolyData(fname, cut);
  }

  // Create polygonal mesh tesselation of closed intersection polygon
  vtkSmartPointer<vtkPolyData> divider = Divider(cut);
  if (debug) {
    static int callId = 0; ++callId;
    char fname[64];
    snprintf(fname, 64, "debug_split_surface_polygon_%d.vtp", callId);
    WritePolyData(fname, divider);
  }
  divider = TesselateDivider(divider, mean_edge_length);

  // Prepare point/cell data before appending polygonal data sets
  vtkCellData * const surfaceCD = surface->GetCellData();
  vtkCellData * const dividerCD = divider->GetCellData();
  vtkCellData * const splitCD   = split->GetCellData();

  arr = surfaceCD->GetArray(SOURCE_ARRAY_NAME);
  int divider_source_index = min(0, static_cast<int>(arr->GetRange(0)[0])) - 1;

  arr = splitCD->GetArray(SOURCE_ARRAY_NAME);
  if (!arr) {
    arr = NewVtkDataArray(VTK_SHORT, split->GetNumberOfCells(), 1, SOURCE_ARRAY_NAME);
    splitCD->AddArray(arr);
  }
  arr->FillComponent(0, 0.);

  arr = dividerCD->GetArray(SOURCE_ARRAY_NAME);
  if (arr == nullptr) {
    arr = NewVtkDataArray(VTK_SHORT, divider->GetNumberOfCells(), 1, SOURCE_ARRAY_NAME);
    dividerCD->AddArray(arr);
  }
  arr->FillComponent(0, static_cast<double>(divider_source_index));

  divider = CalculateNormals(divider, surface->GetPointData()->GetNormals() != nullptr,
                                      surface->GetCellData ()->GetNormals() != nullptr);

  if (debug) {
    static int callId = 0; ++callId;
    char fname[64];
    snprintf(fname, 64, "debug_split_surface_divider_%d.vtp", callId);
    WritePolyData(fname, divider);
  }

  // Merge surface with intersected cells and tesselated divider polygon
  vtkNew<vtkAppendPolyData> appender;
  AddVTKInput(appender, surface);
  AddVTKInput(appender, split);
  // divider added as input below

  vtkNew<vtkCleanPolyData> merger;
  SetVTKConnection(merger, appender);
  merger->ConvertStripsToPolysOff();
  merger->ConvertPolysToLinesOff();
  merger->ConvertLinesToPointsOff();
  merger->PointMergingOn();
  merger->ToleranceIsAbsoluteOn();
  merger->SetAbsoluteTolerance(1e-12);

  if (debug) {
    static int callId = 0; ++callId;
    char fname[64];
    snprintf(fname, 64, "debug_split_surface_interim_%d.vtp", callId);
    merger->Update();
    WritePolyData(fname, merger->GetOutput());
  }

  AddVTKInput(appender, divider);
  merger->Update();
  vtkSmartPointer<vtkPolyData> output = merger->GetOutput();

  // Remove those triangles from resulting mesh originating from the tesselated
  // divider whose three vertices are all on the divider boundary and that have
  // a corresponding cell in the split surface mesh, i.e., a duplicate triangle
  output->BuildLinks();
  vtkNew<vtkIdList> ptIds1, ptIds2;
  for (vtkIdType cellId = 0; cellId < output->GetNumberOfCells(); ++cellId) {
    if (output->GetCellType(cellId) != VTK_EMPTY_CELL) {
      output->GetCellPoints(cellId, ptIds1.GetPointer());
      for (vtkIdType i = 0; i < ptIds1->GetNumberOfIds(); ++i) {
        output->GetPointCells(ptIds1->GetId(i), cellIds.GetPointer());
        for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
          if (cellIds->GetId(j) > cellId && output->GetCellType(cellIds->GetId(j)) != VTK_EMPTY_CELL) {
            output->GetCellPoints(cellIds->GetId(j), ptIds2.GetPointer());
            if (ptIds1->GetNumberOfIds() == ptIds2->GetNumberOfIds()) {
              ptIds2->IntersectWith(ptIds1.GetPointer());
              if (ptIds1->GetNumberOfIds() == ptIds2->GetNumberOfIds()) {
                output->DeleteCell(cellIds->GetId(j));
              }
            }
          }
        }
      }
    }
  }
  output->RemoveDeletedCells();

  return output;
}

// -----------------------------------------------------------------------------
void LabelSurfacePatches(vtkPolyData *surface)
{
  SurfacePatches cc;
  cc.Input(surface);
  cc.Ordering(SurfacePatches::LargestFirst);
  cc.Run();

  vtkDataArray * const patch_labels  = cc.GetLabelsArray();
  vtkDataArray * const source_labels = surface->GetCellData()->GetArray(SOURCE_ARRAY_NAME);

  UnorderedSet<double>      used;
  UnorderedMap<double, int> bins;
  decltype(bins)::iterator  bin;
  double                    patch_label, label;
  int                       count;

  double min_label = source_labels->GetRange(0)[0];
  double max_label = source_labels->GetRange(0)[1];

  for (int i = 0; i < cc.NumberOfPatches(); ++i) {
    bins.clear();
    patch_label = static_cast<double>(i + 1);
    for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
      if (patch_labels->GetComponent(cellId, 0) == patch_label) {
        label = source_labels->GetComponent(cellId, 0);
        if (label != 0.) {
          bin = bins.find(label);
          if (bin == bins.end()) bins[label]  = 1;
          else                   bin->second += 1;
        }
      }
    }
    label = 0.;
    count = 0;
    for (bin = bins.begin(); bin != bins.end(); ++bin) {
      if (bin->second > count) {
        label = bin->first;
        count = bin->second;
      }
    }
    if (label == 0. || used.find(label) != used.end()) {
      if (label < 0.) {
        min_label -= 1.;
        label = min_label;
      } else {
        max_label += 1.;
        label = max_label;
      }
    }
    for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
      if (patch_labels->GetComponent(cellId, 0) == patch_label) {
        source_labels->SetComponent(cellId, 0, label);
      }
    }
    used.insert(label);
  }
}

// -----------------------------------------------------------------------------
void GrowSourceRegion(vtkPolyData *surface, bool ignore_border_edges)
{
  vtkDataArray * const source = surface->GetCellData()->GetArray(SOURCE_ARRAY_NAME);

  Queue<vtkIdType>          active;
  UnorderedMap<double, int> bins;
  decltype(bins)::iterator  bin;
  double                    label;
  int                       count;
  vtkIdType                 cellId, nbrId, npts, *pts;
  vtkNew<vtkIdList>         cellIds;

  for (cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
    if (source->GetComponent(cellId, 0) == 0.) {
      surface->GetCellPoints(cellId, npts, pts);
      for (vtkIdType i = 0; i < npts; ++i) {
        surface->GetCellEdgeNeighbors(cellId, pts[i], pts[(i+1)%npts], cellIds.GetPointer());
        if (ignore_border_edges || cellIds->GetNumberOfIds() == 1) {
          for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
            nbrId = cellIds->GetId(j);
            label = source->GetComponent(nbrId, 0);
            if (label != 0.) {
              active.push(cellId);
              i = npts; // break also cell points loop
              break;
            }
          }
        }
      }
    }
  }
  while (!active.empty()) {
    cellId = active.front();
    active.pop();
    if (source->GetComponent(cellId, 0) == 0.) {
      bins.clear();
      surface->GetCellPoints(cellId, npts, pts);
      for (vtkIdType i = 0; i < npts; ++i) {
        surface->GetCellEdgeNeighbors(cellId, pts[i], pts[(i+1)%npts], cellIds.GetPointer());
        if (ignore_border_edges || cellIds->GetNumberOfIds() == 1) {
          for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
            nbrId = cellIds->GetId(j);
            label = source->GetComponent(nbrId, 0);
            if (label == 0.) {
              active.push(nbrId);
            } else {
              bin = bins.find(label);
              if (bin == bins.end()) bins[label]  = 1;
              else                   bin->second += 1;
            }
          }
        }
      }
      label = 0.;
      count = 0;
      for (bin = bins.begin(); bin != bins.end(); ++bin) {
        if (bin->second > count) {
          label = bin->first;
          count = bin->second;
        }
      }
      source->SetComponent(cellId, 0, label);
    }
  }
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  REQUIRES_POSARGS(0);

  Array<const char *> input_names;
  const char *output_name = nullptr;
  const char *labels_name = nullptr;

  double labels_percentage    = 33.;
  double tolerance            = NaN;
  double snap_tolerance       = 0.;
  int    smooth_boundaries    = -1;
  bool   join_boundaries      = true;
  double min_edge_length      = NaN;
  double max_edge_length      = NaN;
  double edge_length_sd       = 2.;  // SD factor of edge length range about mean
  int    max_remesh_steps     = 3;   // max no. of remeshing steps
  int    max_smooth_steps     = 0;   // (max) no. of smoothing steps
  int    remesh_radius        = 3;   // intersection cell neighborhood connectivity radius
  int    smooth_radius        = 2;   // intersection cell neighborhood connectivity radius
  double smooth_lambda        = NaN;
  double smooth_mu            = NaN;
  bool   fill_source_label    = true;
  int    max_smooth_source    = 0;
  bool   output_largest_comp  = false;
  bool   add_dividers         = false;
  bool   output_point_normals = true;
  bool   output_cell_normals  = true;

  const char *point_source_name = nullptr;
  const char *cell_source_name  = nullptr;

  if (NUM_POSARGS == 1) {
    output_name = POSARG(1);
  } else if (NUM_POSARGS > 1) {
    input_names.reserve(NUM_POSARGS-1);
    for (int i = 1; i < NUM_POSARGS; ++i) {
      input_names.push_back(POSARG(i));
    }
    output_name = POSARG(NUM_POSARGS);
  }

  for (ALL_OPTIONS) {
    if (OPTION("-labels")) {
      labels_name = ARGUMENT;
    }
    else if (OPTION("-labels-percentage")) {
      PARSE_ARGUMENT(labels_percentage);
    }
    else if (OPTION("-i") || OPTION("-input")) {
      do {
        input_names.push_back(ARGUMENT);
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-o") || OPTION("-output")) {
      if (output_name != nullptr) {
        input_names.insert(input_names.begin() + NUM_POSARGS - 1, output_name);
      }
      output_name = ARGUMENT;
    }
    else if (OPTION("-source-array")) {
      point_source_name = cell_source_name = ARGUMENT;
    }
    else if (OPTION("-point-source-array")) {
      point_source_name = ARGUMENT;
    }
    else if (OPTION("-cell-source-array")) {
      cell_source_name = ARGUMENT;
    }
    else if (OPTION("-tolerance") || OPTION("-tol")) {
      PARSE_ARGUMENT(tolerance);
    }
    else if (OPTION("-snap-tolerance") || OPTION("-snap-tol")) {
      PARSE_ARGUMENT(snap_tolerance);
    }
    else if (OPTION("-smooth-boundaries") || OPTION("-boundary-smoothing")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(smooth_boundaries);
      else smooth_boundaries = 3;
    }
    else if (OPTION("-nosmooth-boundaries") || OPTION("-noboundary-smoothing")) {
      smooth_boundaries = 0;
    }
    else if (OPTION("-edge-length")) {
      PARSE_ARGUMENT(min_edge_length);
      if (HAS_ARGUMENT) PARSE_ARGUMENT(max_edge_length);
      else max_edge_length = min_edge_length;
      if (min_edge_length < max_edge_length) {
        swap(min_edge_length, max_edge_length);
      }
    }
    else if (OPTION("-min-edge-length")) {
      PARSE_ARGUMENT(min_edge_length);
    }
    else if (OPTION("-max-edge-length")) {
      PARSE_ARGUMENT(max_edge_length);
    }
    else if (OPTION("-edge-length-sigma") || OPTION("-edge-length-sd")) {
      PARSE_ARGUMENT(edge_length_sd);
    }
    else if (OPTION("-remesh") || OPTION("-remeshing")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(max_remesh_steps);
      else max_remesh_steps = 3;
    }
    else if (OPTION("-remeshing-iterations") || OPTION("-remesh-iterations")) {
      PARSE_ARGUMENT(max_remesh_steps);
    }
    else if (OPTION("-noremeshing") || OPTION("-noremesh")) {
      max_remesh_steps = 0;
    }
    else if (OPTION("-smooth") || OPTION("-smoothing")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(max_smooth_steps);
      else max_smooth_steps = 100;
    }
    else if (OPTION("-smoothing-iterations") || OPTION("-smooth-iterations")) {
      PARSE_ARGUMENT(max_smooth_steps);
    }
    else if (OPTION("-smoothing-lambda") || OPTION("-smooth-lambda")) {
      PARSE_ARGUMENT(smooth_lambda);
    }
    else if (OPTION("-smoothing-mu") || OPTION("-smooth-mu")) {
      PARSE_ARGUMENT(smooth_mu);
    }
    else if (OPTION("-nosmoothing") || OPTION("-nosmooth")) {
      max_smooth_steps = 0;
    }
    else if (OPTION("-neighborhood") || OPTION("-neighbourhood") || OPTION("-radius")) {
      PARSE_ARGUMENT(remesh_radius);
      smooth_radius = remesh_radius;
    }
    else if (OPTION("-remeshing-neighborhood") || OPTION("-remeshing-neighbourhood") || OPTION("-remeshing-radius") ||
             OPTION("-remesh-neighborhood")    || OPTION("-remesh-neighbourhood")    || OPTION("-remesh-radius")) {
      PARSE_ARGUMENT(remesh_radius);
    }
    else if (OPTION("-smoothing-neighborhood") || OPTION("-smoothing-neighbourhood") || OPTION("-smoothing-radius") ||
             OPTION("-smooth-neighborhood")    || OPTION("-smooth-neighbourhood")    || OPTION("-smooth-radius")) {
      PARSE_ARGUMENT(smooth_radius);
    }
    else if (OPTION("-smooth-source") || OPTION("-source-smoothing")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(max_smooth_source);
      else max_smooth_source = 10;
    }
    else if (OPTION("-nosmooth-source") || OPTION("-nosource-smoothing")) {
      max_smooth_source = 0;
    }
    else if (OPTION("-normals")) {
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(output_point_normals);
        output_cell_normals = output_point_normals;
      } else {
        output_cell_normals = output_point_normals = true;
      }
    }
    else if (OPTION("-nonormals")) {
      output_cell_normals = output_point_normals = false;
    }
    else HANDLE_BOOLEAN_OPTION("point-normals", output_point_normals);
    else HANDLE_BOOLEAN_OPTION("cell-normals", output_point_normals);
    else HANDLE_BOOLEAN_OPTION("join", join_boundaries);
    else HANDLE_BOOLEAN_OPTION("largest", output_largest_comp);
    else HANDLE_BOOLEAN_OPTION("dividers", add_dividers);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (input_names.size() < 2) {
    FatalError("At least two -input surfaces required!");
  }
  if (output_name == nullptr) {
    FatalError("No -output surface file name specified!");
  }

  if (join_boundaries) {
    if (smooth_boundaries < 0) smooth_boundaries = 3;
  } else {
    smooth_boundaries = max(smooth_boundaries, 0);
    remesh_radius    = smooth_radius    = 0;
    max_remesh_steps = max_smooth_steps = 0;
  }
  if (!IsNaN(min_edge_length) && min_edge_length < 0. &&
      !IsNaN(max_edge_length) && IsInf(max_edge_length)) {
    remesh_radius    = 0;
    max_remesh_steps = 0;
  }
  if (IsNaN(smooth_lambda)) {
    if (IsNaN(smooth_mu)) smooth_mu = -.75;
    if (-1. < smooth_mu && smooth_mu < 0.) {
      smooth_lambda = abs(smooth_mu) - .01;
    } else {
      smooth_lambda = smooth_mu;
    }
  }
  if (cell_source_name == nullptr && point_source_name == nullptr) {
    fill_source_label = false;
    max_smooth_source = 0;
  }

  // Read segmentation labels
  GreyImage labels_image;
  if (labels_name) {
    InitializeIOLibrary();
    labels_image.Read(labels_name);
  } else {
    FatalError("Can only merge surfaces with corresponding -labels image at the moment!");
  }

  // Read input surfaces
  Array<vtkSmartPointer<vtkPolyData> > surfaces;
  surfaces.resize(input_names.size());
  for (size_t i = 0; i < input_names.size(); ++i) {
    surfaces[i] = ReadPolyData(input_names[i]);
  }

  // Default tolerance
  if (IsNaN(tolerance)) {
    tolerance = max(max(labels_image.XSize(),
                        labels_image.YSize()),
                        labels_image.ZSize());
  }

  // Add cell source arrays
  size_t source_offset = 1;
  vtkSmartPointer<vtkDataArray> cell_source;
  if (point_source_name) {
    cell_source = surfaces[0]->GetPointData()->GetArray(point_source_name);
    if (cell_source && cell_source->GetNumberOfComponents() == 1) {
      source_offset = max(source_offset, static_cast<size_t>(ceil(cell_source->GetRange(0)[1])));
    }
    cell_source = nullptr; // not a cell data array!
  }
  if (cell_source_name) {
    cell_source = surfaces[0]->GetCellData()->GetArray(cell_source_name);
    if (cell_source && cell_source->GetNumberOfComponents() == 1) {
      source_offset = max(source_offset, static_cast<size_t>(ceil(cell_source->GetRange(0)[1])));
    }
  }
  if (cell_source) {
    cell_source->SetName(SOURCE_ARRAY_NAME);
    cell_source = nullptr;
  } else {
    source_offset = 1;
    AddCellSourceArray(surfaces[0], source_offset);
  }
  for (size_t i = 1; i < surfaces.size(); ++i) {
    AddCellSourceArray(surfaces[i], source_offset + i);
  }

  // Merge surfaces
  if (verbose > 0) {
    cout << "Merging surfaces at segmentation boundaries...";
    if (verbose > 1) cout << "\n";
    cout.flush();
  }
  vtkSmartPointer<vtkPolyData> output = surfaces[0];
  Array<vtkSmartPointer<vtkPolyData>> boundaries;
  {
    vtkSmartPointer<vtkPolyData> boundary;
    UnorderedSet<int> output_labels = InsideLabels(labels_image, labels_percentage, output);
    if (output_labels.empty()) {
      FatalError("Could not determine labels corresponding to the inside of the first input surface!");
    }
    for (size_t i = 1; i < surfaces.size(); ++i) {
      UnorderedSet<int> surface_labels = InsideLabels(labels_image, labels_percentage, surfaces[i]);
      if (surface_labels.empty()) {
        FatalError("Could not determine labels corresponding to the inside of the "
                   << (i == 1 ? "1st" : (i == 2 ? "2nd" : (i == 3 ? "3rd" : ((ToString(i) + "th").c_str()))))
                   << " input surface!");
      }
      boundary = LabelBoundary(labels_image, output_labels, surface_labels);
      output = Merge(output, surfaces[i], boundary, tolerance, smooth_boundaries, join_boundaries);
      output_labels.insert(surface_labels.begin(), surface_labels.end());
      boundaries.push_back(boundary);
    }
    surfaces.clear();
  }
  if (verbose > 0) {
    if (verbose > 1) cout << "Merging surfaces at segmentation boundaries...";
    cout << " done" << endl;
  }

  if (join_boundaries) {

    // Extract largest connected component
    if (output_largest_comp) {
      if (verbose > 0) {
        cout << "Extracting largest surface component...";
        cout.flush();
      }
      vtkNew<vtkPolyDataConnectivityFilter> lcc;
      SetVTKInput(lcc, output);
      lcc->SetExtractionModeToLargestRegion();
      lcc->Update();
      output = lcc->GetOutput();
      if (verbose > 0) cout << " done" << endl;
    }

    // Recalculate surface normals and fix vertex order if necessary
    if (verbose > 0) {
      cout << "Calculating surface normals...";
      cout.flush();
    }
    const bool calc_cell_normals = true;
    output = CalculateNormals(output, output_point_normals, calc_cell_normals);
    if (verbose > 0) cout << " done" << endl;

    if (debug > 0) {
      WritePolyData("debug_joined_boundaries.vtp", output);
    }

    // Remesh newly inserted cells which are joining the intersection boundaries
    if (remesh_radius > 0 && max_remesh_steps > 0) {
      if (verbose > 0) {
        cout << "Remeshing surface near intersection boundaries...";
        cout.flush();
      }

      output = Triangulate(output);
      output->BuildLinks();

      vtkSmartPointer<vtkDataArray> mask;
      mask = IntersectionCellMask(output, remesh_radius);

      double min_global_length, max_global_length;
      GetMinMaxEdgeLength(output, min_global_length, max_global_length);
      GetEdgeLengthRange(output, mask, min_edge_length, max_edge_length, edge_length_sd);

      min_global_length -= 1e-12;
      max_global_length += 1e-12;

      const vtkIdType ncells = output->GetNumberOfCells();
      vtkSmartPointer<vtkDataArray> min_edge_length_arr, max_edge_length_arr;
      min_edge_length_arr = NewVtkDataArray(VTK_FLOAT, ncells, 1, MIN_EDGE_LENGTH_ARRAY_NAME);
      max_edge_length_arr = NewVtkDataArray(VTK_FLOAT, ncells, 1, MAX_EDGE_LENGTH_ARRAY_NAME);

      int n = 0;
      for (vtkIdType cellId = 0; cellId < ncells; ++cellId) {
        if (mask->GetComponent(cellId, 0) != 0.) {
          min_edge_length_arr->SetComponent(cellId, 0, min_edge_length);
          max_edge_length_arr->SetComponent(cellId, 0, max_edge_length);
          ++n;
        } else {
          min_edge_length_arr->SetComponent(cellId, 0, min_global_length);
          max_edge_length_arr->SetComponent(cellId, 0, max_global_length);
        }
      }
      mask = nullptr;

      if (n > 0) {
        if (verbose > 1) {
          cout << "\n";
          cout << "  No. of boundary cells = " << n << "\n";
          cout << "  Edge length range     = [" << min_edge_length << ", " << max_edge_length << "]\n";
          if (debug_time > 0) cout << "\n";
          cout.flush();
        }
        int nmelt = 0, ninvs = 0, nsubd = 0;
        for (int iter = 0; iter < max_remesh_steps; ++iter) {
          SurfaceRemeshing remesher;
          remesher.MeltTrianglesOff();
          remesher.MeltNodesOff();
          remesher.InvertTrianglesToIncreaseMinHeightOff();
          remesher.Input(output);
          if (iter == 0) {
            remesher.MinCellEdgeLengthArray(min_edge_length_arr);
            remesher.MaxCellEdgeLengthArray(max_edge_length_arr);
            min_edge_length_arr = nullptr;
            max_edge_length_arr = nullptr;
          }
          remesher.Run();
          if (remesher.NumberOfChanges() == 0) break;
          output = remesher.Output();
          if (debug > 0) {
            char fname[64];
            snprintf(fname, 64, "debug_remeshed_%d.vtp", iter+1);
            WritePolyData(fname, output);
          }
          nmelt += remesher.NumberOfMeltedEdges();
          ninvs += remesher.NumberOfInversions();
          nsubd += remesher.NumberOfBisections();
          nsubd += 2 * remesher.NumberOfTrisections();
          nsubd += 3 * remesher.NumberOfQuadsections();
        }
        if (verbose > 1) {
          if (debug_time > 0) cout << "\n";
          cout << "  No. of edge-meltings  = " << nmelt << "\n";
          cout << "  No. of inversions     = " << ninvs << "\n";
          cout << "  No. of subdivisions   = " << nsubd << "\n";
          cout.flush();
        }
        output->GetPointData()->RemoveArray(SurfaceRemeshing::MIN_EDGE_LENGTH);
        output->GetPointData()->RemoveArray(SurfaceRemeshing::MAX_EDGE_LENGTH);
      }

      if (verbose > 0) {
        if (verbose > 1) {
          cout << "Remeshing surface near intersection boundaries...";
        }
        cout << " done" << endl;
      }
    }

    // Smooth surface in vicinity of joined intersection boundaries to remove
    // small self-intersections introduced by joining and/or remeshing steps
    if (smooth_radius > 0 && max_smooth_steps > 0) {
      if (verbose > 0) {
        cout << "Smoothing surface near intersection boundaries...";
        cout.flush();
      }
      output->BuildLinks();
      MeshSmoothing smoother;
      smoother.Input(output);
      smoother.Mask(IntersectionPointMask(output, smooth_radius));
      smoother.Lambda(smooth_lambda);
      smoother.Mu(smooth_mu);
      smoother.NumberOfIterations(max_smooth_steps);
      if (IsNaN(smooth_mu) || fequal(smooth_lambda, smooth_mu)) {
        smoother.Weighting(MeshSmoothing::Gaussian);
        smoother.AdjacentValuesOnlyOff();
      } else {
        smoother.Weighting(MeshSmoothing::Combinatorial);
        smoother.AdjacentValuesOnlyOn();
      }
      smoother.SmoothPointsOn();
      smoother.Run();
      output->SetPoints(smoother.Output()->GetPoints());
      if (verbose > 0) cout << " done" << endl;
      if (debug > 0) {
        const int i = output->GetPointData()->AddArray(smoother.Mask());
        WritePolyData("debug_smoothed_transition.vtp", output);
        output->GetPointData()->RemoveArray(i);
      }
    }

    // Insert cutting planes in reverse order
    if (add_dividers) {
      if (verbose > 0) {
        cout << "Adding internal division planes...";
        if (verbose > 1) cout << "\n";
        cout.flush();
      }
      const EdgeTable edgeTable(output);
      const double ds = AverageEdgeLength(output->GetPoints(), edgeTable);
      vtkSmartPointer<vtkPolyData> plane, polygon, cut;
      for (int i = static_cast<int>(boundaries.size()-1); i >= 0; --i) {
        string msg;
        if (verbose > 1) {
          const int n = i + 1;
          msg = "  Adding division plane for ";
          if      (n == 1) msg += "1st";
          else if (n == 2) msg += "2nd";
          else if (n == 3) msg += "3rd";
          else             msg += ToString(n) + "th";
          msg += " boundary";
          cout << msg << "..." << endl;
        }
        if (FindCuttingPlane(output, boundaries[i], plane, cut, 2. * tolerance, ds, snap_tolerance)) {
          if (debug > 0) {
            char fname[64];
            snprintf(fname, 64, "debug_cutting_plane_%d.vtp", static_cast<int>(boundaries.size()) - i);
            WritePolyData(fname, plane);
          }
          if (debug > 0) {
            char fname[64];
            snprintf(fname, 64, "debug_cutting_polygon_%d.vtp", static_cast<int>(boundaries.size()) - i);
            WritePolyData(fname, cut);
          }
          output = AddClosedIntersectionDivider(output, cut, snap_tolerance);
          if (debug > 0) {
            char fname[64];
            snprintf(fname, 64, "debug_output+divider_%d.vtp", static_cast<int>(boundaries.size()) - i);
            WritePolyData(fname, output);
          }
          if (!msg.empty()) cout << msg << "... done" << endl;
        } else {
          if (verbose > 0) {
            if (!msg.empty()) cout << msg << "...";
            cout << " failed\n" << endl;
          }
          FatalError("Could not find a closed intersection with finite cutting plane near segmentation boundary!");
        }
      }
      if (verbose > 0) {
        if (verbose > 1) cout << "Adding internal division planes...";
        cout << " done" << endl;
      }
    }

  } // if (join_boundaries)

  // Clean up and rename or remove data source arrays
  output->BuildLinks();

  if (fill_source_label) {
    if (verbose > 0) {
      cout << "Filling cell source labels...";
      cout.flush();
    }
    if (add_dividers) {
      LabelSurfacePatches(output);
    } else {
      GrowSourceRegion(output, false);
      GrowSourceRegion(output, true); // just in case... usually NOP
    }
    if (verbose > 0) cout << " done" << endl;
  }

  if (max_smooth_source > 0) {
    if (verbose > 0) {
      cout << "Smoothing cell source labels...";
      cout.flush();
    }

    UnorderedMap<double, int> bins;
    decltype(bins)::iterator  bin;
    double                    label;
    int                       count;
    unsigned short            ncells;
    vtkIdType                 npts, *pts, *cells;

    vtkSmartPointer<vtkDataArray> cell_source;
    cell_source = output->GetCellData()->GetArray(SOURCE_ARRAY_NAME);

    for (int iter = 0; iter < max_smooth_source; ++iter) {
      bool modified = false;
      for (vtkIdType cellId = 0; cellId < output->GetNumberOfCells(); ++cellId) {
        bins.clear();
        output->GetCellPoints(cellId, npts, pts);
        for (vtkIdType i = 0; i < npts; ++i) {
          output->GetPointCells(pts[i], ncells, cells);
          for (unsigned short j = 0; j < ncells; ++j) {
            if (cells[j] != cellId) {
              label = cell_source->GetComponent(cells[j], 0);
              bin   = bins.find(label);
              if (bin == bins.end()) bins[label]  = 1;
              else                   bin->second += 1;
            }
          }
        }
        label = cell_source->GetComponent(cellId, 0);
        bin   = bins.find(label);
        count = (bin == bins.end() ? 0 : bin->second);
        for (bin = bins.begin(); bin != bins.end(); ++bin) {
          if (bin->second > count) {
            label    = bin->first;
            count    = bin->second;
            modified = true;
          }
        }
        cell_source->SetComponent(cellId, 0, label);
      }
      if (!modified) break;
    }
    if (verbose > 0) cout << " done" << endl;
  }

  if ((output_point_normals && add_dividers) || point_source_name) {
    if (verbose > 0) {
      cout << "Generating point source labels...";
      cout.flush();
    }

    double         label;
    unsigned short ncells;
    vtkIdType      *cells;
    vtkSmartPointer<vtkDataArray> cell_source, point_source;

    cell_source  = output->GetCellData()->GetArray(SOURCE_ARRAY_NAME);
    point_source = NewVtkDataArray(cell_source->GetDataType(),
                                   output->GetNumberOfPoints(), 1,
                                   SOURCE_ARRAY_NAME); // renamed later
    output->GetPointData()->AddArray(point_source);

    for (vtkIdType ptId = 0; ptId < output->GetNumberOfPoints(); ++ptId) {
      output->GetPointCells(ptId, ncells, cells);
      if (ncells == 0) {
        label = 0.;
      } else {
        label = cell_source->GetComponent(cells[0], 0);
        for (unsigned short i = 1; i < ncells; ++i) {
          if (cell_source->GetComponent(cells[i], 0) != label) {
            label = 0.;
            break;
          }
        }
      }
      point_source->SetComponent(ptId, 0, label);
    }

    if (verbose > 0) cout << " done" << endl;
  }

  // Only now rename/remove not requested data arrays
  if (output_point_normals) {
    if (add_dividers) FixBorderPointNormals(output); // **before** renaming the source arrays
  } else {
    output->GetPointData()->SetNormals(nullptr);
  }
  if (!output_cell_normals) {
    output->GetCellData()->SetNormals(nullptr);
  }
  if (cell_source_name) {
    output->GetCellData()->GetArray(SOURCE_ARRAY_NAME)->SetName(cell_source_name);
  } else {
    output->GetCellData()->RemoveArray(SOURCE_ARRAY_NAME);
  }
  if (point_source_name) {
    output->GetPointData()->GetArray(SOURCE_ARRAY_NAME)->SetName(point_source_name);
  }

  // Write output surface mesh
  if (!WritePolyData(output_name, output)) {
    FatalError("Failed to write output surface to " << output_name);
  }

  return 0;
}

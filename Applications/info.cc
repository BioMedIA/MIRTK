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

#include <mirtkBaseImage.h>
#include <mirtkImageReader.h>

#ifdef HAVE_MIRTK_Transformation
#  include <mirtkTransformation.h>
#  include <mirtkAffineTransformation.h>
#endif // HAVE_MIRTK_Transformation

#ifdef HAVE_MIRTK_PointSet
#  include <mirtkEdgeTable.h>
#  include <mirtkTriangle.h>
#  include <mirtkPolyhedron.h>
#  include <mirtkPointSetUtils.h>
#  include <mirtkDataStatistics.h>

#  include <vtkSmartPointer.h>
#  include <vtkPolyData.h>
#  include <vtkDataArray.h>
#  include <vtkCellArray.h>
#  include <vtkPointData.h>
#  include <vtkCellData.h>
#  include <vtkGenericCell.h>
#  include <vtkOctreePointLocator.h>
#  include <vtkMath.h>
#  include <vtkUnsignedCharArray.h>
#  include <vtkFloatArray.h>
#  include <vtkDataSetSurfaceFilter.h>
#  include <vtkCleanPolyData.h>
#  include <vtkPolyDataConnectivityFilter.h>
#  include <vtkIntersectionPolyDataFilter.h>
#endif // HAVE_MIRTK_PointSet

using namespace mirtk;



// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <file> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Prints some useful information about the given input file, which\n";
  cout << "  can be an image, transformation (requires MIRTK Transformation module),\n";
  cout << "  or point set (requires MIRTK Point Set module).";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input   Input image";
  #ifdef HAVE_MIRTK_Transformation
    cout << " or transformation";
  #endif
  #ifdef HAVE_MIRTK_PointSet
    cout << " or point set";
  #endif
  cout << " file.\n";
  cout << "\n";
  cout << "Image options:\n";
  // Note: No two first named options can be identical, otherwise Sphinx
  //       complains that the -attributes point set option below is a duplicate.
  //       (cf. help-rst.py output).
  cout << "  -a -attributes   Print attributes of image.\n";
  #ifdef HAVE_MIRTK_Transformation
    cout << "\n";
    cout << "Transformation options:\n";
    cout << "  -attributes   Print attributes of transformation. (default)\n";
    cout << "  -type-name    Print name of transformation type.\n";
    cout << "  -type-id      Print enum of transformation type.\n";
  #endif // HAVE_MIRTK_Transformation
  #ifdef HAVE_MIRTK_PointSet
  cout << "\n";
  cout << "Point set options:\n";
  cout << "  -edgelength              Report edge length statistics. (default: off)\n";
  cout << "  -self-intersections      Check for self-intersections. (default: off)\n";
  cout << "  -output-surface <file>   Write surface mesh to specified file. (default: none)\n";
  #endif // HAVE_MIRTK_PointSet
  PrintStandardOptions(cout);
  cout << "\n";
  cout.flush();
}

// =============================================================================
// Point set auxiliaries
// =============================================================================
#ifdef HAVE_MIRTK_PointSet

// -----------------------------------------------------------------------------
/// Convert input dataset to polygonal data, i.e., boundary surface
vtkSmartPointer<vtkPolyData> ConvertToPolyData(vtkDataSet *dataset)
{
  vtkSmartPointer<vtkPolyData> polydata = vtkPolyData::SafeDownCast(dataset);
  if (polydata) return polydata;
  vtkSmartPointer<vtkDataSetSurfaceFilter> surface;
  surface = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  SetVTKInput(surface, dataset);
  surface->PassThroughCellIdsOff();
  surface->PassThroughPointIdsOff();
  surface->Update();
  return surface->GetOutput();
}

// -----------------------------------------------------------------------------
/// Clean polydata
vtkSmartPointer<vtkPolyData> CleanPolyData(vtkSmartPointer<vtkPolyData> polydata)
{
  vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->SetAbsoluteTolerance(.0);
  cleaner->PointMergingOn();
  cleaner->ConvertLinesToPointsOn();
  cleaner->ConvertPolysToLinesOn();
  cleaner->ConvertStripsToPolysOn();
  SetVTKInput(cleaner, polydata);
  cleaner->Update();
  return cleaner->GetOutput();
}

// -----------------------------------------------------------------------------
/// Get number of edges
int NumberOfEdges(vtkDataSet *polydata)
{
  EdgeTable edgeTable(polydata);
  return edgeTable.NumberOfEdges();
}

// -----------------------------------------------------------------------------
/// Determine number of connected components
int NumberOfComponents(vtkPolyData *polydata)
{
  vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivity;
  connectivity = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
  connectivity->SetExtractionModeToAllRegions();
  SetVTKInput(connectivity, polydata);
  connectivity->Update();
  return connectivity->GetNumberOfExtractedRegions();
}

// -----------------------------------------------------------------------------
/// Get center point and radius of bounding sphere
inline double GetBoundingSphereRadius(vtkGenericCell *cell, double c[3])
{
  double p[3];
  // Get center of bounding sphere
  c[0] = .0, c[1] = .0, c[2] = .0;
  double *weights = new double[cell->GetNumberOfPoints()];
  int subId = cell->GetParametricCenter(p);
  cell->EvaluateLocation(subId, p, c, weights);
  delete[] weights;
  // Get radius of bounding sphere
  double r = .0;
  vtkPoints *points = cell->GetPoints();
  for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
    points->GetPoint(i, p);
    r = max(r, vtkMath::Distance2BetweenPoints(c, p));
  }
  return sqrt(r);
}

// -----------------------------------------------------------------------------
/// Auxiliary functor which counts the number of self-intersections of a triangulated surface mesh
struct CountTriangleTriangleIntersections
{
  vtkPolyData             *_DataSet;
  vtkAbstractPointLocator *_PointLocator;
  vtkDataArray            *_Mask;
  int                      _NumberOfIntersections;

  CountTriangleTriangleIntersections() : _NumberOfIntersections(0) {}

  CountTriangleTriangleIntersections(const CountTriangleTriangleIntersections &other, split)
  :
    _DataSet(other._DataSet),
    _PointLocator(other._PointLocator),
    _Mask(other._Mask),
    _NumberOfIntersections(0)
  {}

  void join(const CountTriangleTriangleIntersections &other)
  {
    _NumberOfIntersections += other._NumberOfIntersections;
  }

  void operator ()(const blocked_range<vtkIdType> &re)
  {
    vtkSmartPointer<vtkGenericCell> cell1   = vtkSmartPointer<vtkGenericCell>::New();
    vtkSmartPointer<vtkGenericCell> cell2   = vtkSmartPointer<vtkGenericCell>::New();
    vtkSmartPointer<vtkIdList>      ptIds   = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList>      cellIds = vtkSmartPointer<vtkIdList>::New();

    unsigned short ncells;
    vtkIdType *cells, cellId2;
    double r1, c1[3], v1[3], v2[3], v3[3], u1[3], u2[3], u3[3];

    for (vtkIdType cellId1 = re.begin(); cellId1 != re.end(); ++cellId1) {
      _DataSet->GetCell(cellId1, cell1);
      if (cell1->GetNumberOfPoints() != 3) continue;
      _DataSet->GetPoint(cell1->GetPointId(0), v1);
      _DataSet->GetPoint(cell1->GetPointId(1), v2);
      _DataSet->GetPoint(cell1->GetPointId(2), v3);
      r1 = GetBoundingSphereRadius(cell1, c1);
      _PointLocator->FindPointsWithinRadius(3.0 * r1, c1, ptIds);
      cellIds->Reset();
      for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
        _DataSet->GetPointCells(ptIds->GetId(i), ncells, cells);
        for (unsigned short j = 0; j < ncells; ++j) {
          cellIds->InsertUniqueId(cells[j]);
        }
      }
      for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
        cellId2 = cellIds->GetId(i);
        if (cellId2 <= cellId1) continue; // unordered pairs of cells
        _DataSet->GetCell(cellId2, cell2);
        if (cell2->GetNumberOfPoints() != 3) continue;
        if (cell1->GetPointIds()->IsId(cell2->GetPointId(0)) >= 0) continue;
        if (cell1->GetPointIds()->IsId(cell2->GetPointId(1)) >= 0) continue;
        if (cell1->GetPointIds()->IsId(cell2->GetPointId(2)) >= 0) continue;
        _DataSet->GetPoint(cell2->GetPointId(0), u1);
        _DataSet->GetPoint(cell2->GetPointId(1), u2);
        _DataSet->GetPoint(cell2->GetPointId(2), u3);
        if (Triangle::TriangleTriangleIntersection(v1, v2, v3, u1, u2, u3) != 0) {
          _Mask->SetComponent(cellId1, 0, 1.0);
          _Mask->SetComponent(cellId2, 0, 1.0);
          ++_NumberOfIntersections;
        }
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Count number of self-intersections of triangulated surface mesh
int NumberOfTriangleTriangleIntersections(vtkPolyData *polydata)
{
  vtkSmartPointer<vtkDataArray> mask = vtkSmartPointer<vtkUnsignedCharArray>::New();
  mask->SetName("TriangleTriangleIntersection");
  mask->SetNumberOfComponents(1);
  mask->SetNumberOfTuples(polydata->GetNumberOfCells());
  mask->FillComponent(0, .0);
  polydata->GetCellData()->AddArray(mask);
  vtkSmartPointer<vtkOctreePointLocator> octree = vtkSmartPointer<vtkOctreePointLocator>::New();
  octree->SetDataSet(polydata);
  octree->BuildLocator();
  CountTriangleTriangleIntersections count;
  count._DataSet      = polydata;
  count._PointLocator = octree;
  count._Mask         = mask;
  parallel_reduce(blocked_range<vtkIdType>(0, polydata->GetNumberOfCells()), count);
  return count._NumberOfIntersections;
}

// -----------------------------------------------------------------------------
/// Auxiliary functor which counts the number of vertices inside the polyhedron
struct CountVerticesInsidePolyhedron
{
  Polyhedron   *_Polyhedron;
  vtkDataArray *_Mask;
  int           _Num;

  CountVerticesInsidePolyhedron() : _Num(0) {}
  CountVerticesInsidePolyhedron(const CountVerticesInsidePolyhedron &other, split)
  :
    _Polyhedron(other._Polyhedron), _Mask(other._Mask), _Num(0)
  {}

  void join(const CountVerticesInsidePolyhedron &other)
  {
    _Num += other._Num;
  }

  void operator ()(const blocked_range<int> &re)
  {
    double p[3];
    for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
      _Polyhedron->GetPoint(ptId, p);
      if (_Polyhedron->IsInside(p)) {
        _Mask->SetComponent(ptId, 0, 1);
        ++_Num;
      } else {
        _Mask->SetComponent(ptId, 0, 0);
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Point-inside-polyhedron test
int NumberOfVerticesInsidePolyhedron(vtkPolyData *polydata)
{
  vtkSmartPointer<vtkDataArray> mask = vtkSmartPointer<vtkUnsignedCharArray>::New();
  mask->SetName("VertexInsidePolyhedron");
  mask->SetNumberOfComponents(1);
  mask->SetNumberOfTuples(polydata->GetNumberOfPoints());
  polydata->GetPointData()->AddArray(mask);
  Polyhedron polyhedron(polydata);
  CountVerticesInsidePolyhedron count;
  count._Polyhedron = &polyhedron;
  count._Mask       = mask;
  parallel_reduce(blocked_range<int>(0, polyhedron.NumberOfPoints()), count);
  return count._Num;
}

// -----------------------------------------------------------------------------
/// Compute volume of cells
struct ComputeCellVolumes
{
  vtkPointSet  *_PointSet;
  vtkDataArray *_Volume;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double volume;
    vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
    for (vtkIdType cellId = re.begin(); cellId != re.end(); ++cellId) {
      _PointSet->GetCell(cellId, cell);
      volume = ComputeVolume(cell);
      _Volume->SetComponent(cellId, 0, volume);
    }
  }

  static vtkSmartPointer<vtkDataArray> Run(vtkPointSet *pointset)
  {
    vtkSmartPointer<vtkDataArray> volume;
    volume = vtkSmartPointer<vtkFloatArray>::New();
    volume->SetName("Volume");
    volume->SetNumberOfComponents(1);
    volume->SetNumberOfTuples(pointset->GetNumberOfCells());
    pointset->GetCellData()->AddArray(volume);
    ComputeCellVolumes body;
    body._PointSet = pointset;
    body._Volume   = volume;
    parallel_for(blocked_range<vtkIdType>(0, pointset->GetNumberOfCells()), body);
    return volume;
  }
};

// -----------------------------------------------------------------------------
const char *DataObjectTypeString(int type)
{
  switch (type) {
    case VTK_IMAGE_DATA:        return "Image";
    case VTK_POLY_DATA:         return "Polydata";
    case VTK_UNSTRUCTURED_GRID: return "Unstructured grid";
    case VTK_STRUCTURED_GRID:   return "Structured grid";
    case VTK_STRUCTURED_POINTS: return "Structured points";
    default: return "unknown";
  }
}

// -----------------------------------------------------------------------------
const char *CellTypeString(int type)
{
  switch (type) {
    case VTK_EMPTY_CELL:      return "empty cells";
    case VTK_VERTEX:          return "vertices";
    case VTK_POLY_VERTEX:     return "poly-vertices";
    case VTK_LINE:            return "lines";
    case VTK_TRIANGLE:        return "triangles";
    case VTK_TRIANGLE_STRIP:  return "triangle strips";
    case VTK_POLYGON:         return "polygons";
    case VTK_QUAD:            return "quads";
    case VTK_TETRA:           return "tetrahedra";
    case VTK_VOXEL:           return "voxels";
    case VTK_HEXAHEDRON:      return "hexadrons";
    case VTK_HEXAGONAL_PRISM: return "hexagonal prisms";
    case VTK_PYRAMID:         return "pyramids";
    default:                  return "unknown cells";
  }
}

#endif // HAVE_MIRTK_PointSet
// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  EXPECTS_POSARGS(1);

  bool image_info = false, dof_info = false, pointset_info = false;

  // Read input file
  const char *fname = POSARG(1);
  unique_ptr<ImageReader> image_reader;
  #ifdef HAVE_MIRTK_Transformation
    unique_ptr<Transformation> dof;
    if (Transformation::CheckHeader(fname)) {
      dof.reset(Transformation::New(fname));
      dof_info = true;
    }
    else
  #endif // HAVE_MIRTK_Transformation
  {
    image_reader.reset(ImageReader::TryNew(fname));
    image_info = (image_reader != nullptr);
  }
  #ifdef HAVE_MIRTK_PointSet
  vtkSmartPointer<vtkPointSet> pointset;
  vtkSmartPointer<vtkPolyData> polydata;
  if (!image_info && !dof_info) {
    pointset = ReadPointSet(POSARG(1), NULL, false);
    if (pointset) {
      polydata = ConvertToPolyData(pointset);
      polydata->BuildLinks();
      pointset_info = true;
    }
  }
  #endif // HAVE_MIRTK_PointSet

  if (!image_info && !dof_info && !pointset_info) {
    FatalError("Unknown input file type or file format not supported!");
  }
 
  //////////////////////////////////////////////////////////////////////////////
  // Image info
  if (image_info) {

    bool attributes = false;

    for (ALL_OPTIONS) {
      if (OPTION("-a") || OPTION("-attributes") || OPTION("-attr")) attributes = true;
      else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
    }

    unique_ptr<BaseImage> image(image_reader->Run());
    if (!attributes) {
      cout << "Information from ImageReader::Print\n\n";
      image_reader->Print();
      cout << "\nInformation from BaseImage::Print\n\n";
    }
    image->Print();
    if (!attributes) {
      double min, max;
      image->GetMinMaxAsDouble(&min, &max);
      cout << "\nMinimum intensity: " << min << "\n";
      cout << "Maximum intensity: " << max << "\n";
      cout << "\nImage to world matrix\n";
      image->GetImageToWorldMatrix().Print();
      cout << "\nWorld to image matrix\n";
      image->GetWorldToImageMatrix().Print();
      if (!image->GetAffineMatrix().IsIdentity()) {
        cout << "Including affine transformation\n";
        #ifdef HAVE_MIRTK_Transformation
          AffineTransformation transformation;
          transformation.PutMatrix(image->GetAffineMatrix());
          transformation.Print(Indent(1));
        #else // HAVE_MIRTK_Transformation
          image->GetAffineMatrix().Print();
        #endif // HAVE_MIRTK_Transformation
      }
    }

  } // image_info

  //////////////////////////////////////////////////////////////////////////////
  // Transformation info
#ifdef HAVE_MIRTK_Transformation
  if (dof_info) {

    bool attributes = false;
    bool type_name  = false;
    bool type_id    = false;

    for (ALL_OPTIONS) {
      if      (OPTION("-attributes") || OPTION("-attr") || OPTION("-a")) attributes = true;
      else if (OPTION("-type-name")  || OPTION("-type")) type_name  = true;
      else if (OPTION("-type-id")    || OPTION("-id"))   type_id    = true;
      else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
    }
    if (!type_name && !type_id && !attributes) {
      attributes = true;
    }

    if (type_name) {
      cout << dof->NameOfClass();
    }
    if (type_id) {
      if (type_name) cout << " (";
      cout << dof->TypeOfClass();
      if (type_name) cout << ")";
    }
    if (attributes) {
      if (type_name || type_id) cout << "\n";
      dof->Print();
    }

  } // dof_info
#endif // HAVE_MIRTK_Transformation

  //////////////////////////////////////////////////////////////////////////////
  // Point set info
#ifdef HAVE_MIRTK_PointSet
  if (pointset_info) {

    using namespace mirtk::data::statistic;

    bool check_intersections = false;
    bool report_cell_volumes = false;
    bool report_edge_lengths = false;
    bool report_surface_area = false;

    const char *output_pointset_name = NULL;
    const char *output_surface_name  = NULL;

    for (ALL_OPTIONS) {
      if (OPTION("-self-intersections")) check_intersections = true;
      else if (OPTION("-edgelength")) report_edge_lengths = true;
      else if (OPTION("-area")) report_surface_area = true;
      else if (OPTION("-vol") || OPTION("-volume") || OPTION("-volumes")) {
        report_cell_volumes = true;
      }
      else if (OPTION("-o") || OPTION("-out") || OPTION("-output")) {
        output_pointset_name = ARGUMENT;
      }
      else if (OPTION("-s") || OPTION("-output-surface")) {
        output_surface_name = ARGUMENT;
      }
      else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
    }

    check_intersections = (check_intersections || verbose > 1);
    report_cell_volumes = (report_cell_volumes || verbose > 0);
    report_edge_lengths = (report_edge_lengths || verbose > 0);
    report_surface_area = (report_surface_area || verbose > 0);

    // Read input dataset
    vtkSmartPointer<vtkPointSet> pointset = ReadPointSet(POSARG(1));
    vtkSmartPointer<vtkPolyData> polydata = ConvertToPolyData(pointset);
    polydata->BuildLinks();

    // Cell types
    cout << DataObjectTypeString(pointset->GetDataObjectType()) << ":" << endl;
    cout << "  No. of points:             " << pointset->GetNumberOfPoints() << endl;
    if (verbose) {
      cout << "  No. of edges:              " << NumberOfEdges(pointset) << endl;
    }
    cout << "  No. of cells:              " << pointset->GetNumberOfCells() << endl;
    if (pointset->GetNumberOfCells() > 0) {
      OrderedMap<int, int> cellTypeHist;
      OrderedMap<int, int>::iterator it;
      for (vtkIdType cellId = 0; cellId < pointset->GetNumberOfCells(); ++cellId) {
        int type = pointset->GetCellType(cellId);
        it = cellTypeHist.find(type);
        if (it == cellTypeHist.end()) cellTypeHist[type] = 1;
        else ++it->second;
      }
      for (it = cellTypeHist.begin(); it != cellTypeHist.end(); ++it) {
        string s = CellTypeString(it->first);
        cout << "    No. of " << setw(18) << left << (s + ": ") << it->second << endl;
      }
      cout << "  Maximum cell size:         " << pointset->GetMaxCellSize() << endl;
    }

    // Bounding box
    double bounds[6];
    pointset->GetBounds(bounds);

    cout << "  Bounding box:              [" << bounds[0] << ", " << bounds[1]
                                  << "] x [" << bounds[2] << ", " << bounds[3]
                                  << "] x [" << bounds[4] << ", " << bounds[5] << "]" << endl;
    cout << "  Center point:              (" << .5 * (bounds[1] + bounds[0]) << ", "
                                             << .5 * (bounds[3] + bounds[2]) << ", "
                                             << .5 * (bounds[5] + bounds[4]) << ")" << endl;

    // Attributes
    if (pointset->GetPointData()->GetNumberOfArrays() > 0) {
    	cout << "\nPoint data arrays:" << endl;
      for (int i = 0; i < pointset->GetPointData()->GetNumberOfArrays(); ++i) {
        vtkDataArray *scalars = pointset->GetPointData()->GetArray(i);
        printf(" %2d %-24s (dim: %2d)\n", i, scalars->GetName(), scalars->GetNumberOfComponents());
      }
    }
    if (pointset->GetCellData()->GetNumberOfArrays() > 0) {
      cout << "\nCell data arrays:" << endl;
      for (int i = 0; i < pointset->GetCellData()->GetNumberOfArrays(); ++i) {
        vtkDataArray *scalars = pointset->GetCellData()->GetArray(i);
        printf(" %2d %-24s (dim: %2d)\n", i, scalars->GetName(), scalars->GetNumberOfComponents());
      }
    }

    if (pointset->GetNumberOfCells() > 0) {

      // Compute cell volumes
      if (report_cell_volumes && !IsSurfaceMesh(pointset)) {
        vtkSmartPointer<vtkDataArray> volume = ComputeCellVolumes::Run(pointset);
        double min_volume, max_volume, avg_volume, std_volume;
        Extrema::Calculate(min_volume, max_volume, volume);
        NormalDistribution::Calculate(avg_volume, std_volume, volume);
        cout << "\n";
        cout << "  Average cell volume: " << avg_volume << "\n";
        cout << "  Cell volume StDev:   " << std_volume << "\n";
        cout << "  Minimum cell volume: " << min_volume << "\n";
        cout << "  Maximum cell volume: " << max_volume << "\n";
        cout.flush();
      }

      // Discard degenerate cells and unused points
      polydata = CleanPolyData(polydata);
      polydata->SetVerts(NULL);
      polydata->SetLines(NULL);
      polydata = CleanPolyData(polydata);
      polydata->BuildLinks();

      // Compute edge table
      EdgeTable edgeTable(polydata);

      // Euler characteristic / Genus
      const int    noOfVerts = polydata->GetNumberOfPoints();
      const int    noOfFaces = polydata->GetNumberOfCells();
      const int    noOfComps = NumberOfComponents(polydata);
      const int    noOfEdges = edgeTable.NumberOfEdges();
      const int    euler     = noOfVerts - noOfEdges + noOfFaces;
      const double genus     = 0.5 * (2 * noOfComps - euler);

      cout << "\n";
      cout << "Surface mesh:\n";
      cout << "  V " << noOfVerts << "\n";
      cout << "  E " << noOfEdges << "\n";
      cout << "  F " << noOfFaces << "\n";
      cout << "  C " << noOfComps << "\n";
      cout << "\n";
      cout << "  Euler characteristic / Genus (V - E + F = 2C - 2g)\n";
      cout << "    Euler: " << euler << "\n";
      cout << "    Genus: " << genus << "\n";
      cout.flush();

      if (report_edge_lengths) {
        double min_length, max_length, avg_length, std_length;
        EdgeLengthNormalDistribution(polydata->GetPoints(), edgeTable, avg_length, std_length);
        GetMinMaxEdgeLength(polydata->GetPoints(), edgeTable, min_length, max_length);
        cout << "\n";
        cout << "  Average edge length: " << avg_length << "\n";
        cout << "  Edge length StDev:   " << std_length << "\n";
        cout << "  Minimum edge length: " << min_length << "\n";
        cout << "  Maximum edge length: " << max_length << "\n";
        cout.flush();
      }

      if (report_surface_area) {
        double sum_area = Area(polydata, true);
        vtkDataArray *area = polydata->GetCellData()->GetArray("Area");
        double min_area, max_area, avg_area, std_area;
        Extrema::Calculate(min_area, max_area, area);
        NormalDistribution::Calculate(avg_area, std_area, area);
        cout << "\n";
        cout << "  Average cell area:   " << avg_area << "\n";
        cout << "  Cell area StDev:     " << std_area << "\n";
        cout << "  Minimum cell area:   " << min_area << "\n";
        cout << "  Maximum cell area:   " << max_area << "\n";
        cout << "  Total surface area:  " << sum_area << "\n";
        cout.flush();
      }

      // Self-intersections
      if (check_intersections) {
        cout << endl;
        cout << "  No. of triangle/triangle intersections: ", cout.flush();
        cout << NumberOfTriangleTriangleIntersections(polydata) << endl;
        if (verbose > 2) {
          int count = 0;
          vtkDataArray *tritri = polydata->GetCellData()->GetArray("TriangleTriangleIntersection");
          for (vtkIdType cellId = 0; cellId < polydata->GetNumberOfCells(); ++cellId) {
            if (tritri->GetComponent(cellId, 0) != .0) {
              ++count;
              cout << "  " << right << setw(5) << count;
              if      (count == 1) cout << "st";
              else if (count == 2) cout << "nd";
              else if (count == 3) cout << "rd";
              else                 cout << "th";
              cout << " intersected cell: ID = " << cellId;
              double area = ComputeArea(polydata->GetCell(cellId));
              cout << ", area = " << area << "\n";
            }
          }
          cout.flush();
        }
    //    cout << "  No. of vertices inside polyhedron: ", cout.flush();
    //    cout << NumberOfVerticesInsidePolyhedron(polydata) << endl;
      }
    }

    // Write output point set (possibly with additional attributes)
    if (output_pointset_name) WritePointSet(output_pointset_name, pointset);

    // Write output surface mesh (possibly with additional attributes)
    if (output_surface_name && polydata->GetNumberOfCells() > 0) {
      WritePointSet(output_surface_name, polydata);
    }

  } // dof_info
#endif // HAVE_MIRTK_PointSet

  return 0;
}

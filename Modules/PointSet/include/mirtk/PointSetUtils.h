/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#ifndef MIRTK_PointSetUtils_H
#define MIRTK_PointSetUtils_H

#include "mirtk/UnorderedSet.h"
#include "mirtk/Array.h"
#include "mirtk/List.h"
#include "mirtk/Pair.h"

#include "mirtk/Math.h"
#include "mirtk/Vtk.h"
#include "mirtk/VtkMath.h"

class vtkDataSet;
class vtkDataSetAttributes;
class vtkPointSet;
class vtkPolyData;
class vtkImageData;
class vtkImageStencilData;
class vtkDataArray;
class vtkCell;


namespace mirtk {


struct ImageAttributes;
class BaseImage;
class PointSet;
class EdgeTable;
class Vector;

template <typename> struct Vector3D;

/// List of pairs of edge end point IDs
typedef List<Pair<int, int>> EdgeList;

// =============================================================================
// VTK / MIRTK type conversion
// =============================================================================

/// Add points of vtkPointSet to mirtk::PointSet
///
/// \param[out] oset Point set to which the vtkPointSet points are added.
/// \param[in]  iset VTK point set.
void AddPoints(PointSet &oset, vtkPointSet *iset);

// =============================================================================
// Point set domain
// =============================================================================

/// Determine dimension of data set
int Dimension(vtkDataSet *);

/// Determine bounding box of point set
///
/// \param[in] data Data set.
/// \param[in] dx   Desired lattice spacing along x axis.
///                 If non-positive, the average interval between the x
///                 coordinates of the input points is used.
/// \param[in] dy   Desired lattice spacing along y axis.
///                 If non-positive, the average interval between the y
///                 coordinates of the input points is used.
/// \param[in] dz   Desired lattice spacing along z axis.
///                 If non-positive, the average interval between the z
///                 coordinates of the input points is used.
///
/// \returns Attributes of oriented minimum-volume bounding box.
ImageAttributes PointSetDomain(vtkPointSet *data, double dx = -1, double dy = -1, double dz = -1);

/// Determine bounding box of point set
///
/// \param[in] data Data set.
/// \param[in] ds   Desired lattice spacing.
///                 If an entry is non-positive, the average interval between
///                 the coordinates of the input points along this axis is used.
///
/// \returns Attributes of oriented minimum-volume bounding box.
ImageAttributes PointSetDomain(vtkPointSet *data, const Vector3D<double> &ds);

/// Determine bounding box of polydata points
///
/// \param[in] data Data set.
/// \param[in] dx   Desired lattice spacing along x axis.
///                 If non-positive, the average interval between the x
///                 coordinates of the input points is used.
/// \param[in] dy   Desired lattice spacing along y axis.
///                 If non-positive, the average interval between the y
///                 coordinates of the input points is used.
/// \param[in] dz   Desired lattice spacing along z axis.
///                 If non-positive, the average interval between the z
///                 coordinates of the input points is used.
///
/// \returns Attributes of oriented minimum-volume bounding box.
ImageAttributes PolyDataDomain(vtkPolyData *data, double dx = -1, double dy = -1, double dz = -1);

/// Determine bounding box of polydata points
///
/// \param[in] data Data set.
/// \param[in] ds   Desired lattice spacing.
///                 If an entry is non-positive, the average interval between
///                 the coordinates of the input points along this axis is used.
///
/// \returns Attributes of oriented minimum-volume bounding box.
ImageAttributes PolyDataDomain(vtkPolyData *data, const Vector3D<double> &ds);

// =============================================================================
// Point/cell data
// =============================================================================

/// Map string to vtkDataSetAttributes::AttributeTypes enumeration value
///
/// @param type Case insensitive attribute type name (e.g. "scalars", "NORMALS")
///
/// @return VTK data set attribute type ID (e.g., vtkDataSetAttributes::SCALARS)
///         or -1 if string does not name a known attribute type.
int PolyDataAttributeType(const char *type);

/// Get point data array using case insensitive name
///
/// @param[in]  data Point or cell data attributes (cf. vtkPolyData::GetPointData, vtkPolyData::GetCellData).
/// @param[in]  name Case insenitive name of data array.
/// @param[out] loc  Set to array index if not @c NULL.
///
/// @return Pointer to data array or @c NULL if not found.
vtkDataArray *GetArrayByCaseInsensitiveName(vtkDataSetAttributes *data, const char *name, int *loc = NULL);

/// Copy named data array from one dataset attributes to another
///
/// If an array with the given name exists in the destination dataset attributes,
/// it's data is overridden using the vtkDataSetAttributes::DeepCopy function.
/// Otherwise, a new array is added.
///
/// @param dst  Destination dataset attributes.
/// @param src  Source dataset attributes.
/// @param name Case insensitive name of data array.
///
/// @return Index of array in destination or -1 if source array not found.
int DeepCopyArrayUsingCaseInsensitiveName(vtkDataSetAttributes *dst,
                                          vtkDataSetAttributes *src,
                                          const char *name);

/// Determine whether data array name suggests it contains categorical values
///
/// \param[in] name Data array name.
bool IsCategoricalArrayName(const string &name);

// =============================================================================
// Cell attributes
// =============================================================================

/// Compute area of cell
///
/// @return Area of cell or NaN if cell type is not supported.
double ComputeArea(vtkCell *cell);

/// Compute volume of cell
///
/// @return Volume of cell or NaN if cell type is not supported.
double ComputeVolume(vtkCell *cell);

/// Compute orthogonal vectors spanning the tangent plane of a cell
///
/// @param[in]  n  Cell normal vector.
/// @param[out] e1 First  normal vector in tangent plane.
/// @param[out] e2 Second normal vector in tangent plane.
///
/// @returns Whether tangent vectors are valid.
inline bool ComputeTangents(const double n[3], double e1[3], double e2[3])
{
  double v[3] = {n[1], n[2], n[0]};
  vtkMath::Cross(n, v, e1);
  if (vtkMath::Dot(e1, e1) < 1e-6) {
    v[1] *= -1.0;
    vtkMath::Cross(n, v, e1);
    if (vtkMath::Dot(e1, e1) < 1e-6) return false;
  }
  vtkMath::Cross(n, e1, e2);
  vtkMath::Normalize(e1);
  vtkMath::Normalize(e2);
  return true;
}

// -----------------------------------------------------------------------------
/// Compute 4 (8) equally spaced tangent vectors (angular sampling of 45 degrees)
inline bool ComputeTangents(const double n[3], double e1[3], double e2[3], double e3[3], double e4[3])
{
  if (!ComputeTangents(n, e1, e3)) return false;
  const double sqrt2 = sqrt(2.0);
  e2[0] = (e1[0] + e3[0]) / sqrt2;
  e2[1] = (e1[1] + e3[1]) / sqrt2;
  e2[2] = (e1[2] + e3[2]) / sqrt2;
  e4[0] = (e1[0] - e3[0]) / sqrt2;
  e4[1] = (e1[1] - e3[1]) / sqrt2;
  e4[2] = (e1[2] - e3[2]) / sqrt2;
  return true;
}

// =============================================================================
// Basic point set manipulation
// =============================================================================

/// Get boundary surface mesh
///
/// @param[in] dataset     Dataset whose boundary surface is extracted.
/// @param[in] passPtIds   Whether to pass point array with IDs of points in @p dataset.
/// @param[in] passCellIds Whether to pass cell array with IDs of cells in @p dataset.
///
/// @return Boundary surface mesh.
vtkSmartPointer<vtkPolyData> DataSetSurface(vtkSmartPointer<vtkDataSet> dataset,
                                            bool passPtIds   = false,
                                            bool passCellIds = false);

/// Translate point set such that center is at origin
void Center(vtkSmartPointer<vtkPointSet> pointset);

/// Scale point set around center
void Scale(vtkSmartPointer<vtkPointSet> pointset, double);

// =============================================================================
// Surface meshes
// =============================================================================

/// Calculate are of triangle with given edge length
inline double EdgeLengthToTriangleArea(double l)
{
  return 0.25 * sqrt(3.0) * l * l;
}

/// Determine whether a point set is a surface mesh
bool IsSurfaceMesh(vtkDataSet *);

/// Check whether given point set is a triangular mesh
bool IsTriangularMesh(vtkDataSet *);

/// Check whether given point set is a tetrahedral mesh
bool IsTetrahedralMesh(vtkDataSet *);

/// Get IDs of end points of boundary edges
UnorderedSet<int> BoundaryPoints(vtkDataSet *, const EdgeTable * = nullptr);

/// Get list of all boundary edges
EdgeList BoundaryEdges(vtkDataSet *);

/// Get list of all boundary edges
EdgeList BoundaryEdges(vtkDataSet *, const EdgeTable &);

/// Get list of edges with given end point
///
/// \param[in] edges List of edges.
/// \param[in] ptId  Edge end point.
///
/// \returns List of edges, where ptId is always the first entry.
EdgeList GetPointEdges(const EdgeList &edges, int ptId);

/// Get list of edges with given end point and remove them from input list
///
/// \param[in] edges List of edges.
/// \param[in] ptId  Edge end point.
///
/// \returns List of edges, where ptId is always the first entry.
EdgeList PopPointEdges(EdgeList &edges, int ptId);

/// Get connected boundary segments as (closed) line strips
Array<Array<int> > BoundarySegments(vtkDataSet *, const EdgeTable * = nullptr);

/// Number of points
int NumberOfPoints(vtkDataSet *);

/// Number of edges
int NumberOfEdges(vtkDataSet *, const EdgeTable * = nullptr);

/// Number of faces
int NumberOfFaces(vtkDataSet *);

/// Number of empty/deleted cells
int NumberOfEmptyCells(vtkDataSet *);

/// Number of connected components
int NumberOfConnectedComponents(vtkDataSet *);

/// Number of connected boundary components
int NumberOfBoundarySegments(vtkDataSet *, const EdgeTable * = nullptr);

/// Euler characeteristic, i.e., V - E + F
int EulerCharacteristic(vtkDataSet *dataset,
                        const EdgeTable &,
                        int *npoints = nullptr,
                        int *nedges  = nullptr,
                        int *nfaces  = nullptr);

/// Euler characeteristic, i.e., V - E + F
int EulerCharacteristic(vtkDataSet *dataset,
                        int *npoints = nullptr,
                        int *nedges  = nullptr,
                        int *nfaces  = nullptr);

/// Genus of surface mesh
double Genus(vtkDataSet *dataset,
             const EdgeTable &,
             int *npoints = nullptr,
             int *nedges  = nullptr,
             int *nfaces  = nullptr,
             int *nbounds = nullptr,
             int *ncomps  = nullptr,
             int *euler   = nullptr);

/// Genus of surface mesh
double Genus(vtkDataSet *dataset,
             int *npoints = nullptr,
             int *nedges  = nullptr,
             int *nfaces  = nullptr,
             int *nbounds = nullptr,
             int *ncomps  = nullptr,
             int *euler   = nullptr);

/// Area of surface mesh
double Area(vtkPolyData *, bool per_cell = false);

/// Area of surface mesh
inline double Area(vtkSmartPointer<vtkPolyData> surface, bool per_cell = false)
{
  return Area(surface.GetPointer(), per_cell);
}


/// Area of point set surface
double Area(vtkSmartPointer<vtkPointSet>);

/// Compute edge lengths of point set given a precomputed edge table
Vector EdgeLengths(vtkSmartPointer<vtkPoints>, const EdgeTable &);

/// Compute squared edge lengths of point set given a precomputed edge table
Vector SquaredEdgeLengths(vtkSmartPointer<vtkPoints>, const EdgeTable &);

/// Determine average edge length of point set given a precomputed edge table
double AverageEdgeLength(vtkSmartPointer<vtkPoints>, const EdgeTable &);

/// Determine average edge length of point set
double AverageEdgeLength(vtkSmartPointer<vtkPointSet>);

/// Determine median edge length of point set given a precomputed edge table
double MedianEdgeLength(vtkSmartPointer<vtkPoints>, const EdgeTable &);

/// Determine median edge length of point set
double MedianEdgeLength(vtkSmartPointer<vtkPointSet>);

/// Determine average edge length of point set given a precomputed edge table
///
/// This function only considers edges with a length within in the 5th and 95th
/// percentile of all edge lengths. It thus ignores extrem short/long edges.
double RobustAverageEdgeLength(vtkSmartPointer<vtkPoints>, const EdgeTable &);

/// Determine average edge length of point set
double RobustAverageEdgeLength(vtkSmartPointer<vtkPointSet>);

/// Determine minimum and maximum edge length
///
/// \param[in]  points    Points.
/// \param[in]  edgeTable Edge table.
/// \param[out] min       Minimum edge length.
/// \param[out] max       Maximum edge length.
void GetMinMaxEdgeLength(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable, double &min, double &max);

/// Determine minimum and maximum edge length
///
/// \param[in]  pointset Point set.
/// \param[out] min      Minimum edge length.
/// \param[out] max      Maximum edge length.
void GetMinMaxEdgeLength(vtkSmartPointer<vtkPointSet> pointset, double &min, double &max);

/// Determine minimum edge length
///
/// \param[in] points    Points.
/// \param[in] edgeTable Edge table.
///
/// \return Minimum edge length.
double MinEdgeLength(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable);

/// Determine minimum edge length
///
/// \param[in] pointset Point set.
///
/// \return Minimum edge length.
double MinEdgeLength(vtkSmartPointer<vtkPointSet> pointset);

/// Determine maximum edge length
///
/// \param[in] points    Points.
/// \param[in] edgeTable Edge table.
///
/// \return Maximum edge length.
double MaxEdgeLength(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable);

/// Determine maximum edge length
///
/// \param[in] pointset Point set.
///
/// \return Maximum edge length.
double MaxEdgeLength(vtkSmartPointer<vtkPointSet> pointset);

/// Compute statistics of edge lengths
///
/// \param[in]  points    Points.
/// \param[in]  edgeTable Edge table.
/// \param[out] mean      Average edge length.
/// \param[out] sigma     Standard deviation of edge length.
void EdgeLengthNormalDistribution(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable, double &mean, double &sigma);

/// Compute statistics of edge lengths
///
/// \param[in]  pointset Point set.
/// \param[out] mean     Average edge length.
/// \param[out] sigma    Standard deviation of edge length.
void EdgeLengthNormalDistribution(vtkSmartPointer<vtkPointSet> pointset, double &mean, double &sigma);

/// Get approximate volume enclosed by polygonal mesh
double Volume(vtkSmartPointer<vtkPolyData>);

/// Get convex hull of point set
///
/// @param[in] pointset Input point set.
/// @param[in] levels   Parameter of vtkHull::AddRecursiveSpherePlanes.
///
/// @return Convex hull of input point set.
vtkSmartPointer<vtkPolyData> ConvexHull(vtkSmartPointer<vtkPointSet> pointset, int levels = 3);

/// Triangulate surface mesh
vtkSmartPointer<vtkPolyData> Triangulate(vtkSmartPointer<vtkPolyData>);

/// Tetrahedralize the interior of a piecewise linear complex (PLC)
vtkSmartPointer<vtkPointSet> Tetrahedralize(vtkSmartPointer<vtkPointSet>);

/// Instantiate new VTK image mask
///
/// Note that vtkImageData has no implicit orientation. Therefore we just
/// convert the image in voxel coordinates (origin at 0 and voxel size 1x1x1)
/// and instead convert points to voxel coordinates.
vtkSmartPointer<vtkImageData> NewVtkMask(int nx, int ny, int nz);

/// Map point set points to voxel coordinates
vtkSmartPointer<vtkPointSet> WorldToImage(vtkSmartPointer<vtkPointSet> pointset,
                                          const BaseImage             *image);

/// Get inside surface image stencil
vtkSmartPointer<vtkImageStencilData> ImageStencil(vtkSmartPointer<vtkImageData> image,
                                                  vtkSmartPointer<vtkPointSet>  pointset);

/// Convert surface image stencil to binary mask
///
/// \param[out] image   Image data with voxel type IRTK_VOXEL_BINARY.
/// \param[in]  stencil Image stencil.
void ImageStencilToMask(vtkSmartPointer<vtkImageStencilData> stencil,
                        vtkSmartPointer<vtkImageData>        image);


} // namespace mirtk

#endif // MIRTK_PointSetUtils_H

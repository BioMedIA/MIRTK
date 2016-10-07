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

#ifndef MIRTK_RegisteredPointSet_H
#define MIRTK_RegisteredPointSet_H

#include "mirtk/Object.h"

#include "mirtk/Array.h"
#include "mirtk/PointSet.h"
#include "mirtk/Transformation.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/EdgeConnectivity.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkAbstractPointLocator.h"
#include "vtkAbstractCellLocator.h"
#include "vtkIdTypeArray.h"


namespace mirtk {


/**
 * Registered point set
 *
 * A registered point set is aligned with another data set (i.e., an image or
 * another point set) through a given transformation. If no transformation is
 * set (i.e., NULL), the point set is fixed and usually serves as target for
 * the registration. Given the pull-back character of image transformations,
 * the point set which is being transformed is generally the one defined in
 * the space of the target image (i.e., the fixed reference image). We therefore
 * refer to the transformed point set as the "target" instead, whereas the
 * "source" point set is usually not transformed. Eventually, however, both
 * target and source point sets may be transformed, e.g., in case of a
 * symmetric registration.
 *
 * Besides the geometry, i.e., point positions, a point set may also have a
 * topological structure as in case of a tetrahedralized simplicial complex or
 * a boundary surface mesh. Point set distance measures and constraints may be
 * specialized for surface meshes only, such as the surface curvature constraint,
 * for example (CurvatureConstraint). These objective function terms only
 * consider the surface of a given point set which is provided by this class.
 * If the input point set itself is already a vtkPolyData object, the point set
 * and its surface are identical. Otherwise, an instance of this class extracts
 * the surface of the point set using the vtkDataSetSurfaceFilter upon
 * initialization and keeps the point positions and data of the registered point
 * set and its corresponding surface in sync.
 *
 * \note VTK does not make use of the const keyword and none of the member
 *       functions are declared const. Therefore, returning a pointer to a
 *       const VTK object is useless. Thus, this class grants anyone non-const
 *       access to the internal VTK data objects even if these are not to be
 *       modified by the caller.
 */
class RegisteredPointSet : public Object
{
  mirtkObjectMacro(RegisteredPointSet);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Feature scaling function parameters
  struct ScalingFunction
  {
    int    _InputIndex;  ///< Feature point data index in input  dataset
    int    _OutputIndex; ///< Feature point data index in output dataset
    double _Slope;       ///< Slope of linear scaling function
    double _Intercept;   ///< Intercept of linear scaling function
  };

  /// Indices and scaling function parameters of transformed point data
  typedef Array<ScalingFunction> ScalingFunctions;

  /// Adjacency matrix with edge IDs
  typedef mirtk::EdgeTable EdgeTable;

  /// Table of n-connected node neighbors
  typedef EdgeConnectivity NodeNeighbors;

  // ---------------------------------------------------------------------------
  // Attributes

  /// Get (transformed) point set
  vtkPointSet *PointSet() const;

protected:

  /// Untransformed input point set
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPointSet>, InputPointSet);

  /// Untransformed surface of input point set
  mirtkReadOnlyAttributeMacro(vtkSmartPointer<vtkPolyData>, InputSurface);

  /// Whether this point is a surface mesh
  mirtkReadOnlyAttributeMacro(bool, IsSurfaceMesh);

  /// (Cached) Untransformed input points
  /// used by energy terms for gradient calculation
  mirtkReadOnlyAttributeMacro(class PointSet, InputPoints);

  /// (Cached) Untransformed input surface points
  /// used by energy terms for gradient calculation
  class PointSet *_InputSurfacePoints;

  /// Time point of input dataset
  mirtkPublicAttributeMacro(double, InputTime);

  /// Time point of (transformed) dataset
  mirtkPublicAttributeMacro(double, Time);

  /// Current transformation estimate
  mirtkPublicAggregateMacro(const class Transformation, Transformation);

  /// Indices of point data to copy and (optionally) rescale/normalize
  mirtkPublicAttributeMacro(ScalingFunctions, PointDataToCopy);

  /// Whether to copy all point and cell data
  mirtkPublicAttributeMacro(bool, CopyAll);

  /// (Minimum) radius (maximum edge-connectivity) of node neighborhood
  mirtkPublicAttributeMacro(int, NeighborhoodRadius);

  /// Average edge length of input point set
  mirtkReadOnlyAttributeMacro(double, AverageInputEdgeLength);

  /// Average edge length of input surface
  mirtkReadOnlyAttributeMacro(double, AverageInputSurfaceEdgeLength);

  /// Transformed output point set
  mirtkReadOnlyAttributeMacro(vtkSmartPointer<vtkPointSet>, OutputPointSet);

  /// Transformed output surface
  mirtkReadOnlyAttributeMacro(vtkSmartPointer<vtkPolyData>, OutputSurface);

  /// Length of diagonal of input point set bounding box
  mirtkReadOnlyAttributeMacro(double, InputDiameter);

  /// Length of diagonal of point set bounding box
  mirtkReadOnlyAttributeMacro(double, Diameter);

  /// Edge table of point set (computed on demand)
  mutable SharedPtr<const EdgeTable> _EdgeTable;

  /// Edge table of point set surface (computed on demand)
  mutable SharedPtr<const EdgeTable> _SurfaceEdgeTable;

  /// Edge-connectivities / neighborhood of point set nodes (computed on demand)
  mutable NodeNeighbors _NodeNeighbors;

  /// Edge-connectivities / neighborhood of point set surface nodes (computed on demand)
  mutable NodeNeighbors _SurfaceNodeNeighbors;

  /// Input point set surface area (computed on demand)
  mutable double _InputSurfaceArea;

  /// Point set surface area (computed on demand)
  mutable double _SurfaceArea;

  /// Point locator (built on demand)
  vtkSmartPointer<vtkAbstractPointLocator> _PointLocator;

  /// Surface point locator (built on demand)
  vtkSmartPointer<vtkAbstractPointLocator> _SurfacePointLocator;

  /// Cell locator (built on demand)
  vtkSmartPointer<vtkAbstractCellLocator> _CellLocator;

  /// Surface cell locator (built on demand)
  vtkSmartPointer<vtkAbstractCellLocator> _SurfaceCellLocator;

  /// Whether self-update is enabled
  mirtkPublicAttributeMacro(bool, SelfUpdate);

  /// Whether point normals of output surface need to be recomputed (on demand)
  mirtkPublicAttributeMacro(bool, UpdateSurfaceNormals);

  /// Whether face normals of output surface need to be recomputed (on demand)
  mirtkPublicAttributeMacro(bool, UpdateSurfaceFaceNormals);

  /// Domain on which to evaluate transformation if it requires caching
  /// The obtained deformation field is interpolated linearly.
  mirtkPublicAttributeMacro(ImageAttributes, Domain);

  /// Externally pre-computed displacements to use
  mirtkPublicAggregateMacro(GenericImage<double>, ExternalDisplacement);

  /// Cached displacement field evaluated at each lattice point of _Domain
  mirtkComponentMacro(GenericImage<double>, Displacement);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const RegisteredPointSet &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  RegisteredPointSet(vtkPointSet * = NULL, const class Transformation * = NULL);

  /// Copy constructor
  RegisteredPointSet(const RegisteredPointSet &);

  /// Assignment operator
  RegisteredPointSet &operator =(const RegisteredPointSet &);

  /// Destructor
  ~RegisteredPointSet();

  /// Initialize (transformed) point set after input and parameters are set
  ///
  /// \param[in] deep_copy_points  Whether to force a deep copy of the input
  ///                              points. If \c false, the points of the
  ///                              output points is only made when Update
  ///                              transforms the input points.
  void Initialize(bool deep_copy_points = false);

  /// Re-copy input point set/surface points
  void InputPointsChanged();

  /// Called when the output points have been modified
  /// (cf. DeformableSurfaceModel without FFD)
  void PointsChanged();

  /// Pre-initialize edge tables
  ///
  /// The edge tables returned by the Edges and SurfaceEdges functions are
  /// built upon demand when the respective function is called the first time.
  /// If required, a pre-initialization of the edge tables can be requested by
  /// calling this function.
  void BuildEdgeTables();

  /// Pre-initialize node neighborhood tables
  ///
  /// The sets of node neighbors with the given maximum edge-connectivity
  /// (or neighborhood radius) are determined on demand upon the first call
  /// of the Neighbors and SurfaceNeighbors functions. If required,
  /// a pre-initialization of the tables can be requested by calling this
  /// function. Beforehand, the required minimum neighborhood radius must be
  /// set using the NeighborhoodRadius(int) function unless the desired maximum
  /// node connectivity of the neighbors is passed as argument.
  void BuildNeighborhoodTables(int n = -1);

  /// Pre-initialize point/cell locators
  void BuildLocators();

  // ---------------------------------------------------------------------------
  // Input point set

  /// Get number of points
  int NumberOfPoints() const;

  /// Get number of cells
  int NumberOfCells() const;

  /// Number of (input) point set edges
  int NumberOfEdges() const;

  /// Get untransformed input point with specified index
  void GetInputPoint(int, double &, double &, double &) const;

  /// Get untransformed input point with specified index
  void GetInputPoint(int, double *) const;

  /// Get untransformed input point with specified index
  void GetInputPoint(int, Point &) const;

  /// Get untransformed points of input data set
  void GetInputPoints(class PointSet &) const;

  // ---------------------------------------------------------------------------
  // Input point set surface

  /// Whether (input) point set is a polygonal surface mesh
  /// Otherwise, the surface of the point set is extracted by Initialize.
  bool IsSurface() const;

  /// Get number of surface points
  int NumberOfSurfacePoints() const;

  /// Get number of surface cells
  int NumberOfSurfaceCells() const;

  /// Number of (input) surface edges
  int NumberOfSurfaceEdges() const;

  /// Get array which stores for each surface point the input point set point ID
  vtkIdTypeArray *OriginalSurfacePointIds() const;

  /// Get array which stores for each surface cell the input point set cell ID
  vtkIdTypeArray *OriginalSurfaceCellIds() const;

  /// Get untransformed input surface point with specified index
  void GetInputSurfacePoint(int, double &, double &, double &) const;

  /// Get untransformed input surface point with specified index
  void GetInputSurfacePoint(int, double *) const;

  /// Get untransformed input surface point with specified index
  void GetInputSurfacePoint(int, Point &) const;

  /// Get untransformed points of input point set surface
  void GetInputSurfacePoints(class PointSet &) const;

  /// Untransformed points of input point set surface
  const class PointSet &InputSurfacePoints() const;

  /// Area of input point set surface
  double InputSurfaceArea() const;

  // ---------------------------------------------------------------------------
  // Point set access

  /// Implicit conversion to vtkDataSet pointer
  operator vtkDataSet *() const;

  /// Implicit conversion to vtkPointSet pointer
  operator vtkPointSet *() const;

  /// Get points of point set
  vtkPoints *Points() const;

  /// Get edge table of point set mesh
  ///
  /// Not thread-safe unless \c _EdgeTable is already initialized (cf. BuildEdgeTables).
  const EdgeTable *Edges() const;

  /// Get shared pointer to edge table
  ///
  /// Not thread-safe unless \c _EdgeTable is already initialized (cf. BuildEdgeTables).
  SharedPtr<const EdgeTable> SharedEdgeTable() const;

  /// Get edge-connectivity table of point set node neighbors
  ///
  /// Not thread-safe unless \c _NodeNeighbors is already initialized with a
  /// neighborhood radius greater or equal the requested minimum radius
  /// (cf. BuildNeighborhoodTables).
  const NodeNeighbors *Neighbors(int = -1) const;

  /// Get initial point status array if any
  vtkDataArray *InitialStatus() const;

  /// Get point status array if any
  vtkDataArray *Status() const;

  /// Get point locator
  vtkAbstractPointLocator *PointLocator() const;

  /// Get cell locator
  ///
  /// \attention The VTK cell locators are not thread-safe (at least up to VTK 6.3),
  ///            not even after vtkLocator::BuildLocator has been called.
  vtkAbstractCellLocator *CellLocator() const;

  /// Get point with specified index
  void GetPoint(int, double &, double &, double &) const;

  /// Get point with specified index
  void GetPoint(int, double *) const;

  /// Get point with specified index
  void GetPoint(int, Point &) const;

  /// Get points of point set
  void GetPoints(class PointSet &) const;

  // ---------------------------------------------------------------------------
  // Point set surface access

  /// Get output surface
  vtkPolyData *Surface() const;

  /// Get points of point set surface
  vtkPoints *SurfacePoints() const;

  /// Get output surface (point) normals
  ///
  /// Not thread-safe unless \c _UpdateSurfaceNormals is \c false,
  /// i.e., when this function was called by main thread after Update.
  vtkDataArray *SurfaceNormals() const;

  /// Get output surface (cell) normals
  ///
  /// Not thread-safe unless \c _UpdateSurfaceFaceNormals is \c false,
  /// i.e., when this function was called by main thread after Update.
  vtkDataArray *SurfaceFaceNormals() const;

  /// Get edge table of point set surface mesh
  ///
  /// Not thread-safe unless \c _SurfaceEdgeTable (or \c _EdgeTable if input
  /// point set is a surface mesh) is already initialized (cf. BuildEdgeTables).
  const EdgeTable *SurfaceEdges() const;

  /// Get shared pointer to surface edge table
  ///
  /// Not thread-safe unless \c _SurfaceEdgeTable (or \c _EdgeTable if input
  /// point set is a surface mesh) is already initialized (cf. BuildEdgeTables).
  SharedPtr<const EdgeTable> SharedSurfaceEdgeTable() const;

  /// Get edge-connectivity table of point set surface node neighbors
  ///
  /// Not thread-safe unless \c _SurfaceNodeNeighbors is already initialized
  /// with a neighborhood radius greater or equal the requested minimum radius
  /// (cf. BuildNeighborhoodTables).
  const NodeNeighbors *SurfaceNeighbors(int = -1) const;

  /// Get point status array if any
  vtkDataArray *InitialSurfaceStatus() const;

  /// Get point status array if any
  vtkDataArray *SurfaceStatus() const;

  /// Area of point set surface
  double SurfaceArea() const;

  /// Get surface point locator
  vtkAbstractPointLocator *SurfacePointLocator() const;

  /// Get surface cell locator
  ///
  /// \attention The VTK cell locators are not thread-safe (at least up to VTK 6.3),
  ///            not even after vtkLocator::BuildLocator has been called.
  vtkAbstractCellLocator *SurfaceCellLocator() const;

  /// Get point with specified index
  void GetSurfacePoint(int, double &, double &, double &) const;

  /// Get point with specified index
  void GetSurfacePoint(int, double *) const;

  /// Get point with specified index
  void GetSurfacePoint(int, Point &) const;

  /// Get point set surface points
  void GetSurfacePoints(class PointSet &) const;

  // ---------------------------------------------------------------------------
  // Update

  /// Update (transformed) dataset
  ///
  /// This function only updates the output points if the self-update attribute
  /// is enabled and only if a (changing) transformation is set. If the dataset
  /// is not transformed or only mapped by a fixed transformation, this function
  /// does nothing. Use \c force=true to initialize this point set even if it
  /// does not change over the course of the registration.
  ///
  /// \param[in] force Force update in any case.
  void Update(bool force = false);

  // ---------------------------------------------------------------------------
  // Debugging

  /// Default file name extension
  const char *DefaultExtension() const;

  /// Write transformed dataset to file
  void Write(const char *, vtkAbstractArray * = NULL, vtkAbstractArray * = NULL) const;

  /// Write transformed dataset to file
  void Write(const char *, vtkAbstractArray **, int, vtkAbstractArray ** = NULL, int = 0) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Input point set
// =============================================================================

// -----------------------------------------------------------------------------
inline int RegisteredPointSet::NumberOfPoints() const
{
  return const_cast<vtkPointSet *>(_InputPointSet.GetPointer())->GetNumberOfPoints();
}

// -----------------------------------------------------------------------------
inline int RegisteredPointSet::NumberOfCells() const
{
  return const_cast<vtkPointSet *>(_InputPointSet.GetPointer())->GetNumberOfCells();
}

// -----------------------------------------------------------------------------
inline int RegisteredPointSet::NumberOfEdges() const
{
  return Edges()->NumberOfEdges();
}

// -----------------------------------------------------------------------------
inline void RegisteredPointSet::GetInputPoint(int i, double &x, double &y, double &z) const
{
  double p[3];
  _InputPointSet->GetPoints()->GetPoint(i, p);
  x = p[0], y = p[1], z = p[2];
}

// -----------------------------------------------------------------------------
inline void RegisteredPointSet::GetInputPoint(int i, double *p) const
{
  _InputPointSet->GetPoints()->GetPoint(i, p);
}

// -----------------------------------------------------------------------------
inline void RegisteredPointSet::GetInputPoint(int i, Point &pt) const
{
  double p[3];
  _InputPointSet->GetPoints()->GetPoint(i, p);
  pt._x = p[0], pt._y = p[1], pt._z = p[2];
}

// =============================================================================
// Input point set surface
// =============================================================================

// -----------------------------------------------------------------------------
inline bool RegisteredPointSet::IsSurface() const
{
  return _IsSurfaceMesh;
}

// -----------------------------------------------------------------------------
inline int RegisteredPointSet::NumberOfSurfacePoints() const
{
  return const_cast<vtkPolyData *>(_InputSurface.GetPointer())->GetNumberOfPoints();
}

// -----------------------------------------------------------------------------
inline int RegisteredPointSet::NumberOfSurfaceCells() const
{
  return const_cast<vtkPolyData *>(_InputSurface.GetPointer())->GetNumberOfCells();
}

// -----------------------------------------------------------------------------
inline int RegisteredPointSet::NumberOfSurfaceEdges() const
{
  return SurfaceEdges()->NumberOfEdges();
}

// -----------------------------------------------------------------------------
inline void RegisteredPointSet::GetInputSurfacePoint(int i, double &x, double &y, double &z) const
{
  double p[3];
  _InputSurface->GetPoints()->GetPoint(i, p);
  x = p[0], y = p[1], z = p[2];
}

// -----------------------------------------------------------------------------
inline void RegisteredPointSet::GetInputSurfacePoint(int i, double *p) const
{
  _InputSurface->GetPoints()->GetPoint(i, p);
}

// -----------------------------------------------------------------------------
inline void RegisteredPointSet::GetInputSurfacePoint(int i, Point &pt) const
{
  double p[3];
  _InputSurface->GetPoints()->GetPoint(i, p);
  pt._x = p[0], pt._y = p[1], pt._z = p[2];
}

// -----------------------------------------------------------------------------
inline const PointSet &RegisteredPointSet::InputSurfacePoints() const
{
  return *_InputSurfacePoints;
}

// =============================================================================
// Point set
// =============================================================================

// -----------------------------------------------------------------------------
inline RegisteredPointSet::operator vtkDataSet *() const
{
  return _OutputPointSet.GetPointer();
}

// -----------------------------------------------------------------------------
inline RegisteredPointSet::operator vtkPointSet *() const
{
  return _OutputPointSet.GetPointer();
}

// -----------------------------------------------------------------------------
inline vtkPointSet *RegisteredPointSet::PointSet() const
{
  return _OutputPointSet.GetPointer();
}

// -----------------------------------------------------------------------------
inline vtkPoints *RegisteredPointSet::Points() const
{
  return _OutputPointSet->GetPoints();
}

// -----------------------------------------------------------------------------
inline void RegisteredPointSet::GetPoint(int i, double &x, double &y, double &z) const
{
  double p[3];
  _OutputPointSet->GetPoints()->GetPoint(i, p);
  x = p[0], y = p[1], z = p[2];
}

// -----------------------------------------------------------------------------
inline void RegisteredPointSet::GetPoint(int i, double *p) const
{
  _OutputPointSet->GetPoints()->GetPoint(i, p);
}

// -----------------------------------------------------------------------------
inline void RegisteredPointSet::GetPoint(int i, Point &pt) const
{
  double p[3];
  _OutputPointSet->GetPoints()->GetPoint(i, p);
  pt._x = p[0], pt._y = p[1], pt._z = p[2];
}

// =============================================================================
// Point set surface
// =============================================================================

// -----------------------------------------------------------------------------
inline vtkPolyData *RegisteredPointSet::Surface() const
{
  return _OutputSurface.GetPointer();
}

// -----------------------------------------------------------------------------
inline vtkPoints *RegisteredPointSet::SurfacePoints() const
{
  return _OutputSurface->GetPoints();
}

// -----------------------------------------------------------------------------
inline void RegisteredPointSet::GetSurfacePoint(int i, double &x, double &y, double &z) const
{
  double p[3];
  _OutputSurface->GetPoints()->GetPoint(i, p);
  x = p[0], y = p[1], z = p[2];
}

// -----------------------------------------------------------------------------
inline void RegisteredPointSet::GetSurfacePoint(int i, double *p) const
{
  _OutputSurface->GetPoints()->GetPoint(i, p);
}

// -----------------------------------------------------------------------------
inline void RegisteredPointSet::GetSurfacePoint(int i, Point &pt) const
{
  double p[3];
  _OutputSurface->GetPoints()->GetPoint(i, p);
  pt._x = p[0], pt._y = p[1], pt._z = p[2];
}


} // namespace mirtk

#endif // MIRTK_RegisteredPointSet_H

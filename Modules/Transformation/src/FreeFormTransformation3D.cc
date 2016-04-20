/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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

#include "mirtk/FreeFormTransformation3D.h"

#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"

#ifdef HAVE_VTK
#  include "vtkSmartPointer.h"
#  include "vtkPoints.h"
#  include "vtkPolyData.h"
#  include "vtkIdTypeArray.h"
#  include "vtkOctreePointLocator.h"
#endif


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
ImageAttributes FreeFormTransformation3D
::DefaultAttributes(double x1, double y1, double z1,
                    double x2, double y2, double z2,
                    double dx, double dy, double dz,
                    const double *xaxis, const double *yaxis, const double *zaxis)
{
  return FreeFormTransformation::DefaultAttributes(
             x1, y1, z1, .0, x2, y2, z2, .0,
             dx, dy, dz, .0, xaxis, yaxis, zaxis
         );
}

// -----------------------------------------------------------------------------
FreeFormTransformation3D
::FreeFormTransformation3D(CPInterpolator &func, CPInterpolator *func2D)
:
  FreeFormTransformation(func, func2D)
{
}

// -----------------------------------------------------------------------------
FreeFormTransformation3D
::FreeFormTransformation3D(const FreeFormTransformation3D &ffd,
                               CPInterpolator &func, CPInterpolator *func2D)
:
  FreeFormTransformation(ffd, func, func2D)
{
}

// -----------------------------------------------------------------------------
FreeFormTransformation3D::~FreeFormTransformation3D()
{
}

// -----------------------------------------------------------------------------
void FreeFormTransformation3D::Initialize(const ImageAttributes &attr)
{
  if (attr._dt && attr._t != 1) {
    cerr << this->NameOfClass() << "::Initialize: Image attributes must represent a 3D image domain" << endl;
    exit(1);
  }
  FreeFormTransformation::Initialize(attr);
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
/// (Multi-threaded) body of FreeFormTransformation3D::ParametricGradient
///
/// This default implementation is more efficient than the general implementation
/// provided by FreeFormTransformation for FFDs which are parameterized by
/// displacements. It takes advantage of the fact that the derivatives w.r.t.
/// the control point parameters are non-zero only within a local support region
/// centered at the control point position. Therefore, this implementation first
/// iterates over all control points and in the inner loop only over those discrete
/// image positions (target voxels) which are within the support region of the
/// respective control point. Note that in case of FFDs parameterized by velocity
/// fields, each control point may still affect any point trajectory. The local
/// support region only applies to the velocity field, but not the generated
/// displacement field. Such transformations must therefore override the
/// ParametricGradient function and either use the more general base class
/// implementation given by FreeFormTransformation::ParametricGradient or
/// provide their own specialized implementation.
class FreeFormTransformation3DParametricGradientBody
{
public:
  const FreeFormTransformation3D *_FFD;
  const GenericImage<double>     *_Input;
  const WorldCoordsImage         *_Image2World;
  const WorldCoordsImage         *_WorldCoords;
  const Matrix                   *_World2Image;
  double                          _MarginX, _MarginY, _MarginZ;
  double                         *_Output;
  double                          _Weight;

  int _X, _Y, _Z, _N;

  // ---------------------------------------------------------------------------
  /// Calculates the gradient of the similarity term w.r.t. the transformation
  /// parameters for each voxel in the specified image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    int           cp, xdof, ydof, zdof; // indices of control point and corresponding DoFs
    const double *gx, *gy, *gz;         // pointers to input gradient data
    double        jac[3];               // B-spline basis
    FreeFormTransformation3D::CPStatus status;

    // With transformed world coordinates
    if (_WorldCoords) {
      const Matrix &w2i = *_World2Image; // approx. affine world to image space transformation
      double        x[2], y[2], z[2];        // bounding box of control point in world space
      double        x1, y1, z1, x2, y2, z2;  // bounding box of control point in image space
      int           i1, i2, j1, j2, k1, k2;  // bounding box of target voxels affected by control point
      double        px, py, pz;              // affine mapped bounding box corner
      const double *wx, *wy, *wz;            // transformed world coordinates
      int           s2, s3;                  // stride of data/increment of data pointers
      // Loop over control points
      for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
      for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
      for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
        // Get index of the control point
        cp = _FFD->LatticeToIndex(ci, cj, ck);
        // Check if any DoF corresponding to the control point is active
        _FFD->GetStatus(cp, status);
        if (status._x == Passive && status._y == Passive && status._z == Passive) continue;
        // Calculate bounding box of control point in world coordinates
        _FFD->BoundingBox(cp, x[0], y[0], z[0], x[1], y[1], z[1], 1.0 / _FFD->SpeedupFactor());
        // Map bounding box into (affine) image space
        x1 = y1 = z1 = + numeric_limits<double>::infinity();
        x2 = y2 = z2 = - numeric_limits<double>::infinity();
        for (int c = 0; c <= 1; ++c)
        for (int b = 0; b <= 1; ++b)
        for (int a = 0; a <= 1; ++a) {
          Transform(w2i, x[a], y[b], z[c], px, py, pz);
          if (px < x1) x1 = px;
          if (px > x2) x2 = px;
          if (py < y1) y1 = py;
          if (py > y2) y2 = py;
          if (pz < z1) z1 = pz;
          if (pz > z2) z2 = pz;
        }
        // Add approximation error margin and round to nearest voxels
        i1 = iround(x1 - _MarginX);
        i2 = iround(x2 + _MarginX);
        j1 = iround(y1 - _MarginY);
        j2 = iround(y2 + _MarginY);
        k1 = iround(z1 - _MarginZ);
        k2 = iround(z2 + _MarginZ);
        // When both indices are outside in opposite directions,
        // use the full range [0, N[. If they are both outside in
        // the same direction, the condition i1 <= i2 is false which
        // indicates that the bounding box is empty in this case
        i1 = (i1 < 0 ?  0 : (i1 >= _X ? _X     : i1));
        i2 = (i2 < 0 ? -1 : (i2 >= _X ? _X - 1 : i2));
        j1 = (j1 < 0 ?  0 : (j1 >= _Y ? _Y     : j1));
        j2 = (j2 < 0 ? -1 : (j2 >= _Y ? _Y - 1 : j2));
        k1 = (k1 < 0 ?  0 : (k1 >= _Z ? _Z     : k1));
        k2 = (k2 < 0 ? -1 : (k2 >= _Z ? _Z - 1 : k2));
        if (i2 < i1 || j2 < j1 || k2 < k1) continue;
        // Get indices of DoFs corresponding to the control point
        _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
        // Loop over target voxels
        s2 =  _X - (i2 - i1 + 1);
        s3 = (_Y - (j2 - j1 + 1)) * _X;
        wx = _WorldCoords->Data(i1, j1, k1), wy = wx + _N, wz = wy + _N;
        gx = _Input      ->Data(i1, j1, k1), gy = gx + _N, gz = gy + _N;
        for (int k = k1; k <= k2; ++k, wx += s3, wy += s3, wz += s3, gx += s3, gy += s3, gz += s3)
        for (int j = j1; j <= j2; ++j, wx += s2, wy += s2, wz += s2, gx += s2, gy += s2, gz += s2)
        for (int i = i1; i <= i2; ++i, wx +=  1, wy +=  1, wz +=  1, gx +=  1, gy +=  1, gz +=  1) {
          // Check whether reference point is valid
          if (*gx == .0 && *gy == .0 && *gz == .0) continue;
          // Check if world coordinate is in support region of control point
          if (x[0] <= *wx && *wx <= x[1] &&
              y[0] <= *wy && *wy <= y[1] &&
              z[0] <= *wz && *wz <= z[1]) {
            // Convert non-parametric gradient into parametric gradient
            _FFD->JacobianDOFs(jac, ci, cj, ck, *wx, *wy, *wz);
            if (status._x == Active) _Output[xdof] += _Weight * jac[0] * (*gx);
            if (status._y == Active) _Output[ydof] += _Weight * jac[1] * (*gy);
            if (status._z == Active) _Output[zdof] += _Weight * jac[2] * (*gz);
          }
        }
      }
    }
    // With pre-computed voxel coordinates
    else if (_Image2World) {
      const double *wx, *wy, *wz;           // pre-computed voxel coordinates
      int           i1, i2, j1, j2, k1, k2; // bounding box of target voxels affected by control point
      int           s2, s3;                 // stride of data/increment of data pointers
      // Loop over control points
      for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
      for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
      for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
        // Get index of the control point
        cp = _FFD->LatticeToIndex(ci, cj, ck);
        // Check if any DoF corresponding to the control point is active
        _FFD->GetStatus(cp, status);
        if (status._x == Passive && status._y == Passive && status._z == Passive) continue;
        // Calculate bounding box of control point in image coordinates
        if (!_FFD->BoundingBox(_Input, cp, i1, j1, k1, i2, j2, k2, 1.0 / _FFD->SpeedupFactor())) continue;
        // Get indices of DoFs corresponding to the control point
        _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
        // Loop over all voxels in the target (reference) volume
        s2 =  _X - (i2 - i1 + 1);
        s3 = (_Y - (j2 - j1 + 1)) * _X;
        wx = _Image2World->Data(i1, j1, k1), wy = wx + _N, wz = wy + _N;
        gx = _Input      ->Data(i1, j1, k1), gy = gx + _N, gz = gy + _N;
        for (int k = k1; k <= k2; ++k, wx += s3, wy += s3, wz += s3, gx += s3, gy += s3, gz += s3)
        for (int j = j1; j <= j2; ++j, wx += s2, wy += s2, wz += s2, gx += s2, gy += s2, gz += s2)
        for (int i = i1; i <= i2; ++i, wx +=  1, wy +=  1, wz +=  1, gx +=  1, gy +=  1, gz +=  1) {
          // Check whether reference point is valid
          if (*gx == .0 && *gy == .0 && *gz == .0) continue;
          // Convert non-parametric gradient into parametric gradient
          _FFD->JacobianDOFs(jac, ci, cj, ck, *wx, *wy, *wz);
          if (status._x == Active) _Output[xdof] += _Weight * jac[0] * (*gx);
          if (status._y == Active) _Output[ydof] += _Weight * jac[1] * (*gy);
          if (status._z == Active) _Output[zdof] += _Weight * jac[2] * (*gz);
        }
      }
    }
    // Without pre-computed world coordinates
    else {
      int    i1, i2, j1, j2, k1, k2; // bounding box of target voxels affected by control point
      int    s2, s3;                 // stride of data/increment of data pointers
      double x,  y,  z;              // world coordinates of target voxel
      // Loop over control points
      for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
      for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
      for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
        // Get index of the control point
        cp = _FFD->LatticeToIndex(ci, cj, ck);
        // Check if any DoF corresponding to the control point is active
        _FFD->GetStatus(cp, status);
        if (status._x == Passive && status._y == Passive && status._z == Passive) continue;
        // Calculate bounding box of control point in image coordinates
        if (!_FFD->BoundingBox(_Input, cp, i1, j1, k1, i2, j2, k2, 1.0 / _FFD->SpeedupFactor())) continue;
        // Get indices of DoFs corresponding to the control point
        _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
        // Loop over all voxels in the target (reference) volume
        s2 =  _X - (i2 - i1 + 1);
        s3 = (_Y - (j2 - j1 + 1)) * _X;
        gx = _Input->Data(i1, j1, k1), gy = gx + _N, gz = gy + _N;
        for (int k = k1; k <= k2; ++k, gx += s3, gy += s3, gz += s3)
        for (int j = j1; j <= j2; ++j, gx += s2, gy += s2, gz += s2)
        for (int i = i1; i <= i2; ++i, gx +=  1, gy +=  1, gz +=  1) {
          // Check whether reference point is valid
          if (*gx == .0 && *gy == .0 && *gz == .0) continue;
          // Convert voxel coordinates to world coordinates
          x = i, y = j, z = k;
          _Input->ImageToWorld(x, y, z);
          _FFD->JacobianDOFs(jac, ci, cj, ck, x, y, z);
          if (status._x == Active) _Output[xdof] += _Weight * jac[0] * (*gx);
          if (status._y == Active) _Output[ydof] += _Weight * jac[1] * (*gy);
          if (status._z == Active) _Output[zdof] += _Weight * jac[2] * (*gz);
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  void operator ()()
  {
    // Domain of voxel-wise gradient
    _X = _Input->X();
    _Y = _Input->Y();
    _Z = _Input->Z();
    _N = _X * _Y * _Z;

    // Check input
    if (_Image2World && (_Image2World->X() != _X ||
                         _Image2World->Y() != _Y ||
                         _Image2World->Z() != _Z ||
                         _Image2World->T() !=  3)) {
      cerr << "FreeFormTransformation3D::ParametricGradient: Invalid voxel coordinates map" << endl;
      exit(1);
    }
    if (_WorldCoords && (_WorldCoords->X() != _X ||
                         _WorldCoords->Y() != _Y ||
                         _WorldCoords->Z() != _Z ||
                         _WorldCoords->T() !=  3)) {
      cerr << "FreeFormTransformation3D::ParametricGradient: Invalid world coordinates map" << endl;
      exit(1);
    }

    if (_Input->GetXSize() > _FFD->GetXSpacing() ||
        _Input->GetYSize() > _FFD->GetYSpacing() ||
        (_FFD->Z() > 1 && _Input->GetZSize() > _FFD->GetZSpacing())) {
      // FIXME: In this case, the non-parametric input gradient should be
      //        resampled using linear interpolation to obtain a value at
      //        each control point.
      cerr << "Warning: FFD spacing smaller than image resolution!" << endl;
      cerr << "         This may lead to artifacts in the transformation because" << endl;
      cerr << "         not all control points are within the vicinity of a voxel center." << endl;
    }

    Matrix w2i;
    _World2Image = &w2i;
    _MarginX = _MarginY = _MarginZ = .0;

    if (_WorldCoords) {
      double x, y, z;
      int    i, j, k;

      // Approximate world to image map by affine transformation
      PointSet voxel(8), world(voxel.Size());
      int n = 0;
      for (int c = 0; c <= 1; ++c)
      for (int b = 0; b <= 1; ++b)
      for (int a = 0; a <= 1; ++a, ++n) {
        i = a * (_X - 1), j = b * (_Y - 1), k = c * (_Z - 1);
        voxel(n) = Point(i, j, k);
        world(n)._x = _WorldCoords->Get(i, j, k, 0);
        world(n)._y = _WorldCoords->Get(i, j, k, 1);
        world(n)._z = _WorldCoords->Get(i, j, k, 2);
      }
      w2i = ApproximateAffineMatrix(world, voxel);

      // Determine maximum absolute error of approximation
      double i2, j2, k2, dx, dy, dz;
      for (k = 0; k < _Z; ++k)
      for (j = 0; j < _Y; ++j)
      for (i = 0; i < _X; ++i) {
        x = _WorldCoords->Get(i, j, k, 0);
        y = _WorldCoords->Get(i, j, k, 1);
        z = _WorldCoords->Get(i, j, k, 2);
        Transform(w2i, x, y, z, i2, j2, k2);
        dx = abs(i2 - i);
        dy = abs(j2 - j);
        dz = abs(k2 - k);
        if (dx > _MarginX) _MarginX = dx;
        if (dy > _MarginY) _MarginY = dy;
        if (dz > _MarginZ) _MarginZ = dz;
      }
    }

    // Calculate parametric gradient
    blocked_range3d<int> cps(0, _FFD->Z(), 0, _FFD->Y(), 0, _FFD->X());
    parallel_for(cps, *this);
  }

}; // FreeFormTransformation3DParametricGradientBody

#ifdef HAVE_VTK

// -----------------------------------------------------------------------------
class FreeFormTransformation3DPointWiseParametricGradientBody
{
public:

  const FreeFormTransformation3D *_Transformation;
  const PointSet                 *_Point;
  const Vector3D<double>         *_Input;
  double                         *_Output;
  double                          _Weight;

  double _t;  ///< Time corrresponding to input gradient image (in ms)
  double _t0; ///< Second time argument for velocity-based transformations

  // Attention: vtkKdTree/vtkKdTreePointLocator is not thread-safe
  //            (cf. http://www.vtk.org/Bug/view.php?id=15206 )!
  vtkOctreePointLocator *_Locator; ///< Point locator used to speed up computation

  // ---------------------------------------------------------------------------
  /// Default constructor
  FreeFormTransformation3DPointWiseParametricGradientBody()
  :
    _Transformation(NULL),
    _Point         (NULL),
    _Input         (NULL),
    _Output        (NULL),
    _Weight        ( 1.0),
    _t             ( 0.0),
    _t0            (-1.0),
    _Locator       (NULL)
  {}

  // ---------------------------------------------------------------------------
  /// Copy constructor
  FreeFormTransformation3DPointWiseParametricGradientBody(
    const FreeFormTransformation3DPointWiseParametricGradientBody &other
  ) :
    _Transformation(other._Transformation),
    _Point         (other._Point),
    _Input         (other._Input),
    _Output        (other._Output),
    _Weight        (other._Weight),
    _t             (other._t),
    _t0            (other._t0),
    _Locator       (other._Locator)
  {}

  // ---------------------------------------------------------------------------
  /// Destructor
  ~FreeFormTransformation3DPointWiseParametricGradientBody() {}

  // ---------------------------------------------------------------------------
  /// Calculates the gradient of the similarity term w.r.t. the transformation
  /// parameters for each point in 3D.
  void operator ()(const blocked_range<int> &re) const
  {
    double     bounds[6];        // bounding box of control point support
    int        xdof, ydof, zdof; // indices of control point DoFs
    int        ci, cj, ck;       // indices of control point
    double     jac[3];           // Jacobian of transformation w.r.t control point
    FreeFormTransformation3D::CPStatus status;

    // Loop over CPs
    vtkSmartPointer<vtkIdTypeArray> ids = vtkSmartPointer<vtkIdTypeArray>::New();
    for (int cp = re.begin(); cp != re.end(); ++cp) {
      // Check status of DoFs
      _Transformation->GetStatus(cp, status);
      if (status._x == Passive && status._y == Passive && status._z == Passive) continue;
      // Get indices of control point parameters (DoFs)
      _Transformation->IndexToDOFs(cp, xdof, ydof, zdof);
      _Transformation->IndexToLattice(cp, ci, cj, ck);
      // Find points within kernel support region
      _Transformation->BoundingBox(cp, bounds[0], bounds[2], bounds[4],
                                       bounds[1], bounds[3], bounds[5],
                                   1.0 / _Transformation->SpeedupFactor());
      _Locator->FindPointsInArea(bounds, ids);
      // Loop over points within support region
      for (vtkIdType i = 0; i < ids->GetNumberOfTuples(); ++i) {
        const int idx = static_cast<int>(ids->GetTuple1(i));
        const Vector3D<double> &g = _Input[idx];
        // Check whether reference point is valid
        if (g._x != .0 || g._y != .0 || g._z != .0) {
          // Convert non-parametric gradient into parametric gradient
          const Point &p = (*_Point)(idx);
          _Transformation->JacobianDOFs(jac, ci, cj, ck, p._x, p._y, p._z);
          if (status._x == Active) _Output[xdof] += _Weight * jac[0] * g._x;
          if (status._y == Active) _Output[ydof] += _Weight * jac[1] * g._y;
          if (status._z == Active) _Output[zdof] += _Weight * jac[2] * g._z;
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  void operator ()()
  {
    if (_Point->Size() == 0) {
      return;
    }

    vtkSmartPointer<vtkPoints> points;
    points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(_Point->Size());
    for (int idx = 0; idx < _Point->Size(); ++idx) {
      const Point &p = (*_Point)(idx);
      points->SetPoint(idx, p._x, p._y, p._z);
    }

    vtkSmartPointer<vtkPolyData> dataset;
    dataset = vtkSmartPointer<vtkPolyData>::New();
    dataset->SetPoints(points);

    vtkSmartPointer<vtkOctreePointLocator> locator;
    locator = vtkSmartPointer<vtkOctreePointLocator>::New();
    locator->SetDataSet(dataset);
    locator->BuildLocator();

    _Locator = locator;
    blocked_range<int> cps(0, _Transformation->NumberOfCPs());
    parallel_for(cps, *this);
  }

}; // FreeFormTransformation3DPointWiseParametricGradientBody

#endif // HAVE_VTK

// -----------------------------------------------------------------------------
void FreeFormTransformation3D
::ParametricGradient(const GenericImage<double> *in, double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     double, double w) const
{
  MIRTK_START_TIMING();
  FreeFormTransformation3DParametricGradientBody body;
  body._FFD         = this;
  body._Input       = in;
  body._Output      = out;
  body._Weight      = w;
  body._Image2World = i2w;
  body._WorldCoords = wc;
  body();
  MIRTK_DEBUG_TIMING(2, "parametric gradient computation (3D FFD)");
}

// -----------------------------------------------------------------------------
void FreeFormTransformation3D
::ParametricGradient(const PointSet &pos, const Vector3D<double> *in,
                     double *out, double t, double t0, double w) const
{
#ifdef HAVE_VTK
  MIRTK_START_TIMING();
  FreeFormTransformation3DPointWiseParametricGradientBody body;
  body._Transformation = this;
  body._Point          = &pos;
  body._Input          = in;
  body._Output         = out;
  body._Weight         = w;
  body._t              = t;
  body._t0             = t0;
  body();
  MIRTK_DEBUG_TIMING(2, "point-wise parametric gradient computation (3D FFD)");
#else
  FreeFormTransformation::ParametricGradient(pos, in, out, t, t0, w);
#endif
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
Cifstream &FreeFormTransformation3D::ReadDOFs(Cifstream &from, TransformationType format)
{
  ImageAttributes attr;

  // Read no of control points
  from.ReadAsInt(&attr._x, 1);
  from.ReadAsInt(&attr._y, 1);
  from.ReadAsInt(&attr._z, 1);
  
  // Read orientation of bounding box
  from.ReadAsDouble(attr._xaxis, 3);
  from.ReadAsDouble(attr._yaxis, 3);

  if (format == TRANSFORMATION_BSPLINE_FFD_3D_v1 || format == TRANSFORMATION_LINEAR_FFD_3D_v1) {
    // Deprecated older formats which did not store full orientation matrix
    attr._zaxis[0] = attr._xaxis[1] * attr._yaxis[2] - attr._xaxis[2] * attr._yaxis[1];
    attr._zaxis[1] = attr._xaxis[2] * attr._yaxis[0] - attr._xaxis[0] * attr._yaxis[2];
    attr._zaxis[2] = attr._xaxis[0] * attr._yaxis[1] - attr._xaxis[1] * attr._yaxis[0];
  } else {
    from.ReadAsDouble(attr._zaxis, 3);
  }

  // Read spacing of bounding box
  from.ReadAsDouble(&attr._dx, 1);
  from.ReadAsDouble(&attr._dy, 1);
  from.ReadAsDouble(&attr._dz, 1);

  // Read origin of bounding box
  from.ReadAsDouble(&attr._xorigin, 1);
  from.ReadAsDouble(&attr._yorigin, 1);
  from.ReadAsDouble(&attr._zorigin, 1);

  // Read extrapolation mode
  if (static_cast<int>(format) >= 70) {
    unsigned int mode;
    from.ReadAsUInt(&mode, 1);
    if (mode >= Extrapolation_Last) {
      cerr << this->NameOfClass() << "::ReadDOFs: Invalid extrapolation mode: " << mode << endl;
      cerr << "  Using default extrapolation mode instead." << endl;
    } else {
      this->ExtrapolationMode(static_cast<enum ExtrapolationMode>(mode));
    }
  }

  // Initialize free-form transformation
  this->Initialize(attr);

  // Read control point data
  if (static_cast<int>(format) < 50) {
    // Older transformations stored control points in different order, i.e.,
    // inner loop over z and outer loop over x instead of the other way round.
    DOFValue *param = new DOFValue[_NumberOfDOFs];
    from.ReadAsDouble(param, _NumberOfDOFs);
    int dof = 0;
    for (int i = 0; i < attr._x; ++i)
    for (int j = 0; j < attr._y; ++j)
    for (int k = 0; k < attr._z; ++k) {
      _CPImage(i, j, k) = CPValue(param[dof], param[dof + 1], param[dof + 2]);
      dof += 3;
    }
  } else {
    from.ReadAsDouble(_Param, _NumberOfDOFs);
  }

  // Read control point status
  if (static_cast<int>(format) < 50) {
    // Older transformations stored status as "xxx...yyy...zzz..."
    // even though the parameters were stored as "xyzxyzxyz..."
    DOFStatus *status = new DOFStatus[_NumberOfDOFs];
    from.ReadAsInt((int *)status, _NumberOfDOFs);
    for (int dof = 0, cp = 0; dof < _NumberOfDOFs; dof += 3, cp += 1) {
      _Status[dof    ] = status[cp                    ];
      _Status[dof + 1] = status[cp +     NumberOfCPs()];
      _Status[dof + 2] = status[cp + 2 * NumberOfCPs()];
    }
    delete[] status;
  } else {
    from.ReadAsInt((int *)_Status, _NumberOfDOFs);
  }

  return from;
}

// -----------------------------------------------------------------------------
Cofstream &FreeFormTransformation3D::WriteDOFs(Cofstream &to) const
{
  // Write no of control points
  to.WriteAsInt(&_attr._x, 1);
  to.WriteAsInt(&_attr._y, 1);
  to.WriteAsInt(&_attr._z, 1);
  
  // Write orientation of bounding box
  to.WriteAsDouble(_attr._xaxis, 3);
  to.WriteAsDouble(_attr._yaxis, 3);
  to.WriteAsDouble(_attr._zaxis, 3);
  
  // Write spacing of bounding box
  to.WriteAsDouble(&_attr._dx, 1);
  to.WriteAsDouble(&_attr._dy, 1);
  to.WriteAsDouble(&_attr._dz, 1);
  
  // Write origin of bounding box
  to.WriteAsDouble(&_attr._xorigin, 1);
  to.WriteAsDouble(&_attr._yorigin, 1);
  to.WriteAsDouble(&_attr._zorigin, 1);

  // Write extrapolation mode
  const unsigned int emode = _ExtrapolationMode;
  to.WriteAsUInt(&emode, 1);

  return WriteCPs(to);
}


} // namespace mirtk

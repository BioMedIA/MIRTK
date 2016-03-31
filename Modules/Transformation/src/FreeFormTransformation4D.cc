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

#include "mirtk/FreeFormTransformation4D.h"

#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
FreeFormTransformation4D::FreeFormTransformation4D(CPInterpolator &func)
:
  FreeFormTransformation(func)
{
}

// -----------------------------------------------------------------------------
FreeFormTransformation4D
::FreeFormTransformation4D(const FreeFormTransformation4D &ffd,
                           CPInterpolator                 &func)
:
  FreeFormTransformation(ffd, func)
{
}

// -----------------------------------------------------------------------------
FreeFormTransformation4D::~FreeFormTransformation4D()
{
}

// -----------------------------------------------------------------------------
void FreeFormTransformation4D::Initialize(const ImageAttributes &attr)
{
  if (attr._dt == 0) {
    cerr << this->NameOfClass() << "::Initialize: Image attributes must represent a 4D image domain" << endl;
    exit(1);
  }
  FreeFormTransformation::Initialize(attr);
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
/// (Multi-threaded) body of FreeFormTransformation4D::ParametricGradient
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
class FreeFormTransformation4DParametricGradientBody
{
public:
  const FreeFormTransformation4D *_FFD;
  const GenericImage<double>     *_Input;
  const WorldCoordsImage         *_Image2World;
  const WorldCoordsImage         *_WorldCoords;
  double                         *_Output;
  double                          _Weight;

  int _X, _Y, _N;

  double _t;  ///< Time corrresponding to input gradient image (in ms)
  double _t0; ///< Second time argument for velocity-based transformations.
  int    _cl; ///< Current temporal control point index

  // ---------------------------------------------------------------------------
  void operator ()(const blocked_range3d<int> &re) const
  {
    double        t1, t2;               // temporal bounding box (in ms)
    int           cp, xdof, ydof, zdof; // indices of DoFs corresponding to control point
    const double *gx, *gy, *gz;         // pointers to input gradient data
    Matrix    jac;                  // Jacobian of transformation w.r.t control point
    FreeFormTransformation4D::CPStatus status;

    // With (transformed) world coordinates
    if (_WorldCoords) {
      const double *wx, *wy, *wz;
      double        x1, y1, z1, x2, y2, z2;
      // Loop over control points
      for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
      for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
      for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
        // Compute DoFs corresponding to the control point
        cp = _FFD->LatticeToIndex(ci, cj, ck, _cl);
        // Check if any DoF corresponding to the control point is active
        _FFD->GetStatus(cp, status);
        if (status._x == Passive && status._y == Passive && status._z == Passive) continue;
        // Calculate temporal bounding box and check if spatial frame is in it
        _FFD->BoundingBox(cp, t1, t2);
        if (t1 > _t || t2 < _t) continue;
        // Calculate spatial bounding box of control point in world coordinates
        _FFD->BoundingBox(cp, x1, y1, z1, x2, y2, z2, 1.0 / _FFD->SpeedupFactor());
        // Get indices of DoFs corresponding to the control point
        _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
        // Loop over target voxels
        wx = _WorldCoords->Data(), wy = wx + _N, wz = wy + _N;
        gx = _Input      ->Data(), gy = gx + _N, gz = gy + _N;
        for (int i = 0; i < _N; ++i, ++wx, ++wy, ++wz, ++gx, ++gy, ++gz) {
          // Check whether reference point is valid
          if (*gx == .0 && *gy == .0 && *gz == .0) continue;
          // Check if world coordinate is in support region of control point
          if (x1 <= *wx && *wx <= x2 &&
              y1 <= *wy && *wy <= y2 &&
              z1 <= *wz && *wz <= z2) {
            // Convert non-parametric gradient into parametric gradient
            _FFD->JacobianDOFs(jac, ci, cj, ck, _cl, *wx, *wy, *wz, _t, _t0);
            if (status._x == Active) _Output[xdof] += _Weight * (jac(0, 0) * (*gx) + jac(1, 0) * (*gy) + jac(2, 0) * (*gz));
            if (status._y == Active) _Output[ydof] += _Weight * (jac(0, 1) * (*gx) + jac(1, 1) * (*gy) + jac(2, 1) * (*gz));
            if (status._z == Active) _Output[zdof] += _Weight * (jac(0, 2) * (*gx) + jac(1, 2) * (*gy) + jac(2, 2) * (*gz));
          }
        }
      }
    }
    // With pre-computed voxel coordinates
    else if (_Image2World) {
      int           i1, i2, j1, j2, k1, k2; // spatial bounding box (in voxels)
      int           s2, s3;                 // stride of data/increment of data pointers
      const double *wx, *wy, *wz;           // pre-computed voxel coordinates
      // Loop over control points
      for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
      for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
      for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
        // Compute DoFs corresponding to the control point
        cp = _FFD->LatticeToIndex(ci, cj, ck, _cl);
        // Check if any DoF corresponding to the control point is active
        _FFD->GetStatus(cp, status);
        if (status._x == Passive && status._y == Passive && status._z == Passive) continue;
        // Calculate temporal bounding box and check if spatial frame is in it
        _FFD->BoundingBox(cp, t1, t2);
        if (t1 > _t || t2 < _t) continue;
        // Calculate spatial bounding box of control point in image coordinates
        if (!_FFD->BoundingBox(_Input, cp, i1, j1, k1, i2, j2, k2, 1.0 / _FFD->SpeedupFactor())) continue;
        _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
        // Loop over all voxels in the target (reference) volume
        //s1=1
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
          _FFD->JacobianDOFs(jac, ci, cj, ck, _cl, *wx, *wy, *wz, _t, _t0);
          if (status._x == Active) _Output[xdof] += _Weight * (jac(0, 0) * (*gx) + jac(1, 0) * (*gy) + jac(2, 0) * (*gz));
          if (status._y == Active) _Output[ydof] += _Weight * (jac(0, 1) * (*gx) + jac(1, 1) * (*gy) + jac(2, 1) * (*gz));
          if (status._z == Active) _Output[zdof] += _Weight * (jac(0, 2) * (*gx) + jac(1, 2) * (*gy) + jac(2, 2) * (*gz));
        }
      }
    }
    // Without pre-computed voxel coordinates
    else {
      int    i1, i2, j1, j2, k1, k2; // spatial bounding box (in voxels)
      int    s2, s3;                 // stride of data/increment of data pointers
      double x,  y,  z;              // world coordinates of target voxel
      // Loop over control points
      for (int ck = re.pages().begin(); ck != re.pages().end(); ++ck)
      for (int cj = re.rows ().begin(); cj != re.rows ().end(); ++cj)
      for (int ci = re.cols ().begin(); ci != re.cols ().end(); ++ci) {
        // Compute DoFs corresponding to the control point
        cp = _FFD->LatticeToIndex(ci, cj, ck, _cl);
        // Check if any DoF corresponding to the control point is active
        _FFD->GetStatus(cp, status);
        if (status._x == Passive && status._y == Passive && status._z == Passive) continue;
        // Calculate temporal bounding box and check if spatial frame is in it
        _FFD->BoundingBox(cp, t1, t2);
        if (t1 > _t || t2 < _t) continue;
        // Calculate spatial bounding box of control point in image coordinates
        if (!_FFD->BoundingBox(_Input, cp, i1, j1, k1, i2, j2, k2, 1.0 / _FFD->SpeedupFactor())) continue;
        _FFD->IndexToDOFs(cp, xdof, ydof, zdof);
        // Loop over all voxels in the target (reference) volume
        //s1=1
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
          _FFD->JacobianDOFs(jac, ci, cj, ck, _cl, x, y, z, _t, _t0);
          if (status._x == Active) _Output[xdof] += _Weight * (jac(0, 0) * (*gx) + jac(1, 0) * (*gy) + jac(2, 0) * (*gz));
          if (status._y == Active) _Output[ydof] += _Weight * (jac(0, 1) * (*gx) + jac(1, 1) * (*gy) + jac(2, 1) * (*gz));
          if (status._z == Active) _Output[zdof] += _Weight * (jac(0, 2) * (*gx) + jac(1, 2) * (*gy) + jac(2, 2) * (*gz));
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  void operator ()()
  {
    // Domain of per voxel gradient
    _X           = _Input->X();
    _Y           = _Input->Y();
    const int _Z = _Input->Z();
    _N           = _X * _Y * _Z;

    // Check input
    if (_Image2World && (_Image2World->X() != _X ||
                         _Image2World->Y() != _Y ||
                         _Image2World->Z() != _Z ||
                         _Image2World->T() !=  3)) {
      cerr << "FreeFormTransformation4D::ParametricGradient: Invalid voxel coordinates map" << endl;
      exit(1);
    }
    if (_WorldCoords && (_WorldCoords->X() != _X ||
                         _WorldCoords->Y() != _Y ||
                         _WorldCoords->Z() != _Z ||
                         _WorldCoords->T() !=  3)) {
      cerr << "FreeFormTransformation4D::ParametricGradient: Invalid world coordinates map" << endl;
      exit(1);
    }

    if (_Input->GetXSize() > _FFD->GetXSpacing() ||
        _Input->GetYSize() > _FFD->GetYSpacing() ||
        _Input->GetZSize() > _FFD->GetZSpacing()) {
      // FIXME: In this case, the non-parametric input gradient should be
      //        resampled using linear interpolation to obtain a value at
      //        each control point.
      cerr << "Warning: FFD spacing smaller than image resolution!" << endl;
      cerr << "         This may lead to artifacts in the transformation because" << endl;
      cerr << "         not all control points are within the vicinity of a voxel center." << endl;
    }

    // Calculate parametric gradient
    for (_cl = 0; _cl < _FFD->T(); _cl++) {
      blocked_range3d<int> cps(0, _FFD->Z(), 0, _FFD->Y(), 0, _FFD->X());
      parallel_for(cps, FreeFormTransformation4DParametricGradientBody(*this));
    }
  }

}; // FreeFormTransformation4DParametricGradientBody

// -----------------------------------------------------------------------------
void FreeFormTransformation4D
::ParametricGradient(const GenericImage<double> *in, double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     double t0, double w) const
{
  MIRTK_START_TIMING();
  FreeFormTransformation4DParametricGradientBody body;
  body._FFD         = this;
  body._Input       = in;
  body._Output      = out;
  body._Weight      = w;
  body._Image2World = i2w;
  body._WorldCoords = wc;
  body._t           = in->ImageToTime(.0);
  body._t0          = t0;
  body();
  MIRTK_DEBUG_TIMING(2, "parametric gradient computation (4D FFD)");
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
Cifstream &FreeFormTransformation4D::ReadDOFs(Cifstream &from, TransformationType format)
{
  ImageAttributes attr;

  // Read no of control points
  from.ReadAsInt(&attr._x, 1);
  from.ReadAsInt(&attr._y, 1);
  from.ReadAsInt(&attr._z, 1);
  from.ReadAsInt(&attr._t, 1);

  // Read orientation of bounding box
  from.ReadAsDouble(attr._xaxis, 3);
  from.ReadAsDouble(attr._yaxis, 3);
  from.ReadAsDouble(attr._zaxis, 3);

  // Read spacing of bounding box
  from.ReadAsDouble(&attr._dx, 1);
  from.ReadAsDouble(&attr._dy, 1);
  from.ReadAsDouble(&attr._dz, 1);
  from.ReadAsDouble(&attr._dt, 1);

  // Read origin of bounding box
  from.ReadAsDouble(&attr._xorigin, 1);
  from.ReadAsDouble(&attr._yorigin, 1);
  from.ReadAsDouble(&attr._zorigin, 1);
  from.ReadAsDouble(&attr._torigin, 1);

  double tmax;
  from.ReadAsDouble(&tmax, 1);

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
  if (static_cast<int>(format) <= 23) {
    // Older transformations stored control points in different order, i.e.,
    // inner loop over t and outer loop over x instead of the other way round.
    DOFValue *param = new DOFValue[_NumberOfDOFs];
    from.ReadAsDouble(param, _NumberOfDOFs);
    int dof = 0;
    for (int i = 0; i < attr._x; ++i)
    for (int j = 0; j < attr._y; ++j)
    for (int k = 0; k < attr._z; ++k)
    for (int l = 0; l < attr._t; ++l, dof += 3) {
      _CPImage(i, j, k, l) = CPValue(param[dof], param[dof + 1], param[dof + 2]);
    }
  } else {
    from.ReadAsDouble(_Param, _NumberOfDOFs);
  }

  // Read control point status
  if (static_cast<int>(format) <= 23) {
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
Cofstream &FreeFormTransformation4D::WriteDOFs(Cofstream &to) const
{
  // Write no of control points
  to.WriteAsInt(&_attr._x, 1);
  to.WriteAsInt(&_attr._y, 1);
  to.WriteAsInt(&_attr._z, 1);
  to.WriteAsInt(&_attr._t, 1);

  // Write orientation of bounding box
  to.WriteAsDouble(_attr._xaxis, 3);
  to.WriteAsDouble(_attr._yaxis, 3);
  to.WriteAsDouble(_attr._zaxis, 3);

  // Write spacing of bounding box
  to.WriteAsDouble(&_attr._dx, 1);
  to.WriteAsDouble(&_attr._dy, 1);
  to.WriteAsDouble(&_attr._dz, 1);
  to.WriteAsDouble(&_attr._dt, 1);

  // Write origin of bounding box
  to.WriteAsDouble(&_attr._xorigin, 1);
  to.WriteAsDouble(&_attr._yorigin, 1);
  to.WriteAsDouble(&_attr._zorigin, 1);

  // Write time domain
  double torigin = this->LatticeToTime(0);
  double tmax    = this->LatticeToTime(_attr._t - 1);
  to.WriteAsDouble(&torigin, 1);
  to.WriteAsDouble(&tmax,    1);

  // Write extrapolation mode
  const unsigned int emode = _ExtrapolationMode;
  to.WriteAsUInt(&emode, 1);

  return WriteCPs(to);
}


} // namespace mirtk

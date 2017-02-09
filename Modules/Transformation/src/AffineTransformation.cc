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

#include "mirtk/AffineTransformation.h"

#include "mirtk/Math.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void AffineTransformation::UpdateShearingTangent()
{
  _tansxy = tan(_Param[SXY] * rad_per_deg);
  _tansxz = tan(_Param[SXZ] * rad_per_deg);
  _tansyz = tan(_Param[SYZ] * rad_per_deg);
}

// -----------------------------------------------------------------------------
AffineTransformation::AffineTransformation(int ndofs)
:
  SimilarityTransformation(ndofs)
{
  _Param[SY] = _Param[SZ] = _Param[SX]; // SX set by SimilarityTransformation
  UpdateShearingTangent();
}

// -----------------------------------------------------------------------------
AffineTransformation::AffineTransformation(const RigidTransformation &t, int ndofs)
:
  SimilarityTransformation(t, ndofs)
{
  _Param[SY] = _Param[SZ] = _Param[SX]; // SX set by SimilarityTransformation
  UpdateShearingTangent();
}

// -----------------------------------------------------------------------------
AffineTransformation::AffineTransformation(const SimilarityTransformation &t, int ndofs)
:
  SimilarityTransformation(t, ndofs)
{
  _Param[SY] = _Param[SZ] = _Param[SX]; // SX set by SimilarityTransformation
  UpdateShearingTangent();
}

// -----------------------------------------------------------------------------
AffineTransformation::AffineTransformation(const AffineTransformation &t, int ndofs)
:
  SimilarityTransformation(t, ndofs),
  _tansxy(t._tansxy),
  _tansxz(t._tansxz),
  _tansyz(t._tansyz)
{
}

// -----------------------------------------------------------------------------
AffineTransformation::AffineTransformation()
:
  SimilarityTransformation(12)
{
  _Param[SY] = _Param[SZ] = _Param[SX]; // SX set by SimilarityTransformation
  UpdateShearingTangent();
}

// -----------------------------------------------------------------------------
AffineTransformation::AffineTransformation(const RigidTransformation &t)
:
  SimilarityTransformation(t, 12)
{
  _Param[SY] = _Param[SZ] = _Param[SX]; // SX set by SimilarityTransformation
  UpdateShearingTangent();
}

// -----------------------------------------------------------------------------
AffineTransformation::AffineTransformation(const SimilarityTransformation &t)
:
  SimilarityTransformation(t, 12)
{
  _Param[SY] = _Param[SZ] = _Param[SX]; // SX set by SimilarityTransformation
  UpdateShearingTangent();
}

// -----------------------------------------------------------------------------
AffineTransformation::AffineTransformation(const AffineTransformation &t)
:
  SimilarityTransformation(t, 12),
  _tansxy(t._tansxy),
  _tansxz(t._tansxz),
  _tansyz(t._tansyz)
{
}

// -----------------------------------------------------------------------------
AffineTransformation::~AffineTransformation()
{
}

// =============================================================================
// Approximation
// =============================================================================

// -----------------------------------------------------------------------------
void AffineTransformation
::ApproximateDOFs(const double *x,  const double *y,  const double *z, const double *t,
                  const double *dx, const double *dy, const double *dz, int no)
{
  // Solve linear system to find 12 DoFs affine transformation which minimizes
  // the mean squared error if none of the parameters is passive/disabled
  if (this->NumberOfActiveDOFs() == 12) {
    PointSet target(no), source(no);
    for (int n = 0; n < no; ++n) {
      target(n) = Point(x[n],         y[n],         z[n]);
      source(n) = Point(x[n] + dx[n], y[n] + dy[n], z[n] + dz[n]);
    }
    this->PutMatrix(ApproximateAffineMatrix(target, source));
    return;
  }

  // Initialize transformation using closed form solution for rigid parameters
  // if neither of these parameters is passive/disabled
  if (this->AllowTranslations()) {
    if (this->AllowRotations()) {
      // Use temporary rigid transformation instead of subclass implementation
      // RigidTransformation::ApproximateDOFs such that
      // RigidTransformation::UpdateDOFs is used instead of
      // AffineTransformation::UpdateDOFs as the latter may result in
      // non-zero non-rigid parameters due to numerical errors
      RigidTransformation rigid;
      rigid.ApproximateDOFs(x, y, z, t, dx, dy, dz, no);
      this->CopyFrom(&rigid);
      if (this->NumberOfActiveDOFs() == 6) return;
    } else {
      this->Reset(); // Not done by HomogeneousTransformation::ApproximateDOFs
    }
  }

  // Find (active) transformation parameters which minimize the mean squared error
  HomogeneousTransformation::ApproximateDOFs(x, y, z, t, dx, dy, dz, no);
}

// ===========================================================================
// Parameters (non-DoFs)
// ===========================================================================

// ---------------------------------------------------------------------------
bool AffineTransformation::Set(const char *param, const char *value)
{
  bool active;
  // Status of translation parameters
  if (strcmp(param, "Allow translations") == 0) {
    if (!FromString(value, active)) return false;
    PutStatus(TX, active ? Active : Passive);
    PutStatus(TY, active ? Active : Passive);
    PutStatus(TZ, active ? Active : Passive);
    return true;
  } else if (strcmp(param, "Allow translation in X") == 0) {
    if (!FromString(value, active)) return false;
    PutStatus(TX, active ? Active : Passive);
    return true;
  } else if (strcmp(param, "Allow translation in Y") == 0) {
    if (!FromString(value, active)) return false;
    PutStatus(TY, active ? Active : Passive);
    return true;
  } else if (strcmp(param, "Allow translation in Z") == 0) {
    if (!FromString(value, active)) return false;
    PutStatus(TZ, active ? Active : Passive);
    return true;
  // Status of rotation parameters
  } else if (strcmp(param, "Allow rotations") == 0) {
    if (!FromString(value, active)) return false;
    PutStatus(RX, active ? Active : Passive);
    PutStatus(RY, active ? Active : Passive);
    PutStatus(RZ, active ? Active : Passive);
    return true;
  } else if (strcmp(param, "Allow rotation in X") == 0) {
    if (!FromString(value, active)) return false;
    PutStatus(RX, active ? Active : Passive);
    return true;
  } else if (strcmp(param, "Allow rotation in Y") == 0) {
    if (!FromString(value, active)) return false;
    PutStatus(RY, active ? Active : Passive);
    return true;
  } else if (strcmp(param, "Allow rotation in Z") == 0) {
    if (!FromString(value, active)) return false;
    PutStatus(RZ, active ? Active : Passive);
    return true;
  // Status of scaling parameters
  } else if (strcmp(param, "Allow scaling") == 0) {
    if (!FromString(value, active)) return false;
    PutStatus(SX, active ? Active : Passive);
    PutStatus(SY, active ? Active : Passive);
    PutStatus(SZ, active ? Active : Passive);
    return true;
  } else if (strcmp(param, "Allow scaling in X") == 0) {
    if (!FromString(value, active)) return false;
    PutStatus(SX, active ? Active : Passive);
    return true;
  } else if (strcmp(param, "Allow scaling in Y") == 0) {
    if (!FromString(value, active)) return false;
    PutStatus(SY, active ? Active : Passive);
    return true;
  } else if (strcmp(param, "Allow scaling in Z") == 0) {
    if (!FromString(value, active)) return false;
    PutStatus(SZ, active ? Active : Passive);
    return true;
  // Status of shearing parameters
  } else if (strcmp(param, "Allow shearing") == 0) {
    if (!FromString(value, active)) return false;
    PutStatus(SXY, active ? Active : Passive);
    PutStatus(SXZ, active ? Active : Passive);
    PutStatus(SYZ, active ? Active : Passive);
    return true;
  } else if (strcmp(param, "Allow shearing of XY") == 0 ||
             strcmp(param, "Allow shearing of YX") == 0) {
    if (!FromString(value, active)) return false;
    PutStatus(SXY, active ? Active : Passive);
    return true;
  } else if (strcmp(param, "Allow shearing of XZ") == 0 ||
             strcmp(param, "Allow shearing of ZX") == 0) {
    if (!FromString(value, active)) return false;
    PutStatus(SXZ, active ? Active : Passive);
    return true;
  } else if (strcmp(param, "Allow shearing of YZ") == 0 ||
             strcmp(param, "Allow shearing of ZY") == 0) {
    if (!FromString(value, active)) return false;
    PutStatus(SYZ, active ? Active : Passive);
    return true;
  }
  return false;
}

// ---------------------------------------------------------------------------
ParameterList AffineTransformation::Parameter() const
{
  ParameterList params;
  // Status of translation parameters
  if (GetStatus(TX) == GetStatus(TY) && GetStatus(TY) == GetStatus(TZ)) {
    Insert(params, "Allow translations",     (GetStatus(TX) == Active ? "Yes" : "No"));
  } else {
    Insert(params, "Allow translation in X", (GetStatus(TX) == Active ? "Yes" : "No"));
    Insert(params, "Allow translation in Y", (GetStatus(TY) == Active ? "Yes" : "No"));
    Insert(params, "Allow translation in Z", (GetStatus(TZ) == Active ? "Yes" : "No"));
  }
  // Status of rotation parameters
  if (GetStatus(RX) == GetStatus(RY) && GetStatus(RY) == GetStatus(RZ)) {
    Insert(params, "Allow rotations",     (GetStatus(RX) == Active ? "Yes" : "No"));
  } else {
    Insert(params, "Allow rotation in X", (GetStatus(RX) == Active ? "Yes" : "No"));
    Insert(params, "Allow rotation in Y", (GetStatus(RY) == Active ? "Yes" : "No"));
    Insert(params, "Allow rotation in Z", (GetStatus(RZ) == Active ? "Yes" : "No"));
  }
  // Status of scaling parameters
  if (GetStatus(SX) == GetStatus(SY) && GetStatus(SY) == GetStatus(SZ)) {
    Insert(params, "Allow scaling",      (GetStatus(SX) == Active ? "Yes" : "No"));
  } else {
    Insert(params, "Allow scaling in X", (GetStatus(SX) == Active ? "Yes" : "No"));
    Insert(params, "Allow scaling in Y", (GetStatus(SY) == Active ? "Yes" : "No"));
    Insert(params, "Allow scaling in Z", (GetStatus(SZ) == Active ? "Yes" : "No"));
  }
  // Status of scaling parameters
  if (GetStatus(SXY) == GetStatus(SXZ) && GetStatus(SXZ) == GetStatus(SYZ)) {
    Insert(params, "Allow shearing",       (GetStatus(SXY) == Active ? "Yes" : "No"));
  } else {
    Insert(params, "Allow shearing of XY", (GetStatus(SXY) == Active ? "Yes" : "No"));
    Insert(params, "Allow shearing of XZ", (GetStatus(SXZ) == Active ? "Yes" : "No"));
    Insert(params, "Allow shearing of YZ", (GetStatus(SYZ) == Active ? "Yes" : "No"));
  }
  return params;
}

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
Matrix AffineTransformation::DOFs2Matrix(const double *param)
{
  Matrix matrix(4, 4);
  AffineParametersToMatrix(param[TX], param[TY], param[TZ],
                           param[RX]  * rad_per_deg,
                           param[RY]  * rad_per_deg,
                           param[RZ]  * rad_per_deg,
                           param[SX]  / 100.0,
                           param[SY]  / 100.0,
                           param[SZ]  / 100.0,
                           param[SXY] * rad_per_deg,
                           param[SXZ] * rad_per_deg,
                           param[SYZ] * rad_per_deg,
                           matrix);
  return matrix;
}

// -----------------------------------------------------------------------------
void AffineTransformation::Matrix2DOFs(const Matrix &matrix, double *param)
{
  // Get affine parameters
  MatrixToAffineParameters(matrix, param[TX],  param[TY],  param[TZ],
                                   param[RX],  param[RY],  param[RZ],
                                   param[SX],  param[SY],  param[SZ],
                                   param[SXY], param[SXZ], param[SYZ]);
  
  // Convert angles to degrees
  param[RX]  *= deg_per_rad;
  param[RY]  *= deg_per_rad;
  param[RZ]  *= deg_per_rad;
  param[SXY] *= deg_per_rad;
  param[SXZ] *= deg_per_rad;
  param[SYZ] *= deg_per_rad;

  // Convert scales to percentages
  param[SX] *= 100;
  param[SY] *= 100;
  param[SZ] *= 100;
}

// ---------------------------------------------------------------------------
void AffineTransformation::UpdateMatrix()
{
  // Convert angles to radians
  const double rx  = _Param[RX]  * rad_per_deg;
  const double ry  = _Param[RY]  * rad_per_deg;
  const double rz  = _Param[RZ]  * rad_per_deg;
  const double sxy = _Param[SXY] * rad_per_deg;
  const double sxz = _Param[SXZ] * rad_per_deg;
  const double syz = _Param[SYZ] * rad_per_deg;

  // Update sines, cosines, and tans
  _cosrx  = cos(rx);
  _cosry  = cos(ry);
  _cosrz  = cos(rz);
  _sinrx  = sin(rx);
  _sinry  = sin(ry);
  _sinrz  = sin(rz);
  _tansxy = tan(sxy);
  _tansxz = tan(sxz);
  _tansyz = tan(syz);

  // Update transformation matrix
  AffineParametersToMatrix(_Param[TX], _Param[TY], _Param[TZ], rx, ry, rz,
                           _Param[SX] / 100.0, _Param[SY] / 100.0, _Param[SZ] / 100.0,
                           sxy, sxz, syz, _matrix);
  _inverse = _matrix.Inverse();
}

// ---------------------------------------------------------------------------
void AffineTransformation::UpdateDOFs()
{
  // Get affine parameters
  MatrixToAffineParameters(_matrix, _Param[TX],  _Param[TY],  _Param[TZ],
                                    _Param[RX],  _Param[RY],  _Param[RZ],
                                    _Param[SX],  _Param[SY],  _Param[SZ],
                                    _Param[SXY], _Param[SXZ], _Param[SYZ]);

  // Reset disabled/passive parameters
  if (_Status[TX ] == Passive) _Param[TX ] = .0;
  if (_Status[TY ] == Passive) _Param[TY ] = .0;
  if (_Status[TZ ] == Passive) _Param[TZ ] = .0;
  if (_Status[RX ] == Passive) _Param[RX ] = .0;
  if (_Status[RY ] == Passive) _Param[RY ] = .0;
  if (_Status[RZ ] == Passive) _Param[RZ ] = .0;
  if (_Status[SX ] == Passive) _Param[SX ] = 1.0;
  if (_Status[SY ] == Passive) _Param[SY ] = 1.0;
  if (_Status[SZ ] == Passive) _Param[SZ ] = 1.0;
  if (_Status[SXY] == Passive) _Param[SXY] = .0;
  if (_Status[SYZ] == Passive) _Param[SYZ] = .0;
  if (_Status[SXZ] == Passive) _Param[SXZ] = .0;

  // Update sines, cosines, and tans
  _cosrx  = cos(_Param[RX]);
  _cosry  = cos(_Param[RY]);
  _cosrz  = cos(_Param[RZ]);
  _sinrx  = sin(_Param[RX]);
  _sinry  = sin(_Param[RY]);
  _sinrz  = sin(_Param[RZ]);
  _tansxy = tan(_Param[SXY]);
  _tansxz = tan(_Param[SXZ]);
  _tansyz = tan(_Param[SYZ]);

  // Convert angles to degrees
  _Param[RX]  *= deg_per_rad;
  _Param[RY]  *= deg_per_rad;
  _Param[RZ]  *= deg_per_rad;
  _Param[SXY] *= deg_per_rad;
  _Param[SXZ] *= deg_per_rad;
  _Param[SYZ] *= deg_per_rad;

  // Convert scales to percentages
  _Param[SX] *= 100.0;
  _Param[SY] *= 100.0;
  _Param[SZ] *= 100.0;
}

// -----------------------------------------------------------------------------
bool AffineTransformation::CopyFrom(const Transformation *other)
{
  if (strcmp(other->NameOfClass(), "SimilarityTransformation") == 0) {
    const SimilarityTransformation *sim;
    sim = dynamic_cast<const SimilarityTransformation *>(other);
    this->Reset();
    if (_Status[TX] == Active) _Param[TX] = sim->GetTranslationX();
    if (_Status[TY] == Active) _Param[TY] = sim->GetTranslationY();
    if (_Status[TZ] == Active) _Param[TZ] = sim->GetTranslationZ();
    if (_Status[RX] == Active) _Param[RX] = sim->GetRotationX();
    if (_Status[RY] == Active) _Param[RY] = sim->GetRotationY();
    if (_Status[RZ] == Active) _Param[RZ] = sim->GetRotationZ();
    if (_Status[SX] == Active) _Param[SX] = sim->GetScale();
    if (_Status[SY] == Active) _Param[SY] = sim->GetScale();
    if (_Status[SZ] == Active) _Param[SZ] = sim->GetScale();
    this->Update(MATRIX);
    return true;
  } else {
    return HomogeneousTransformation::CopyFrom(other);
  }
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
void AffineTransformation::JacobianDOFs(double jac[3], int dof, double x, double y, double z, double, double) const
{
  // T(x) = R W S x + t
  // where
  //        R is rotation matrix
  //        W is shear matrix
  //        S is scale matrix
  //        t is translation vector

  double r[3][3] = {{.0, .0, .0}, {.0, .0, .0}, {.0, .0, .0}};

  // Derivative w.r.t translation or initialization of those rotation matrix
  // elements, which are needed below for the computation of the derivative
  // w.r.t the scale or shear parameter
  switch (dof) {
    case TX:
      jac[0] = 1.0;
      jac[1] = jac[2] = .0;
      return;
    case TY:
      jac[1] = 1.0;
      jac[0] = jac[2] = .0;
      return;
    case TZ:
      jac[2] = 1.0;
      jac[0] = jac[1] = .0;
      return;
    case SX: case SXY: case SXZ:
      r[0][0] = _cosry * _cosrz;
      r[1][0] = _sinrx * _sinry * _cosrz - _cosrx * _sinrz;
      r[2][0] = _cosrx * _sinry * _cosrz + _sinrx * _sinrz;
      if (dof != SX) break;
    case SY: case SYZ: // SX continued, i.e., do not insert a break before
      r[0][1] = _cosry * _sinrz;
      r[1][1] = _sinrx * _sinry * _sinrz + _cosrx * _cosrz;
      r[2][1] = _cosrx * _sinry * _sinrz - _sinrx * _cosrz;
      if (dof == SYZ) break;
    case SZ: // SX|SY continued, i.e., do not insert a break before
      r[0][2] = -_sinry;
      r[1][2] = _sinrx * _cosry;
      r[2][2] = _cosrx * _cosry;
      break;
  }

  // Derivative w.r.t scale or shear parameter
  switch (dof) {
    case SX: {
      jac[0] = r[0][0] * x / 100.0;
      jac[1] = r[1][0] * x / 100.0;
      jac[2] = r[2][0] * x / 100.0;
    } break;
    case SY: {
      jac[0] = (r[0][0] * _tansxy + r[0][1]) * y / 100.0;
      jac[1] = (r[1][0] * _tansxy + r[1][1]) * y / 100.0;
      jac[2] = (r[2][0] * _tansxy + r[2][1]) * y / 100.0;
    } break;
    case SZ: {
      jac[0] = (r[0][0] * _tansxz + r[0][1] * _tansyz + r[0][2]) * z / 100.0;
      jac[1] = (r[1][0] * _tansxz + r[1][1] * _tansyz + r[1][2]) * z / 100.0;
      jac[2] = (r[2][0] * _tansxz + r[2][1] * _tansyz + r[2][2]) * z / 100.0;
    } break;
    case SXY: {
      const double dsxy = (rad_per_deg / pow(cos(_Param[SXY] * rad_per_deg), 2.0)) * y * _Param[SY] / 100.0;
      jac[0] = r[0][0] * dsxy;
      jac[1] = r[1][0] * dsxy;
      jac[2] = r[2][0] * dsxy;
    } break;
    case SYZ: {
      const double dsyz = (rad_per_deg / pow(cos(_Param[SYZ] * rad_per_deg), 2.0)) * z * _Param[SZ] / 100.0;
      jac[0] = r[0][1] * dsyz;
      jac[1] = r[1][1] * dsyz;
      jac[2] = r[2][1] * dsyz;
    } break;
    case SXZ: {
      const double dsxz = (rad_per_deg / pow(cos(_Param[SXZ] * rad_per_deg), 2.0)) * z * _Param[SZ] / 100.0;
      jac[0] = r[0][0] * dsxz;
      jac[1] = r[1][0] * dsxz;
      jac[2] = r[2][0] * dsxz;
    } break;
    // Derivate w.r.t rotation parameter
    default: {
      // apply scaling
      x *= _Param[SX] / 100.0;
      y *= _Param[SY] / 100.0;
      z *= _Param[SZ] / 100.0;
      // apply shearing
      x += _tansxy * y + _tansxz * z;
      y += _tansyz * z;
      // compute derivative w.r.t rigid parameter
      RigidTransformation::JacobianDOFs(jac, dof, x, y, z);
    }
  }
}

// -----------------------------------------------------------------------------
void AffineTransformation::DeriveJacobianWrtDOF(Matrix &dJdp, int dof, double x, double y, double z, double, double) const
{
  dJdp.Initialize(3, 3);

  double r[3][3] = {{.0, .0, .0}, {.0, .0, .0}, {.0, .0, .0}};

  // Derivative w.r.t translation or initialization of those rotation matrix
  // elements, which are needed below for the computation of the derivative
  // w.r.t the scale or shear parameter
  switch (dof) {
    case SX: case SXY: case SXZ:
      r[0][0] = _cosry * _cosrz;
      r[1][0] = _sinrx * _sinry * _cosrz - _cosrx * _sinrz;
      r[2][0] = _cosrx * _sinry * _cosrz + _sinrx * _sinrz;
      if (dof != SX) break;
    case SY: case SYZ: // SX continued, i.e., do not insert a break before
      r[0][1] = _cosry * _sinrz;
      r[1][1] = _sinrx * _sinry * _sinrz + _cosrx * _cosrz;
      r[2][1] = _cosrx * _sinry * _sinrz - _sinrx * _cosrz;
      if (dof == SYZ) break;
    case SZ: // SX|SY continued, i.e., do not insert a break before
      r[0][2] = -_sinry;
      r[1][2] = _sinrx * _cosry;
      r[2][2] = _cosrx * _cosry;
      break;
  }

  // Derivative w.r.t scale or shear parameter
  switch (dof) {
    case SX: {
      dJdp(0, 0) = r[0][0] / 100.0;
      dJdp(1, 0) = r[1][0] / 100.0;
      dJdp(2, 0) = r[2][0] / 100.0;
    } break;
    case SY: {
      dJdp(0, 1) = (r[0][0] * _tansxy + r[0][1]) / 100.0;
      dJdp(1, 1) = (r[1][0] * _tansxy + r[1][1]) / 100.0;
      dJdp(2, 1) = (r[2][0] * _tansxy + r[2][1]) / 100.0;
    } break;
    case SZ: {
      dJdp(0, 2) = (r[0][0] * _tansxz + r[0][1] * _tansyz + r[0][2]) / 100.0;
      dJdp(1, 2) = (r[1][0] * _tansxz + r[1][1] * _tansyz + r[1][2]) / 100.0;
      dJdp(2, 2) = (r[2][0] * _tansxz + r[2][1] * _tansyz + r[2][2]) / 100.0;
    } break;
    case SXY: {
      const double dsxy = (rad_per_deg / pow(cos(_Param[SXY] * rad_per_deg), 2.0)) * _Param[SY] / 100.0;
      dJdp(0, 1) = r[0][0] * dsxy;
      dJdp(1, 1) = r[1][0] * dsxy;
      dJdp(2, 1) = r[2][0] * dsxy;
    } break;
    case SYZ: {
      const double dsyz = (rad_per_deg / pow(cos(_Param[SYZ] * rad_per_deg), 2.0)) * _Param[SZ] / 100.0;
      dJdp(0, 2) = r[0][1] * dsyz;
      dJdp(1, 2) = r[1][1] * dsyz;
      dJdp(2, 2) = r[2][1] * dsyz;
    } break;
    case SXZ: {
      const double dsxz = (rad_per_deg / pow(cos(_Param[SXZ] * rad_per_deg), 2.0)) * _Param[SZ] / 100.0;
      dJdp(0, 2) = r[0][0] * dsxz;
      dJdp(1, 2) = r[1][0] * dsxz;
      dJdp(2, 2) = r[2][0] * dsxz;
    } break;
    // Derivate w.r.t rotation parameter
    default: {
      // compute derivative w.r.t rigid parameter
      RigidTransformation::DeriveJacobianWrtDOF(dJdp, dof, x, y, z);

      //apply W * S (product of shearing and scaling matrices)
      Matrix ws;
      ws.Initialize(3, 3);

      ws(0, 0) = _Param[SX] / 100.0;
      ws(1, 0) = .0;
      ws(2, 0) = .0;
      ws(0, 1) = _tansxy * _Param[SY] / 100.0;
      ws(1, 1) = _Param[SY] / 100.0;
      ws(2, 1) = .0;
      ws(0, 2) = _tansxz * _Param[SZ] / 100.0;
      ws(1, 2) = _tansyz * _Param[SZ] / 100.0;
      ws(2, 2) = _Param[SZ] / 100.0;

      dJdp *= ws;
    }
  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void AffineTransformation::Print(ostream &os, Indent indent) const
{
  RigidTransformation::Print(os, indent);

  os.setf(ios::right);
  os.setf(ios::fixed);
  streamsize previous_precision = os.precision(4);

  if (_Status[SX] == Active || !fequal(_Param[SX], 100.0, 1e-4) ||
      _Status[SY] == Active || !fequal(_Param[SY], 100.0, 1e-4) ||
      _Status[SZ] == Active || !fequal(_Param[SZ], 100.0, 1e-4)) {
    os << indent;
    if (_Status[SX]  == Active) os << "sx  = "  << setw(8) << _Param[SX]  << "  ";
    if (_Status[SY]  == Active) os << "sy  = "  << setw(8) << _Param[SY]  << "  ";
    if (_Status[SZ]  == Active) os << "sz  = "  << setw(8) << _Param[SZ]  << "  ";
    os << endl;
  }
  if (_Status[SXY] == Active || !fequal(_Param[SXY], .0, 1e-4) ||
      _Status[SXZ] == Active || !fequal(_Param[SXZ], .0, 1e-4) ||
      _Status[SYZ] == Active || !fequal(_Param[SYZ], .0, 1e-4)) {
    os << indent;
    if (_Status[SXY] == Active) os << "sxy = " << setw(8) << _Param[SXY] << "  ";
    if (_Status[SYZ] == Active) os << "syz = " << setw(8) << _Param[SYZ] << "  ";
    if (_Status[SXZ] == Active) os << "sxz = " << setw(8) << _Param[SXZ] << "  ";
    os << endl;
  }

  os.precision(previous_precision);
  os.unsetf(ios::right);
  os.unsetf(ios::fixed);
}

// -----------------------------------------------------------------------------
bool AffineTransformation::CanRead(TransformationType format) const
{
  switch (format) {
    case TRANSFORMATION_RIGID:
    case TRANSFORMATION_SIMILARITY:
    case TRANSFORMATION_AFFINE:
      return true;
    default:
      return false;
  }
}

// -----------------------------------------------------------------------------
Cofstream &AffineTransformation::Write(Cofstream &to) const
{
  unsigned int ndofs = this->NumberOfDOFs();

  // Write rigid parameters only if additional affine components not used
  if ((abs(_Param[SX] - 100.0) < 0.0001) &&
      (abs(_Param[SY] - 100.0) < 0.0001) &&
      (abs(_Param[SZ] - 100.0) < 0.0001) &&
      (abs(_Param[SXY])        < 0.0001) &&
      (abs(_Param[SXZ])        < 0.0001) &&
      (abs(_Param[SYZ])        < 0.0001)) {
    ndofs = 6;
  }

  // Write magic no. for transformations
  unsigned int magic_no = TRANSFORMATION_MAGIC;
  to.WriteAsUInt(&magic_no, 1);

  // Write transformation type
  unsigned int trans_type = (ndofs == 6) ? TRANSFORMATION_RIGID
                                         : TRANSFORMATION_AFFINE;
  to.WriteAsUInt(&trans_type, 1);

  // Write number of rigid parameters
  to.WriteAsUInt(&ndofs, 1);

  // Write rigid transformation parameters
  to.WriteAsDouble(_Param, ndofs);

  return to;
}

// -----------------------------------------------------------------------------
Cifstream &AffineTransformation::ReadDOFs(Cifstream &from, TransformationType)
{
  // Read number of parameters
  unsigned int ndofs;
  from.ReadAsUInt(&ndofs, 1);
  if (ndofs != 6 && ndofs != 7 && ndofs != 9 && ndofs != 12) {
    cerr << this->NameOfClass() << "::Read: Invalid no. of parameters: " << ndofs << endl;
    exit(1);
  }

  // Read transformation parameters
  from.ReadAsDouble(_Param, ndofs);
  if (ndofs == 6)  _Param[SX]  = _Param[SY]  = _Param[SZ]  = 100;
  if (ndofs == 7)                _Param[SY]  = _Param[SZ]  = _Param[SX];
  if (ndofs != 12) _Param[SXY] = _Param[SXZ] = _Param[SYZ] = 0;

  // Update transformation matrix
  this->Update(MATRIX);

  return from;
}


} // namespace mirtk

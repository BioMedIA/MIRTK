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

#include "mirtk/InverseAffineTransformation.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
InverseAffineTransformation::InverseAffineTransformation(AffineTransformation *t)
:
  _Transformation(NULL)
{
  _TransformationObserver.Bind(
    ModifiedEvent,
    MakeDelegate(this, &InverseAffineTransformation::OnTransformationChanged)
  );
  Transformation(t);
}

// -----------------------------------------------------------------------------
InverseAffineTransformation::~InverseAffineTransformation()
{
  Transformation(NULL);
}

// -----------------------------------------------------------------------------
void InverseAffineTransformation::Transformation(AffineTransformation *t)
{
  if (_Transformation) _Transformation->DeleteObserver(_TransformationObserver);
  _Transformation = t;
  if (_Transformation) _Transformation->AddObserver(_TransformationObserver);
  this->OnTransformationChanged();
}

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool InverseAffineTransformation::HasSameDOFsAs(const class Transformation *t) const
{
  return HaveSameDOFs(_Transformation, t);
}

// -----------------------------------------------------------------------------
void InverseAffineTransformation::OnTransformationChanged()
{
  // Do nothing if no decorated transformation set
  if (!_Transformation) return;
  // Copy status of transformation parameters
  const int ndofs = _Transformation->NumberOfDOFs();
  if (ndofs <= 12) {
    _Status[TX] = _Transformation->GetStatus(TX);
    _Status[TY] = _Transformation->GetStatus(TY);
    _Status[TZ] = _Transformation->GetStatus(TZ);
    _Status[RX] = _Transformation->GetStatus(RX);
    _Status[RY] = _Transformation->GetStatus(RY);
    _Status[RZ] = _Transformation->GetStatus(RZ);
  }
  if (ndofs == 7) {
    _Status[SX] = _Status[SY] = _Status[SZ] = _Transformation->GetStatus(SG);
  } else if (ndofs == 9 || ndofs == 12) {
    _Status[SX] = _Transformation->GetStatus(SX);
    _Status[SY] = _Transformation->GetStatus(SY);
    _Status[SZ] = _Transformation->GetStatus(SZ);
  }
  if (ndofs == 12) {
    _Status[SXY] = _Transformation->GetStatus(SXY);
    _Status[SXZ] = _Transformation->GetStatus(SXZ);
    _Status[SYZ] = _Transformation->GetStatus(SYZ);
  }
  // Put matrix of inverse transformation
  this->PutMatrix(_Transformation->GetInverseMatrix());
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
// Multiply 3x3 matrix and 3x1 vector
static inline void mul(const double m[3][3], double v[3])
{
  const double x = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
  const double y = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
  const double z = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];
  v[0] = x, v[1] = y, v[2] = z;
}

// -----------------------------------------------------------------------------
void InverseAffineTransformation::JacobianDOFs(double jac[3], int dof, double x, double y, double z, double, double) const
{
  // T(x) = inv(S) inv(W) inv(R) (x - t)
  // where
  //        R is rotation matrix
  //        W is shear matrix
  //        S is scale matrix
  //        t is translation vector

  double m[3][3];

  const double *param  = _Transformation->_Param;
  const double &cosrx  = _Transformation->_cosrx;
  const double &cosry  = _Transformation->_cosry;
  const double &cosrz  = _Transformation->_cosrz;
  const double &sinrx  = _Transformation->_sinrx;
  const double &sinry  = _Transformation->_sinry;
  const double &sinrz  = _Transformation->_sinrz;
  const double &tansxy = _Transformation->_tansxy;
  const double &tansxz = _Transformation->_tansxz;
  const double &tansyz = _Transformation->_tansyz;

  // inverse translation
  if (dof == TX) {
    jac[0] = -1.0;
    jac[1] = .0;
    jac[2] = .0;
  } else if (dof == TY) {
    jac[0] = .0;
    jac[1] = -1.0;
    jac[2] = .0;
  } else if (dof == TZ) {
    jac[0] = .0;
    jac[1] = .0;
    jac[2] = -1.0;
  } else {
    jac[0] = x - param[TX];
    jac[1] = y - param[TY];
    jac[2] = z - param[TZ];
  }

  // inverse rotation
  if (dof == RX) {
    m[0][0] = .0;
    m[0][1] = sinrx * sinrz + cosrx * cosrz * sinry;
    m[0][2] = cosrx * sinrz - cosrz * sinrx * sinry;
    m[1][0] = .0;
    m[1][1] = cosrx * sinry * sinrz - cosrz * sinrx;
    m[1][2] = - cosrx * cosrz - sinrx * sinry * sinrz;
    m[2][0] = .0;
    m[2][1] = cosrx * cosry;
    m[2][2] = - cosry * sinrx;
    jac[0] *= rad_per_deg, jac[1] *= rad_per_deg, jac[2] *= rad_per_deg;
  } else if (dof == RY) {
    m[0][0] = - cosrz * sinry;
    m[0][1] = cosry * cosrz * sinrx;
    m[0][2] = cosrx * cosry * cosrz;
    m[1][0] = - sinry * sinrz;
    m[1][1] = cosry * sinrx * sinrz;
    m[1][2] = cosrx * cosry * sinrz;
    m[2][0] = - cosry;
    m[2][1] = - sinrx * sinry;
    m[2][2] = - cosrx * sinry;
    jac[0] *= rad_per_deg, jac[1] *= rad_per_deg, jac[2] *= rad_per_deg;
  } else if (dof == RZ) {
    m[0][0] = - cosry * sinrz;
    m[0][1] = - cosrx * cosrz - sinrx * sinry * sinrz;
    m[0][2] = cosrz * sinrx - cosrx * sinry * sinrz;
    m[1][0] = cosry * cosrz;
    m[1][1] = cosrz * sinrx * sinry - cosrx * sinrz;
    m[1][2] = sinrx * sinrz + cosrx * cosrz * sinry;
    m[2][0] = .0;
    m[2][1] = .0;
    m[2][2] = .0;
    jac[0] *= rad_per_deg, jac[1] *= rad_per_deg, jac[2] *= rad_per_deg;
  } else {
    m[0][0] = cosry * cosrz;
    m[0][1] = cosrz * sinrx * sinry - cosrx * sinrz;
    m[0][2] = sinrx * sinrz + cosrx * cosrz * sinry;
    m[1][0] = cosry * sinrz;
    m[1][1] = cosrx * cosrz + sinrx * sinry * sinrz;
    m[1][2] = cosrx * sinry * sinrz - cosrz * sinrx;
    m[2][0] = - sinry;
    m[2][1] = cosry * sinrx;
    m[2][2] = cosrx * cosry;
  }

  mul(m, jac);

  // inverse shearing
  memset(m, 0, sizeof(m));
  if (dof == SXY) {
    const double derivative = rad_per_deg / pow(cos(param[SXY] * rad_per_deg), 2);
    m[0][1] = - derivative;
    m[0][2] = tansyz * derivative;
  } else if (dof == SXZ) {
    const double derivative = rad_per_deg / pow(cos(param[SXZ] * rad_per_deg), 2);
    m[0][2] = - derivative;
  } else if (dof == SYZ) {
    const double derivative = rad_per_deg / pow(cos(param[SYZ] * rad_per_deg), 2);
    m[0][2] = tansxy * derivative;
    m[1][2] = - derivative;
  } else {
    m[0][0] = m[1][1] = m[2][2] = 1.0;
    m[0][1] = - tansxy;
    m[0][2] = tansxy * tansyz - tansxz;
    m[1][2] = - tansyz;
  }

  mul(m, jac);

  // inverse scaling
  memset(m, 0, sizeof(m));
  if (dof == SX) {
    m[0][0] = - 100.0 / pow(param[SX], 2);
  } else if (dof == SY) {
    m[1][1] = - 100.0 / pow(param[SY], 2);
  } else if (dof == SZ) {
    m[2][2] = - 100.0 / pow(param[SZ], 2);
  } else {
    m[0][0] = 100.0 / param[SX];
    m[1][1] = 100.0 / param[SY];
    m[2][2] = 100.0 / param[SZ];
  }

  mul(m, jac);
}


} // namespace mirtk

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

#include "mirtk/SimilarityTransformation.h"

#include "mirtk/Math.h"
#include "mirtk/AffineTransformation.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
SimilarityTransformation::SimilarityTransformation(int ndofs)
:
  RigidTransformation(ndofs)
{
  _Param[SG] = 100;
}

// -----------------------------------------------------------------------------
SimilarityTransformation::SimilarityTransformation(const RigidTransformation &t, int ndofs)
:
  RigidTransformation(t, ndofs)
{
  _Param[SG] = 100;
}

// -----------------------------------------------------------------------------
SimilarityTransformation::SimilarityTransformation(const SimilarityTransformation &t, int ndofs)
:
  RigidTransformation(t, ndofs)
{
}

// -----------------------------------------------------------------------------
SimilarityTransformation::SimilarityTransformation()
:
  RigidTransformation(7)
{
  _Param[SG] = 100;
}

// -----------------------------------------------------------------------------
SimilarityTransformation::SimilarityTransformation(const RigidTransformation &t)
:
  RigidTransformation(t, 7)
{
  _Param[SG] = 100;
}

// -----------------------------------------------------------------------------
SimilarityTransformation::SimilarityTransformation(const SimilarityTransformation &t)
:
  RigidTransformation(t, 7)
{
}

// -----------------------------------------------------------------------------
SimilarityTransformation::~SimilarityTransformation()
{
}

// =============================================================================
// Approximation
// =============================================================================

// -----------------------------------------------------------------------------
void SimilarityTransformation
::ApproximateDOFs(const double *x,  const double *y,  const double *z, const double *t,
                  const double *dx, const double *dy, const double *dz, int no)
{
  // Initialize transformation using closed form solution for rigid parameters
  RigidTransformation::ApproximateDOFs(x, y, z, t, dx, dy, dz, no);

  // Iteratively refine transformation
  HomogeneousTransformation::ApproximateDOFs(x, y, z, t, dx, dy, dz, no);
}

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
Matrix SimilarityTransformation::DOFs2Matrix(const double *param)
{
  // Get the rigid components
  Matrix matrix = RigidTransformation::DOFs2Matrix(param);

  // Update similarity transformation: Add global scaling
  Matrix scale(4, 4);
  scale.Ident();
  scale(0, 0) = param[SG] / 100.0;
  scale(1, 1) = param[SG] / 100.0;
  scale(2, 2) = param[SG] / 100.0;
  matrix *= scale;

  return matrix;
}

// -----------------------------------------------------------------------------
void SimilarityTransformation::UpdateMatrix()
{
  // Update rigid transformation
  RigidTransformation::UpdateMatrix();

  // Update similarity transformation: Add global scaling
  Matrix scale(4, 4);
  scale.Ident();
  scale(0, 0) = _Param[SG] / 100.0;
  scale(1, 1) = _Param[SG] / 100.0;
  scale(2, 2) = _Param[SG] / 100.0;
  _matrix *= scale;

  // Update inverse transformation
  _inverse = _matrix.Inverse();
}

// -----------------------------------------------------------------------------
void SimilarityTransformation::UpdateDOFs()
{
  // Use current matrix to extract the transformation parameters based on
  // 7 DOF model, assuming _sx == _sy == _sz and no shearing, see:
  //    http://www.acm.org/pubs/tog/GraphicsGems/gemsii/unmatrix.c
  //    https://github.com/erich666/GraphicsGems/blob/master/gemsii/unmatrix.c
  // (from Graphics Gems II)

  // Use _matrix to evaluate the similarity transformation parameters.
  // It is assumed that there is no perspective transformation
  int i;
  const double TOL = 0.000001;

  if (abs(_matrix(3,3) - 1.0) > TOL) {
    cerr << "SimilarityTransformation::UpdateParameter" << endl;
    cerr << "Value at _matrix(3,3) must equal 1." << endl;
    exit(1);
  }

  if (abs(_matrix.Det()) < TOL) {
    cerr << "SimilarityTransformation::UpdateParameter" << endl;
    cerr << "Matrix singular (or very close to singular!)." << endl;
    exit(1);
  }

  // First Part Of Graphics Gems Code Ignored Because It Relates To
  // Perspective Transformation
  if (abs(_matrix(3, 0)) > TOL ||
      abs(_matrix(3, 1)) > TOL ||
      abs(_matrix(3, 2)) > TOL ) {
    cerr << "SimilarityTransformation::UpdateParameter" << endl;
    cerr << "Matrix contains perspective components." << endl;
    exit(1);
  }

  Matrix copy(4, 4);
  copy = _matrix;

  // Get scale y manipulating the columns of the upper left 3x3 sub-matrix
  Vector col_0, col_1, col_2;
  col_0.Initialize(3);
  col_1.Initialize(3);
  col_2.Initialize(3);
  for (i = 0; i < 3; ++i) {
    col_0(i) = copy(i, 0);
    col_1(i) = copy(i, 1);
    col_2(i) = copy(i, 2);
  }

  // Compute X scale factor and normalize first col
  _Param[SG] = col_0.Norm();
  col_0     /= _Param[SG];
  col_1     /= _Param[SG];
  col_2     /= _Param[SG];

  // At this point, the columns should be orthonormal.  Check for a coordinate
  // system flip.  If the determinant is -1, then negate the matrix and the
  // scaling factor
  Vector col_1_x_col_2;
  col_1_x_col_2.Initialize(3);
  col_1_x_col_2 = col_1.CrossProduct(col_2);

  if (col_0.ScalarProduct(col_1_x_col_2) < 0) {
    _Param[SG] *= -1;
    col_0      *= -1;
    col_1      *= -1;
    col_2      *= -1;
  }

  // Convert scales to percentages
  _Param[SG] *= 100;

  // Put the rotation matrix components into the upper left 3x3 submatrix
  for (i = 0; i < 3; ++i) {
    copy(i, 0) = col_0(i);
    copy(i, 1) = col_1(i);
    copy(i, 2) = col_2(i);
  }

  // Now get the rigid transformation parameters
  RigidTransformation::Matrix2DOFs(copy, _Param);
  UpdateRotationSineCosine();
}

// ---------------------------------------------------------------------------
bool SimilarityTransformation::CopyFrom(const Transformation *other)
{
  const AffineTransformation *aff = NULL;
  if ((aff = dynamic_cast<const AffineTransformation *>(other))) {
    this->Reset();
    if (_Status[TX] == Active) _Param[TX] = aff->GetTranslationX();
    if (_Status[TY] == Active) _Param[TY] = aff->GetTranslationY();
    if (_Status[TZ] == Active) _Param[TZ] = aff->GetTranslationZ();
    if (_Status[RX] == Active) _Param[RX] = aff->GetRotationX();
    if (_Status[RY] == Active) _Param[RY] = aff->GetRotationY();
    if (_Status[RZ] == Active) _Param[RZ] = aff->GetRotationZ();
    if (_Status[SG] == Active) _Param[SG] = (aff->GetScaleX() + aff->GetScaleY() + aff->GetScaleZ()) / 3.0;
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
void SimilarityTransformation::JacobianDOFs(double jac[3], int dof, double x, double y, double z, double, double) const
{
  if (dof == SG) {
    jac[0] = (  (_cosry * _cosrz                           ) * x
              + (_cosry * _sinrz                           ) * y
              + (-_sinry                                   ) * z) / 100.0;
    jac[1] = (  (_sinrx * _sinry * _cosrz - _cosrx * _sinrz) * x
              + (_sinrx * _sinry * _sinrz + _cosrx * _cosrz) * y
              + (_sinrx * _cosry                           ) * z) / 100.0;
    jac[2] = (  (_cosrx * _sinry * _cosrz + _sinrx * _sinrz) * x
              + (_cosrx * _sinry * _sinrz - _sinrx * _cosrz) * y
              + (_cosrx * _cosry                           ) * z) / 100.0;
  } else {
    const double s = _Param[SG] / 100.0;
    RigidTransformation::JacobianDOFs(jac, dof, x * s, y * s, z * s);
  }
}

// -----------------------------------------------------------------------------
void SimilarityTransformation::DeriveJacobianWrtDOF(Matrix &dJdp, int dof, double x, double y, double z, double, double) const
{
  if (dof == SG) {
    dJdp(0, 0) = (_cosry * _cosrz) / 100.0;
    dJdp(1, 0) = (_sinrx * _sinry * _cosrz - _cosrx * _sinrz) / 100.0;
    dJdp(2, 0) = (_cosrx * _sinry * _cosrz + _sinrx * _sinrz) / 100.0;
    dJdp(0, 1) = (_cosry * _sinrz) / 100.0;
    dJdp(1, 1) = (_sinrx * _sinry * _sinrz + _cosrx * _cosrz) / 100.0;
    dJdp(2, 1) = (_cosrx * _sinry * _sinrz - _sinrx * _cosrz) / 100.0;
    dJdp(0, 2) = (-_sinry) / 100.0;
    dJdp(1, 2) = (_sinrx * _cosry) / 100.0;
    dJdp(2, 2) = (_cosrx * _cosry) / 100.0;
  } else {
    RigidTransformation::DeriveJacobianWrtDOF(dJdp, dof, x, y, z);
    dJdp *= (_Param[SG] / 100.0);
  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void SimilarityTransformation::Print(ostream &os, Indent indent) const
{
  RigidTransformation::Print(os, indent);

  os.setf(ios::right);
  os.setf(ios::fixed);
  streamsize previous_precision = os.precision(4);

  if (_Status[SX] == Active || !fequal(_Param[SX], 100.0, 1e-4)) {
    os << indent << "s   = " << setw(8) << _Param[SG] << endl;
  }

  os.precision(previous_precision);
  os.unsetf(ios::right);
  os.unsetf(ios::fixed);
}

// -----------------------------------------------------------------------------
bool SimilarityTransformation::CanRead(TransformationType format) const
{
  switch (format) {
    case TRANSFORMATION_RIGID:
    case TRANSFORMATION_SIMILARITY:
      return true;
    default:
      return false;
  }
}

// -----------------------------------------------------------------------------
Cofstream &SimilarityTransformation::Write(Cofstream &to) const
{
  unsigned int ndofs = 7;
  if (abs(_Param[SG] - 100.0) < 0.0001) ndofs = 6;

  // Write magic no. for transformations
  unsigned int magic_no = TRANSFORMATION_MAGIC;
  to.WriteAsUInt(&magic_no, 1);

  // Write transformation type
  unsigned int trans_type = (ndofs == 6) ? TRANSFORMATION_RIGID
                                         : TRANSFORMATION_SIMILARITY;
  to.WriteAsUInt(&trans_type, 1);

  // Write number of parameters
  to.WriteAsUInt(&ndofs, 1);

  // Write transformation parameters
  for (int dof = 0; dof < static_cast<int>(ndofs); ++dof) {
    double data = _Param[dof];
    to.WriteAsDouble(&data, 1);
  }

  return to;
}

// -----------------------------------------------------------------------------
Cifstream &SimilarityTransformation::ReadDOFs(Cifstream &from, TransformationType)
{
  // Read number of parameters
  unsigned int ndofs;
  from.ReadAsUInt(&ndofs, 1);
  if (ndofs < 6 || ndofs > 7) {
    cerr << this->NameOfClass() << "::Read: Invalid no. of parameters: " << ndofs << endl;
    exit(1);
  }

  // Read parameters
  double x;
  for (int i = 0; i < static_cast<int>(ndofs); ++i) {
    from.ReadAsDouble(&x, 1);
    _Param[i] = x;
  }
  if (ndofs == 6) _Param[SG] = 100;

  // Update matrix
  this->Update(MATRIX);

  return from;
}


} // namespace mirtk

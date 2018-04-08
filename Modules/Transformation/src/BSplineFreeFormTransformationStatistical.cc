/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2015 Stefan Pszczolkowski Parraguez
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

#include "mirtk/BSplineFreeFormTransformationStatistical.h"

#include "mirtk/Memory.h"
#include "mirtk/Vector3D.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
BSplineFreeFormTransformationStatistical
::BSplineFreeFormTransformationStatistical()
{
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformationStatistical
::BSplineFreeFormTransformationStatistical(const ImageAttributes    &attr,
                                           CPStatus              ****status,
                                           const Matrix             &mat,
                                           const Vector             &vec)
:
  _BasisVectors(mat),
  _MeanVector  (vec)
{
  // Initialize DOFs
  InitializeDOFs(mat.Cols());

  // Initialize control points with dofs flag set to false, i.e.,
  // specifying that in case of this transformation, the CPs are *not* the DOFs
  InitializeCPs(attr, false);

  // Copy control point status
  for (int k = 0; k < attr._z; ++k)
  for (int j = 0; j < attr._y; ++j)
  for (int i = 0; i < attr._x; ++i) {
    _CPStatus[0][k][j][i] = status[0][k][j][i];
  }
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformationStatistical
::BSplineFreeFormTransformationStatistical(const BSplineFreeFormTransformationStatistical &ffd)
:
  _BasisVectors(ffd._BasisVectors),
  _MeanVector  (ffd._MeanVector)
{
  // Initialize DOFs
  InitializeDOFs(ffd._BasisVectors.Cols());
  const ImageAttributes &attr = ffd.Attributes();
  if (attr._x > 0 && attr._y > 0 && attr._z > 0) {
    // Initialize control points with dofs flag set to false, i.e.,
    // specifying that in case of this transformation, the CPs are *not* the DOFs
    InitializeCPs(attr, false);
  }
}

// -----------------------------------------------------------------------------
BSplineFreeFormTransformationStatistical::~BSplineFreeFormTransformationStatistical()
{
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationStatistical::Initialize(const ImageAttributes &)
{
  // For the moment, initialize the trasformation as one with no DOFs
  // until the statistical deformation model is read from disk
  InitializeDOFs(0);
}

// =============================================================================
// Approximation/Interpolation
// =============================================================================

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationStatistical
::ApproximateDOFs(const double *x,  const double *y,  const double *z,  const double *t,
                  const double *dx, const double *dy, const double *dz, int no)
{
  // Approximate control point displacements
  BSplineFreeFormTransformation3D::ApproximateDOFs(x, y, z, t, dx, dy, dz, no);
  // Update transformation parameters
  UpdateDOFs();
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationStatistical
::ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                          const double *, const double *, const double *,
                          int, double *, double) const
{
  cerr << this->NameOfClass() << "::ApproximateDOFsGradient: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationStatistical
::Interpolate(const double *dx, const double *dy, const double *dz)
{
  // Interpolate control point displacements
  BSplineFreeFormTransformation3D::Interpolate(dx, dy, dz);
  // Update transformation parameters
  UpdateDOFs();
}

// =============================================================================
// Lattice
// =============================================================================

// -----------------------------------------------------------------------------
bool BSplineFreeFormTransformationStatistical::CropPadPassiveCPs(int, int, int, int, bool)
{
  // Do nothing as the FFD lattice attributes depend on the statistical
  // deformation model and may not be modified once such model was loaded
  return true;
}

// =============================================================================
// Updating
// =============================================================================

// -----------------------------------------------------------------------------
struct BSplineFreeFormTransformationStatisticalCPUpdate
{
  BSplineFreeFormTransformationStatistical          *_FFD;
  BSplineFreeFormTransformationStatistical::CPValue *_Data;

  const Vector *_Mean;
  const Matrix *_Bases;
  double       *_Input;
  int           _NumberOfDOFs;

  void operator()(const blocked_range<int> &re) const
  {
    int    dof_x, dof_y, dof_z;
    double val_x, val_y, val_z;

    for (int cp = re.begin(); cp != re.end(); ++cp) {
      _FFD->IndexToDOFs(cp, dof_x, dof_y, dof_z);
      val_x = _Mean->Get(dof_x);
      val_y = _Mean->Get(dof_y);
      val_z = _Mean->Get(dof_z);
      for(int q = 0; q < _NumberOfDOFs; q++) {
        val_x += _Bases->Get(dof_x, q) * _Input[q];
        val_y += _Bases->Get(dof_y, q) * _Input[q];
        val_z += _Bases->Get(dof_z, q) * _Input[q];
      }
      _Data[cp]._x = val_x;
      _Data[cp]._y = val_y;
      _Data[cp]._z = val_z;
    }
  }
};

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationStatistical::UpdateCPs()
{
  BSplineFreeFormTransformationStatisticalCPUpdate updater;

  updater._FFD          = this;
  updater._Data         = _CPImage.Data();
  updater._Mean         = &_MeanVector;
  updater._Bases        = &_BasisVectors;
  updater._NumberOfDOFs = _NumberOfDOFs;
  updater._Input        = _Param;

  blocked_range<int> range(0, this->NumberOfCPs());
  parallel_for(range, updater);
}

// -----------------------------------------------------------------------------
struct BSplineFreeFormTransformationStatisticalDOFUpdate
{
  BSplineFreeFormTransformationStatistical          *_FFD;
  BSplineFreeFormTransformationStatistical::CPValue *_Data;

  const Vector *_Mean;
  const Matrix *_Bases;
  double       *_Output;
  int           _NumberOfCPs;

  void operator()(const blocked_range<int> &re) const
  {
    int dof_x, dof_y, dof_z;
    for (int q = re.begin(); q != re.end(); ++q) {
      for (int cp = 0; cp < _NumberOfCPs; cp++) {
        _FFD->IndexToDOFs(cp, dof_x, dof_y, dof_z);
        _Output[q] += _Bases->Get(dof_x, q) * (_Data[cp]._x - _Mean->Get(dof_x)) +
                      _Bases->Get(dof_y, q) * (_Data[cp]._y - _Mean->Get(dof_y)) +
                      _Bases->Get(dof_z, q) * (_Data[cp]._z - _Mean->Get(dof_z));
      }
    }
  }
};

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationStatistical::UpdateDOFs()
{
  // Clear DOFs first
  memset(_Param, 0, _NumberOfDOFs * sizeof(DOFValue));

  BSplineFreeFormTransformationStatisticalDOFUpdate updater;
  updater._FFD          = this;
  updater._Data         = _CPImage.Data();
  updater._Mean         = &_MeanVector;
  updater._Bases        = &_BasisVectors;
  updater._NumberOfCPs  = this->NumberOfCPs();
  updater._Output       = _Param;

  blocked_range<int> range(0, this->NumberOfDOFs());
  parallel_for(range, updater);
}

// =============================================================================
// Parameters (non-DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool BSplineFreeFormTransformationStatistical::Set(const char *name, const char *value)
{
  if (strcmp(name, "Statistical deformation model file") == 0) {
    ReadSDM(value);
    return true;
  }
  return Transformation::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList BSplineFreeFormTransformationStatistical::Parameter() const
{
  ParameterList params = BSplineFreeFormTransformation3D::Parameter();
  if (!_ModelFile.empty()) {
    Insert(params, "Statistical deformation model file", _ModelFile);
  }
  return params;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
struct BSplineFreeFormTransformationStatisticalCPGradientProjection
{
  const Matrix *_Bases;
  int           _NumberOfCPs;
  double       *_Input;
  double       *_Output;

  void operator()(const blocked_range<int> &re) const
  {
    for (int q = re.begin(); q != re.end(); ++q) {
      for (int p = 0; p != 3 * _NumberOfCPs; ++p) {
        _Output[q] += _Input[p] * _Bases->Get(p, q);
      }
    }
  }
};

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationStatistical
::ParametricGradient(const GenericImage<double> *in,  double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     double t0, double w) const
{
  // Compute parametric gradient w.r.t control point coefficients
  double *tmp_gradient = CAllocate<double>(3 * this->NumberOfCPs());
  memset(tmp_gradient, 0, 3 * this->NumberOfCPs() * sizeof(double));
  BSplineFreeFormTransformation3D::ParametricGradient(in, tmp_gradient, i2w, wc, t0, w);

  // Apply chain rule to obtain gradient w.r.t statistical parameters
  BSplineFreeFormTransformationStatisticalCPGradientProjection proj;
  proj._Bases       = &_BasisVectors;
  proj._NumberOfCPs = this->NumberOfCPs();
  proj._Input       = tmp_gradient;
  proj._Output      = out;

  blocked_range<int> range(0, this->NumberOfDOFs());
  parallel_for(range, proj);

  Deallocate(tmp_gradient);
}

// =============================================================================
// Properties
// =============================================================================

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationStatistical
::BendingEnergyGradient(double *gradient, double w, bool incl_passive, bool wrt_world, bool use_spacing) const
{
  // Compute bending gradient w.r.t control point coefficients
  double *tmp_gradient = CAllocate<double>(3 * this->NumberOfCPs());
  memset(tmp_gradient, 0, 3 * this->NumberOfCPs() * sizeof(double));
  BSplineFreeFormTransformation3D::BendingEnergyGradient(tmp_gradient, w, incl_passive, wrt_world, use_spacing);

  // Apply chain rule to obtain gradient w.r.t statistical parameters
  BSplineFreeFormTransformationStatisticalCPGradientProjection proj;
  proj._Bases       = &_BasisVectors;
  proj._NumberOfCPs = this->NumberOfCPs();
  proj._Input       = tmp_gradient;
  proj._Output      = gradient;

  blocked_range<int> range(0, this->NumberOfDOFs());
  parallel_for(range, proj);

  Deallocate(tmp_gradient);
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationStatistical::Print(ostream &os, Indent indent) const
{
  os << indent << "3D Statistical B-spline FFD:" << endl;
  FreeFormTransformation3D::Print(os, indent + 1);
}

// -----------------------------------------------------------------------------
Cofstream &BSplineFreeFormTransformationStatistical::Write(Cofstream &to) const
{
  // Write magic no. for transformations
  unsigned int magic_no = TRANSFORMATION_MAGIC;
  to.WriteAsUInt(&magic_no, 1);

  // Write transformation type (write it as a 3D Bspline FFD)
  unsigned int trans_type = TRANSFORMATION_BSPLINE_FFD_3D_v3;
  to.WriteAsUInt(&trans_type, 1);

  return this->WriteDOFs(to);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationStatistical::ReadSDM(const char *file)
{
  char type[8], version[5];

  Cifstream from;
  from.Open(file);

  // Read header with file format version information
  from.ReadAsChar(type,    8);
  from.ReadAsChar(version, 5);

  if (strncmp(type, "SDM", 8) != 0) {
    cerr << "The file '" << file << "' is not a valid statistical deformation model!" << endl;
    exit(1);
  }

  // Read lattice attributes
  ImageAttributes attr;

  from.ReadAsInt(&attr._x, 1);
  from.ReadAsInt(&attr._y, 1);
  from.ReadAsInt(&attr._z, 1);
  attr._t = 1;

  from.ReadAsDouble(attr._xaxis, 3);
  from.ReadAsDouble(attr._yaxis, 3);
  from.ReadAsDouble(attr._zaxis, 3);

  from.ReadAsDouble(&attr._dx, 1);
  from.ReadAsDouble(&attr._dy, 1);
  from.ReadAsDouble(&attr._dz, 1);
  attr._dt = 0;

  from.ReadAsDouble(&attr._xorigin, 1);
  from.ReadAsDouble(&attr._yorigin, 1);
  from.ReadAsDouble(&attr._zorigin, 1);

  // Initialize control points whith dofs flag in false, i.e.,
  // specifying that in our case CPs are *not* DOFs
  InitializeCPs(attr, false);

  // Read status
  from.ReadAsInt(reinterpret_cast<int *>(_CPStatus[0][0][0]),
                 3 * attr._x * attr._y * attr._z);

  // Read matrix and vector
  from >> _BasisVectors >> _MeanVector;

  if (_BasisVectors.Rows() != _MeanVector.Rows()) {
    cerr << "Invalid statistical deformation model! "
                 "Basis vectors and mean vector must have the same number of elements." << endl;
    exit(1);
  }

  if (_MeanVector.Rows() != 3 * attr._x * attr._y * attr._z) {
    cerr << "Incompatible statistical deformation model! "
                 "The number of vector elements must be equal to 3 times the number of control points." << endl;
    exit(1);
  }

  // Reinitialize DOFs now that we know how many of them we need
  InitializeDOFs(_BasisVectors.Cols());

  //Initialize the control points in order to be consistent with the newly read SDM
  UpdateCPs();

  // Initialize interpolator
  InitializeInterpolator();

  // Remember model file name
  ModelFile(file);
}

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationStatistical::WriteSDM(const char *file)
{
  Cofstream to;
  to.Open(file);

  // Write header with file format version information
  to.WriteAsChar("SDM", 8);
  to.WriteAsChar("v1.0",    5);

  // Write lattice attributes
  to.WriteAsInt(&_attr._x, 1);
  to.WriteAsInt(&_attr._y, 1);
  to.WriteAsInt(&_attr._z, 1);

  to.WriteAsDouble(_attr._xaxis, 3);
  to.WriteAsDouble(_attr._yaxis, 3);
  to.WriteAsDouble(_attr._zaxis, 3);

  to.WriteAsDouble(&_attr._dx, 1);
  to.WriteAsDouble(&_attr._dy, 1);
  to.WriteAsDouble(&_attr._dz, 1);

  to.WriteAsDouble(&_attr._xorigin, 1);
  to.WriteAsDouble(&_attr._yorigin, 1);
  to.WriteAsDouble(&_attr._zorigin, 1);

  to.WriteAsInt(reinterpret_cast<const int *>(_CPStatus[0][0][0]),
                3 * _attr._x * _attr._y * _attr._z);

  // Write statistical deformation model data
  to << _BasisVectors << _MeanVector;
}

// -----------------------------------------------------------------------------
Cifstream &BSplineFreeFormTransformationStatistical
::ReadDOFs(Cifstream &from, TransformationType format)
{
  from = BSplineFreeFormTransformation3D::ReadDOFs(from, format);
  UpdateDOFs();
  return from;
}

// =============================================================================
// Others
// =============================================================================

// -----------------------------------------------------------------------------
void BSplineFreeFormTransformationStatistical::Verify()
{
  if (_NumberOfDOFs == 0) {
    cerr << "No statistical deformation model was provided! "
                 "Make sure 'Statistical deformation model file' is set on the parameter file" << endl;
    exit(1);
  }

  if (_BasisVectors.Norm() == 0) {
    cerr << "Invalid statistical deformation model: "
                 "All basis verctors have zero norm!" << endl;
    exit(1);
  }
}


} // namespace mirtk

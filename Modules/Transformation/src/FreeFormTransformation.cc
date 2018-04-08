/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2017 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
 * Copyright 2013-2017 Andreas Schuh
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

#include "mirtk/FreeFormTransformation.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// Auxiliary macro used to initialize references to often needed
// control point image attributes in initialization list of constructor
#define _References()                                                          \
  _attr  (_CPImage.Attributes()),                                              \
  _x     (_CPImage.Attributes()._x),                                           \
  _y     (_CPImage.Attributes()._y),                                           \
  _z     (_CPImage.Attributes()._z),                                           \
  _t     (_CPImage.Attributes()._t),                                           \
  _dx    (_CPImage.Attributes()._dx),                                          \
  _dy    (_CPImage.Attributes()._dy),                                          \
  _dz    (_CPImage.Attributes()._dz),                                          \
  _dt    (_CPImage.Attributes()._dt),                                          \
  _matW2L(_CPImage.GetWorldToImageMatrix()),                                   \
  _matL2W(_CPImage.GetImageToWorldMatrix())

// -----------------------------------------------------------------------------
FreeFormTransformation
::FreeFormTransformation(CPInterpolator &func, CPInterpolator *func2D)
:
  _ExtrapolationMode(Extrapolation_Const),
  _SpeedupFactor(1.0),
  _CPValue (NULL),
  _CPFunc  (&func),
  _CPFunc2D(func2D),
  _CPStatus(NULL),
  _References()
{
}

// -----------------------------------------------------------------------------
FreeFormTransformation
::FreeFormTransformation(const FreeFormTransformation &ffd,
                         CPInterpolator &func, CPInterpolator *func2D)
:
  Transformation(ffd),
  _ExtrapolationMode(ffd._ExtrapolationMode),
  _SpeedupFactor(ffd._SpeedupFactor),
  _CPValue (NULL),
  _CPFunc  (&func),
  _CPFunc2D(func2D),
  _CPStatus(NULL),
  _References()
{
  if (_NumberOfDOFs > 0) {
    _CPImage.Initialize(ffd.Attributes(), reinterpret_cast<CPValue  *>(_Param));
    Allocate(_CPStatus, _x, _y, _z, _t,   reinterpret_cast<CPStatus *>(_Status));
  }
}

// -----------------------------------------------------------------------------
FreeFormTransformation::~FreeFormTransformation()
{
  Deallocate(_CPStatus, _Status);
  Delete(_CPValue);
}


#undef _References // only used/needed for constructors

// -----------------------------------------------------------------------------
ImageAttributes
FreeFormTransformation
::DefaultAttributes(const ImageAttributes &attr,
                    double dx, double dy, double dz, double dt)
{
  // Set spacing equal to voxel size if not specified (i.e., negative)
  if (dx < .0) dx = attr._dx;
  if (dy < .0) dy = attr._dy;
  if (dz < .0) dz = attr._dz;
  if (dt < .0) dt = attr._dt;

  // Orthogonalize image grid if necessary
  const ImageAttributes grid = OrthogonalFieldOfView(attr);

  // Set lattice attributes to image attributes
  ImageAttributes lattice = grid;

  // Adjust lattice dimensions
  lattice._x = ((dx > .0 && grid._x > 1) ? iround((grid._x - 1) * grid._dx / dx) + 1 : 1);
  lattice._y = ((dy > .0 && grid._y > 1) ? iround((grid._y - 1) * grid._dy / dy) + 1 : 1);
  lattice._z = ((dz > .0 && grid._z > 1) ? iround((grid._z - 1) * grid._dz / dz) + 1 : 1);
  lattice._t = ((dt > .0 && grid._t > 1) ? iround((grid._t - 1) * grid._dt / dt) + 1 : 1);

  // Update lattice spacing
  lattice._dx = ((lattice._x > 1) ? dx : .0);
  lattice._dy = ((lattice._y > 1) ? dy : .0);
  lattice._dz = ((lattice._z > 1) ? dz : .0);
  lattice._dt = ((lattice._t > 1) ? dt : .0);

  return lattice;
}

// -----------------------------------------------------------------------------
ImageAttributes
FreeFormTransformation
::DefaultAttributes(double x1, double y1, double z1, double t1,
                    double x2, double y2, double z2, double t2,
                    double dx, double dy, double dz, double dt,
                    const double *xaxis,
                    const double *yaxis,
                    const double *zaxis)
{
  ImageAttributes attr;

  // Lattice origin
  attr._xorigin = (x2 + x1) / 2.0;
  attr._yorigin = (y2 + y1) / 2.0;
  attr._zorigin = (z2 + z1) / 2.0;
  attr._torigin = t1;

  // Lattice orientation
  memcpy(attr._xaxis, xaxis, 3 * sizeof(double));
  memcpy(attr._yaxis, yaxis, 3 * sizeof(double));
  memcpy(attr._zaxis, zaxis, 3 * sizeof(double));

  // FOV in lattice orientation
  double a, b, c;
  a  = x1 * xaxis[0] + y1 * xaxis[1] + z1 * xaxis[2];
  b  = x1 * yaxis[0] + y1 * yaxis[1] + z1 * yaxis[2];
  c  = x1 * zaxis[0] + y1 * zaxis[1] + z1 * zaxis[2];
  x1 = a, y1 = b, z1 = c;
  a  = x2 * xaxis[0] + y2 * xaxis[1] + z2 * xaxis[2];
  b  = x2 * yaxis[0] + y2 * yaxis[1] + z2 * yaxis[2];
  c  = x2 * zaxis[0] + y2 * zaxis[1] + z2 * zaxis[2];
  x2 = a, y2 = b, z2 = c;

  // Lattice dimensions
  attr._x = ((dx > .0) ? iround((x2 - x1) / dx) + 1 : 1);
  attr._y = ((dy > .0) ? iround((y2 - y1) / dy) + 1 : 1);
  attr._z = ((dz > .0) ? iround((z2 - z1) / dz) + 1 : 1);
  attr._t = ((dt > .0) ? iround((t2 - t1) / dt) + 1 : 1);

  // Lattice spacing
  attr._dx = ((attr._x > 1) ? dx : .0);
  attr._dy = ((attr._y > 1) ? dy : .0);
  attr._dz = ((attr._z > 1) ? dz : .0);
  attr._dt = ((attr._t > 1) ? dt : .0);

  return attr;
}

// -----------------------------------------------------------------------------
void FreeFormTransformation::InitializeInterpolator()
{
  // Attention: Returns NULL for _ExtrapolationMode == Extrapolation_None!
  Delete(_CPValue);
  _CPValue = CPExtrapolator::New(_ExtrapolationMode, &_CPImage);
  _CPFunc->Input(&_CPImage);
  _CPFunc->Extrapolator(_CPValue);
  _CPFunc->Initialize(true); // also initializes _CPValue
  if (_CPFunc2D) {
    _CPFunc2D->Input(&_CPImage);
    _CPFunc2D->Extrapolator(_CPValue);
    _CPFunc2D->Initialize(true);
  }
}

// -----------------------------------------------------------------------------
void FreeFormTransformation::InitializeCPs(const ImageAttributes &attr, bool dofs)
{
  if (attr._x > 1000 || attr._y > 1000 || attr._z > 1000) {
    // Otherwise we may not only run out of memory, but also encounter an
    // integer overflow when computing the number of control points...
    cerr << "FreeFormTransformation::InitializeCPs: Number of control points may not exceed 1000 in each dimension!" << endl;
    exit(1);
  }
  Deallocate(_CPStatus, _Status);
  const int ncps = attr.NumberOfPoints();
  if (ncps > 0) {
    // If NULL, allocate memory for control points which differs from
    // the memory referenced by _Param and _Status of the base class
    CPValue  *param  = NULL;
    CPStatus *status = NULL;
    // If control points are directly the parameters of the transformation...
    if (dofs) {
      // Assert size of types
      if (sizeof(CPValue) != 3 * sizeof(DOFValue)) {
        cerr << "FreeFormTransformation::InitializeCPs: Assertion \"sizeof(CPValue) == 3 * sizeof(DOFValue)\" failed!" << endl;
        cerr << "    Make sure that the vector type used fulfills this assertion such that FFDs can wrap the linear" << endl;
        cerr << "    memory used to store the transformation parameters in an image instance with 3D vector voxel type." << endl;
        exit(1);
      }
      // ...allocate memory for parameters and their respective status
      InitializeDOFs(3 * ncps);
      // ...and reinterpret this memory as 3D(+t) images
      param  = reinterpret_cast<CPValue  *>(_Param);
      status = reinterpret_cast<CPStatus *>(_Status);
    }
    // Allocate 4D array for convenient access to parameters of each control point
    _CPImage.Initialize(attr, param);
    Allocate(_CPStatus, attr._x, attr._y, attr._z, attr._t, status);
    InitializeInterpolator();
    InitializeStatus();
  } else {
    Delete(_CPValue);
    _CPImage.Clear();
    InitializeDOFs(0);
  }
}

// -----------------------------------------------------------------------------
void FreeFormTransformation::InitializeCPs(const FreeFormTransformation &ffd, bool dofs)
{
  if (ffd.NumberOfCPs() > 0) {
    // If NULL, allocate memory for control points which differs from
    // the memory referenced by _Param and _Status of the base class
    CPValue  *param  = NULL;
    CPStatus *status = NULL;
    // If control points are directly the parameters of the transformation...
    if (dofs) {
      // ...copy transformation parameters and their respective status
      InitializeDOFs(ffd);
      // ...and reinterpret this memory as 3D(+t) images
      param  = reinterpret_cast<CPValue  *>(_Param);
      status = reinterpret_cast<CPStatus *>(_Status);
    }
    // Allocate 4D array for convenient access to parameters of each control point
    _CPImage.Initialize(ffd.Attributes(), param);
    Allocate(_CPStatus, _x, _y, _z, _t,   status);
    InitializeInterpolator();
    InitializeStatus();
  } else {
    Delete(_CPValue);
    _CPImage.Clear();
    Deallocate(_CPStatus, _Status);
    InitializeDOFs(0);
  }
}

// -----------------------------------------------------------------------------
void FreeFormTransformation::InitializeStatus()
{
  // Set status of unused control point dimensions to Passive
  // TODO: Consider removing this commented code or uncomment it again once
  //       the control point coefficients are in lattice units/orientation.
  //       A 2D FFD in a plane not parallel to the xy plane has non-zero
  //       (i.e., active) z coefficients in world space!
//  const CPStatus status(((_x > 1) ? Active : Passive),
//                        ((_y > 1) ? Active : Passive),
//                        ((_z > 1) ? Active : Passive));
  const CPStatus status = Active;
  for (int cp = 0; cp < this->NumberOfCPs(); ++cp) PutStatus(cp, status);
}

// -----------------------------------------------------------------------------
void FreeFormTransformation::Initialize(const ImageAttributes &attr)
{
  InitializeCPs(attr);
  this->Changed(true);
}

// -----------------------------------------------------------------------------
void FreeFormTransformation::Initialize(const ImageAttributes &attr,
                                        double dx, double dy, double dz, double dt)
{
  this->Initialize(DefaultAttributes(attr, dx, dy, dz, dt));
}

// -----------------------------------------------------------------------------
void FreeFormTransformation::Initialize(const CPImage &image, bool displacement)
{
  // Initialize free-form transformation
  this->Initialize(image.Attributes());
  // Copy control point coefficients
  _CPImage.CopyFrom(image);
  // Convert control point deformation to displacement
  if (!displacement) {
    double x, y, z;
    CPValue *data = _CPImage.Data();
    for (int l = 0; l < _t; ++l)
    for (int k = 0; k < _z; ++k)
    for (int j = 0; j < _y; ++j)
    for (int i = 0; i < _x; ++i) {
      x = i, y = j, z = k;
      this->LatticeToWorld(x, y, z);
      data->_x -= x;
      data->_y -= y;
      data->_z -= z;
      ++data;
    }
  }
}

// -----------------------------------------------------------------------------
void FreeFormTransformation
::Initialize(const GenericImage<double> &image, bool displacement)
{
  if (image.T() < 2 || image.T() > 3) {
    cerr << "FreeFormTransformation::Initialize:"
            " Input image must be 2D/3D vector field with _t == 2 or 3)" << endl;
    exit(1);
  }
  // Initialize control points
  ImageAttributes attr = image.Attributes();
  attr._t  = 1;
  attr._dt = .0;
  this->Initialize(attr);
  // Copy vector field
  CPValue *data = _CPImage.Data();
  for (int k = 0; k < image.GetZ(); ++k)
  for (int j = 0; j < image.GetY(); ++j)
  for (int i = 0; i < image.GetX(); ++i) {
    data->_x = image(i, j, k, 0);
    data->_y = image(i, j, k, 1);
    if (image.GetT() >= 3) data->_z = image(i, j, k, 2);
    else                   data->_z = .0;
    ++data;
  }
  // Convert control point deformation to displacement
  if (!displacement) {
    double x, y, z;
    CPValue *data = _CPImage.Data();
    for (int l = 0; l < _t; ++l)
    for (int k = 0; k < _z; ++k)
    for (int j = 0; j < _y; ++j)
    for (int i = 0; i < _x; ++i) {
      x = i, y = j, z = k;
      this->LatticeToWorld(x, y, z);
      data->_x -= x;
      data->_y -= y;
      data->_z -= z;
      ++data;
    }
  }
}

// =============================================================================
// Parameters (non-DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
void FreeFormTransformation::ExtrapolationMode(enum ExtrapolationMode mode)
{
  // Set extrapolation mode
  _ExtrapolationMode = mode;
  // Update control point function, if _CPValue == NULL the FFD is uninitialized yet
  if (_CPValue && _CPValue->ExtrapolationMode() != _ExtrapolationMode) {
    Delete(_CPValue);
    _CPValue = CPExtrapolator::New(_ExtrapolationMode, &_CPImage);
    if (_CPValue) {
      _CPValue->Input(&_CPImage);
      _CPValue->Initialize();
    }
    _CPFunc->Extrapolator(_CPValue);
    if (_CPFunc2D) _CPFunc2D->Extrapolator(_CPValue);
    // Note: No need to call _CPFunc(2D)->Initialize() again
  }
}

// -----------------------------------------------------------------------------
bool FreeFormTransformation::Set(const char *name, const char *value)
{
  // Extrapolation mode
  if (strcmp(name, "Transformation extrapolation") == 0) {
    enum ExtrapolationMode mode;
    if (!FromString(value, mode)) return false;
    this->ExtrapolationMode(mode);
    this->Changed(true);
    return true;
  }
  // Speedup factor for gradient computation
  if (strcmp(name, "Speedup factor") == 0) {
    double d;
    if (!FromString(value, d) || d < 1.0) return false;
    _SpeedupFactor = d;
    this->Changed(true);
    return true;
  }
  // Control point spacing
  if (strncmp(name, "Control point spacing", 21) == 0) {
    double dx, dy, dz, dt, ds;
    if (!FromString(value, ds) || ds <= .0) return false;
    if      (strstr(name, "[m]")  != NULL) ds *= 1.0e+3;
    else if (strstr(name, "[mu]") != NULL) ds *= 1.0e-3;
    _CPImage.GetPixelSize(dx, dy, dz, dt);
    if (strncmp(name + 21, " in ", 4) == 0) {
      if      (name[25] == 'X') dx = ds;
      else if (name[25] == 'Y') dy = ds;
      else if (name[25] == 'Z') dz = ds;
      else if (name[25] == 'T') dt = ds;
    } else {
      dx = dy = dz = dt = ds;
    }
    if (_x == 1) dx = .0;
    if (_y == 1) dy = .0;
    if (_z == 1) dz = .0;
    if (_t == 1) dt = .0;
    _CPImage.PutPixelSize(dx, dy, dz, dt);
    this->Changed(true);
    return true;
  }
  if (strcmp(name, "Control point xyz units") == 0) {
    double s;
    if      (strcmp(value, "Meter")      == 0 || strcmp(value, "m")  == 0) s = 1.0e+3;
    else if (strcmp(value, "Millimeter") == 0 || strcmp(value, "mm") == 0) s = 1.0;
    else if (strcmp(value, "Micrometer") == 0 || strcmp(value, "mu") == 0) s = 1.0e-3;
    else return false;
    if (s != 1.0) {
      double ox, oy, oz, dx, dy, dz;
      _CPImage.GetOrigin   (ox, oy, oz);
      _CPImage.PutOrigin   (ox * s, oy * s, oz * s);
      _CPImage.GetPixelSize(dx, dy, dz);
      _CPImage.PutPixelSize(dx * s, dy * s, dz * s);
      _CPImage *= s;
    }
    this->Changed(true);
    return true;
  }
  return Transformation::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList FreeFormTransformation::Parameter() const
{
  ParameterList params = Transformation::Parameter();
  Insert(params, "Transformation extrapolation", ToString(_ExtrapolationMode));
  Insert(params, "Speedup factor",               ToString(_SpeedupFactor));
  return params;
}

// =============================================================================
// Approximation/Interpolation
// =============================================================================

// -----------------------------------------------------------------------------
ImageAttributes FreeFormTransformation
::ApproximationDomain(const ImageAttributes &domain, const Transformation *)
{
  return domain;
}

// -----------------------------------------------------------------------------
double FreeFormTransformation::EvaluateRMSError(const Transformation *dof) const
{
  return EvaluateRMSError(_attr, dof);
}

// -----------------------------------------------------------------------------
double FreeFormTransformation::Approximate(const Transformation *t, int niter, double max_error)
{
  return this->Approximate(_attr, t, niter, max_error);
}

// -----------------------------------------------------------------------------
double FreeFormTransformation
::Approximate(const ImageAttributes &domain, double *dx, double *dy, double *dz,
              int niter, double max_error)
{
  const int no = domain.NumberOfPoints();
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Compute world coordinates of lattice points
  double *x = Allocate<double>(no);
  double *y = Allocate<double>(no);
  double *z = Allocate<double>(no);
  double *t = Allocate<double>(no);
  domain.LatticeToWorld(x, y, z, t);

  // Copy original displacements
  double *tx = CAllocate<double>(no);
  double *ty = CAllocate<double>(no);
  double *tz = CAllocate<double>(no);
  memcpy(tx, dx, no * sizeof(double));
  memcpy(ty, dy, no * sizeof(double));
  memcpy(tz, dz, no * sizeof(double));

  // Evaluate error of approximation and residual displacements
  if (this->RequiresCachingOfDisplacements()) {
    error = EvaluateRMSError(domain, dx, dy, dz);
  } else {
    error = EvaluateRMSError(x, y, z, t, dx, dy, dz, no);
  }

  // Repeat approximation n times or until error drops below a threshold
  DOFValue *param = Allocate<DOFValue>(_NumberOfDOFs);
  for (int iter = 0; iter < niter && error > max_error; ++iter) {

    // Copy current parameters
    this->Get(param);

    // Approximate residual displacements by new parameters
    this->ApproximateDOFs(x, y, z, t, dx, dy, dz, no);

    // Add previous parameters
    this->Add(param);

    // Evaluate error of approximation and residual displacements
    memcpy(dx, tx, no * sizeof(double));
    memcpy(dy, ty, no * sizeof(double));
    memcpy(dz, tz, no * sizeof(double));
    if (this->RequiresCachingOfDisplacements()) {
      error = EvaluateRMSError(domain, dx, dy, dz);
    } else {
      error = EvaluateRMSError(x, y, z, t, dx, dy, dz, no);
    }
  }

  // Free memory
  Deallocate(param);
  Deallocate(tx);
  Deallocate(ty);
  Deallocate(tz);
  Deallocate(x);
  Deallocate(y);
  Deallocate(z);
  Deallocate(t);

  return error;
}

// -----------------------------------------------------------------------------
double FreeFormTransformation
::Approximate(const double *x,  const double *y,  const double *z,
              double       *dx, double       *dy, double       *dz, int no,
              int niter, double max_error)
{
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Compute (fixed) time coordinates
  const double t0 = .0;
  double *t = CAllocate<double>(no, &t0);

  // Copy original displacements
  double *tx = Allocate<double>(no);
  double *ty = Allocate<double>(no);
  double *tz = Allocate<double>(no);
  memcpy(tx, dx, no * sizeof(double));
  memcpy(ty, dy, no * sizeof(double));
  memcpy(tz, dz, no * sizeof(double));

  // Evaluate error of approximation and residual displacements
  error = EvaluateRMSError(x, y, z, t0, dx, dy, dz, no);

  // Repeat approximation n times or until error drops below a threshold
  DOFValue *param = Allocate<DOFValue>(_NumberOfDOFs);
  for (int iter = 0; iter < niter && error > max_error; ++iter) {

    // Copy current parameters
    this->Get(param);

    // Approximate residual displacements by new parameters
    this->ApproximateDOFs(x, y, z, t, dx, dy, dz, no);

    // Add previous parameters
    this->Add(param);

    // Evaluate error of approximation and residual displacements
    memcpy(dx, tx, no * sizeof(double));
    memcpy(dy, ty, no * sizeof(double));
    memcpy(dz, tz, no * sizeof(double));
    error = EvaluateRMSError(x, y, z, t0, dx, dy, dz, no);
  }

  // Free memory
  Deallocate(param);
  Deallocate(tx);
  Deallocate(ty);
  Deallocate(tz);
  Deallocate(t);

  return error;
}

// -----------------------------------------------------------------------------
double FreeFormTransformation
::Approximate(const double *x,  const double *y,  const double *z,  const double *t,
              double       *dx, double       *dy, double       *dz, int no,
              int niter, double max_error)
{
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Copy original displacements
  double *tx = Allocate<double>(no);
  double *ty = Allocate<double>(no);
  double *tz = Allocate<double>(no);
  memcpy(tx, dx, no * sizeof(double));
  memcpy(ty, dy, no * sizeof(double));
  memcpy(tz, dz, no * sizeof(double));

  // Evaluate error of approximation and residual displacements
  error = EvaluateRMSError(x, y, z, t, dx, dy, dz, no);

  // Repeat approximation n times or until error drops below a threshold
  DOFValue *param = Allocate<DOFValue>(_NumberOfDOFs);
  for (int iter = 0; iter < niter && error > max_error; ++iter) {

    // Copy current parameters
    this->Get(param);

    // Approximate residual displacements by new parameters
    this->ApproximateDOFs(x, y, z, t, dx, dy, dz, no);

    // Add previous parameters
    this->Add(param);

    // Evaluate error of approximation and residual displacements
    memcpy(dx, tx, no * sizeof(double));
    memcpy(dy, ty, no * sizeof(double));
    memcpy(dz, tz, no * sizeof(double));
    error = EvaluateRMSError(x, y, z, t, dx, dy, dz, no);
  }

  // Free memory
  Deallocate(param);
  Deallocate(tx);
  Deallocate(ty);
  Deallocate(tz);

  return error;
}

// -----------------------------------------------------------------------------
double FreeFormTransformation::ApproximateAsNew(const Transformation *t, int niter, double max_error)
{
  return this->ApproximateAsNew(_attr, t, niter, max_error);
}

// -----------------------------------------------------------------------------
double FreeFormTransformation
::ApproximateAsNew(const ImageAttributes &domain, double *dx, double *dy, double *dz, int niter, double max_error)
{
  bool ignore_time = (AreEqual(domain._dt, 0.) || domain._t == 1) && (AreEqual(_dt, 0.) || _t == 1);
  if ((!ignore_time && domain == _attr) || (ignore_time && domain.EqualInSpace(_attr))) {
    this->Interpolate(dx, dy, dz);
    return .0; // No approximation error at lattice/control points
  }
  return Transformation::ApproximateAsNew(domain, dx, dy, dz, niter, max_error);
}

// =============================================================================
// Lattice
// =============================================================================

// -----------------------------------------------------------------------------
bool FreeFormTransformation::CropPadPassiveCPs(int margin, bool reset)
{
  return this->CropPadPassiveCPs(margin, margin, margin, margin, reset);
}

// -----------------------------------------------------------------------------
bool FreeFormTransformation::CropPadPassiveCPs(int mx, int my, int mz, int mt, bool reset)
{
  ImageAttributes attr = _attr;
  // Determine lower bound along x axis: i1
  int i1 = _x;
  for (int l = 0; l < _t; ++l)
  for (int k = 0; k < _z; ++k)
  for (int j = 0; j < _y; ++j)
  for (int i = 0; i < _x; ++i) {
    if (IsActive(i, j, k, l)) {
      if (i < i1) i1 = i;
      break;
    }
  }
  // Determine upper bound along x axis: i2
  int i2 = -1;
  for (int l = 0;      l <  _t; ++l)
  for (int k = 0;      k <  _z; ++k)
  for (int j = 0;      j <  _y; ++j)
  for (int i = _x - 1; i >= i1; --i) {
    if (IsActive(i, j, k, l)) {
      if (i > i2) i2 = i;
      break;
    }
  }
  // Determine lower bound along y axis: j1
  int j1 = _y;
  for (int l = 0;  l <  _t; ++l)
  for (int k = 0;  k <  _z; ++k)
  for (int i = i1; i <= i2; ++i)
  for (int j = 0;  j <  _y; ++j) {
    if (IsActive(i, j, k, l)) {
      if (j < j1) j1 = j;
      break;
    }
  }
  // Determine upper bound along y axis: j2
  int j2 = -1;
  for (int l = 0;      l <  _t; ++l)
  for (int k = 0;      k <  _z; ++k)
  for (int i = i1;     i <= i2; ++i)
  for (int j = _y - 1; j >= j1; --j) {
    if (IsActive(i, j, k, l)) {
      if (j > j2) j2 = j;
      break;
    }
  }
  // Determine lower bound along z axis: k1
  int k1 = _z;
  for (int l = 0;  l <  _t; ++l)
  for (int j = j1; j <= j2; ++j)
  for (int i = i1; i <= i2; ++i)
  for (int k = 0;  k <  _z; ++k) {
    if (IsActive(i, j, k, l)) {
      if (k < k1) k1 = k;
      break;
    }
  }
  // Determine upper bound along z axis: k2
  int k2 = -1;
  for (int l = 0;      l <  _t; ++l)
  for (int j = j1;     j <= j2; ++j)
  for (int i = i1;     i <= i2; ++i)
  for (int k = _z - 1; k >= k1; --k) {
    if (IsActive(i, j, k, l)) {
      if (k > k2) k2 = k;
       break;
    }
  }
  // Determine lower bound along t axis: l1
  int l1 = _t;
  for (int k = k1; k <= k2; ++k)
  for (int j = j1; j <= j2; ++j)
  for (int i = i1; i <= i2; ++i)
  for (int l = 0;  l <  _t; ++l) {
    if (IsActive(i, j, k, l)) {
      if (l < l1) l1 = l;
      break;
    }
  }
  // Determine upper bound along t axis: l2
  int l2 = -1;
  for (int k = k1;     k <= k2; ++k)
  for (int j = j1;     j <= j2; ++j)
  for (int i = i1;     i <= i2; ++i)
  for (int l = _t - 1; l >= l1; --l) {
    if (IsActive(i, j, k, l)) {
      if (l > l2) l2 = l;
      break;
    }
  }
  // Do nothing if all control points are passive, but report it
  if (i1 > i2 || j1 > j2 || k1 > k2 || l1 > l2) return false;
  // Convert upper index bounds to margin widths
  i2 = (_x - 1) - i2;
  j2 = (_y - 1) - j2;
  k2 = (_z - 1) - k2;
  l2 = (_t - 1) - l2;
  // Leave a margin of passive control points with specified width
  // Note: Negative value gives the number of control points to add.
  if (_x > 1) i1 -= mx, i2 -= mx;
  if (_y > 1) j1 -= my, j2 -= my;
  if (_z > 1) k1 -= mz, k2 -= mz;
  if (_t > 1) l1 -= mt, l2 -= mt;
  // Do nothing, if nothing to be done
  if (i1 == 0 && i2 == 0 && j1 == 0 && j2 == 0 &&
      k1 == 0 && k2 == 0 && l1 == 0 && l2 == 0) return true;
  // Adjust control point lattice
  attr._x -= i1 + i2;
  attr._y -= j1 + j2;
  attr._z -= k1 + k2;
  attr._t -= l1 + l2;
  attr._xorigin = 0.5 * ((_x - 1) + (i1 - i2));
  attr._yorigin = 0.5 * ((_y - 1) + (j1 - j2));
  attr._zorigin = 0.5 * ((_z - 1) + (k1 - k2));
  this->LatticeToWorld(attr._xorigin, attr._yorigin, attr._zorigin);
  attr._torigin = this->LatticeToTime(l1);
  // Convert upper margin widths to index bounds
  i2 = (_x - 1) - i2;
  j2 = (_y - 1) - j2;
  k2 = (_z - 1) - k2;
  l2 = (_t - 1) - l2;
  // Copy remaining control points and pad lattice where needed
  const int ncps = attr.NumberOfPoints();
  CPValue  *param       = Allocate<CPValue >(ncps);
  CPStatus *status      = Allocate<CPStatus>(ncps);
  CPValue  *param_iter  = param;
  CPStatus *status_iter = status;
  for (int l = l1; l <= l2; ++l)
  for (int k = k1; k <= k2; ++k)
  for (int j = j1; j <= j2; ++j)
  for (int i = i1; i <= i2; ++i, ++param_iter, ++status_iter) {
    if (0 <= i && i < _x &&
        0 <= j && j < _y &&
        0 <= k && k < _z &&
        0 <= l && l < _t) {
      (*status_iter) = _CPStatus[l][k][j][i];
      (*param_iter)  = _CPImage(i, j, k, l);
    } else {
      // Padded control point to extend margin
      (*status_iter) = CPStatus(Passive);
      (*param_iter)  = CPValue (.0);
    }
  }
  // Initialize new control point lattice
  this->Initialize(attr);
  param_iter  = param;
  status_iter = status;
  for (int l = 0; l < _t; ++l)
  for (int k = 0; k < _z; ++k)
  for (int j = 0; j < _y; ++j)
  for (int i = 0; i < _x; ++i, ++param_iter, ++status_iter) {
    _CPStatus[l][k][j][i] = (*status_iter);
    if (// Copy all control point values by default
        !reset ||
        // Always if control point status is not passive...
        (status_iter->_x != Passive) ||
        (status_iter->_y != Passive) ||
        (status_iter->_z != Passive) ||
        // ...or control point is not within the boundary margin
        i >= mx || i < _x - mx ||
        j >= my || j < _y - my ||
        k >= mz || k < _z - mz ||
        l >= mt || l < _t - mt) {
      _CPImage(i, j, k, l) = (*param_iter);
    }
  }
  // Free temporary allocated memory
  Deallocate(param);
  Deallocate(status);
  return true;
}

// =============================================================================
// Bounding box
// =============================================================================

// -----------------------------------------------------------------------------
void FreeFormTransformation::PutBoundingBox(double x1, double y1, double z1,
                                            double x2, double y2, double z2)
{
  const ImageAttributes &attr = Attributes();
  double a, b, c;

  // Update lattice origin
  _CPImage.PutOrigin((x2 + x1) / 2.0, (y2 + y1) / 2.0, (z2 + z1) / 2.0);

  // FOV in lattice orientation
  a = x1 * attr._xaxis[0] + y1 * attr._xaxis[1] + z1 * attr._xaxis[2];
  b = x1 * attr._yaxis[0] + y1 * attr._yaxis[1] + z1 * attr._yaxis[2];
  c = x1 * attr._zaxis[0] + y1 * attr._zaxis[1] + z1 * attr._zaxis[2];
  x1 = a, y1 = b, z1 = c;

  a = x2 * attr._xaxis[0] + y2 * attr._xaxis[1] + z2 * attr._xaxis[2];
  b = x2 * attr._yaxis[0] + y2 * attr._yaxis[1] + z2 * attr._yaxis[2];
  c = x2 * attr._zaxis[0] + y2 * attr._zaxis[1] + z2 * attr._zaxis[2];
  x2 = a, y2 = b, z2 = c;

  // Update lattice spacing
  _CPImage.PutPixelSize(((x2 > x1) ? ((x2 - x1) / (attr._x - 1)) : 1),
                        ((y2 > y1) ? ((y2 - y1) / (attr._y - 1)) : 1),
                        ((z2 > z1) ? ((z2 - z1) / (attr._z - 1)) : 1));

  this->Changed(true);
}

// -----------------------------------------------------------------------------
void FreeFormTransformation::PutBoundingBox(const Point &p1, const Point &p2)
{
  PutBoundingBox(p1._x, p1._y, p1._z, p2._x, p2._y, p2._z);
}

// ---------------------------------------------------------------------------
void FreeFormTransformation::PutBoundingBox(double t1, double t2)
{
  _CPImage.PutTOrigin((t2 + t1) / 2.0);
  _CPImage.PutTSize  ((t2 > t2) ? ((t2 - t1) / (_t - 1)) : 1);
  this->Changed(true);
}

// ---------------------------------------------------------------------------
void FreeFormTransformation::PutBoundingBox(double x1, double y1, double z1, double t1,
                                            double x2, double y2, double z2, double t2)
{
  PutBoundingBox(x1, y1, z1, x2, y2, z2);
  PutBoundingBox(t1, t2);
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
/// (Multi-threaded) body of FreeFormTransformation::ParametricGradient
///
/// This default implementation uses the FreeFormTransformation::JacobianDOFs
/// overload which computes the derivatives of the transformation w.r.t. all
/// transformation parameters at once for a given spatial location (target voxel).
/// It provides better performance in particular for transformations which are
/// parameterized by non-stationary velocity fields (3D+t) in which case the
/// derivatives along each temporal trajectory are computed at once and therefore
/// this trajectory does not have to be re-computed for each and every control
/// point. For classic FFDs, which are parameterized by displacements, less
/// general but more efficient implementations of the ParametricGradient function
/// are provided by the FreeFormTransformation3D and
/// FreeFormTransformation4D subclasses, which therefore override this
/// FreeFormTransformation::ParametricGradient implementation. If a subclass
/// of either of these subclasses is parameterized by velocities, it must therefore
/// override the ParametericGradient function again, for example as follows,
/// unless a more specialized gradient computation is provided by this subclass.
///
/// \code
/// void BSplineFreeFormTransformationTD::ParametricGradient(...)
/// {
///   FreeFormTransformation::ParametricGradient(...);
/// }
/// \endcode
class FreeFormTransformationParametricGradientBody
{
public:
  const FreeFormTransformation *_FFD;
  const GenericImage<double>   *_Input;
  const WorldCoordsImage       *_WorldCoords;
  double                       *_Output;
  double                        _Weight;

  double _t;  ///< Time corrresponding to input gradient image (in ms)
  double _t0; ///< Second time argument for velocity-based transformations

private:

  int    _X;  ///< Number of voxels along x axis
  int    _Y;  ///< Number of voxels along y axis

public:

  // ---------------------------------------------------------------------------
  /// Default constructor
  FreeFormTransformationParametricGradientBody()
  :
    _FFD        (NULL),
    _Input      (NULL),
    _WorldCoords(NULL),
    _Output     (NULL),
    _Weight     ( 1.0),
    _t          ( 0.0),
    _t0         (-1.0),
    _X          (0),
    _Y          (0)
  {}

  // ---------------------------------------------------------------------------
  /// Split constructor
  FreeFormTransformationParametricGradientBody(const FreeFormTransformationParametricGradientBody &other, split)
  :
    _FFD        (other._FFD),
    _Input      (other._Input),
    _WorldCoords(other._WorldCoords),
    _Output     (NULL),
    _Weight     (other._Weight),
    _t          (other._t),
    _t0         (other._t0),
    _X          (other._X),
    _Y          (other._Y)
  {
    const int ndofs = _FFD->NumberOfDOFs();
    _Output         = new double[ndofs];
    memset(_Output, 0, ndofs * sizeof(double));
  }

  // ---------------------------------------------------------------------------
  /// Destructor
  ~FreeFormTransformationParametricGradientBody() {}

  // ---------------------------------------------------------------------------
  /// Join results of right-hand body with this body
  void join(FreeFormTransformationParametricGradientBody &rhs)
  {
    const int ndofs = _FFD->NumberOfDOFs();
    for (int dof = 0; dof < ndofs; dof++) _Output[dof] += rhs._Output[dof];
    delete[] rhs._Output;
    rhs._Output = NULL;
  }

  // ---------------------------------------------------------------------------
  /// Calculates the gradient of the similarity term w.r.t. the transformation
  /// parameters for each voxel in the specified image region
  void operator ()(const blocked_range3d<int> &re)
  {
    const int i1 = re.cols ().begin();
    const int j1 = re.rows ().begin();
    const int k1 = re.pages().begin();
    const int i2 = re.cols ().end();
    const int j2 = re.rows ().end();
    const int k2 = re.pages().end();

    TransformationJacobian                      jac;
    TransformationJacobian::ConstColumnIterator it;

    //s1=1
    const int s2 =  _X - (i2 - i1);
    const int s3 = (_Y - (j2 - j1)) * _X;

    // Non-parametric/voxelwise gradient
    const double *gx = _Input->Data(i1, j1, k1, 0);
    const double *gy = _Input->Data(i1, j1, k1, 1);
    const double *gz = _Input->Data(i1, j1, k1, 2);

    // With (transformed) pre-computed world coordinates
    if (_WorldCoords) {
      const double *wx = _WorldCoords->Data(i1, j1, k1, 0);
      const double *wy = _WorldCoords->Data(i1, j1, k1, 1);
      const double *wz = _WorldCoords->Data(i1, j1, k1, 2);
      for (int k = k1; k < k2; ++k, wx += s3, wy += s3, wz += s3, gx += s3, gy += s3, gz += s3)
      for (int j = j1; j < j2; ++j, wx += s2, wy += s2, wz += s2, gx += s2, gy += s2, gz += s2)
      for (int i = i1; i < i2; ++i, wx +=  1, wy +=  1, wz +=  1, gx +=  1, gy +=  1, gz +=  1) {
        if (*gx || *gy || *gz) {
          // Calculate derivatives of transformation w.r.t. the parameters
          _FFD->JacobianDOFs(jac, *wx, *wy, *wz, _t, _t0);
          // Apply chain rule to obtain similarity gradient w.r.t. the transformation parameters
          for (it = jac.Begin(); it != jac.End(); ++it) {
            _Output[it->first] += _Weight * (it->second._x * (*gx) +
                                             it->second._y * (*gy) +
                                             it->second._z * (*gz));
          }
        }
      }
    // Without pre-computed world coordinates
    } else {
      double x, y, z;
      for (int k = k1; k < k2; ++k, gx += s3, gy += s3, gz += s3)
      for (int j = j1; j < j2; ++j, gx += s2, gy += s2, gz += s2)
      for (int i = i1; i < i2; ++i, gx +=  1, gy +=  1, gz +=  1) {
        if (*gx || *gy || *gz) {
          // Convert voxel to world coordinates
          x = i, y = j, z = k;
          _Input->ImageToWorld(x, y, z);
          // Calculate derivatives of transformation w.r.t. the parameters
          _FFD->JacobianDOFs(jac, x, y, z, _t, _t0);
          // Apply chain rule to obtain similarity gradient w.r.t. the transformation parameters
          for (it = jac.Begin(); it != jac.End(); ++it) {
            _Output[it->first] += _Weight * (it->second._x * (*gx) +
                                             it->second._y * (*gy) +
                                             it->second._z * (*gz));
          }
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  void operator ()()
  {
    // Initialize members to often accessed data
    _X           = _Input->GetX();
    _Y           = _Input->GetY();
    const int _Z = _Input->GetZ();

    // Check input
    if (_WorldCoords && (_WorldCoords->X() != _X ||
                         _WorldCoords->Y() != _Y ||
                         _WorldCoords->Z() != _Z ||
                         _WorldCoords->T() !=  3)) {
      cerr << "FreeFormTransformation::ParametricGradient: Invalid world coordinates map" << endl;
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

    // Calculate parametric gradient
    blocked_range3d<int> voxels(0, _Z, 0, _Y, 0, _X);
    parallel_reduce(voxels, *this);
  }

}; // FreeFormTransformationParametricGradientBody

// -----------------------------------------------------------------------------
class FreeFormTransformationPointWiseParametricGradientBody
{
public:

  const FreeFormTransformation *_FFD;
  const PointSet               *_PointSet;
  const Vector3D<double>       *_Input;
  double                       *_Output;
  double                        _Weight;

  double _t;  ///< Time corrresponding to input gradient image (in ms)
  double _t0; ///< Second time argument for velocity-based transformations

  // ---------------------------------------------------------------------------
  /// Default constructor
  FreeFormTransformationPointWiseParametricGradientBody()
  :
    _FFD     (NULL),
    _PointSet(NULL),
    _Input   (NULL),
    _Output  (NULL),
    _Weight  ( 1.0),
    _t       ( 0.0),
    _t0      (-1.0)
  {}

  // ---------------------------------------------------------------------------
  /// Split constructor
  FreeFormTransformationPointWiseParametricGradientBody(const FreeFormTransformationPointWiseParametricGradientBody &other, split)
  :
    _FFD     (other._FFD),
    _PointSet(other._PointSet),
    _Input   (other._Input),
    _Output  (NULL),
    _Weight  (other._Weight),
    _t       (other._t),
    _t0      (other._t0)
  {
    const int ndofs = _FFD->NumberOfDOFs();
    _Output         = new double[ndofs];
    memset(_Output, 0, ndofs * sizeof(double));
  }

  // ---------------------------------------------------------------------------
  /// Destructor
  ~FreeFormTransformationPointWiseParametricGradientBody() {}

  // ---------------------------------------------------------------------------
  /// Join results of right-hand body with this body
  void join(FreeFormTransformationPointWiseParametricGradientBody &rhs)
  {
    const int ndofs = _FFD->NumberOfDOFs();
    for (int dof = 0; dof < ndofs; ++dof) _Output[dof] += rhs._Output[dof];
    delete[] rhs._Output;
    rhs._Output = NULL;
  }

  // ---------------------------------------------------------------------------
  /// Calculates the gradient of the similarity term w.r.t. the transformation
  /// parameters for the specified points in 3D.
  void operator ()(const blocked_range<int> &re)
  {
    TransformationJacobian                      jac;
    TransformationJacobian::ConstColumnIterator it;
    Point                                       p;

    for (int i = re.begin(); i != re.end(); ++i) {
      const Vector3D<double> &g = _Input[i];
      // Check whether reference point is valid
      if (g._x != .0 || g._y != .0 || g._z != .0) {
        // Calculate derivatives of transformation w.r.t. the parameters
        _PointSet->GetPoint(i, p);
        _FFD->JacobianDOFs(jac, p._x, p._y, p._z, _t, _t0);
        // Apply chain rule to obtain gradient w.r.t. the transformation parameters
        for (it = jac.Begin(); it != jac.End(); ++it) {
          _Output[it->first] += _Weight * (it->second._x * g._x +
                                           it->second._y * g._y +
                                           it->second._z * g._z);
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  void operator ()()
  {
    blocked_range<int> pts(0, _PointSet->Size());
    parallel_reduce(pts, *this);
  }

}; // FreeFormTransformationPointWiseParametricGradientBody

// -----------------------------------------------------------------------------
void FreeFormTransformation
::ParametricGradient(const GenericImage<double> *in, double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     double t0, double w) const
{
  MIRTK_START_TIMING();
  FreeFormTransformationParametricGradientBody body;
  body._FFD         = this;
  body._Input       = in;
  body._Output      = out;
  body._Weight      = w;
  body._WorldCoords = (wc ? wc : i2w);
  body._t           = in->GetTOrigin();
  body._t0          = t0;
  body();
  MIRTK_DEBUG_TIMING(2, "parametric gradient computation (FFD)");
}

// -----------------------------------------------------------------------------
void FreeFormTransformation
::ParametricGradient(const PointSet &pos, const Vector3D<double> *in,
                     double *out, double t, double t0, double w) const
{
  MIRTK_START_TIMING();
  FreeFormTransformationPointWiseParametricGradientBody body;
  body._FFD      = this;
  body._PointSet = &pos;
  body._Input    = in;
  body._Output   = out;
  body._Weight   = w;
  body._t        = t;
  body._t0       = t0;
  body();
  MIRTK_DEBUG_TIMING(2, "point-wise parametric gradient computation (FFD)");
}

// -----------------------------------------------------------------------------
void FreeFormTransformation
::FFDJacobianDetDerivative(double [3], const Matrix &, int, double, double, double, double, double, bool, bool) const
{
  Throw(ERR_NotImplemented, __FUNCTION__, "Not implemented");
}

// =============================================================================
// Properties
// =============================================================================

// -----------------------------------------------------------------------------
double FreeFormTransformation
::BendingEnergy(double x, double y, double z, double t, double t0, bool wrt_world) const
{
  if (!wrt_world) {
    cerr << "FreeFormTransformation::BendingEnergy: Always uses derivatives w.r.t. world coordinates" << endl;
    exit(1);
  }
  // Calculate 2nd order derivatives
  Matrix hessian[3];
  this->LocalHessian(hessian, x, y, z, t, t0);
  // Calculate bending energy
  return Bending3D(hessian);
}

// -----------------------------------------------------------------------------
double FreeFormTransformation
::BendingEnergy(bool incl_passive, bool wrt_world) const
{
  int    nactive = 0;
  double bending = .0;
  double x, y, z, t;

  for (int l = 0; l < _t; ++l)
  for (int k = 0; k < _z; ++k)
  for (int j = 0; j < _y; ++j)
  for (int i = 0; i < _x; ++i) {
    if (incl_passive || this->IsActive(i, j, k, l)) {
      x = i, y = j, z = k;
      this->LatticeToWorld(x, y, z);
      t = this->LatticeToTime(l);
      bending += this->BendingEnergy(x, y, z, t, 1.0, wrt_world);
      ++nactive;
    }
  }

  if (nactive > 0) bending /= nactive;
  return bending;
}

// -----------------------------------------------------------------------------
double FreeFormTransformation
::BendingEnergy(const ImageAttributes &attr, double t0, bool wrt_world) const
{
  const int N = attr.NumberOfPoints();
  const int L = attr._dt ? attr._t : 1;
  if (N == 0) return .0;

  double x, y, z, t, bending = .0;
  for (int l = 0; l < L; ++l) {
    t = attr.LatticeToTime(l);
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i) {
      x = i, y = j, z = k;
      attr.LatticeToWorld(x, y, z);
      bending += this->BendingEnergy(x, y, z, t, t0, wrt_world);
    }
  }

  return bending / N;
}

// -----------------------------------------------------------------------------
void FreeFormTransformation::BendingEnergyGradient(double *, double, bool, bool, bool) const
{
  cerr << this->NameOfClass() << "::BendingEnergyGradient: Not implemented" << endl;
  exit(1);
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void FreeFormTransformation::Print(ostream &os, Indent indent) const
{
  // Print no. of transformation parameters
  os << indent << "Number of DOFs:          " << this->NumberOfDOFs() << endl;
  os << indent << "Number of CPs (active):  " << this->NumberOfActiveCPs()<< endl;
  os << indent << "Number of CPs (passive): " << this->NumberOfPassiveCPs()<< endl;
  os << indent << "Extrapolation mode:      " << ToString(_ExtrapolationMode) << endl;
  // Print lattice attributes
  os << indent << "Control point lattice:" << endl;
  _attr.Print(os, indent + 1);
}

// -----------------------------------------------------------------------------
Cofstream &FreeFormTransformation::WriteCPs(Cofstream &to) const
{
  // Note: this->NumberOfDOFs() may differ for specialized subclasses!
  const int num = 3 * this->NumberOfCPs();
  to.WriteAsDouble(reinterpret_cast<const double *>(_CPImage.Data()),    num);
  to.WriteAsInt   (reinterpret_cast<const int    *>(_CPStatus[0][0][0]), num);
  return to;
}


} // namespace mirtk

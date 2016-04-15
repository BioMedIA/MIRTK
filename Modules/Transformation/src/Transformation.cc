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

#include "mirtk/Transformation.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/VoxelFunction.h"

#include "mirtk/AdaptiveLineSearch.h"
#include "mirtk/ConjugateGradientDescent.h"
#include "mirtk/TransformationApproximationError.h"

#include "mirtk/Transformations.h"

#include "TransformationUtils.h"


namespace mirtk {


// =============================================================================
// Declaration of inverse transformation (cf. TransformationInverse.cc)
// =============================================================================

bool EvaluateGlobalInverse(const Transformation *, double &, double &, double &, double, double);
bool EvaluateLocalInverse (const Transformation *, double &, double &, double &, double, double);
bool EvaluateInverse      (const Transformation *, double &, double &, double &, double, double);

// =============================================================================
// Factory methods
// =============================================================================

// -----------------------------------------------------------------------------
Transformation *Transformation::New(TransformationType type)
{
  switch (type) {
    // -------------------------------------------------------------------------
    // linear transformation
    case TRANSFORMATION_HOMOGENEOUS:
      return new HomogeneousTransformation;
    case TRANSFORMATION_RIGID:
      return new RigidTransformation;
    case TRANSFORMATION_SIMILARITY:
      return new SimilarityTransformation;
    case TRANSFORMATION_AFFINE:
      return new AffineTransformation;
    // -------------------------------------------------------------------------
    // linear FFD
    case TRANSFORMATION_LINEAR_FFD_2D_v1:
      // TODO: Implement separate LinearFreeFormTransformation2D
      //       instead of using the 3D version which deals with both.
    case TRANSFORMATION_LINEAR_FFD_3D_v1:
    case TRANSFORMATION_LINEAR_FFD_3D_v2:
    case TRANSFORMATION_LINEAR_FFD_3D_v3:
    case TRANSFORMATION_LINEAR_FFD_3D_v4:
      return new LinearFreeFormTransformation3D;
    case TRANSFORMATION_LINEAR_FFD_4D_v1:
    case TRANSFORMATION_LINEAR_FFD_4D_v2:
    case TRANSFORMATION_LINEAR_FFD_4D_v3:
      return new LinearFreeFormTransformation4D;
//    case TRANSFORMATION_LINEAR_FFD_SV_v1:
//      return new LinearFreeFormTransformationSV;
    case TRANSFORMATION_LINEAR_FFD_TD_v1:
    case TRANSFORMATION_LINEAR_FFD_TD_v2:
    case TRANSFORMATION_LINEAR_FFD_TD_v3:
      return new LinearFreeFormTransformationTD;
    // -------------------------------------------------------------------------
    // B-spline FFD
    case TRANSFORMATION_BSPLINE_FFD_2D_v1:
      // TODO: Implement separate BSplineFreeFormTransformation2D
      //       instead of using the 3D version which deals with both.
    case TRANSFORMATION_BSPLINE_FFD_3D_v1:
    case TRANSFORMATION_BSPLINE_FFD_3D_v2:
    case TRANSFORMATION_BSPLINE_FFD_3D_v3:
    case TRANSFORMATION_BSPLINE_FFD_3D_v4:
      return new BSplineFreeFormTransformation3D;
    case TRANSFORMATION_BSPLINE_FFD_4D_v1:
    case TRANSFORMATION_BSPLINE_FFD_4D_v2:
    case TRANSFORMATION_BSPLINE_FFD_4D_v3:
      return new BSplineFreeFormTransformation4D;
    case TRANSFORMATION_BSPLINE_FFD_SV_v1:
    case TRANSFORMATION_BSPLINE_FFD_SV_v2:
    case TRANSFORMATION_BSPLINE_FFD_SV_v3:
    case TRANSFORMATION_BSPLINE_FFD_SV_v4:
    case TRANSFORMATION_BSPLINE_FFD_SV_v5:
    case TRANSFORMATION_BSPLINE_FFD_SV_v6:
    case TRANSFORMATION_BSPLINE_FFD_SV_v7:
    case TRANSFORMATION_BSPLINE_FFD_SV_v8:
      return new BSplineFreeFormTransformationSV;
    case TRANSFORMATION_BSPLINE_FFD_TD_v1:
    case TRANSFORMATION_BSPLINE_FFD_TD_v2:
    case TRANSFORMATION_BSPLINE_FFD_TD_v3:
    case TRANSFORMATION_BSPLINE_FFD_TD_v4:
      return new BSplineFreeFormTransformationTD;
    case TRANSFORMATION_BSPLINE_FFD_STATISTICAL:
      return new BSplineFreeFormTransformationStatistical;
    // -------------------------------------------------------------------------
    // composite transformation
    case TRANSFORMATION_MFFD:
      return new MultiLevelFreeFormTransformation;
    case TRANSFORMATION_MFFD_SV:
      return new MultiLevelStationaryVelocityTransformation;
    case TRANSFORMATION_FLUID_v1:
    case TRANSFORMATION_FLUID_v2:
      return new FluidFreeFormTransformation;
    // -------------------------------------------------------------------------
    default:
      cerr << "Transformation::New: Unknown transformation type: " << type << endl;
      return NULL;
  }
}

// -----------------------------------------------------------------------------
template <class TransformationType>
inline Transformation *CopyOf(const Transformation *t)
{
  return new TransformationType(*static_cast<const TransformationType *>(t));
}

// -----------------------------------------------------------------------------
Transformation *Transformation::New(const Transformation *t)
{
  switch (t->TypeOfClass()) {
    // -------------------------------------------------------------------------
    // linear transformation
    case TRANSFORMATION_HOMOGENEOUS:
      return CopyOf<HomogeneousTransformation>(t);
    case TRANSFORMATION_RIGID:
      return CopyOf<RigidTransformation>(t);
    case TRANSFORMATION_SIMILARITY:
      return CopyOf<SimilarityTransformation>(t);
    case TRANSFORMATION_AFFINE:
      return CopyOf<AffineTransformation>(t);
    // -------------------------------------------------------------------------
    // linear FFD
    case TRANSFORMATION_LINEAR_FFD_2D_v1:
      // TODO: return CopyOf<LinearFreeFormTransformation2D>(t);
    case TRANSFORMATION_LINEAR_FFD_3D_v1:
    case TRANSFORMATION_LINEAR_FFD_3D_v2:
    case TRANSFORMATION_LINEAR_FFD_3D_v3:
    case TRANSFORMATION_LINEAR_FFD_3D_v4:
      return CopyOf<LinearFreeFormTransformation3D>(t);
    case TRANSFORMATION_LINEAR_FFD_4D_v1:
    case TRANSFORMATION_LINEAR_FFD_4D_v2:
    case TRANSFORMATION_LINEAR_FFD_4D_v3:
      return CopyOf<LinearFreeFormTransformation4D>(t);
//    case TRANSFORMATION_LINEAR_FFD_SV_v1:
//      return CopyOf<LinearFreeFormTransformationSV>(t);
    case TRANSFORMATION_LINEAR_FFD_TD_v1:
    case TRANSFORMATION_LINEAR_FFD_TD_v2:
    case TRANSFORMATION_LINEAR_FFD_TD_v3:
      return CopyOf<LinearFreeFormTransformationTD>(t);
    // -------------------------------------------------------------------------
    // B-spline FFD
    case TRANSFORMATION_BSPLINE_FFD_2D_v1:
      // TODO: return CopyOf<BSplineFreeFormTransformation2D>(t);
    case TRANSFORMATION_BSPLINE_FFD_3D_v1:
    case TRANSFORMATION_BSPLINE_FFD_3D_v2:
    case TRANSFORMATION_BSPLINE_FFD_3D_v3:
    case TRANSFORMATION_BSPLINE_FFD_3D_v4:
      return CopyOf<BSplineFreeFormTransformation3D>(t);
    case TRANSFORMATION_BSPLINE_FFD_4D_v1:
    case TRANSFORMATION_BSPLINE_FFD_4D_v2:
    case TRANSFORMATION_BSPLINE_FFD_4D_v3:
      return CopyOf<BSplineFreeFormTransformation4D>(t);
    case TRANSFORMATION_BSPLINE_FFD_SV_v1:
    case TRANSFORMATION_BSPLINE_FFD_SV_v2:
    case TRANSFORMATION_BSPLINE_FFD_SV_v3:
    case TRANSFORMATION_BSPLINE_FFD_SV_v4:
    case TRANSFORMATION_BSPLINE_FFD_SV_v5:
    case TRANSFORMATION_BSPLINE_FFD_SV_v6:
    case TRANSFORMATION_BSPLINE_FFD_SV_v7:
    case TRANSFORMATION_BSPLINE_FFD_SV_v8:
      return CopyOf<BSplineFreeFormTransformationSV>(t);
    case TRANSFORMATION_BSPLINE_FFD_TD_v1:
    case TRANSFORMATION_BSPLINE_FFD_TD_v2:
    case TRANSFORMATION_BSPLINE_FFD_TD_v3:
    case TRANSFORMATION_BSPLINE_FFD_TD_v4:
      return CopyOf<BSplineFreeFormTransformationTD>(t);
    case TRANSFORMATION_BSPLINE_FFD_STATISTICAL:
      return CopyOf<BSplineFreeFormTransformationStatistical>(t);
    // -------------------------------------------------------------------------
    // composite transformation
    case TRANSFORMATION_MFFD:
      return CopyOf<MultiLevelFreeFormTransformation>(t);
    case TRANSFORMATION_MFFD_SV:
      return CopyOf<MultiLevelStationaryVelocityTransformation>(t);
    case TRANSFORMATION_FLUID_v1:
    case TRANSFORMATION_FLUID_v2:
      return CopyOf<FluidFreeFormTransformation>(t);
    // -------------------------------------------------------------------------
    default:
      cerr << "Transformation::New: Failed to copy transformation of type: " << t->TypeOfClass() << endl;
      return NULL;
  }
}

// -----------------------------------------------------------------------------
Transformation *Transformation::New(const char *name)
{
  Transformation *t = NULL;

  Cifstream from;
  from.Open(name);
  unsigned int magic_no;
  from.ReadAsUInt(&magic_no, 1);

  if (magic_no != TRANSFORMATION_MAGIC) {
    from.Close();
    cerr << "Transformation::New: Not a transformation file: " << name << endl;
    exit(1);
  }

  unsigned int trans_type;
  from.ReadAsUInt(&trans_type, 1);
  from.Close();

  switch (trans_type) {
    // -------------------------------------------------------------------------
    // linear transformation
    case TRANSFORMATION_HOMOGENEOUS:
      t = new HomogeneousTransformation;
      break;
    case TRANSFORMATION_RIGID:
      t = new RigidTransformation;
      break;
    case TRANSFORMATION_SIMILARITY:
      t = new SimilarityTransformation;
      break;
    case TRANSFORMATION_AFFINE:
      t = new AffineTransformation;
      break;
    // -------------------------------------------------------------------------
    // linear FFD
    case TRANSFORMATION_LINEAR_FFD_2D_v1:
// TODO: Implement LinearFreeFormTransformation2D and uncomment the following.
//      t = new LinearFreeFormTransformation2D;
//      break;
    case TRANSFORMATION_LINEAR_FFD_3D_v1:
    case TRANSFORMATION_LINEAR_FFD_3D_v2:
    case TRANSFORMATION_LINEAR_FFD_3D_v3:
    case TRANSFORMATION_LINEAR_FFD_3D_v4:
      t = new LinearFreeFormTransformation3D;
      break;
    case TRANSFORMATION_LINEAR_FFD_4D_v1:
    case TRANSFORMATION_LINEAR_FFD_4D_v2:
    case TRANSFORMATION_LINEAR_FFD_4D_v3:
      t = new LinearFreeFormTransformation4D;
      break;
// TODO: Implement LinearFreeFormTransformationSV and uncomment the following.
//    case TRANSFORMATION_LINEAR_FFD_SV_v1:
//      t = new LinearFreeFormTransformationSV;
//      break;
    case TRANSFORMATION_LINEAR_FFD_TD_v1:
    case TRANSFORMATION_LINEAR_FFD_TD_v2:
    case TRANSFORMATION_LINEAR_FFD_TD_v3:
      t = new LinearFreeFormTransformationTD;
      break;
    // -------------------------------------------------------------------------
    // B-spline FFD
    case TRANSFORMATION_BSPLINE_FFD_2D_v1:
// TODO: Implement BSplineFreeFormTransformation2D and uncomment the following.
//      t = new BSplineFreeFormTransformation2D;
//      break;
    case TRANSFORMATION_BSPLINE_FFD_3D_v1:
    case TRANSFORMATION_BSPLINE_FFD_3D_v2:
    case TRANSFORMATION_BSPLINE_FFD_3D_v3:
    case TRANSFORMATION_BSPLINE_FFD_3D_v4:
      t = new BSplineFreeFormTransformation3D;
      break;
    case TRANSFORMATION_BSPLINE_FFD_4D_v1:
    case TRANSFORMATION_BSPLINE_FFD_4D_v2:
    case TRANSFORMATION_BSPLINE_FFD_4D_v3:
      t = new BSplineFreeFormTransformation4D;
      break;
    case TRANSFORMATION_BSPLINE_FFD_SV_v1:
    case TRANSFORMATION_BSPLINE_FFD_SV_v2:
    case TRANSFORMATION_BSPLINE_FFD_SV_v3:
    case TRANSFORMATION_BSPLINE_FFD_SV_v4:
    case TRANSFORMATION_BSPLINE_FFD_SV_v5:
    case TRANSFORMATION_BSPLINE_FFD_SV_v6:
    case TRANSFORMATION_BSPLINE_FFD_SV_v7:
    case TRANSFORMATION_BSPLINE_FFD_SV_v8:
      t = new BSplineFreeFormTransformationSV;
      break;
    case TRANSFORMATION_BSPLINE_FFD_TD_v1:
    case TRANSFORMATION_BSPLINE_FFD_TD_v2:
    case TRANSFORMATION_BSPLINE_FFD_TD_v3:
    case TRANSFORMATION_BSPLINE_FFD_TD_v4:
      t = new BSplineFreeFormTransformationTD;
      break;
    // -------------------------------------------------------------------------
    // composite transformation
    case TRANSFORMATION_MFFD:
      t = new MultiLevelFreeFormTransformation;
      break;
    case TRANSFORMATION_MFFD_SV:
      t = new MultiLevelStationaryVelocityTransformation;
      break;
    case TRANSFORMATION_FLUID_v1:
    case TRANSFORMATION_FLUID_v2:
      t = new FluidFreeFormTransformation;
      break;
    // -------------------------------------------------------------------------
    default:
      cerr << "Transformation::New: Transformation file " << name << " has unknown format: " << trans_type << endl;
      exit(1);
  }

  t->Read(name);
  return t;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
Transformation::Transformation(int ndofs)
:
  _NumberOfDOFs(0),
  _Param       (NULL),
  _Status      (NULL)
{
  InitializeDOFs(ndofs);
}

// -----------------------------------------------------------------------------
Transformation::Transformation(const Transformation &t)
:
  _NumberOfDOFs(0),
  _Param       (NULL),
  _Status      (NULL)
{
  InitializeDOFs(t);
}

// -----------------------------------------------------------------------------
Transformation::Transformation(const Transformation &t, int ndofs)
:
  _NumberOfDOFs(0),
  _Param       (NULL),
  _Status      (NULL)
{
  InitializeDOFs(t, ndofs);
}

// -----------------------------------------------------------------------------
void Transformation::InitializeDOFs(int ndofs)
{
  if (_NumberOfDOFs != ndofs) {
    Deallocate(_Param);
    Deallocate(_Status);
    _NumberOfDOFs = ndofs;
    if (_NumberOfDOFs > 0) {
      CAllocate(_Param,  _NumberOfDOFs, .0);
      CAllocate(_Status, _NumberOfDOFs, Active);
    }
  }
}

// -----------------------------------------------------------------------------
void Transformation::InitializeDOFs(const Transformation &t, int ndofs)
{
  if (ndofs < 0) ndofs = t.NumberOfDOFs();
  InitializeDOFs(ndofs);
  if (ndofs > 0) {
    if (t.NumberOfDOFs() < ndofs) ndofs = t.NumberOfDOFs();
    memcpy(_Param,  t._Param,  ndofs * sizeof(DOFValue));
    memcpy(_Status, t._Status, ndofs * sizeof(DOFStatus));
  }
}

// -----------------------------------------------------------------------------
Transformation::~Transformation()
{
  Deallocate(_Param);
  Deallocate(_Status);
}

// =============================================================================
// Others
// =============================================================================

// -----------------------------------------------------------------------------
TransformationType Transformation::TypeOfClass(const char *clsname)
{
  if      (strcmp(clsname, "HomogeneousTransformation")                  == 0) return TRANSFORMATION_HOMOGENEOUS;
  else if (strcmp(clsname, "RigidTransformation")                        == 0) return TRANSFORMATION_RIGID;
  else if (strcmp(clsname, "SimilarityTransformation")                   == 0) return TRANSFORMATION_SIMILARITY;
  else if (strcmp(clsname, "AffineTransformation")                       == 0) return TRANSFORMATION_AFFINE;
  else if (strcmp(clsname, "LinearFreeFormTransformation")               == 0) return TRANSFORMATION_LINEAR_FFD_3D;
  else if (strcmp(clsname, "LinearFreeFormTransformation3D")             == 0) return TRANSFORMATION_LINEAR_FFD_3D;
  else if (strcmp(clsname, "LinearFreeFormTransformation4D")             == 0) return TRANSFORMATION_LINEAR_FFD_4D;
  else if (strcmp(clsname, "LinearFreeFormTransformationTD")             == 0) return TRANSFORMATION_LINEAR_FFD_TD;
  else if (strcmp(clsname, "BSplineFreeFormTransformation")              == 0) return TRANSFORMATION_BSPLINE_FFD_3D;
  else if (strcmp(clsname, "BSplineFreeFormTransformation3D")            == 0) return TRANSFORMATION_BSPLINE_FFD_3D;
  else if (strcmp(clsname, "BSplineFreeFormTransformation4D")            == 0) return TRANSFORMATION_BSPLINE_FFD_4D;
  else if (strcmp(clsname, "BSplineFreeFormTransformationTD")            == 0) return TRANSFORMATION_BSPLINE_FFD_TD;
  else if (strcmp(clsname, "BSplineFreeFormTransformationSV")            == 0) return TRANSFORMATION_BSPLINE_FFD_SV;
  else if (strcmp(clsname, "BSplineFreeFormTransformationPeriodic")      == 0) return TRANSFORMATION_PERIODIC;
  else if (strcmp(clsname, "BSplineFreeFormTransformationStatistical")   == 0) return TRANSFORMATION_BSPLINE_FFD_STATISTICAL;
  else if (strcmp(clsname, "MultiLevelFreeFormTransformation")           == 0) return TRANSFORMATION_MFFD;
  else if (strcmp(clsname, "MultiLevelStationaryVelocityTransformation") == 0) return TRANSFORMATION_MFFD_SV;
  else if (strcmp(clsname, "FluidFreeFormTransformation")                == 0) return TRANSFORMATION_FLUID;
  else if (strcmp(clsname, "LatticeFreeFormTransformation")              == 0) return TRANSFORMATION_LATTICE_FFD;
  else if (strcmp(clsname, "MultiFrameLatticeFreeFormTransformation")    == 0) return TRANSFORMATION_MULTI_FRAME_LATTICE_FFD;

  return TRANSFORMATION_UNKNOWN;
}

// -----------------------------------------------------------------------------
void Transformation::Verify()
{
}

// =============================================================================
// Parameters (non-DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool Transformation::Set(const char *, const char *)
{
  return false;
}

// -----------------------------------------------------------------------------
ParameterList Transformation::Parameter() const
{
  return ParameterList();
}

// =============================================================================
// Parameters (DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool Transformation::CopyFrom(const Transformation *other)
{
  if (other && strcmp(this->NameOfClass(), other->NameOfClass()) == 0
            && _NumberOfDOFs == other->NumberOfDOFs()) {
    this->Reset();
    for (int dof = 0; dof < _NumberOfDOFs; ++dof) {
      if (_Status[dof] == Active) {
        _Param[dof] = other->_Param[dof];
      }
    }
    this->Changed(true);
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
void Transformation::Reset()
{
  memset(_Param, 0, _NumberOfDOFs * sizeof(DOFValue));
  this->Changed(true);
}

// -----------------------------------------------------------------------------
double Transformation::EvaluateRMSError(const ImageAttributes &domain,
                                        const Transformation  *other) const
{
  const int n = domain.NumberOfPoints();
  if (n <= 0) return .0;

  double error = .0, dx, dy, dz, t;

  if (this->RequiresCachingOfDisplacements() && other->RequiresCachingOfDisplacements()) {

    GenericImage<double> disp1(domain, 3);
    GenericImage<double> disp2(domain, 3);

    for (int l = 0; l < domain._t; ++l) {
      t = domain.LatticeToTime(l);
      disp1 = .0;
      disp2 = .0;
      disp1.PutTOrigin(t);
      disp2.PutTOrigin(t);
      other->Displacement(disp1);
      this ->Displacement(disp2);
      for (int k = 0; k < domain._z; ++k)
      for (int j = 0; j < domain._y; ++j)
      for (int i = 0; i < domain._x; ++i) {
        dx = disp1(i, j, k, 0) - disp2(i, j, k, 0);
        dy = disp1(i, j, k, 1) - disp2(i, j, k, 1);
        dz = disp1(i, j, k, 2) - disp2(i, j, k, 2);
        error += dx*dx + dy*dy + dz*dz;
      }
    }

  } else if (this->RequiresCachingOfDisplacements()) {

    GenericImage<double> disp(domain, 3);

    for (int l = 0; l < domain._t; ++l) {
      t = domain.LatticeToTime(l);
      disp = .0;
      disp.PutTOrigin(t);
      this->Displacement(disp);
      for (int k = 0; k < domain._z; ++k)
      for (int j = 0; j < domain._y; ++j)
      for (int i = 0; i < domain._x; ++i) {
        dx = i, dy = j, dz = k;
        domain.LatticeToWorld(dx, dy, dz);
        other->Displacement  (dx, dy, dz);
        dx -= disp(i, j, k, 0);
        dy -= disp(i, j, k, 1);
        dz -= disp(i, j, k, 2);
        error += dx*dx + dy*dy + dz*dz;
      }
    }

  } else if (other->RequiresCachingOfDisplacements()) {

    GenericImage<double> disp(domain, 3);

    for (int l = 0; l < domain._t; ++l) {
      t = domain.LatticeToTime(l);
      disp = .0;
      disp.PutTOrigin(t);
      other->Displacement(disp);
      for (int k = 0; k < domain._z; ++k)
      for (int j = 0; j < domain._y; ++j)
      for (int i = 0; i < domain._x; ++i) {
        dx = i, dy = j, dz = k;
        domain.LatticeToWorld(dx, dy, dz);
        this->Displacement   (dx, dy, dz);
        dx -= disp(i, j, k, 0);
        dy -= disp(i, j, k, 1);
        dz -= disp(i, j, k, 2);
        error += dx*dx + dy*dy + dz*dz;
      }
    }

  } else {

    double x1, y1, z1, x2, y2, z2;
    for (int l = 0; l < domain._t; ++l) {
      t = domain.LatticeToTime(l);
      for (int k = 0; k < domain._z; ++k)
      for (int j = 0; j < domain._y; ++j)
      for (int i = 0; i < domain._x; ++i) {
        x1 = i, y1 = j, z1 = k;
        domain.LatticeToWorld(x1, y1, z1);
        x2 = x1, y2 = y1, z2 = z1;
        other->Transform(x1, y1, z1, t);
        this ->Transform(x2, y2, z2, t);
        dx = x1 - x2;
        dy = y1 - y2;
        dz = z1 - z2;
        error += dx*dx + dy*dy + dz*dz;
      }
    }
    
  }
  
  return sqrt(error / n);
}

// -----------------------------------------------------------------------------
double Transformation::EvaluateRMSError(const ImageAttributes &domain,
                                        double *dx, double *dy) const
{
  const int no = domain.NumberOfPoints();
  if (no <= 0) return .0;

  double error = .0, t;

  if (this->RequiresCachingOfDisplacements()) {

    GenericImage<double> disp(domain, 2);

    for (int l = 0; l < domain._t; ++l) {
      t = domain.LatticeToTime(l);
      disp = .0;
      this->Displacement(disp, t);
      for (int k = 0; k < domain._z; ++k)
      for (int j = 0; j < domain._y; ++j)
      for (int i = 0; i < domain._x; ++i) {
        const int idx = domain.LatticeToIndex(i, j, k, l);
        dx[idx] -= disp(i, j, k, 0);
        dy[idx] -= disp(i, j, k, 1);
        error += dx[idx] * dx[idx] + dy[idx] * dy[idx];
      }
    }

  } else {

    double tx, ty, tz;
    for (int l = 0; l < domain._t; ++l) {
      t = domain.LatticeToTime(l);
      for (int k = 0; k < domain._z; ++k)
      for (int j = 0; j < domain._y; ++j)
      for (int i = 0; i < domain._x; ++i) {
        const int idx = domain.LatticeToIndex(i, j, k, l);
        tx = i, ty = j, tz = k;
        domain.LatticeToWorld(tx, ty, tz);
        this->Displacement(tx, ty, tz, t);
        dx[idx] -= tx;
        dy[idx] -= ty;
        error += dx[idx] * dx[idx] + dy[idx] * dy[idx];
      }
    }

  }

  return sqrt(error / no);
}

// -----------------------------------------------------------------------------
double Transformation::EvaluateRMSError(const ImageAttributes &domain,
                                        double *dx, double *dy, double *dz) const
{
  const int no = domain.NumberOfPoints();
  if (no <= 0) return .0;

  double error = .0, t;

  if (this->RequiresCachingOfDisplacements()) {

    GenericImage<double> disp(domain, 3);

    for (int l = 0; l < domain._t; ++l) {
      t = domain.LatticeToTime(l);
      disp = .0;
      this->Displacement(disp, t);
      for (int k = 0; k < domain._z; ++k)
      for (int j = 0; j < domain._y; ++j)
      for (int i = 0; i < domain._x; ++i) {
        const int idx = domain.LatticeToIndex(i, j, k, l);
        dx[idx] -= disp(i, j, k, 0);
        dy[idx] -= disp(i, j, k, 1);
        dz[idx] -= disp(i, j, k, 2);
        error += pow(dx[idx], 2) + pow(dy[idx], 2) + pow(dz[idx], 2);
      }
    }

  } else {

    double tx, ty, tz;
    for (int l = 0; l < domain._t; ++l) {
      t = domain.LatticeToTime(l);
      for (int k = 0; k < domain._z; ++k)
      for (int j = 0; j < domain._y; ++j)
      for (int i = 0; i < domain._x; ++i) {
        const int idx = domain.LatticeToIndex(i, j, k, l);
        tx = i, ty = j, tz = k;
        domain.LatticeToWorld(tx, ty, tz);
        this->Displacement(tx, ty, tz, t);
        dx[idx] -= tx;
        dy[idx] -= ty;
        dz[idx] -= tz;
        error += pow(dx[idx], 2) + pow(dy[idx], 2) + pow(dz[idx], 2);
      }
    }

  }

  return sqrt(error / no);
}

// -----------------------------------------------------------------------------
double Transformation
::EvaluateRMSError(const double *x,  const double *y,  const double *z, double t,
                   double       *dx, double       *dy, double       *dz, int no) const
{
  if (no <= 0) return .0;
  double error = .0, tx, ty, tz;
  for (int idx = 0; idx < no; ++idx) {
    tx = x[idx], ty = y[idx], tz = z[idx];
    this->Displacement(tx, ty, tz, t);
    dx[idx] -= tx;
    dy[idx] -= ty;
    dz[idx] -= tz;
    error += pow(dx[idx], 2) + pow(dy[idx], 2) + pow(dz[idx], 2);
  }
  return sqrt(error / no);
}

// -----------------------------------------------------------------------------
double Transformation
::EvaluateRMSError(const double *x,  const double *y,  const double *z,  const double *t,
                   double       *dx, double       *dy, double       *dz, int no) const
{
  if (no <= 0) return .0;

  TransformationUtils::SubDisplacements sub(this, x, y, z, t, dx, dy, dz);
  parallel_for(blocked_range<int>(0, no), sub);

  double error = .0;
  for (int idx = 0; idx < no; ++idx) {
    error += pow(dx[idx], 2) + pow(dy[idx], 2) + pow(dz[idx], 2);
  }
  return sqrt(error / no);
}

// -----------------------------------------------------------------------------
double Transformation::Approximate(const ImageAttributes &domain,
                                       const Transformation  *other,
                                       int niter, double max_error)
{
  const int no = domain.NumberOfPoints();
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Evaluate displacements of transformation at lattice points
  double *dx = CAllocate<double>(no);
  double *dy = CAllocate<double>(no);
  double *dz = CAllocate<double>(no);
  other->Displacement(domain, dx, dy, dz);

  // Approximate displacements
  error = this->Approximate(domain, dx, dy, dz, niter, max_error);

  // Free memory
  Deallocate(dx);
  Deallocate(dy);
  Deallocate(dz);

  return error;
}

// -----------------------------------------------------------------------------
double Transformation::Approximate(GenericImage<double> &disp, int niter, double max_error)
{
  if (disp.GetTOrigin() == .0 && disp.T() == 3) {
    return this->Approximate(disp.Attributes(), disp.Data(0, 0, 0, 0),
                                                disp.Data(0, 0, 0, 1),
                                                disp.Data(0, 0, 0, 2),
                             niter, max_error);
  } else {
    cerr << this->NameOfClass() << "::Approximate: Displacement field must have _t = 3 and _dt = 0" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
double Transformation::Approximate(const ImageAttributes &domain,
                                       double *dx, double *dy, double *dz,
                                       int niter, double max_error)
{
  const int no = domain.NumberOfPoints();
  if (no    <= 0) return .0;
  if (niter <  1) return numeric_limits<double>::infinity();
  // Note: Iterative approximation depends on how the parameters of each
  //       approximation iteration are to be combined with the previous
  //       parameters (i.e., addition vs. composition).
  if (niter > 1) {
    cerr << "WARNING: Transformation::Approximate not implemented by " << this->NameOfClass() << endl;
    cerr << "         Using ApproximateAsNew with one iteration instead." << endl;
  }
  // Attention: niter must 1 otherwise this may lead to endless recursion!
  return this->ApproximateAsNew(domain, dx, dy, dz);
}

// -----------------------------------------------------------------------------
double Transformation
::Approximate(const double *x,  const double *y,  const double *z,
              double       *dx, double       *dy, double       *dz, int no,
              int niter, double max_error)
{
  if (no    <= 0) return .0;
  if (niter <  1) return numeric_limits<double>::infinity();
  // Note: Iterative approximation depends on how the parameters of each
  //       approximation iteration are to be combined with the previous
  //       parameters (i.e., addition vs. composition).
  if (niter > 1) {
    cerr << "WARNING: Transformation::Approximate not implemented by " << this->NameOfClass() << endl;
    cerr << "         Using ApproximateAsNew with one iteration instead." << endl;
  }
  // Attention: niter must 1 otherwise this may lead to endless recursion!
  return this->ApproximateAsNew(x, y, z, dx, dy, dz, no);
}

// -----------------------------------------------------------------------------
double Transformation
::Approximate(const double *x,  const double *y,  const double *z,  const double *t,
              double       *dx, double       *dy, double       *dz, int no,
              int niter, double max_error)
{
  if (no    <= 0) return .0;
  if (niter <  1) return numeric_limits<double>::infinity();
  // Note: Iterative approximation depends on how the parameters of each
  //       approximation iteration are to be combined with the previous
  //       parameters (i.e., addition vs. composition).
  if (niter > 1) {
    cerr << "WARNING: Transformation::Approximate not implemented by " << this->NameOfClass() << endl;
    cerr << "         Using ApproximateAsNew with one iteration instead." << endl;
  }
  // Attention: niter must 1 otherwise this may lead to endless recursion!
  return this->ApproximateAsNew(x, y, z, t, dx, dy, dz, no);
}

// -----------------------------------------------------------------------------
double Transformation::ApproximateAsNew(const ImageAttributes &domain,
                                        const Transformation  *other,
                                        int niter, double max_error)
{
  const int no = domain.NumberOfPoints();
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Evaluate displacements of transformation at lattice points
  double *dx = CAllocate<double>(no);
  double *dy = CAllocate<double>(no);
  double *dz = CAllocate<double>(no);
  other->Displacement(domain, dx, dy, dz);

  // Approximate displacements
  error = this->ApproximateAsNew(domain, dx, dy, dz, niter, max_error);

  // Free memory
  Deallocate(dx);
  Deallocate(dy);
  Deallocate(dz);
  
  return error;
}

// -----------------------------------------------------------------------------
double Transformation::ApproximateAsNew(GenericImage<double> &disp,
                                        int niter, double max_error)
{
  if (disp.GetTSize() == .0 && disp.T() == 3) {
    return this->ApproximateAsNew(disp.Attributes(), disp.Data(0, 0, 0, 0),
                                                     disp.Data(0, 0, 0, 1),
                                                     disp.Data(0, 0, 0, 2),
                                  niter, max_error);
  } else {
    cerr << this->NameOfClass() << "::ApproximateAsNew: Displacement field must have _t = 3 and _dt = 0" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
double Transformation::ApproximateAsNew(const ImageAttributes &domain,
                                        double *dx, double *dy, double *dz,
                                        int niter, double max_error)
{
  // Reset transformation
  this->Reset();

  // Check input arguments
  const int no = domain.NumberOfPoints();
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Either approximate displacements iteratively
  if (niter > 1) {

    error = this->Approximate(domain, dx, dy, dz, niter, max_error);

  // or directly without intermediate copies (requires less memory)
  } else {

    // Compute world coordinates of lattice points
    double *x = Allocate<double>(no);
    double *y = Allocate<double>(no);
    double *z = Allocate<double>(no);
    double *t = Allocate<double>(no);
    domain.LatticeToWorld(x, y, z, t);

    // Approximate displacements
    this->ApproximateDOFs(x, y, z, t, dx, dy, dz, no);

    // Evaluate error of approximation
    if (this->RequiresCachingOfDisplacements()) {
      error = EvaluateRMSError(domain,     dx, dy, dz);
    } else {
      error = EvaluateRMSError(x, y, z, t, dx, dy, dz, no);
    }

    // Free memory
    Deallocate(x);
    Deallocate(y);
    Deallocate(z);
    Deallocate(t);
  }

  return error;
}

// -----------------------------------------------------------------------------
double Transformation
::ApproximateAsNew(const double *x,  const double *y,  const double *z,
                   double       *dx, double       *dy, double       *dz, int no,
                   int niter, double max_error)
{
  // Reset transformation
  this->Reset();

  // Check input arguments
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Either approximate displacements iteratively
  if (niter > 1) {

    error = this->Approximate(x, y, z, dx, dy, dz, no, niter, max_error);

  // or directly without intermediate copies (requires less memory)
  } else {

    // Allocate fixed time coordinates
    const double t0 = .0;
    double       *t = CAllocate<double>(no, &t0);

    // Approximate displacements
    this->ApproximateDOFs(x, y, z, t, dx, dy, dz, no);

    // Evaluate error of approximation
    error = EvaluateRMSError(x, y, z, .0, dx, dy, dz, no);

    // Free memory
    Deallocate(t);
  }

  return error;
}

// -----------------------------------------------------------------------------
double Transformation
::ApproximateAsNew(const double *x,  const double *y,  const double *z,  const double *t,
                   double       *dx, double       *dy, double       *dz, int no,
                   int niter, double max_error)
{
  // Reset transformation
  this->Reset();

  // Check input arguments
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Either approximate displacements iteratively
  if (niter > 1) {

    error = this->Approximate(x, y, z, t, dx, dy, dz, no, niter, max_error);

  // or directly without intermediate copies (requires less memory)
  } else {

    // Approximate displacements
    this->ApproximateDOFs(x, y, z, t, dx, dy, dz, no);

    // Evaluate error of approximation
    error = EvaluateRMSError(x, y, z, t, dx, dy, dz, no);

  }

  return error;
}

// -----------------------------------------------------------------------------
void Transformation
::ApproximateGradient(const ImageAttributes &domain,
                      const double *dx, const double *dy, const double *dz,
                      double *gradient, double w) const
{
  const int no = domain.NumberOfPoints();
  if (no <= 0) return;

  // Compute world coordinates of lattice points
  double *x = Allocate<double>(no);
  double *y = Allocate<double>(no);
  double *z = Allocate<double>(no);
  double *t = Allocate<double>(no);
  domain.LatticeToWorld(x, y, z, t);

  // Add gradient of approximation error w.r.t. DoFs
  this->ApproximateDOFsGradient(x, y, z, t, dx, dy, dz, no, gradient, w);

  // Free memory
  Deallocate(x);
  Deallocate(y);
  Deallocate(z);
  Deallocate(t);
}

// -----------------------------------------------------------------------------
void Transformation
::ApproximateGradient(const double *x,  const double *y,  const double *z,
                      const double *dx, const double *dy, const double *dz, int no,
                      double *gradient, double w) const
{
  if (no <= 0) return;
  const double t0 = .0;
  double       *t = CAllocate<double>(no, &t0);
  this->ApproximateGradient(x, y, z, t, dx, dy, dz, no, gradient, w);
  Deallocate(t);
}

// -----------------------------------------------------------------------------
void Transformation
::ApproximateGradient(const double *x,  const double *y,  const double *z,  const double *t,
                      const double *dx, const double *dy, const double *dz, int no,
                      double *gradient, double w) const
{
  // Add gradient of approximation error w.r.t. DoFs
  this->ApproximateDOFsGradient(x, y, z, t, dx, dy, dz, no, gradient, w);
}

// -----------------------------------------------------------------------------
void Transformation
::ApproximateDOFs(const double *x,  const double *y,  const double *z,  const double *t,
                  const double *dx, const double *dy, const double *dz, int no)
{
  MIRTK_START_TIMING();

  // Reset transformation
  this->Reset();

  // Mean squared error function
  TransformationApproximationError error(this, x, y, z, t, dx, dy, dz, no);

  // Optimization method
  AdaptiveLineSearch linesearch;
  linesearch.MaxRejectedStreak(0);
  linesearch.StrictStepLengthRange(false);
  linesearch.ReusePreviousStepLength(true);
  linesearch.NumberOfIterations(20);
  linesearch.MinStepLength(1e-3);
  linesearch.MaxStepLength(1.0);

  ConjugateGradientDescent optimizer;
  optimizer.Function(&error);
  optimizer.LineSearch(&linesearch);
  optimizer.NumberOfSteps(100);
  optimizer.Epsilon(1e-6);
  optimizer.Delta(1e-12);

  // Find transformation parameters which minimize the approximation error
  optimizer.Run();

  MIRTK_DEBUG_TIMING(5, "Transformation::ApproximateDOFs");
}

// -----------------------------------------------------------------------------
void Transformation
::ApproximateDOFsGradient(const double *x,  const double *y,  const double *z, const double *t,
                          const double *dx, const double *dy, const double *dz, int no,
                          double *gradient, double w) const
{
  // Mean squared error function
  TransformationApproximationError error(const_cast<Transformation *>(this),
                                         x, y, z, t, dx, dy, dz, no);

  // Evaluate and add gradient
  double *g = new double[_NumberOfDOFs];
  error.Update(true);
  error.Gradient(g);
  for (int dof = 0; dof < _NumberOfDOFs; ++dof) gradient[dof] += w * g[dof];
  delete[] g;
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
bool Transformation::CanModifyDisplacement(int) const
{
  // The default implementation simply recomputes the entire displacement field,
  // but is not able to efficiently update the given displacement vectors
  return false;
}

// -----------------------------------------------------------------------------
void Transformation::DisplacementAfterDOFChange(int dof, double dv, GenericImage<double> &dx,
                                                double t, double t0, const WorldCoordsImage *i2w) const
{
  dx = .0;
  const double v = this->Get(dof);
  const_cast<Transformation *>(this)->Put(dof, v + dv);
  this->Displacement(dx, t, t0, i2w);
  const_cast<Transformation *>(this)->Put(dof, v);
}

// -----------------------------------------------------------------------------
void Transformation::Transform(int no, double *x, double *y, double *z, double t, double t0) const
{
  TransformationUtils::TransformPoints transform(this, x, y, z, t, t0);
  parallel_for(blocked_range<int>(0, no), transform);
}

// -----------------------------------------------------------------------------
void Transformation::Transform(int no, double *x, double *y, double *z, const double *t, double t0) const
{
  TransformationUtils::TransformPoints transform(this, x, y, z, t, t0);
  parallel_for(blocked_range<int>(0, no), transform);
}

// -----------------------------------------------------------------------------
void Transformation::Transform(WorldCoordsImage &coords, double t0) const
{
  TransformationUtils::TransformWorldCoords transform(this, coords, t0);
  parallel_for(blocked_range<int>(0, coords.NumberOfSpatialVoxels()), transform);
}

// -----------------------------------------------------------------------------
/// Voxel function used to convert transformation to dense displacement field
class TransformationToDisplacementField : public VoxelFunction
{
  const ImageAttributes &_Domain;
  const Transformation  *_Transformation;
  double                 _TargetTime;

public:

  TransformationToDisplacementField(const ImageAttributes &domain,
                                    const Transformation  *transformation,
                                    double                 t0 = 1.0)
  :
    _Domain        (domain),
    _Transformation(transformation),
    _TargetTime    (t0)
  {}

  template <class TReal>
  void operator ()(int i, int j, int k, int l, TReal *dx, TReal *dy, TReal *dz)
  {
    // Transform point into world coordinates
    double x = i, y = j, z = k;
    _Domain.LatticeToWorld(x, y, z);
    double t = _Domain.LatticeToTime(l);
    // Apply current displacement
    x += static_cast<double>(*dx);
    y += static_cast<double>(*dy);
    z += static_cast<double>(*dz);
    // Calculate displacement
    _Transformation->Displacement(x, y, z, t, _TargetTime);
    // Update displacement
    *dx += static_cast<TReal>(x);
    *dy += static_cast<TReal>(y);
    *dz += static_cast<TReal>(z);
  }
};

// -----------------------------------------------------------------------------
void Transformation::Displacement(const ImageAttributes &domain,
                                  double *dx, double *dy, double *dz) const
{
  const int no = domain.NumberOfPoints();
  if (no == 0) return;

  if (this->RequiresCachingOfDisplacements()) {

    int idx;
    GenericImage<double> disp(domain, 3);
    for (int l = 0; l < domain._t; ++l) {
      for (int k = 0; k < domain._z; ++k)
      for (int j = 0; j < domain._y; ++j)
      for (int i = 0; i < domain._x; ++i) {
        idx = domain.LatticeToIndex(i, j, k, l);
        disp(i, j, k, 0) = dx[idx];
        disp(i, j, k, 1) = dy[idx];
        disp(i, j, k, 2) = dz[idx];
      }
      this->Displacement(disp, domain.LatticeToTime(l));
      for (int k = 0; k < domain._z; ++k)
      for (int j = 0; j < domain._y; ++j)
      for (int i = 0; i < domain._x; ++i) {
        idx = domain.LatticeToIndex(i, j, k, l);
        dx[idx] = disp(i, j, k, 0);
        dy[idx] = disp(i, j, k, 1);
        dz[idx] = disp(i, j, k, 2);
      }
    }

  } else {

    GenericImage<double>              xdisp(domain, dx);
    GenericImage<double>              ydisp(domain, dy);
    GenericImage<double>              zdisp(domain, dz);
    TransformationToDisplacementField eval (domain, this);
    ParallelForEachVoxel(domain, xdisp, ydisp, zdisp, eval);

  }
}

// -----------------------------------------------------------------------------
/// Voxel function used to convert transformation to dense 3D displacement field
struct TransformationToDisplacementImage : public VoxelFunction
{
  const Transformation *_Transformation;
  BaseImage            *_Displacement;
  double                _TargetTime;
  double                _SourceTime;

  void Initialize()
  {
    _x = 0;
    _y = _Displacement->X() * _Displacement->Y() * _Displacement->Z();
    if (_Displacement->T() == 2) _z = 0;       // 2D vectors only
    else                         _z = _y + _y; // 2 * number of voxels
  }

  template <class TReal>
  void operator ()(int i, int j, int k, int, TReal *disp)
  {
    // Transform point into world coordinates
    double x = i, y = j, z = (_z != 0 ? k : .0);
    _Displacement->ImageToWorld(x, y, z);
    // Apply current displacement
    x += static_cast<double>(disp[_x]);
    y += static_cast<double>(disp[_y]);
    if (_z != 0) z += static_cast<double>(disp[_z]);
    // Calculate displacement
    _Transformation->Displacement(x, y, z, _SourceTime, _TargetTime);
    // Update displacement
    disp[_x] += static_cast<TReal>(x);
    disp[_y] += static_cast<TReal>(y);
    if (_z != 0) disp[_z] += static_cast<TReal>(z);
  }

  template <class TCoord, class TReal>
  void operator ()(int, int, int k, int, TCoord *wc, TReal *disp)
  {
    // Transform point into world coordinates
    double x = static_cast<double>(wc[_x]);
    double y = static_cast<double>(wc[_y]);
    double z = (_z != 0 ? static_cast<double>(wc[_z]) : .0);
    // Apply current displacement
    x += static_cast<double>(disp[_x]);
    y += static_cast<double>(disp[_y]);
    if (_z != 0) z += static_cast<double>(disp[_z]);
    // Calculate displacement
    _Transformation->Displacement(x, y, z, _SourceTime, _TargetTime);
    // Update displacement
    disp[_x] += static_cast<TReal>(x);
    disp[_y] += static_cast<TReal>(y);
    if (_z != 0) disp[_z] += static_cast<TReal>(z);
  }

private:
  int _x, _y, _z;
};

// -----------------------------------------------------------------------------
void Transformation::Displacement(GenericImage<double> &disp, double t0, const WorldCoordsImage *i2w) const
{
  disp = .0;
  this->Displacement(disp, disp.GetTOrigin(), t0, i2w);
}

// -----------------------------------------------------------------------------
void Transformation::Displacement(GenericImage<float> &disp, double t0, const WorldCoordsImage *i2w) const
{
  disp = .0f;
  this->Displacement(disp, disp.GetTOrigin(), t0, i2w);
}

// -----------------------------------------------------------------------------
void Transformation::Displacement(GenericImage<double> &disp, double t, double t0, const WorldCoordsImage *i2w) const
{
  if (disp.T() < 2 || disp.T() > 3) {
    cerr << "Transformation::Displacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }
  if (i2w) {
    if (i2w->T() != disp.T()) {
      cerr << "Transformation::Displacement: Coordinate map must have as many vector components (_t) as the displacement field" << endl;
      exit(1);
    }
    if (i2w->X() != disp.X() || i2w->Y() != disp.Y() || i2w->Z() != disp.Z()) {
      cerr << "Transformation::Displacement: Coordinate map must have the same size as the input/output image" << endl;
      exit(1);
    }
  }

  TransformationToDisplacementImage vf;
  vf._Displacement   = &disp;
  vf._Transformation = this;
  vf._TargetTime     = t0;
  vf._SourceTime     = t;
  vf.Initialize();

  if (i2w) ParallelForEachVoxel(disp.GetImageAttributes(), *i2w, disp, vf);
  else     ParallelForEachVoxel(disp.GetImageAttributes(),       disp, vf);
}

// -----------------------------------------------------------------------------
void Transformation::Displacement(GenericImage<float> &disp, double t, double t0, const WorldCoordsImage *i2w) const
{
  if (disp.T() < 2 || disp.T() > 3) {
    cerr << "Transformation::Displacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }
  if (i2w) {
    if (i2w->T() < 2 || i2w->T() > 3) {
      cerr << "Transformation::Displacement: Coordinate map must have either 2 or 3 vector components (_t)" << endl;
      exit(1);
    }
    if (i2w->X() != disp.X() || i2w->Y() != disp.Y() || i2w->Z() != disp.Z()) {
      cerr << "Transformation::Displacement: Coordinate map must have the same size as the input/output image" << endl;
      exit(1);
    }
  }

  TransformationToDisplacementImage vf;
  vf._Displacement   = &disp;
  vf._Transformation = this;
  vf._TargetTime     = t0;
  vf._SourceTime     = t;
  vf.Initialize();

  if (i2w) ParallelForEachVoxel(disp.GetImageAttributes(), *i2w, disp, vf);
  else     ParallelForEachVoxel(disp.GetImageAttributes(),       disp, vf);
}

// -----------------------------------------------------------------------------
void Transformation::GlobalInverse(double &x, double &y, double &z, double t, double t0) const
{
  EvaluateGlobalInverse(this, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
bool Transformation::LocalInverse(double &x, double &y, double &z, double t, double t0) const
{
  return EvaluateLocalInverse(this, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
bool Transformation::Inverse(double &x, double &y, double &z, double t, double t0) const
{
  return EvaluateInverse(this, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
/// Voxel function used to convert inverse transformation to dense displacement field
class TransformationToInverseDisplacementField : public VoxelReduction
{
  const ImageAttributes &_Domain;
  const Transformation  *_Transformation;
  double                 _TargetTime;

  mirtkReadOnlyAttributeMacro(int, NumberOfSingularPoints);

public:

  TransformationToInverseDisplacementField(const ImageAttributes &domain,
                                           const Transformation  *transformation,
                                           double                 t0 = 1.0)
  :
    _Domain        (domain),
    _Transformation(transformation),
    _TargetTime    (t0),
    _NumberOfSingularPoints(0)
  {}

  void split(const TransformationToInverseDisplacementField &)
  {
    _NumberOfSingularPoints = 0;
  }

  void join(const TransformationToInverseDisplacementField &other)
  {
    _NumberOfSingularPoints += other._NumberOfSingularPoints;
  }

  template <class TReal>
  void operator ()(int i, int j, int k, int l, TReal *dx, TReal *dy, TReal *dz)
  {
    // Transform point into world coordinates
    double x = i, y = j, z = k;
    _Domain.LatticeToWorld(x, y, z);
    double t = _Domain.LatticeToTime(l);
    // Apply current displacement
    x += static_cast<double>(*dx);
    y += static_cast<double>(*dy);
    z += static_cast<double>(*dz);
    // Calculate inverse displacement
    if (!_Transformation->InverseDisplacement(x, y, z, t, _TargetTime)) {
      ++_NumberOfSingularPoints;
    }
    // Update displacement
    *dx += static_cast<TReal>(x);
    *dy += static_cast<TReal>(y);
    *dz += static_cast<TReal>(z);
  }
};

// -----------------------------------------------------------------------------
int Transformation::InverseDisplacement(const ImageAttributes &domain,
                                        double *dx, double *dy, double *dz) const
{
  const int no = domain.NumberOfPoints();
  if (no == 0) return 0;

  int ninv = 0;

  if (this->RequiresCachingOfDisplacements()) {

    int idx;
    GenericImage<double> disp(domain, 3);
    for (int l = 0; l < domain._t; ++l) {
      for (int k = 0; k < domain._z; ++k)
      for (int j = 0; j < domain._y; ++j)
      for (int i = 0; i < domain._x; ++i) {
        idx = domain.LatticeToIndex(i, j, k, l);
        disp(i, j, k, 0) = dx[idx];
        disp(i, j, k, 1) = dy[idx];
        disp(i, j, k, 2) = dz[idx];
      }
      ninv += this->InverseDisplacement(disp, domain.LatticeToTime(l));
      for (int k = 0; k < domain._z; ++k)
      for (int j = 0; j < domain._y; ++j)
      for (int i = 0; i < domain._x; ++i) {
        idx = domain.LatticeToIndex(i, j, k, l);
        dx[idx] = disp(i, j, k, 0);
        dy[idx] = disp(i, j, k, 1);
        dz[idx] = disp(i, j, k, 2);
      }
    }

  } else {

    GenericImage<double> xdisp(domain, dx);
    GenericImage<double> ydisp(domain, dy);
    GenericImage<double> zdisp(domain, dz);
    TransformationToInverseDisplacementField eval(domain, this);
    ParallelForEachVoxel(domain, xdisp, ydisp, zdisp, eval);
    ninv = eval.NumberOfSingularPoints();

  }

  return ninv;
}

// -----------------------------------------------------------------------------
/// Voxel function used to convert inverse transformation to dense 3D displacement field
struct TransformationToInverseDisplacementImage : public VoxelReduction
{
  const Transformation *_Transformation;
  BaseImage            *_Displacement;
  double                _TargetTime;
  double                _SourceTime;
  int                   _NumberOfSingularPoints;

  void Initialize()
  {
    _x = 0;
    _y = _Displacement->X() * _Displacement->Y() * _Displacement->Z();
    if (_Displacement->T() == 2) _z = 0;       // 2D vectors only
    else                         _z = _y + _y; // 2 * number of voxels
    _NumberOfSingularPoints = 0;
  }

  void split(const TransformationToInverseDisplacementImage &)
  {
    _NumberOfSingularPoints = 0;
  }

  void join(const TransformationToInverseDisplacementImage &other)
  {
    _NumberOfSingularPoints += other._NumberOfSingularPoints;
  }

  template <class TReal>
  void operator ()(int i, int j, int k, int, TReal *disp)
  {
    // Transform point into world coordinates
    double x = i, y = j, z = (_z != 0 ? k : .0);
    _Displacement->ImageToWorld(x, y, z);
    // Apply current displacement
    x += static_cast<double>(disp[_x]);
    y += static_cast<double>(disp[_y]);
    if (_z != 0) z += static_cast<double>(disp[_z]);
    // Calculate inverse displacement
    if (!_Transformation->InverseDisplacement(x, y, z, _SourceTime, _TargetTime)) {
      ++_NumberOfSingularPoints;
    }
    // Update displacement
    disp[_x] += static_cast<TReal>(x);
    disp[_y] += static_cast<TReal>(y);
    if (_z != 0) disp[_z] += static_cast<TReal>(z);
  }

  template <class TCoord, class TReal>
  void operator ()(int, int, int k, int, TCoord *wc, TReal *disp)
  {
    // Transform point into world coordinates
    double x = static_cast<double>(wc[_x]);
    double y = static_cast<double>(wc[_y]);
    double z = (_z != 0 ? static_cast<double>(wc[_z]) : .0);
    // Apply current displacement
    x += static_cast<double>(disp[_x]);
    y += static_cast<double>(disp[_y]);
    if (_z != 0) z += static_cast<double>(disp[_z]);
    // Calculate inverse displacement
    if (!_Transformation->InverseDisplacement(x, y, z, _SourceTime, _TargetTime)) {
      ++_NumberOfSingularPoints;
    }
    // Update displacement
    disp[_x] += static_cast<TReal>(x);
    disp[_y] += static_cast<TReal>(y);
    if (_z != 0) disp[_z] += static_cast<TReal>(z);
  }

private:
  int _x, _y, _z;
};

// -----------------------------------------------------------------------------
int Transformation::InverseDisplacement(GenericImage<double> &disp, double t0, const WorldCoordsImage *i2w) const
{
  disp = .0;
  return this->InverseDisplacement(disp, disp.GetTOrigin(), t0, i2w);
}

// -----------------------------------------------------------------------------
int Transformation::InverseDisplacement(GenericImage<float> &disp, double t0, const WorldCoordsImage *i2w) const
{
  disp = .0f;
  return this->InverseDisplacement(disp, disp.GetTOrigin(), t0, i2w);
}

// -----------------------------------------------------------------------------
int Transformation::InverseDisplacement(GenericImage<double> &disp, double t, double t0, const WorldCoordsImage *i2w) const
{
  if (disp.T() < 2 || disp.T() > 3) {
    cerr << "Transformation::InverseDisplacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }
  if (i2w) {
    if (i2w->T() < 2 || i2w->T() > 3) {
      cerr << "Transformation::InverseDisplacement: Coordinate map must have either 2 or 3 vector components (_t)" << endl;
      exit(1);
    }
    if (i2w->X() != disp.X() || i2w->Y() != disp.Y() || i2w->Z() != disp.Z()) {
      cerr << "Transformation::InverseDisplacement: Coordinate map must have the same size as the input/output image" << endl;
      exit(1);
    }
  }

  TransformationToInverseDisplacementImage vf;
  vf._Displacement   = &disp;
  vf._Transformation = this;
  vf._TargetTime     = t0;
  vf._SourceTime     = t;
  vf.Initialize();

  if (i2w) ParallelForEachVoxel(disp.GetImageAttributes(), *i2w, disp, vf);
  else     ParallelForEachVoxel(disp.GetImageAttributes(),       disp, vf);

  return vf._NumberOfSingularPoints;
}

// -----------------------------------------------------------------------------
int Transformation::InverseDisplacement(GenericImage<float> &disp, double t, double t0, const WorldCoordsImage *i2w) const
{
  if (disp.T() < 2 || disp.T() > 3) {
    cerr << "Transformation::InverseDisplacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }
  if (i2w) {
    if (i2w->T() < 2 || i2w->T() > 3) {
      cerr << "Transformation::InverseDisplacement: Coordinate map must have either 2 or 3 vector components (_t)" << endl;
      exit(1);
    }
    if (i2w->X() != disp.X() || i2w->Y() != disp.Y() || i2w->Z() != disp.Z()) {
      cerr << "Transformation::InverseDisplacement: Coordinate map must have the same size as the input/output image" << endl;
      exit(1);
    }
  }

  TransformationToInverseDisplacementImage vf;
  vf._Displacement   = &disp;
  vf._Transformation = this;
  vf._TargetTime     = t0;
  vf._SourceTime     = t;
  vf.Initialize();

  if (i2w) ParallelForEachVoxel(disp.GetImageAttributes(), *i2w, disp, vf);
  else     ParallelForEachVoxel(disp.GetImageAttributes(),       disp, vf);

  return vf._NumberOfSingularPoints;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
class TransformationParametricGradientBody
{
public:
  const Transformation       *_Transformation;
  const GenericImage<double> *_Input;
  const WorldCoordsImage     *_WorldCoords;
  double                     *_Output;
  double                      _Weight;

  double _t;            ///< Time corrresponding to input gradient image (in ms)
  double _t0;           ///< Second time argument for velocity-based transformations

private:

  int    _X;            ///< Number of voxels along x axis
  int    _Y;            ///< Number of voxels along y axis
  int    _NumberOfDOFs; ///< Number of transformation parameters

public:

  // ---------------------------------------------------------------------------
  /// Default constructor
  TransformationParametricGradientBody()
  :
    _Transformation(NULL),
    _Input         (NULL),
    _WorldCoords   (NULL),
    _Output        (NULL),
    _Weight        ( 1.0),
    _t             ( 0.0),
    _t0            (-1.0),
    _X             (0),
    _Y             (0),
    _NumberOfDOFs  (0)
  {}

  // ---------------------------------------------------------------------------
  /// Split constructor
  TransformationParametricGradientBody(TransformationParametricGradientBody &rhs, split)
  :
    _Transformation(rhs._Transformation),
    _Input         (rhs._Input),
    _WorldCoords   (rhs._WorldCoords),
    _Weight        (rhs._Weight),
    _t             (rhs._t),
    _t0            (rhs._t0),
    _X             (rhs._X),
    _Y             (rhs._Y),
    _NumberOfDOFs  (rhs._NumberOfDOFs)
  {
    _Output = new double[_NumberOfDOFs];
    memset(_Output, 0, _NumberOfDOFs * sizeof(double));
  }

  // ---------------------------------------------------------------------------
  /// Destructor
  ~TransformationParametricGradientBody() {}

  // ---------------------------------------------------------------------------
  /// Join results of right-hand body with this body
  void join(TransformationParametricGradientBody &rhs)
  {
    for (int dof = 0; dof < _NumberOfDOFs; ++dof) {
      _Output[dof] += rhs._Output[dof];
    }
    delete[] rhs._Output;
    rhs._Output = NULL;
  }

  // ---------------------------------------------------------------------------
  /// Calculates the gradient of the similarity term w.r.t. the transformation
  /// parameters for each discrete voxel of the specified image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int i1 = re.cols ().begin();
    const int j1 = re.rows ().begin();
    const int k1 = re.pages().begin();
    const int i2 = re.cols ().end();
    const int j2 = re.rows ().end();
    const int k2 = re.pages().end();

    double jac[3]; // Jacobian of transformation w.r.t DoF

    // s1=1
    const int s2 =  _X - (i2 - i1);
    const int s3 = (_Y - (j2 - j1)) * _X;

    // Non-parametric/voxelwise gradient
    const double *gx = _Input->Data(i1, j1, k1, 0);
    const double *gy = _Input->Data(i1, j1, k1, 1);
    const double *gz = _Input->Data(i1, j1, k1, 2);

    // With (transformed) pre-computed world coordinates
    if (_WorldCoords) {
      // Loop over image voxels
      const double *wx = _WorldCoords->Data(i1, j1, k1, 0);
      const double *wy = _WorldCoords->Data(i1, j1, k1, 1);
      const double *wz = _WorldCoords->Data(i1, j1, k1, 2);
      for (int k = k1; k < k2; ++k, wx += s3, wy += s3, wz += s3, gx += s3, gy += s3, gz += s3)
      for (int j = j1; j < j2; ++j, wx += s2, wy += s2, wz += s2, gx += s2, gy += s2, gz += s2)
      for (int i = i1; i < i2; ++i, wx +=  1, wy +=  1, wz +=  1, gx +=  1, gy +=  1, gz +=  1) {
        // Check whether reference point is valid
        if (*gx || *gy || *gz) {
          // Loop over DoFs
          for (int dof = 0; dof < _NumberOfDOFs; ++dof) {
            // Check status of DoF
            if (_Transformation->GetStatus(dof) == Active) {
              // Convert non-parametric gradient into parametric gradient
              _Transformation->JacobianDOFs(jac, dof, *wx, *wy, *wz, _t, _t0);
              _Output[dof] += _Weight * (jac[0] * (*gx) + jac[1] * (*gy) + jac[2] * (*gz));
            }
          }
        }
      }
    // Without pre-computed world coordinates
    } else {
      double x, y, z;
      // Loop over image voxels
      for (int k = k1; k < k2; ++k, gx += s3, gy += s3, gz += s3)
      for (int j = j1; j < j2; ++j, gx += s2, gy += s2, gz += s2)
      for (int i = i1; i < i2; ++i, gx +=  1, gy +=  1, gz +=  1) {
        // Check whether reference point is valid
        if (*gx || *gy || *gz) {
          // Convert voxel to world coordinates
          x = i, y = j, z = k;
          _Input->ImageToWorld(x, y, z);
          // Loop over DoFs
          for (int dof = 0; dof < _NumberOfDOFs; ++dof) {
            // Check status of DoF
            if (_Transformation->GetStatus(dof) == Active) {
              // Convert non-parametric gradient into parametric gradient
              _Transformation->JacobianDOFs(jac, dof, x, y, z, _t, _t0);
              _Output[dof] += _Weight * (jac[0] * (*gx) + jac[1] * (*gy) + jac[2] * (*gz));
            }
          }
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  void operator ()()
  {
    // Initialize members to often accessed data
    _X            = _Input->X();
    _Y            = _Input->Y();
    const int _Z  = _Input->Z();
    _NumberOfDOFs = _Transformation->NumberOfDOFs();

    // Check input
    if (_WorldCoords && (_WorldCoords->X() != _X ||
                         _WorldCoords->Y() != _Y ||
                         _WorldCoords->Z() != _Z ||
                         _WorldCoords->T() !=  3)) {
      cerr << "Transformation::ParametricGradient: Invalid world coordinates map" << endl;
      exit(1);
    }

    // Calculate parametric gradient
    //
    // FIXME: In case of affine registration with NMI as image similarity measure,
    //        the parallel reduction seems to introduce numeric differences in
    //        the summation of the individual gradient contributions. This
    //        results in different results depending on the split of the image
    //        domain into smaller pieces.
    //        https://gitlab.doc.ic.ac.uk/BioMedIA/MIRTK/issues/31
    blocked_range3d<int> voxels(0, _Z, 0, _Y, 0, _X);
//    parallel_reduce(voxels, *this);
    this->operator()(voxels);
  }

}; // TransformationParametricGradientBody

// -----------------------------------------------------------------------------
class TransformationPointWiseParametricGradientBody
{
public:

  const Transformation       *_Transformation;
  const PointSet             *_Point;
  const Vector3D<double>     *_Input;
  double                     *_Output;
  double                      _Weight;

  double _t;  ///< Time corrresponding to input gradient image (in ms)
  double _t0; ///< Second time argument for velocity-based transformations

  // ---------------------------------------------------------------------------
  /// Default constructor
  TransformationPointWiseParametricGradientBody()
  :
    _Transformation(NULL),
    _Point         (NULL),
    _Input         (NULL),
    _Output        (NULL),
    _Weight        ( 1.0),
    _t             ( 0.0),
    _t0            (-1.0)
  {}

  // ---------------------------------------------------------------------------
  /// Destructor
  ~TransformationPointWiseParametricGradientBody() {}

  // ---------------------------------------------------------------------------
  /// Calculates the gradient of the similarity term w.r.t. the transformation
  /// parameters for each point in the specified index range
  void operator ()(const blocked_range<int> &re) const
  {
    double jac[3]; // Jacobian of transformation w.r.t DoF

    // Loop over DoFs
    for (int dof = re.begin(); dof != re.end(); ++dof) {
      // Check status of DoF
      if (_Transformation->GetStatus(dof) == Active) {
        // Loop over points
        const Vector3D<double> *g = _Input;
        for (int idx = 0; idx < _Point->Size(); ++idx, ++g) {
          // Check whether reference point is valid
          if (g->_x != .0 || g->_y != .0 || g->_z != .0) {
            // Gradient vector position
            const Point &p = (*_Point)(idx);
            // Convert non-parametric gradient into parametric gradient
            _Transformation->JacobianDOFs(jac, dof, p._x, p._y, p._z, _t, _t0);
            _Output[dof] += _Weight * (jac[0] * (g->_x) + jac[1] * (g->_y) + jac[2] * (g->_z));
          }
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  void operator ()()
  {
    blocked_range<int> dofs(0, _Transformation->NumberOfDOFs());
    parallel_for(dofs, *this);
  }

}; // TransformationPointWiseParametricGradientBody

// -----------------------------------------------------------------------------
void Transformation::ParametricGradient(const GenericImage<double> *in, double *out,
                                        const WorldCoordsImage *i2w,
                                        const WorldCoordsImage *wc,
                                        double t0, double w) const
{
  MIRTK_START_TIMING();
  TransformationParametricGradientBody body;
  body._Transformation = this;
  body._Input          = in;
  body._Output         = out;
  body._Weight         = w;
  body._WorldCoords    = (wc ? wc : i2w);
  body._t              = in->ImageToTime(.0);
  body._t0             = t0;
  body();
  MIRTK_DEBUG_TIMING(2, "parametric gradient computation");
}

// -----------------------------------------------------------------------------
void Transformation::ParametricGradient(const PointSet &pos, const Vector3D<double> *in,
                                        double *out, double t, double t0, double w) const
{
  MIRTK_START_TIMING();
  TransformationPointWiseParametricGradientBody body;
  body._Transformation = this;
  body._Point          = &pos;
  body._Input          = in;
  body._Output         = out;
  body._Weight         = w;
  body._t              = t;
  body._t0             = t0;
  body();
  MIRTK_DEBUG_TIMING(2, "point-wise parametric gradient computation");
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
bool Transformation::CheckHeader(const char *name)
{
  Cifstream from(name);
  unsigned int magic_no;
  from.ReadAsUInt(&magic_no, 1);
  from.Close();
  if (magic_no != TRANSFORMATION_MAGIC) {
    swap32((char *)&magic_no, (char *)&magic_no, 1);
  }
  return magic_no == TRANSFORMATION_MAGIC;
}

// -----------------------------------------------------------------------------
void Transformation::Read(const char *name)
{
  Cifstream from(name);
  unsigned int magic_no;
  from.ReadAsUInt(&magic_no, 1);
  if (magic_no != TRANSFORMATION_MAGIC) {
    swap32((char *)&magic_no, (char *)&magic_no, 1);
    if (magic_no == TRANSFORMATION_MAGIC) {
      from.Swapped(!from.Swapped());
    } else {
      cerr << "Transformation::Read: Not a (M)IRTK transformation file" << endl;
      exit(1);
    }
  }
  from.Seek(0);
  Read(from);
  from.Close();
}

// -----------------------------------------------------------------------------
void Transformation::Write(const char *name) const
{
  Cofstream to(name);
  Write(to);
  to.Close();
}

// -----------------------------------------------------------------------------
bool Transformation::CanRead(TransformationType type) const
{
  return (type == this->TypeOfClass());
}

// -----------------------------------------------------------------------------
Cifstream &Transformation::Read(Cifstream &from)
{
  // Read magic no. for transformations
  unsigned int magic_no;
  from.ReadAsUInt(&magic_no, 1);
  if (magic_no != TRANSFORMATION_MAGIC) {
    swap32((char *)&magic_no, (char *)&magic_no, 1);
    if (magic_no == TRANSFORMATION_MAGIC) {
      from.Swapped(!from.Swapped());
    } else {
      cerr << this->NameOfClass() << "::Read: Not a (M)IRTK transformation file" << endl;
      exit(1);
    }
  }

  // Read transformation type
  unsigned int trans_type;
  from.ReadAsUInt(&trans_type, 1);
  TransformationType format = static_cast<TransformationType>(trans_type);
  if (!this->CanRead(format)) {
    cerr << this->NameOfClass() << "::Read: Not a valid transformation type" << endl;
    exit(1);
  }

  // Read transformation parameters
  return this->ReadDOFs(from, format);
}

// -----------------------------------------------------------------------------
Cofstream &Transformation::Write(Cofstream &to) const
{
  // Write magic no. for transformations
  const unsigned int magic_no = TRANSFORMATION_MAGIC;
  to.WriteAsUInt(&magic_no, 1);

  // Write transformation type
  const unsigned int trans_type = this->TypeOfClass();
  to.WriteAsUInt(&trans_type, 1);

  // Write transformation parameters
  return this->WriteDOFs(to);
}

// -----------------------------------------------------------------------------
Cifstream &Transformation::ReadDOFs(Cifstream &from, TransformationType)
{
  // Read number of parameters
  unsigned int ndofs;
  from.ReadAsUInt(&ndofs, 1);
  if (ndofs != static_cast<unsigned int>(_NumberOfDOFs)) {
    cerr << this->NameOfClass() << "::Read: Invalid no. of parameters: " << ndofs << endl;
    exit(1);
  }

  // Read transformation parameters
  from.ReadAsDouble(this->_Param, _NumberOfDOFs);

  // TODO: Read status of transformation parameters ?
  //
  //       If stored in file as well (as is the case for FFDs already),
  //       new linear transformation type IDs have to be added because
  //       the file format of these transformations would change.

  return from;
}

// -----------------------------------------------------------------------------
Cofstream &Transformation::WriteDOFs(Cofstream &to) const
{
  // Write number of parameters
  unsigned int ndofs = _NumberOfDOFs;
  to.WriteAsUInt(&ndofs, 1);

  // Write transformation parameters
  to.WriteAsDouble(_Param, _NumberOfDOFs);

  // TODO: Write status of transformation parameters ?
  //
  //       If stored in file as well (as is the case for FFDs already),
  //       new linear transformation type IDs have to be added because
  //       the file format of these transformations would change.

  return to;
}

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
bool IsTransformation(const char *name)
{
  Cifstream from(name);
  unsigned int magic_no;
  from.ReadAsUInt(&magic_no, 1);
  from.Close();
  if (magic_no != TRANSFORMATION_MAGIC) {
    swap32((char *)&magic_no, (char *)&magic_no, 1);
  }
  return (magic_no == TRANSFORMATION_MAGIC);
}


} // namespace mirtk

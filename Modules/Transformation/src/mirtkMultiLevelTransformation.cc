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

#include <mirtkMultiLevelTransformation.h>

#include <mirtkDeallocate.h>
#include <mirtkHomogeneousTransformation.h>
#include <mirtkFreeFormTransformation.h>
#include <mirtkParallel.h>
#include <mirtkProfiling.h>


namespace mirtk {


// =============================================================================
// Declaration of inverse transformation (cf. TransformationInverse.cc)
// =============================================================================

// -----------------------------------------------------------------------------
bool EvaluateInverse(const MultiLevelTransformation *, int, int,
                     double &, double &, double &, double, double);

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
MultiLevelTransformation::MultiLevelTransformation()
:
  _NumberOfLevels(0)
{
  for (int l = 0; l < MAX_TRANS; ++l) {
    _LocalTransformation      [l] = NULL;
    _LocalTransformationStatus[l] = Passive;
  }
}

// -----------------------------------------------------------------------------
MultiLevelTransformation::MultiLevelTransformation(const RigidTransformation &t)
:
  _GlobalTransformation(t),
  _NumberOfLevels(0)
{
  for (int l = 0; l < MAX_TRANS; ++l) {
    _LocalTransformation      [l] = NULL;
    _LocalTransformationStatus[l] = Passive;
  }
}

// -----------------------------------------------------------------------------
MultiLevelTransformation::MultiLevelTransformation(const AffineTransformation &t)
:
  _GlobalTransformation(t),
  _NumberOfLevels(0)
{
  for (int l = 0; l < MAX_TRANS; ++l) {
    _LocalTransformation      [l] = NULL;
    _LocalTransformationStatus[l] = Passive;
  }
}

// -----------------------------------------------------------------------------
MultiLevelTransformation::MultiLevelTransformation(const MultiLevelTransformation &t)
:
  Transformation(t),
  _NumberOfLevels(t._NumberOfLevels)
{
  for (int l = 0; l < _NumberOfLevels; ++l) {
    _LocalTransformation[l] = dynamic_cast<FreeFormTransformation *>(Transformation::New(t._LocalTransformation[l]));
    if (_LocalTransformation[l] == NULL) {
      cerr << "MultiLevelTransformation::MultiLevelTransformation: Failed to copy local transformation at level ";
      cerr << l << " and of type " << t._LocalTransformation[l]->NameOfClass() << endl;
      exit(1);
    }
    _LocalTransformationStatus[l] = t._LocalTransformationStatus[l];
  }
}

// -----------------------------------------------------------------------------
MultiLevelTransformation::~MultiLevelTransformation()
{
  for (int l = 0; l < MAX_TRANS; ++l) Delete(_LocalTransformation[l]);
}

// =============================================================================
// Approximation
// =============================================================================

// -----------------------------------------------------------------------------
double MultiLevelTransformation
::Approximate(const ImageAttributes &, double *, double *, double *,
              int, double)
{
  cerr << this->NameOfClass() << "::Approximate: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
double MultiLevelTransformation
::Approximate(const double *, const double *, const double *,
              double       *, double       *, double       *, int,
              int, double)
{
  cerr << this->NameOfClass() << "::Approximate: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
double MultiLevelTransformation
::Approximate(const double *, const double *, const double *, const double *,
              double       *, double       *, double       *, int,
              int, double)
{
  cerr << this->NameOfClass() << "::Approximate: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
double MultiLevelTransformation
::ApproximateAsNew(const ImageAttributes &, double *, double *, double *,
                   int, double)
{
  cerr << this->NameOfClass() << "::ApproximateAsNew: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
double MultiLevelTransformation
::ApproximateAsNew(const double *, const double *, const double *,
                   double       *, double       *, double       *, int,
                   int, double)
{
  cerr << this->NameOfClass() << "::ApproximateAsNew: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
double MultiLevelTransformation
::ApproximateAsNew(const double *, const double *, const double *, const double *,
                   double       *, double       *, double       *, int,
                   int, double)
{
  cerr << this->NameOfClass() << "::ApproximateAsNew: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void MultiLevelTransformation
::ApproximateDOFs(const double *, const double *, const double *, const double *,
                  const double *, const double *, const double *, int)
{
  cerr << this->NameOfClass() << "::ApproximateDOFs: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void MultiLevelTransformation
::ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                          const double *, const double *, const double *, int,
                          double *, double) const
{
  cerr << this->NameOfClass() << "::ApproximateDOFsGradient: Not implemented" << endl;
  exit(1);
}

// =============================================================================
// Transformation parameters (DOFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool MultiLevelTransformation::CopyFrom(const Transformation *other)
{
  const HomogeneousTransformation *lin  = NULL;
  const MultiLevelTransformation  *mffd = NULL;
  const FreeFormTransformation    *affd = NULL;

  (lin  = dynamic_cast<const HomogeneousTransformation *>(other)) ||
  (mffd = dynamic_cast<const MultiLevelTransformation  *>(other)) ||
  (affd = dynamic_cast<const FreeFormTransformation    *>(other));

  if (mffd && mffd->NumberOfLevels() == 0) {
    lin  = mffd->GetGlobalTransformation();
    mffd = NULL;
  }

  if (lin) {
    this->Reset();
    _GlobalTransformation.CopyFrom(lin);
    return true;
  }
  if (affd) {
    if (_NumberOfLevels == 0) return false;
    for (int i = 0; i < _NumberOfLevels; ++i) {
      if (_LocalTransformation[i]->CopyFrom(affd)) {
        _GlobalTransformation.Reset();
        for (int j = 0; j < _NumberOfLevels; ++j) {
          if (i != j) _LocalTransformation[j]->Reset();
        }
        return true;
      }
    }
    return false;
  }
  if (mffd && mffd->NumberOfLevels() == _NumberOfLevels
           && strcmp(this->NameOfClass(), mffd->NameOfClass()) == 0) {
    if (!_GlobalTransformation.CopyFrom(mffd->GetGlobalTransformation())) {
      return false;
    }
    for (int i = 0; i < _NumberOfLevels; ++i) {
      if (!_LocalTransformation[i]->CopyFrom(mffd->GetLocalTransformation(i))) {
        this->Reset();
        return false;
      }
    }
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
double MultiLevelTransformation::DOFGradientNorm(const double *gradient) const
{
  double norm, max = .0;
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    if (!this->LocalTransformationIsActive(l)) continue;
    const FreeFormTransformation *ffd = this->GetLocalTransformation(l);
    norm = ffd->DOFGradientNorm(gradient);
    if (norm > max) max = norm;
  }
  return max;
}

// -----------------------------------------------------------------------------
int MultiLevelTransformation::NumberOfDOFs() const
{
  int ndofs = 0;
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    if (!this->LocalTransformationIsActive(l)) continue;
    ndofs += this->GetLocalTransformation(l)->NumberOfDOFs();
  }
  return ndofs;
}

// -----------------------------------------------------------------------------
void MultiLevelTransformation::Put(int idx, double value)
{
  FreeFormTransformation *ffd;
  int                         dof;
  DOFIndexToLocalTransformation(this, idx, ffd, dof);
  ffd->Put(dof, value);
  this->Changed(true);
}

// -----------------------------------------------------------------------------
double MultiLevelTransformation::Get(int idx) const
{
  const FreeFormTransformation *ffd;
  int                               dof;
  DOFIndexToLocalTransformation(this, idx, ffd, dof);
  return ffd->Get(dof);
}

// -----------------------------------------------------------------------------
void MultiLevelTransformation::Put(const DOFValue *x)
{
  FreeFormTransformation *ffd;
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    if (!this->LocalTransformationIsActive(l)) continue;
    ffd = this->GetLocalTransformation(l);
    ffd->Put(x);
    x += ffd->NumberOfDOFs();
  }
  this->Changed(true);
}

// -----------------------------------------------------------------------------
void MultiLevelTransformation::Add(const DOFValue *dx)
{
  FreeFormTransformation *ffd;
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    if (!this->LocalTransformationIsActive(l)) continue;
    ffd = this->GetLocalTransformation(l);
    ffd->Add(dx);
    dx += ffd->NumberOfDOFs();
  }
  this->Changed(true);
}

// -----------------------------------------------------------------------------
double MultiLevelTransformation::Update(const DOFValue *dx)
{
  double max_delta = .0;
  FreeFormTransformation *ffd;
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    if (!this->LocalTransformationIsActive(l)) continue;
    ffd = this->GetLocalTransformation(l);
    max_delta = max(max_delta, ffd->Update(dx));
    dx += ffd->NumberOfDOFs();
  }
  if (max_delta > .0) this->Changed(true);
  return max_delta;
}

// -----------------------------------------------------------------------------
void MultiLevelTransformation::Get(DOFValue *x) const
{
  const FreeFormTransformation *ffd;
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    if (!this->LocalTransformationIsActive(l)) continue;
    ffd = this->GetLocalTransformation(l);
    ffd->Get(x);
    x += ffd->NumberOfDOFs();
  }
}

// -----------------------------------------------------------------------------
void MultiLevelTransformation::PutStatus(int idx, DOFStatus status)
{
  FreeFormTransformation *ffd;
  int                         dof;
  DOFIndexToLocalTransformation(this, idx, ffd, dof);
  ffd->PutStatus(dof, status);
  this->Changed(true);
}

// -----------------------------------------------------------------------------
MultiLevelTransformation::DOFStatus MultiLevelTransformation::GetStatus(int idx) const
{
  const FreeFormTransformation *ffd;
  int                           dof;
  DOFIndexToLocalTransformation(this, idx, ffd, dof);
  return ffd->GetStatus(dof);
}

// -----------------------------------------------------------------------------
bool MultiLevelTransformation::HasSameDOFsAs(const Transformation *t) const
{
  const MultiLevelTransformation *mffd;
  mffd = dynamic_cast<const MultiLevelTransformation *>(t);

  if (mffd) {
    if (this->NumberOfLevels() != mffd->NumberOfLevels()) return false;
    for (int l = 0; l < this->NumberOfLevels(); ++l) {
      if (this->LocalTransformationStatus(l) != mffd->LocalTransformationStatus(l)) return false;
      const FreeFormTransformation *ffd1 = this->GetLocalTransformation(l);
      const FreeFormTransformation *ffd2 = mffd->GetLocalTransformation(l);
      if (!ffd1->HasSameDOFsAs(ffd2)) return false;
    }
    return true;
  } else {
    return (this->NumberOfLevels() == 1) &&
            this->LocalTransformationIsActive(0) &&
            this->GetLocalTransformation(0)->HasSameDOFsAs(t);
  }
}

// -----------------------------------------------------------------------------
bool MultiLevelTransformation::IsIdentity() const
{
  if (!_GlobalTransformation.IsIdentity()) return false;
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    if (!this->GetLocalTransformation(l)->IsIdentity()) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
void MultiLevelTransformation::Reset()
{
  _GlobalTransformation.Reset();
  for (int l = 0; l < _NumberOfLevels; ++l) {
    _LocalTransformation[l]->Reset();
  }
}

// -----------------------------------------------------------------------------
void MultiLevelTransformation::Clear()
{
  _GlobalTransformation.Reset();
  for (int l = 0; l < _NumberOfLevels; ++l) {
    Delete(_LocalTransformation[l]);
  }
  _NumberOfLevels = 0;
}

// =============================================================================
// Parameters (non-DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool MultiLevelTransformation::Set(const char *name, const char *value)
{
  bool known = this->GetGlobalTransformation()->Set(name, value);
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    known = this->GetLocalTransformation(l)->Set(name, value) || known;
  }
  return known;
}

// -----------------------------------------------------------------------------
ParameterList MultiLevelTransformation::Parameter() const
{
  ParameterList params = this->GetGlobalTransformation()->Parameter();
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    Insert(params, this->GetLocalTransformation(l)->Parameter());
  }
  return params;
}

// =============================================================================
// Levels
// =============================================================================

// -----------------------------------------------------------------------------
int MultiLevelTransformation::NumberOfCPs(bool excl_passive) const
{
  int n = 0;
  if (excl_passive) {
    for (int l = 0; l < _NumberOfLevels; ++l) {
      if (_LocalTransformationStatus[l] == Active) {
        n += _LocalTransformation[l]->NumberOfCPs();
      }
    }
  } else {
    for (int l = 0; l < _NumberOfLevels; ++l) {
      n += _LocalTransformation[l]->NumberOfCPs();
    }
  }
  return n;
}

// -----------------------------------------------------------------------------
int MultiLevelTransformation::NumberOfActiveCPs() const
{
  int n = 0;
  for (int l = 0; l < _NumberOfLevels; ++l) {
    if (_LocalTransformationStatus[l] == Active) {
      n += _LocalTransformation[l]->NumberOfActiveCPs();
    }
  }
  return n;
}

// -----------------------------------------------------------------------------
void MultiLevelTransformation::CheckTransformation(FreeFormTransformation *transformation) const
{
  if (transformation == NULL) {
    cerr << this->NameOfClass() << "::CheckTransformation: NULL pointer given" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
FreeFormTransformation *MultiLevelTransformation::PutLocalTransformation(FreeFormTransformation *transformation, int i)
{
  this->CheckTransformation(transformation);
  if (i < 0 || i >= _NumberOfLevels) {
    cerr << this->NameOfClass() << "::PutLocalTransformation: No such transformation: " << i << endl;
    exit(1);
  }
  FreeFormTransformation * const old = _LocalTransformation[i];
  // Copy status of transformation if it is part of stack at different
  // location. This is the case when this method is used to rearrange the
  // local transformations within the stack.
  for (int l = 0; l < _NumberOfLevels; ++l) {
    if (transformation == _LocalTransformation[l]) {
      _LocalTransformationStatus[i] = _LocalTransformationStatus[l];
    }
  }
  // Set new transformation
  _LocalTransformation[i] = transformation;
  return old;
}

// -----------------------------------------------------------------------------
void MultiLevelTransformation::PushLocalTransformation(FreeFormTransformation *transformation)
{
  this->CheckTransformation(transformation);
  if (_NumberOfLevels == MAX_TRANS) {
    cerr << "MultiLevelTransformation::PushLocalTransformation: Stack overflow" << endl;
    exit(1);
  }
  // Change status of transformation that is currently on top of stack
  // to passive if no more than one transformation was active before
  int nactive = 0;
  for (int l = 0; l < _NumberOfLevels; ++l) {
    if (_LocalTransformationStatus[l] == Active) ++nactive;
  }
  if (nactive == 1) {
    _LocalTransformationStatus[_NumberOfLevels-1] = Passive;
  }
  // Add transformation at top of stack
  _LocalTransformation      [_NumberOfLevels] = transformation;
  _LocalTransformationStatus[_NumberOfLevels] = Active;
  _NumberOfLevels++;
}

// -----------------------------------------------------------------------------
void MultiLevelTransformation::InsertLocalTransformation(FreeFormTransformation *transformation, int pos)
{
  this->CheckTransformation(transformation);
  if (_NumberOfLevels == MAX_TRANS) {
    cerr << "MultiLevelTransformation::InsertLocalTransformation: Stack overflow" << endl;
    exit(1);
  }
  _NumberOfLevels++;
  for (int l = pos; l < _NumberOfLevels; ++l) {
    _LocalTransformation      [l+1] = _LocalTransformation      [l];
    _LocalTransformationStatus[l+1] = _LocalTransformationStatus[l];
  }
  _LocalTransformation      [pos] = transformation;
  _LocalTransformationStatus[pos] = Passive;
}

// -----------------------------------------------------------------------------
FreeFormTransformation *MultiLevelTransformation::PopLocalTransformation()
{
  FreeFormTransformation *localTransformation = NULL;
  if (_NumberOfLevels > 0) {
    // Change status of transformation that is now on top of stack
    // to active if no more than one transformation was active before
    int nactive = 0;
    for (int l = 0; l < _NumberOfLevels; ++l) {
      if (_LocalTransformationStatus[l] == Active) ++nactive;
    }
    if (_NumberOfLevels > 1 && nactive == 1) {
      _LocalTransformationStatus[_NumberOfLevels-2] = Active;
    }
    // Remove transformation at top of stack
    _NumberOfLevels--;
    localTransformation = _LocalTransformation[_NumberOfLevels];
    _LocalTransformation      [_NumberOfLevels] = NULL;
    _LocalTransformationStatus[_NumberOfLevels] = Passive;
  }
  return localTransformation;
}

// -----------------------------------------------------------------------------
FreeFormTransformation *MultiLevelTransformation::RemoveLocalTransformation(int pos)
{
  FreeFormTransformation *localTransformation;
  if (0 <= pos && pos < _NumberOfLevels) {
    localTransformation = _LocalTransformation[pos];
    for (int l = pos; l < _NumberOfLevels; ++l) {
      _LocalTransformation      [l] = _LocalTransformation      [l+1];
      _LocalTransformationStatus[l] = _LocalTransformationStatus[l+1];
    }
    _NumberOfLevels--;
    _LocalTransformation      [_NumberOfLevels] = NULL;
    _LocalTransformationStatus[_NumberOfLevels] = Passive;
  } else {
    cerr << "MultiLevelTransformation::RemoveLocalTransformation: No such transformation: " << pos << endl;
    exit(1);
  }
  return localTransformation;
}

// -----------------------------------------------------------------------------
void MultiLevelTransformation::CombineLocalTransformation()
{
  cerr << this->NameOfClass() << "::CombineLocalTransformation: Not implemented for this transformation" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void MultiLevelTransformation::MergeGlobalIntoLocalDisplacement()
{
  cerr << this->NameOfClass() << "::MergeGlobalIntoLocalDisplacement: Not implemented for this transformation" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void MultiLevelTransformation::InterpolateGlobalDisplacement(FreeFormTransformation *f)
{
  double x, y, z;

  const int no = f->NumberOfDOFs() / 3;

  double *dx = new double[no];
  double *dy = new double[no];
  double *dz = new double[no];

  int idx = 0;
  for (int l = 0; l < f->T(); ++l)
  for (int k = 0; k < f->Z(); ++k)
  for (int j = 0; j < f->Y(); ++j)
  for (int i = 0; i < f->X(); ++i) {
    x = i, y = j, z = k;
    f->LatticeToWorld(x, y, z);
    _GlobalTransformation.Displacement(x, y, z);
    dx[idx] = x;
    dy[idx] = y;
    dz[idx] = z;
    idx++;
  }

  f->Interpolate(dx, dy, dz);

  delete[] dx;
  delete[] dy;
  delete[] dz;

  if (f->X() < 4 || f->Y() < 4 || f->Z() < 4) {
    cerr << "MultiLevelTransformation::InterpolateGlobalDisplacement: ";
    cerr << "Very small lattice for interpolation. Result likely to be inaccurate." << endl;
    return;
  }

  double totVol =  f->X()    * f->GetXSpacing() +  f->Y()    * f->GetYSpacing() +  f->Z()    * f->GetZSpacing();
  double effVol = (f->X()-4) * f->GetXSpacing() + (f->Y()-4) * f->GetYSpacing() + (f->Z()-4) * f->GetZSpacing();
/*
  Not sure if the temporal lattice should be considered here as the global transformation is only defined in 3D.
  Probably this method is best only applied to 3D MFFD's, but not 3D+t MFFD's.
  -as12321

  if (f->GetT() > 1) {
    totVol +=  f->T()    * f->GetTSpacing();
    effVol += (f->T()-4) * f->GetTSpacing();
  }
*/
  cout << "MultiLevelTransformation::InterpolateGlobalDisplacement: ";
  cout << "Accurate interpolation of affine transformation over ";
  printf("% .1f %% of lattice volume\n", 100.0 * effVol / totVol);
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
bool MultiLevelTransformation
::Inverse(int m, int n, double &x, double &y, double &z, double t, double t0) const
{
  return EvaluateInverse(this, m, n, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
class MultiLevelTransformationToDisplacementField
{
public:

  const MultiLevelTransformation *_Transformation;
  GenericImage<VoxelType>        *_Output;
  double                          _TargetTime;
  double                          _SourceTime;
  int                             _M, _N;

  // ---------------------------------------------------------------------------
  /// Evaluate transformation at displacement field voxels in 2D
  void operator() (const blocked_range2d<int>& r) const
  {
    double x,  y, z;
    double dx, dy;

    for (int j = r.rows().begin(); j != r.rows().end(); ++j)
    for (int i = r.cols().begin(); i != r.cols().end(); ++i) {
      // Transform point into world coordinates
      x = i, y = j, z = 0;
      _Output->ImageToWorld(x, y, z);
      // Apply current displacement
      dx = _Output->GetAsDouble(i, j, 0, 0);
      dy = _Output->GetAsDouble(i, j, 0, 1);
      x += dx, y += dy;
      // Calculate displacement
      _Transformation->Displacement(_M, _N, x, y, z, _SourceTime, _TargetTime);
      // Update displacement
      _Output->PutAsDouble(i, j, 0, 0, x + dx);
      _Output->PutAsDouble(i, j, 0, 1, y + dy);
    }
  }

  // ---------------------------------------------------------------------------
  /// Evaluate transformation at displacement field voxels in 3D
  void operator() (const blocked_range3d<int>& r) const
  {
    double x,  y,  z;
    double dx, dy, dz;

    for (int k = r.pages().begin(); k != r.pages().end(); ++k)
    for (int j = r.rows ().begin(); j != r.rows ().end(); ++j)
    for (int i = r.cols ().begin(); i != r.cols ().end(); ++i) {
      // Transform point into world coordinates
      x = i, y = j, z = k;
      _Output->ImageToWorld(x, y, z);
      // Apply current displacement
      dx = _Output->GetAsDouble(i, j, k, 0);
      dy = _Output->GetAsDouble(i, j, k, 1);
      dz = _Output->GetAsDouble(i, j, k, 2);
      x += dx, y += dy, z += dz;
      // Calculate displacement
      _Transformation->Displacement(_M, _N, x, y, z, _SourceTime, _TargetTime);
      // Update displacement
      _Output->PutAsDouble(i, j, k, 0, x + dx);
      _Output->PutAsDouble(i, j, k, 1, y + dy);
      _Output->PutAsDouble(i, j, k, 2, z + dz);
    }
  }

};

// -----------------------------------------------------------------------------
template <class VoxelType>
class MultiLevelTransformationToInverseDisplacementField
{
public:

  const MultiLevelTransformation *_Transformation;
  GenericImage<VoxelType>        *_Output;
  double                          _TargetTime;
  double                          _SourceTime;
  int                             _M, _N;
  int                             _NumberOfSingularPoints;

  // ---------------------------------------------------------------------------
  /// Default constructor
  MultiLevelTransformationToInverseDisplacementField()
  :
    _Transformation(NULL),
    _Output(NULL),
    _TargetTime(.0),
    _SourceTime(.0),
    _M(0), _N(0),
    _NumberOfSingularPoints(0)
  {}

  // ---------------------------------------------------------------------------
  /// Split constructor
  MultiLevelTransformationToInverseDisplacementField(
    const MultiLevelTransformationToInverseDisplacementField &other, split
  ) :
    _Transformation(other._Transformation),
    _Output(other._Output),
    _TargetTime(other._TargetTime),
    _SourceTime(other._SourceTime),
    _M(other._M), _N(other._N),
    _NumberOfSingularPoints(0)
  {}

  // ---------------------------------------------------------------------------
  /// Join results
  void join(const MultiLevelTransformationToInverseDisplacementField &other)
  {
    _NumberOfSingularPoints += other._NumberOfSingularPoints;
  }

  // ---------------------------------------------------------------------------
  /// Evaluate transformation at displacement field voxels in 2D
  void operator() (const blocked_range2d<int>& r)
  {
    double x,  y, z;
    double dx, dy;

    for (int j = r.rows().begin(); j != r.rows().end(); ++j)
    for (int i = r.cols().begin(); i != r.cols().end(); ++i) {
      // Transform point into world coordinates
      x = i, y = j, z = 0;
      _Output->ImageToWorld(x, y, z);
      // Apply current displacement
      dx = _Output->GetAsDouble(i, j, 0, 0);
      dy = _Output->GetAsDouble(i, j, 0, 1);
      x += dx, y += dy;
      // Calculate inverse displacement
      if (!_Transformation->InverseDisplacement(_M, _N, x, y, z, _SourceTime, _TargetTime)) {
        ++_NumberOfSingularPoints;
      }
      // Update displacement
      _Output->PutAsDouble(i, j, 0, 0, x + dx);
      _Output->PutAsDouble(i, j, 0, 1, y + dy);
    }
  }

  // ---------------------------------------------------------------------------
  /// Evaluate transformation at displacement field voxels in 3D
  void operator() (const blocked_range3d<int>& r)
  {
    double x,  y,  z;
    double dx, dy, dz;

    for (int k = r.pages().begin(); k != r.pages().end(); ++k)
    for (int j = r.rows ().begin(); j != r.rows ().end(); ++j)
    for (int i = r.cols ().begin(); i != r.cols ().end(); ++i) {
      // Transform point into world coordinates
      x = i, y = j, z = k;
      _Output->ImageToWorld(x, y, z);
      // Apply current displacement
      dx = _Output->GetAsDouble(i, j, k, 0);
      dy = _Output->GetAsDouble(i, j, k, 1);
      dz = _Output->GetAsDouble(i, j, k, 2);
      x += dx, y += dy, z += dz;
      // Calculate inverse displacement
      if (!_Transformation->InverseDisplacement(_M, _N, x, y, z, _SourceTime, _TargetTime)) {
        ++_NumberOfSingularPoints;
      }
      // Update displacement
      _Output->PutAsDouble(i, j, k, 0, x + dx);
      _Output->PutAsDouble(i, j, k, 1, y + dy);
      _Output->PutAsDouble(i, j, k, 2, z + dz);
    }
  }

};

// -----------------------------------------------------------------------------
void MultiLevelTransformation::Displacement(int m, int n, GenericImage<double> &disp, double t, double t0, const WorldCoordsImage *) const
{
  if (disp.T() < 2 || disp.T() > 3) {
    cerr << "MultiLevelTransformation::Displacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

  MultiLevelTransformationToDisplacementField<double> eval;
  eval._Transformation = this;
  eval._Output         = &disp;
  eval._TargetTime     = t0;
  eval._SourceTime     = t;
  eval._M              = m;
  eval._N              = n;

  if (disp.Z() <= 1) {
    blocked_range2d<int> voxels(0, disp.Y(), 0, disp.X());
    parallel_for(voxels, eval);
  } else {
    blocked_range3d<int> voxels(0, disp.Z(), 0, disp.Y(), 0, disp.X());
    parallel_for(voxels, eval);
  }
}

// -----------------------------------------------------------------------------
void MultiLevelTransformation::Displacement(int m, int n, GenericImage<float> &disp, double t, double t0, const WorldCoordsImage *) const
{
  if (disp.T() < 2 || disp.T() > 3) {
    cerr << "MultiLevelTransformation::Displacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

  MultiLevelTransformationToDisplacementField<float> eval;
  eval._Transformation = this;
  eval._Output         = &disp;
  eval._TargetTime     = t0;
  eval._SourceTime     = t;
  eval._M              = m;
  eval._N              = n;

  if (disp.Z() <= 1) {
    blocked_range2d<int> voxels(0, disp.Y(), 0, disp.X());
    parallel_for(voxels, eval);
  } else {
    blocked_range3d<int> voxels(0, disp.Z(), 0, disp.Y(), 0, disp.X());
    parallel_for(voxels, eval);
  }
}

// -----------------------------------------------------------------------------
int MultiLevelTransformation::InverseDisplacement(int m, int n, GenericImage<double> &disp, double t, double t0, const WorldCoordsImage *) const
{
  if (disp.T() < 2 || disp.T() > 3) {
    cerr << "MultiLevelTransformation::InverseDisplacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

  MultiLevelTransformationToInverseDisplacementField<double> eval;
  eval._Transformation = this;
  eval._Output         = &disp;
  eval._TargetTime     = t0;
  eval._SourceTime     = t;
  eval._M              = m;
  eval._N              = n;

  if (disp.Z() <= 1) {
    blocked_range2d<int> voxels(0, disp.Y(), 0, disp.X());
    parallel_reduce(voxels, eval);
  } else {
    blocked_range3d<int> voxels(0, disp.Z(), 0, disp.Y(), 0, disp.X());
    parallel_reduce(voxels, eval);
  }

  return eval._NumberOfSingularPoints;
}

// -----------------------------------------------------------------------------
int MultiLevelTransformation::InverseDisplacement(int m, int n, GenericImage<float> &disp, double t, double t0, const WorldCoordsImage *) const
{
  if (disp.T() < 2 || disp.T() > 3) {
    cerr << "MultiLevelTransformation::InverseDisplacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

  MultiLevelTransformationToInverseDisplacementField<float> eval;
  eval._Transformation = this;
  eval._Output         = &disp;
  eval._TargetTime     = t0;
  eval._SourceTime     = t;
  eval._M              = m;
  eval._N              = n;

  if (disp.Z() <= 1) {
    blocked_range2d<int> voxels(0, disp.Y(), 0, disp.X());
    parallel_reduce(voxels, eval);
  } else {
    blocked_range3d<int> voxels(0, disp.Z(), 0, disp.Y(), 0, disp.X());
    parallel_reduce(voxels, eval);
  }

  return eval._NumberOfSingularPoints;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
void MultiLevelTransformation
::JacobianDOFs(double jac[3], int idx, double x, double y, double z, double t, double t0) const
{
  const FreeFormTransformation *ffd;
  int                           dof;
  DOFIndexToLocalTransformation(this, idx, ffd, dof);
  ffd->JacobianDOFs(jac, dof, x, y, z, t, t0);
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void MultiLevelTransformation::Print(Indent indent) const
{
  Indent subindent(indent + 1);

  // Print global transformation
  cout << indent << "Global transformation:" << endl;
  this->GetGlobalTransformation()->Print(subindent);

  // Print local transformations
  for (int l = 0; l < _NumberOfLevels; ++l) {
    cout << indent << "Local transformation";
    cout << (this->LocalTransformationIsActive(l) ? " (active)" : " (passive)");
    cout << ":" << endl;
    this->GetLocalTransformation(l)->Print(subindent);
  }
}

// -----------------------------------------------------------------------------
Cifstream &MultiLevelTransformation::ReadDOFs(Cifstream &from, TransformationType)
{
  unsigned int magic_no, trans_type;

  // Delete old local transformations
  for (int l = 0; l < _NumberOfLevels; ++l) {
    Delete(_LocalTransformation[l]);
  }

  // Read number of local transformations
  from.ReadAsInt(&_NumberOfLevels, 1);

  // Read global transformation
  this->GetGlobalTransformation()->Read(from);

  // Read local transformations
  for (int l = 0; l < _NumberOfLevels; ++l) {

    // Remember current file location
    int offset = from.Tell();

    // Read magic no. for transformations
    from.ReadAsUInt(&magic_no, 1);
    if (magic_no != TRANSFORMATION_MAGIC) {
      cerr << this->NameOfClass() << "::Read: Not a valid transformation found at file offset " << offset << endl;
      exit(1);
    }

    // Read transformation type
    from.ReadAsUInt(&trans_type, 1);

    // Instantiate new local transformation
    Transformation *ffd = Transformation::New(static_cast<TransformationType>(trans_type));
    _LocalTransformation[l] = dynamic_cast<FreeFormTransformation *>(ffd);
    if (_LocalTransformation[l] == NULL) {
      delete ffd;
      cerr << this->NameOfClass() << "::Read: Not a valid FFD (ID " << trans_type << ") found at file offset " << offset << endl;
      exit(1);
    }

    // Read local transformation
    from.Seek(offset);
    _LocalTransformation[l]->Read(from);
  }

  return from;
}

// -----------------------------------------------------------------------------
Cofstream &MultiLevelTransformation::WriteDOFs(Cofstream &to) const
{
  // Write number of levels
  to.WriteAsInt(&_NumberOfLevels, 1);

  // Write global transformation
  _GlobalTransformation.Write(to);

  // Write local transformations
  for (int l = 0; l < _NumberOfLevels; ++l) _LocalTransformation[l]->Write(to);

  return to;
}


} // namespace mirtk

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

#ifndef MIRTK_FreeFormTransformation4D_H
#define MIRTK_FreeFormTransformation4D_H

#include "mirtk/FreeFormTransformation.h"

#include <cstdlib>
#include <iostream>


namespace mirtk {


/**
 * Base class for 4D free-form transformations.
 */
class FreeFormTransformation4D : public FreeFormTransformation
{
  mirtkAbstractTransformationMacro(FreeFormTransformation4D);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  FreeFormTransformation4D(CPInterpolator &);

  /// Copy Constructor
  FreeFormTransformation4D(const FreeFormTransformation4D &, CPInterpolator &);

public:

  /// Destructor
  virtual ~FreeFormTransformation4D();

  // Import other Initialize overloades from base class
  using FreeFormTransformation::Initialize;

  /// Initialize free-form transformation
  virtual void Initialize(const ImageAttributes &);

  // ---------------------------------------------------------------------------
  // Derivatives
  using FreeFormTransformation::JacobianDOFs;
  using FreeFormTransformation::ParametricGradient;

  /// Calculates the Jacobian of the transformation w.r.t a control point
  virtual void JacobianDOFs(Matrix &, int, int, int, int,
                            double, double, double, double, double = NaN) const;

  /// Calculates the Jacobian of the transformation w.r.t a control point
  virtual void JacobianDOFs(Matrix &, int, double, double, double, double, double = NaN) const;

  /// Calculates the Jacobian of the transformation w.r.t a control point
  virtual void JacobianDOFs(double [3], int, int, int, int, double, double, double, double, double = NaN) const;

  /// Calculates the Jacobian of the transformation w.r.t a transformation parameter
  virtual void JacobianDOFs(double [3], int, double, double, double, double, double = NaN) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation
  ///
  /// The temporal coordinate t used for the computation of the Jacobian of the transformation
  /// w.r.t the transformation parameters for all spatial (x, y, z) voxel coordinates in world
  /// units, is assumed to correspond to the temporal origin of the given gradient image.
  /// For 4D transformations parameterized by velocities, a second time for the upper integration
  /// bound can be provided as last argument to this method. This last argument is ignored by
  /// transformations parameterized by displacements.
  virtual void ParametricGradient(const GenericImage<double> *, double *,
                                  const WorldCoordsImage *,
                                  const WorldCoordsImage *,
                                  double = NaN, double = 1) const;

  // ---------------------------------------------------------------------------
  // I/O

protected:

  /// Reads transformation parameters from a file stream
  virtual Cifstream &ReadDOFs(Cifstream &, TransformationType);

  /// Writes transformation parameters to a file stream
  virtual Cofstream &WriteDOFs(Cofstream &) const;

public:

  // ---------------------------------------------------------------------------
  // Deprecated

  /// \deprecated Use overloaded PutStatus method instead.
  void PutStatusCP(int, int, int, int, DOFStatus, DOFStatus, DOFStatus);

  /// \deprecated Use overloaded GetStatus method instead.
  void GetStatusCP(int, int, int, int, DOFStatus &, DOFStatus &, DOFStatus &) const;

  /// \deprecated Use overloaded BoundingBox method instead.
  ///             Note, however, that the new function swaps coordinates to
  ///             enforce p1[i] <= p2[i]. This is not done by this function.
  void BoundingBoxCP(int, Point &, Point &, double = 1) const;

  /// \deprecated Use overloaded BoundingBox method instead.
  void BoundingBoxImage(const GreyImage *, int, int &, int &, int &,
                                                    int &, int &, int &, double = 1) const;
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
inline void FreeFormTransformation4D
::JacobianDOFs(Matrix &, int, int, int, int, double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::JacobianDOFs: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void FreeFormTransformation4D
::JacobianDOFs(Matrix &jac, int cp, double x, double y, double z, double t, double t0) const
{
  int ci, cj, ck, cl;
  this->IndexToLattice(cp, ci, cj, ck, cl);
  return this->JacobianDOFs(jac, ci, cj, ck, cl, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
inline void FreeFormTransformation4D
::JacobianDOFs(double jac[3], int ci, int cj, int ck, int cl, double x, double y, double z, double t, double t0) const
{
  Matrix tmp(3, 3);
  this->JacobianDOFs(tmp, ci, cj, ck, cl, x, y, z, t, t0);
  if (tmp(0, 0) != .0 || tmp(0, 2) != .0 || tmp(1, 0) != .0 || tmp(1, 2) != .0 || tmp(2, 0) != .0 || tmp(2, 1) != .0) {
    cerr << "FreeFormTransformation4D::JacobianDOFs: Jacobian is full 3x3 matrix" << endl;
    exit(1);
  }
  jac[0] = tmp(0, 0);
  jac[1] = tmp(1, 1);
  jac[2] = tmp(2, 2);
}

// -----------------------------------------------------------------------------
inline void FreeFormTransformation4D::JacobianDOFs(double jac[3], int dof, double x, double y, double z, double t, double t0) const
{
  int ci, cj, ck, cl;
  this->IndexToLattice(dof / 3, ci, cj, ck, cl);
  this->JacobianDOFs(jac, ci, cj, ck, cl, x, y, z, t, t0);
  const int c = dof % 3;
  for (int i = 0; i < 3; ++i) {
    if (i != c) jac[i] = .0;
  }
}

// =============================================================================
// Deprecated
// =============================================================================

// -----------------------------------------------------------------------------
inline void FreeFormTransformation4D::PutStatusCP(int i, int j, int k, int l,
                                                      DOFStatus sx, DOFStatus sy, DOFStatus sz)
{
  PutStatus(i, j, k, l, sx, sy, sz);
}

// -----------------------------------------------------------------------------
inline void FreeFormTransformation4D::GetStatusCP(int i, int j, int k, int l,
                                                      DOFStatus &sx, DOFStatus &sy, DOFStatus &sz) const
{
  GetStatus(i, j, k, l, sx, sy, sz);
}

// -----------------------------------------------------------------------------
inline void FreeFormTransformation4D::BoundingBoxCP(int cp, Point &p1,
                                                                Point &p2, double fraction) const
{
  int i, j, k;
  IndexToLattice(cp, i, j, k);
  p1 = Point(i - 2 * fraction, j - 2 * fraction, k - 2 * fraction);
  p2 = Point(i + 2 * fraction, j + 2 * fraction, k + 2 * fraction);
  this->LatticeToWorld(p1);
  this->LatticeToWorld(p2);
}

// -----------------------------------------------------------------------------
inline void
FreeFormTransformation4D
::BoundingBoxImage(const GreyImage *image, int cp, int &i1, int &j1, int &k1,
                                                   int &i2, int &j2, int &k2, double fraction) const
{
  // Calculate bounding box in world coordinates
  Point p1, p2;
  this->BoundingBoxCP(cp, p1, p2, fraction);

  // Transform world coordinates to image coordinates
  image->WorldToImage(p1._x, p1._y, p1._z);
  image->WorldToImage(p2._x, p2._y, p2._z);

  // Calculate bounding box in image coordinates
  i1 = (p1._x < 0) ? 0 : int(p1._x)+1;
  j1 = (p1._y < 0) ? 0 : int(p1._y)+1;
  k1 = (p1._z < 0) ? 0 : int(p1._z)+1;
  i2 = (int(p2._x) >= image->GetX()) ? image->GetX()-1 : int(p2._x);
  j2 = (int(p2._y) >= image->GetY()) ? image->GetY()-1 : int(p2._y);
  k2 = (int(p2._z) >= image->GetZ()) ? image->GetZ()-1 : int(p2._z);
}


} // namespace mirtk

#endif // MIRTK_FreeFormTransformation4D_H

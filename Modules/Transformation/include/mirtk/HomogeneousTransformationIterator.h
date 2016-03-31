/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_HomogeneousTransformationIterator_H
#define MIRTK_HomogeneousTransformationIterator_H

#include "mirtk/Point.h"
#include "mirtk/BaseImage.h"
#include "mirtk/HomogeneousTransformation.h"

#include <cstdlib>
#include <iostream>


namespace mirtk {


/**
 * Class for iterator for homogeneous matrix transformations.
 *
 * This class implements a fast access iterator 3D for homogeneous
 * matrix transformations.
 *
 * NOTE: This class has NO copy constructor
 */
class HomogeneousTransformationIterator : public Point
{
  mirtkObjectMacro(HomogeneousTransformationIterator);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Pointer to transformation
  mirtkPublicAggregateMacro(const HomogeneousTransformation, Transformation);

  /// Current x,y,z position in x-direction
  double _xx, _xy, _xz;

  /// Current x,y,z position in y-direction
  double _yx, _yy, _yz;

  /// Current x,y,z position in z-direction
  double _zx, _zy, _zz;

  /// Current x,y,z offsets in x-direction
  double _xdx, _xdy, _xdz;

  /// Current x,y,z offsets in y-direction
  double _ydx, _ydy, _ydz;

  /// Current x,y,z offsets in z-direction
  double _zdx, _zdy, _zdz;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  HomogeneousTransformationIterator(const HomogeneousTransformation * = NULL);

  // ---------------------------------------------------------------------------
  // Iteration

  /// Initialize iterator. This function initializes the transformation
  /// iterator at point x, y, z with voxel offsets 1, 1, 1
  void Initialize(const BaseImage *target, const BaseImage *source,
                  double x = .0, double y = .0, double z = .0,
                  bool inv = false);

  /// Advance iterator in x-direction
  void NextX();

  /// Advance iterator in x-direction by a certain amount
  void NextX(double);

  /// Advance iterator in y-direction
  void NextY();

  /// Advance iterator in y-direction by a certain amount
  void NextY(double);

  /// Advance iterator in z-direction
  void NextZ();

  /// Advance iterator in z-direction by a certain amount
  void NextZ(double);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline HomogeneousTransformationIterator
::HomogeneousTransformationIterator(const HomogeneousTransformation *transformation)
:
  _Transformation(transformation)
{
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformationIterator
::Initialize(const BaseImage *target, const BaseImage *source, double x, double y, double z, bool inv)
{
  if (_Transformation == NULL) {
    cout << "HomogeneousTransformationIterator::Initialize(): Transformation has not been set" << endl;
    exit(1);
  }

  // Image to image transformation matrix
  Matrix matrix = _Transformation->GetMatrix();
  if (inv) matrix.Invert();
  matrix = source->GetWorldToImageMatrix() * matrix * target->GetImageToWorldMatrix();

  // Calculate starting point
  _x = matrix(0, 0) * x + matrix(0, 1) * y + matrix(0, 2) * z + matrix(0, 3);
  _y = matrix(1, 0) * x + matrix(1, 1) * y + matrix(1, 2) * z + matrix(1, 3);
  _z = matrix(2, 0) * x + matrix(2, 1) * y + matrix(2, 2) * z + matrix(2, 3);

  _xx = _yx = _zx = _x;
  _xy = _yy = _zy = _y;
  _xz = _yz = _zz = _z;

  // Calculate offsets
  _zdx = matrix(0, 2);
  _zdy = matrix(1, 2);
  _zdz = matrix(2, 2);
  _ydx = matrix(0, 1);
  _ydy = matrix(1, 1);
  _ydz = matrix(2, 1);
  _xdx = matrix(0, 0);
  _xdy = matrix(1, 0);
  _xdz = matrix(2, 0);
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformationIterator::NextX()
{
  _x = (_xx += _xdx);
  _y = (_xy += _xdy);
  _z = (_xz += _xdz);
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformationIterator::NextX(double offset)
{
  _x = (_xx += _xdx * offset);
  _y = (_xy += _xdy * offset);
  _z = (_xz += _xdz * offset);
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformationIterator::NextY()
{
  _x = _xx = (_yx += _ydx);
  _y = _xy = (_yy += _ydy);
  _z = _xz = (_yz += _ydz);
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformationIterator::NextY(double offset)
{
  _x = _xx = (_yx + _ydx * offset);
  _y = _xy = (_yy + _ydy * offset);
  _z = _xz = (_yz + _ydz * offset);
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformationIterator::NextZ()
{
  _x = _xx = _yx = (_zx += _zdx);
  _y = _xy = _yy = (_zy += _zdy);
  _z = _xz = _yz = (_zz += _zdz);
}

// -----------------------------------------------------------------------------
inline void HomogeneousTransformationIterator::NextZ(double offset)
{
  _x = _xx = _yx = (_zx += _zdx * offset);
  _y = _xy = _yy = (_zy += _zdy * offset);
  _z = _xz = _yz = (_zz += _zdz * offset);
}


} // namespace mirtk

#endif // MIRTK_HomogeneousTransformationIterator_H

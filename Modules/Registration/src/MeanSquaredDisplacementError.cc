/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2017 Imperial College London
 * Copyright 2017 Andreas Schuh
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

#include "mirtk/MeanSquaredDisplacementError.h"

#include "mirtk/Math.h"
#include "mirtk/ObjectFactory.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(MeanSquaredDisplacementError);


// -----------------------------------------------------------------------------
MeanSquaredDisplacementError::MeanSquaredDisplacementError(const char *name, double w)
:
  DataFidelity(name, w),
  _TargetTransformation(nullptr),
  _ExternalDisplacement(nullptr)
{
}

// -----------------------------------------------------------------------------
void MeanSquaredDisplacementError::Initialize()
{
  // Initialize base class
  DataFidelity::Initialize();
  if (_TargetTransformation == nullptr) return;

  // Check required parameters
  if (!_Domain) {
    if (_ExternalDisplacement) {
      _Domain = _ExternalDisplacement->Attributes();
    } else {
      Throw(ERR_LogicError, __FUNCTION__, "Finite regular grid domain not set!");
    }
  }

  // Evaluate target displacement field
  _TargetDisplacement.Initialize(_Domain, 3);
  _TargetTransformation->Displacement(_TargetDisplacement);
  if (_ExternalDisplacement) {
    if (!_ExternalDisplacement->Attributes().EqualInSpace(_Domain)) {
      Throw(ERR_InvalidArgument, __FUNCTION__, "Externally managed displacement field must have same spatial attributes as domain on which error is evaluated!");
    }
    if (_ExternalDisplacement->T() != 3) {
      Throw(ERR_InvalidArgument, __FUNCTION__, "Externally managed displacement field must have 3 components/channels!");
    }
  }

  // Allocate memory for cached displacement fields
  _NonParametricGradient.Initialize(_Domain, 3);
  if (Transformation()->RequiresCachingOfDisplacements() && !_ExternalDisplacement) {
    _CurrentDisplacement.Initialize(_Domain, 3);
  }
}

// -----------------------------------------------------------------------------
void MeanSquaredDisplacementError::Update(bool gradient)
{
  // Update base class
  DataFidelity::Update(gradient);
  if (_TargetTransformation == nullptr) return;

  // Update transformation displacements
  if (Transformation()->RequiresCachingOfDisplacements() && !_ExternalDisplacement) {
    Transformation()->Displacement(_CurrentDisplacement);
  }
}

// -----------------------------------------------------------------------------
double MeanSquaredDisplacementError::Evaluate()
{
  if (_TargetTransformation == nullptr) {
    return 0.;
  }

  double error = 0.;
  if (_ExternalDisplacement || !_CurrentDisplacement.IsEmpty()) {

    const DisplacementImageType *currentdisplacement = _ExternalDisplacement;
    if (currentdisplacement == nullptr) {
      currentdisplacement = &_CurrentDisplacement;
    }

    if (currentdisplacement->X() != _TargetDisplacement.X() ||
        currentdisplacement->Y() != _TargetDisplacement.Y() ||
        currentdisplacement->Z() != _TargetDisplacement.Z()) {
      Throw(ERR_LogicError, __FUNCTION__, "Displacement fields have differing size!");
    }

    Point q1, q2;
    for (int k = 0; k < _TargetDisplacement.Z(); ++k)
    for (int j = 0; j < _TargetDisplacement.Y(); ++j)
    for (int i = 0; i < _TargetDisplacement.X(); ++i) {
      q1._x = currentdisplacement->Get(i, j, k, 0);
      q1._y = currentdisplacement->Get(i, j, k, 1);
      q1._z = currentdisplacement->Get(i, j, k, 2);
      q2._x = _TargetDisplacement(i, j, k, 0);
      q2._y = _TargetDisplacement(i, j, k, 1);
      q2._z = _TargetDisplacement(i, j, k, 2);
      error += q1.Distance(q2);
    }

  } else {

    Point p, q1, q2;
    for (int k = 0; k < _TargetDisplacement.Z(); ++k)
    for (int j = 0; j < _TargetDisplacement.Y(); ++j)
    for (int i = 0; i < _TargetDisplacement.X(); ++i) {
      p = Point(i, j, k);
      _TargetDisplacement.ImageToWorld(p);
      q1 = p, Transformation()->Transform(q1);
      q2._x = p._x + _TargetDisplacement(i, j, k, 0);
      q2._y = p._y + _TargetDisplacement(i, j, k, 1);
      q2._z = p._z + _TargetDisplacement(i, j, k, 2);
      error += q1.Distance(q2);
    }

  }
  error /= _Domain.NumberOfSpatialPoints();
  return error;
}

// -----------------------------------------------------------------------------
void MeanSquaredDisplacementError::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_TargetTransformation == nullptr) return;

  if (_ExternalDisplacement || !_CurrentDisplacement.IsEmpty()) {

    const DisplacementImageType *currentdisplacement = _ExternalDisplacement;
    if (currentdisplacement == nullptr) {
      currentdisplacement = &_CurrentDisplacement;
    }

    if (currentdisplacement->X() != _TargetDisplacement.X() ||
        currentdisplacement->Y() != _TargetDisplacement.Y() ||
        currentdisplacement->Z() != _TargetDisplacement.Z()) {
      Throw(ERR_LogicError, __FUNCTION__, "Displacement fields have differing size!");
    }

    Point q1, q2;
    for (int k = 0; k < _TargetDisplacement.Z(); ++k)
    for (int j = 0; j < _TargetDisplacement.Y(); ++j)
    for (int i = 0; i < _TargetDisplacement.X(); ++i) {
      q1._x = currentdisplacement->Get(i, j, k, 0);
      q1._y = currentdisplacement->Get(i, j, k, 1);
      q1._z = currentdisplacement->Get(i, j, k, 2);
      q2._x = _TargetDisplacement(i, j, k, 0);
      q2._y = _TargetDisplacement(i, j, k, 1);
      q2._z = _TargetDisplacement(i, j, k, 2);
      _NonParametricGradient(i, j, k, 0) = 2. * static_cast<GradientType>(q1._x - q2._x);
      _NonParametricGradient(i, j, k, 1) = 2. * static_cast<GradientType>(q1._y - q2._y);
      _NonParametricGradient(i, j, k, 2) = 2. * static_cast<GradientType>(q1._z - q2._z);
    }

  } else {

    Point p, q1, q2;
    for (int k = 0; k < _TargetDisplacement.Z(); ++k)
    for (int j = 0; j < _TargetDisplacement.Y(); ++j)
    for (int i = 0; i < _TargetDisplacement.X(); ++i) {
      p = Point(i, j, k);
      _TargetDisplacement.ImageToWorld(p);
      q1 = p, Transformation()->Transform(q1);
      q2._x = p._x + _TargetDisplacement(i, j, k, 0);
      q2._y = p._y + _TargetDisplacement(i, j, k, 1);
      q2._z = p._z + _TargetDisplacement(i, j, k, 2);
      _NonParametricGradient(i, j, k, 0) = 2. * static_cast<GradientType>(q1._x - q2._x);
      _NonParametricGradient(i, j, k, 1) = 2. * static_cast<GradientType>(q1._y - q2._y);
      _NonParametricGradient(i, j, k, 2) = 2. * static_cast<GradientType>(q1._z - q2._z);
    }

  }
  _NonParametricGradient /= _Domain.NumberOfSpatialPoints();

  const double t0 = _Domain._torigin;
  _NonParametricGradient.PutTOrigin(t0); // TODO: Set actual source time?
  Transformation()->ParametricGradient(&_NonParametricGradient, gradient, nullptr, t0, weight);
}


} // namespace mirtk

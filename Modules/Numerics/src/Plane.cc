/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#include "mirtk/Plane.h"

#include "mirtk/Vector.h"
#include "mirtk/Matrix.h"


namespace mirtk {


// -----------------------------------------------------------------------------
void Plane::UpdateOrigin()
{
  _Origin._x = _Offset * _Normal[0];
  _Origin._y = _Offset * _Normal[1];
  _Origin._z = _Offset * _Normal[2];
}

// -----------------------------------------------------------------------------
void Plane::UpdateTangents()
{
  Vector3 v(_Normal[1], _Normal[2], _Normal[0]);
  _Tangent1 = _Normal.Cross(v);
  if (_Tangent1.SquaredLength() < 1e-6) {
    v[1] *= -1.0;
    _Tangent1 = _Normal.Cross(v);
    if (_Tangent1.SquaredLength() < 1e-6) {
      cerr << "Plane::UpdateTangents: Degenerate normal vector!" << endl;
      exit(1);
    }
  }
  _Tangent2 = _Normal.Cross(_Tangent1);
  _Tangent1.Normalize();
  _Tangent2.Normalize();
}

// -----------------------------------------------------------------------------
double Plane::Fit(const PointSet &points)
{
  if (points.Size() == 1) {
    _Normal = Vector3(.0, .0, 1.0);
    _Offset = _Normal.Dot(Vector3(points(0)._x, points(0)._y, points(0)._z));
    UpdateTangents();
    return .0;
  }

  // Build design matrix
  Matrix x(points.Size(), 3);
  for (int i = 0; i < x.Rows(); ++i) {
    const Point &p = points(i);
    x(i, 0) = p._x;
    x(i, 1) = p._y;
    x(i, 2) = p._z;
  }

  // Normalize design matrix
  _Origin._x = x.ColMean(0);
  _Origin._y = x.ColMean(1);
  _Origin._z = x.ColMean(2);
  x.AddToCol(0, -_Origin._x);
  x.AddToCol(1, -_Origin._y);
  x.AddToCol(2, -_Origin._z);
  x /= sqrt(x.Rows() - 1);

  // Perform PCA using SVD
  Matrix u, v;
  Vector s;

  x.SVD(u, s, v);

  // Get plane normal and orthonormal tangent vectors
  int c = 2;
  if (s(0) < s(c)) c = 0;
  if (s(1) < s(c)) c = 1;

  _Normal[0] = v(0, c);
  _Normal[1] = v(1, c);
  _Normal[2] = v(2, c);

  c = (c == 2 ? 0 : c + 1);
  _Tangent1[0] = v(0, c);
  _Tangent1[1] = v(1, c);
  _Tangent1[2] = v(2, c);

  c = (c == 2 ? 0 : c + 1);
  _Tangent2[0] = v(0, c);
  _Tangent2[1] = v(1, c);
  _Tangent2[2] = v(2, c);

  Vector3 cross = _Tangent2.Cross(_Normal);
  if (_Tangent1.Dot(cross) < .0) _Normal *= -1.0;

  // Compute distance of plane to world origin
  _Offset = (_Origin._x * _Normal[0] + _Origin._y * _Normal[1] + _Origin._z * _Normal[2]);

  // Evaluate RMS error
  double d, error = .0;
  for (int i = 0; i < points.Size(); ++i) {
    d = Distance(points(i));
    error += d * d;
  }
  return sqrt(error / points.Size());
}

// -----------------------------------------------------------------------------
ostream &Plane::Print(ostream &os, Indent indent) const
{
  os << indent << "<n, x> = b, where n = ["
     << _Normal[0] << ", " << _Normal[1] << ", " << _Normal[2]
     << "] and b = " << _Offset << " with local origin at x0 = ["
     << _Origin._x << ", " << _Origin._y << ", " << _Origin._z << "]\n";
  return os;
}


} // namespace mirtk

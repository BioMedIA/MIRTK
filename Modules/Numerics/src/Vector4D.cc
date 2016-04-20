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

#include "mirtk/Vector4D.h"

namespace mirtk {


// -----------------------------------------------------------------------------
template <typename T> Vector4D<T> Vector4D<T>::operator /(const Vector4D<T> &v) const
{
  const T zero(0);
  Vector4D<T> val(zero);
  if (v._x != zero) val._x = _x / v._x;
  if (v._y != zero) val._y = _y / v._y;
  if (v._z != zero) val._z = _z / v._z;
  if (v._t != zero) val._t = _t / v._t;
  return val;
}

// -----------------------------------------------------------------------------
template <typename T> Vector4D<T> &Vector4D<T>::operator /=(const Vector4D<T> &v)
{
  const T zero(0);
  if (v._x != zero) _x /= v._x;
  if (v._y != zero) _y /= v._y;
  if (v._z != zero) _z /= v._z;
  if (v._t != zero) _t /= v._t;
  return *this;
}

// -----------------------------------------------------------------------------
template struct Vector4D<char>;
template struct Vector4D<short>;
template struct Vector4D<float>;
template struct Vector4D<double>;


} // namespace mirtk

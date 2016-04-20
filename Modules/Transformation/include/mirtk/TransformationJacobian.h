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

#ifndef MIRTK_TransformationJacobian_H
#define MIRTK_TransformationJacobian_H

#include "mirtk/Object.h"

#include "mirtk/Math.h"
#include "mirtk/Matrix.h"
#include "mirtk/Vector3D.h"
#include "mirtk/OrderedMap.h"


namespace mirtk {


/**
 * Sparse matrix for the transformation Jacobian of derivatives w.r.t the parameters.
 *
 * This matrix type only stores the non-zero columns of the Jacobian matrix of
 * the transformation, which contains the derivatives of the transformation
 * w.r.t the transformation parameters. The full Jacobian matrix has dimension
 * 3xN, where N is the number of transformation parameters and the number of
 * rows correspond to the deformation in each spatial dimension (T_x, T_y, T_z).
 */
class TransformationJacobian : public Object
{
  mirtkObjectMacro(TransformationJacobian);

public:
  typedef Vector3D<double>                   ColumnType;
  typedef OrderedMap<int, ColumnType>        SparseMatrixType;
  typedef SparseMatrixType::iterator         ColumnIterator;
  typedef SparseMatrixType::const_iterator   ConstColumnIterator;
  typedef Matrix                             DenseMatrixType;

protected:

  /// Non-zero columns of transformation Jacobian
  SparseMatrixType _Columns;

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Constructor
  TransformationJacobian();

  /// Destructor  
  ~TransformationJacobian();

  /// Remove all non-zero columns
  void Clear();

  // ---------------------------------------------------------------------------
  // Element access

  /// Get number of non-zero columns
  int NumberOfNonZeroColumns() const;

  /// Get iterator to first non-zero column
  ColumnIterator Begin();

  /// Get iterator to first non-zero column
  ConstColumnIterator Begin() const;

  /// Get iterator to position after last non-zero column
  ColumnIterator End();

  /// Get iterator to position after last non-zero column
  ConstColumnIterator End() const;

  /// Get iterator to i-th non-zero column
  ColumnIterator GetNonZeroColumn(int);

  /// Get iterator to i-th non-zero column
  ConstColumnIterator GetNonZeroColumn(int) const;

  /// Get i-th non-zero column vector
  ColumnType &operator [](int);

  /// Get i-th non-zero column vector
  const ColumnType &operator [](int) const;

  /// Get i-th column vector
  ColumnType &operator ()(int);

  /// Get column index of the i-th non-zero column.
  int ColumnIndex(int) const;

  /// Get i-th non-zero column vector.
  ColumnType &ColumnVector(int);

  /// Get i-th non-zero column vector.
  const ColumnType &ColumnVector(int) const;

  /// Get i-th column vector, inserts new zero column if necessary.
  ColumnType &Column(int);

  /// Get iterator to i-th column.
  ColumnIterator Find(int);

  /// Get iterator to i-th column.
  ConstColumnIterator Find(int) const;

  // ---------------------------------------------------------------------------
  // Operators
  
  /// Add transformation Jacobian to this Jacobian matrix
  TransformationJacobian &operator +=(const TransformationJacobian &);

  /// Multiply this transformation Jacobian by a scalar
  TransformationJacobian &operator *=(const double);

  /// Pre-multiply (!) this transformation Jacobian with the given 3x3 matrix
  TransformationJacobian &operator *=(const Matrix &);

  // ---------------------------------------------------------------------------
  // Mathematical operations

  /// Add scaled transformation Jacobian to this Jacobian matrix
  TransformationJacobian &add(const TransformationJacobian &, double);

};


////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
inline TransformationJacobian::TransformationJacobian()
{
}

// -----------------------------------------------------------------------------
inline TransformationJacobian::~TransformationJacobian()
{
}

// =============================================================================
// Element access
// =============================================================================

// -----------------------------------------------------------------------------
inline void TransformationJacobian::Clear()
{
  return _Columns.clear();
}

// -----------------------------------------------------------------------------
inline int TransformationJacobian::NumberOfNonZeroColumns() const
{
  return static_cast<int>(_Columns.size());
}

// -----------------------------------------------------------------------------
inline TransformationJacobian::ColumnIterator TransformationJacobian::Begin()
{
  return _Columns.begin();
}

// -----------------------------------------------------------------------------
inline TransformationJacobian::ConstColumnIterator TransformationJacobian::Begin() const
{
  return _Columns.begin();
}

// -----------------------------------------------------------------------------
inline TransformationJacobian::ColumnIterator TransformationJacobian::End()
{
  return _Columns.end();
}

// -----------------------------------------------------------------------------
inline TransformationJacobian::ConstColumnIterator TransformationJacobian::End() const
{
  return _Columns.end();
}

// -----------------------------------------------------------------------------
inline TransformationJacobian::ColumnIterator TransformationJacobian::GetNonZeroColumn(int i)
{
  ColumnIterator it = _Columns.begin();
  for (int j = 0; j < i; j++) {
    ++it;
    if (it == _Columns.end()) {
      cerr << "TransformationJacobian::GetNonZeroColumn: Index is out of bounds: " << i << endl;
      exit(1);
    }
  }
  return it;
}

// -----------------------------------------------------------------------------
inline TransformationJacobian::ConstColumnIterator TransformationJacobian::GetNonZeroColumn(int i) const
{
  ConstColumnIterator it = _Columns.begin();
  for (int j = 0; j < i; j++) {
    ++it;
    if (it == _Columns.end()) {
      cerr << "TransformationJacobian::GetNonZeroColumn: Index is out of bounds: " << i << endl;
      exit(1);
    }
  }
  return it;
}

// -----------------------------------------------------------------------------
inline int TransformationJacobian::ColumnIndex(int i) const
{
  return GetNonZeroColumn(i)->first;
}

// -----------------------------------------------------------------------------
inline TransformationJacobian::ColumnType &TransformationJacobian::ColumnVector(int i)
{
  return GetNonZeroColumn(i)->second;
}

// -----------------------------------------------------------------------------
inline const TransformationJacobian::ColumnType &TransformationJacobian::ColumnVector(int i) const
{
  return GetNonZeroColumn(i)->second;
}

// -----------------------------------------------------------------------------
inline TransformationJacobian::ColumnType &TransformationJacobian::Column(int c)
{
  return _Columns[c];
}

// -----------------------------------------------------------------------------
inline TransformationJacobian::ColumnType &TransformationJacobian::operator [](int i)
{
  return GetNonZeroColumn(i)->second;
}

// -----------------------------------------------------------------------------
inline const TransformationJacobian::ColumnType &TransformationJacobian::operator [](int i) const
{
  return GetNonZeroColumn(i)->second;
}

// -----------------------------------------------------------------------------
inline TransformationJacobian::ColumnType &TransformationJacobian::operator ()(int i)
{
  return Column(i);
}

// -----------------------------------------------------------------------------
inline TransformationJacobian::ColumnIterator TransformationJacobian::Find(int c)
{
  return _Columns.find(c);
}

// -----------------------------------------------------------------------------
inline TransformationJacobian::ConstColumnIterator TransformationJacobian::Find(int c) const
{
  return _Columns.find(c);
}

// =============================================================================
// Unary operators
// =============================================================================

// -----------------------------------------------------------------------------
inline TransformationJacobian &TransformationJacobian::operator +=(const TransformationJacobian &b)
{
  for (ConstColumnIterator it = b.Begin(); it != b.End(); ++it) {
    _Columns[it->first] += it->second;
  }
  return *this;
}

// -----------------------------------------------------------------------------
inline TransformationJacobian &TransformationJacobian::operator *=(const double s)
{
  for (ColumnIterator it = Begin(); it != End(); ++it) it->second *= s;
  return *this;
}

// -----------------------------------------------------------------------------
inline TransformationJacobian &TransformationJacobian::operator *=(const Matrix &a)
{
  // Note: Read this operator as: a * (*this)!
  ColumnType v;
  for (ColumnIterator it = Begin(); it != End(); ++it) {
    v._x = a(0, 0) * it->second._x + a(0, 1) * it->second._y + a(0, 2) * it->second._z;
    v._y = a(1, 0) * it->second._x + a(1, 1) * it->second._y + a(1, 2) * it->second._z;
    v._z = a(2, 0) * it->second._x + a(2, 1) * it->second._y + a(2, 2) * it->second._z;
    it->second = v;
  }
  return *this;
}

// =============================================================================
// Binary operators
// =============================================================================

// -----------------------------------------------------------------------------
/// Calculate column-by-column sum of transformation Jacobian
inline TransformationJacobian operator +(TransformationJacobian &a, TransformationJacobian &b)
{
  TransformationJacobian c = a;
  c += b;
  return c;
}

// -----------------------------------------------------------------------------
/// Multiply transformation Jacobian and scalar
inline TransformationJacobian operator *(TransformationJacobian &a, double s)
{
  TransformationJacobian b = a;
  b *= s;
  return b;
}

// -----------------------------------------------------------------------------
/// Multiply transformation Jacobian and scalar
inline TransformationJacobian operator *(double s, TransformationJacobian &a)
{
  return a * s;
}

// -----------------------------------------------------------------------------
/// Calculate product of 3x3 matrix and transformation Jacobian
inline TransformationJacobian operator *(Matrix &a, TransformationJacobian &b)
{
  TransformationJacobian c = b;
  c *= a;
  return c;
}

// =============================================================================
// Mathematical operations
// =============================================================================

// -----------------------------------------------------------------------------
inline TransformationJacobian &TransformationJacobian::add(const TransformationJacobian &b, double s)
{
  for (ConstColumnIterator it = b.Begin(); it != b.End(); ++it) {
    _Columns[it->first] += it->second * s;
  }
  return *this;
}

// =============================================================================
// Debugging
// =============================================================================

inline bool has_nan(const TransformationJacobian &a)
{
  for (TransformationJacobian::ConstColumnIterator it = a.Begin(); it != a.End(); ++it) {
    if (IsNaN(it->second._x) || IsNaN(it->second._y) || IsNaN(it->second._z)) {
      cerr << "TransformationJacobian::has_nan: Found NaN in column " << it->first << endl;
      return true;
    }
  }
  return false;
}


} // namespace mirtk

#endif // MIRTK_TransformationJacobian_H

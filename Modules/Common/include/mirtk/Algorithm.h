/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#ifndef MIRTK_Algorithm_H
#define MIRTK_Algorithm_H

#include "mirtk/Array.h"

#include <algorithm>
#include <iterator>


namespace mirtk {


using std::find;
using std::binary_search;
using std::sort;
using std::stable_sort;
using std::partial_sort;
using std::nth_element;
using std::transform;
using std::reverse;
using std::shuffle;

using std::copy;
using std::back_inserter;
using std::front_inserter;
using std::inserter;


/// Compare functor used to sort array indices increasing element attribute
template <class AttributeType>
class CompareIndicesOfArrayByIncreasingValue
{
  const Array<AttributeType> &_Values;

public:

  CompareIndicesOfArrayByIncreasingValue(const Array<AttributeType> &values)
  :
    _Values(values)
  {}

  template <class IndexType>
  bool operator ()(IndexType i, IndexType j) const
  {
    return _Values[i] < _Values[j];
  }
};

/// Compare functor used to sort array indices by decreasing element attribute
template <class AttributeType>
class CompareIndicesOfArrayByDecreasingValue
{
  const Array<AttributeType> &_Values;

public:

  CompareIndicesOfArrayByDecreasingValue(const Array<AttributeType> &values)
  :
    _Values(values)
  {}

  template <class IndexType>
  bool operator ()(IndexType i, IndexType j) const
  {
    return _Values[i] > _Values[j];
  }
};

/// Sort values in array
template <class T>
void Sort(Array<T> &values)
{
  sort(values.begin(), values.end());
}

/// Sort values in array using custom comparator
template <class T, class Compare>
void Sort(Array<T> &values, Compare comp)
{
  sort(values.begin(), values.end(), comp);
}

/// Sort values in array while preserving order of equal entries
template <class T>
void StableSort(Array<T> &values)
{
  stable_sort(values.begin(), values.end());
}

/// Sort values in array while preserving order of equal entries
template <class T, class Compare>
void StableSort(Array<T> &values, Compare comp)
{
  stable_sort(values.begin(), values.end(), comp);
}

/// Sort first n values in array
template <class T>
void PartialSort(Array<T> &values, int n)
{
  partial_sort(values.begin(), values.begin() + n + 1, values.end());
}

/// Sort values in array using custom comparator
template <class T, class Compare>
void PartialSort(Array<T> &values, int n, Compare comp)
{
  partial_sort(values.begin(), values.begin() + n + 1, values.end(), comp);
}

/// Sort values of array such that values before the n-th element are
/// smaller than this element and elements after are greater or equal
template <class T>
T &NthElement(Array<T> &values, int n)
{
  nth_element(values.begin(), values.begin() + n, values.end());
  return values[n];
}

/// Sort values of array such that values before the n-th element are
/// smaller than this element and elements after are greater or equal
/// according to a custom comparator functor
template <class T, class Compare>
T &NthElement(Array<T> &values, int n, Compare comp)
{
  nth_element(values.begin(), values.begin() + n, values.end(), comp);
  return values[n];
}

/// Get permutation of array indices corresponding to sorted order of values
template <class T>
Array<int> IncreasingOrder(const Array<T> &values)
{
  Array<int> indices(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    indices[i] = static_cast<int>(i);
  }
  Sort(indices, CompareIndicesOfArrayByIncreasingValue<T>(values));
  return indices;
}

/// Get permutation of array indices corresponding to sorted order of values
template <class T>
Array<int> DecreasingOrder(const Array<T> &values)
{
  Array<int> indices(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    indices[i] = static_cast<int>(i);
  }
  Sort(indices, CompareIndicesOfArrayByDecreasingValue<T>(values));
  return indices;
}

/// Get permutation of array values given a permutation of array indices
///
/// \sa IncreasingOrder, DecreasingOrder
template <class T>
Array<T> Permutation(const Array<int> &order, const Array<T> &values)
{
  Array<T> permuted_values;
  permuted_values.reserve(values.size());
  for (auto i : order) permuted_values.push_back(values[i]);
  return permuted_values;
}

/// Find value in unsorted Array
template <class T>
typename Array<T>::iterator
Find(Array<T> &values, const T &value)
{
  return find(values.begin(), values.end(), value);
}

/// Find value in unsorted Array
template <class T>
typename Array<T>::const_iterator
Find(const Array<T> &values, const T &value)
{
  return find(values.begin(), values.end(), value);
}

/// Find value in unsorted Array
template <class T>
int FindIndex(const Array<T> &values, const T &value)
{
  auto it = find(values.begin(), values.end(), value);
  if (it == values.end()) return -1;
  return static_cast<int>(std::distance(values.begin(), it));
}

/// Find value in sorted Array using binary search
template <class T>
typename Array<T>::iterator
BinarySearch(Array<T> &values, const T &value)
{
  return binary_search(values.begin(), values.end(), value);
}

/// Find value in sorted Array using binary search
template <class T>
typename Array<T>::const_iterator
BinarySearch(const Array<T> &values, const T &value)
{
  return binary_search(values.begin(), values.end(), value);
}

/// Find value in sorted Array using binary search with custom comparator
template <class T, class Compare>
typename Array<T>::iterator
BinarySearch(Array<T> &values, const T &value, Compare comp)
{
  return binary_search(values.begin(), values.end(), value, comp);
}

/// Find value in sorted Array using binary search with custom comparator
template <class T, class Compare>
typename Array<T>::const_iterator
BinarySearch(const Array<T> &values, const T &value, Compare comp)
{
  return binary_search(values.begin(), values.end(), value, comp);
}


} // namespace mirtk

#endif // MIRTK_Algorithm_H

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

#ifndef MIRTK_Algorithm_H
#define MIRTK_Algorithm_H

#include "mirtk/Array.h"

#include <algorithm>


namespace mirtk {


using std::sort;
using std::partial_sort;
using std::transform;
using std::reverse;
using std::shuffle;


/**
 * Compare functor used to sort array indices based on some element attribute
 *
 * http://stackoverflow.com/questions/3909272/sorting-two-corresponding-arrays
 */
template <class AttributeType>
class SortIndicesOfArray
{
  const Array<AttributeType> &_Values;

public:
  SortIndicesOfArray(const Array<AttributeType> &values) : _Values(values) {}

  template <class IndexType>
  bool operator ()(IndexType i, IndexType j) const {
    return _Values[i] < _Values[j];
  }
};


} // namespace mirtk

#endif // MIRTK_Algorithm_H

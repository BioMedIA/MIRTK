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

/**
 * \file  mirtk/Utils.h
 * \brief Collection of a few very basic functions.
 */

#ifndef MIRTK_Utils_H
#define MIRTK_Utils_H

#include "mirtk/OrderedSet.h"


namespace mirtk {


// -----------------------------------------------------------------------------
/// Get average interval between consecutive values of an ordered set
inline double AverageInterval(const OrderedSet<double> &values)
{
  double avg = .0;
  if (values.size() > 1) {
    OrderedSet<double>::const_iterator j = values.begin();
    OrderedSet<double>::const_iterator i = j++;
    for (; j != values.end(); ++i, ++j) avg += (*j) - (*i);
    avg /= (values.size() - 1);
  }
  return avg;
}


} // namespace mirtk

#endif // MIRTK_Utils_H

/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2016 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
 * Copyright 2016      Andreas Schuh
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

#include <istream>
#include <ostream>

#include "mirtk/Point.h"


namespace mirtk {


// -----------------------------------------------------------------------------
ostream &operator <<(ostream &o, const Point &p)
{
  return o << p._x <<  " " << p._y << " " << p._z;
}

// -----------------------------------------------------------------------------
istream &operator >>(istream &i, Point &p)
{
  char c = i.peek();
  while ((c == ',') || (c == ' ') || (c == '\t') || (c == '\n')) {
    i.get();
    if (!i.eof()) {
      c = i.peek();
    } else {
      break;
    }
  }
  return i >> p._x >> p._y >> p._z;
}


} // namespace mirtk

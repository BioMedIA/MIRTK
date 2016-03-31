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

#ifndef MIRTK_Indent_H
#define MIRTK_Indent_H

#include "mirtk/Stream.h"


namespace mirtk {


/**
 * Auxiliary class for output indentation.
 */
class Indent
{
protected:

  /// Number of indentation levels
  int _N;

  /// Number of whitespace characters per indentation level
  int _NumberOfSpaces;

public:

  /// Constructor
  Indent(int n = 0, int s = 2) : _N(n), _NumberOfSpaces(s) {}

  /// Copy constructor
  Indent(const Indent &o) : _N(o._N), _NumberOfSpaces(o._NumberOfSpaces) {}

  /// Assignment operator
  Indent &operator= (int n)
  {
    _N = n;
    return *this;
  }

  /// Assignment operator
  Indent &operator= (const Indent &indent)
  {
    _N              = indent._N;
    _NumberOfSpaces = indent._NumberOfSpaces;
    return *this;
  }

  /// Pre-decrement operator
  Indent &operator-- ()
  {
    _N--;
    return *this;
  }

  /// Pre-increment operator
  Indent &operator++ ()
  {
    _N++;
    return *this;
  }

  /// Post-decrement operator
  Indent operator-- (int)
  {
    Indent pre(*this);
    --(*this);
    return pre;
  }

  /// Post-increment operator
  Indent operator++ (int)
  {
    Indent pre(*this);
    ++(*this);
    return pre;
  }

  /// Add indentation to this indentation
  Indent &operator+= (const Indent &indent)
  {
    _N += indent._N;
    return *this;
  }

  /// Add two indentations
  Indent operator+ (const Indent &indent) const
  {
    return Indent(_N + indent._N);
  }

  /// Subtract indentation from this indentation
  Indent &operator-= (const Indent &indent)
  {
    _N -= indent._N;
    return *this;
  }

  /// Subtract two indentations from another
  Indent operator- (const Indent &indent) const
  {
    return Indent(_N + indent._N);
  }

  /// Get indentation level
  int Level() const
  {
    return _N;
  }

  /// Set number of indentation whitespace characters per level
  void SpacesPerLevel(int n)
  {
    _NumberOfSpaces = n;
  }

  /// Get number of indentation whitespace characters per level
  int SpacesPerLevel() const
  {
    return _NumberOfSpaces;
  }

  /// Number of space characters
  int Spaces() const
  {
    return _N * _NumberOfSpaces;
  }
};

// ---------------------------------------------------------------------------
// Streaming operator
inline ostream &operator<< (ostream &os, const Indent &indent)
{
  for (int i = 0; i < indent.Spaces(); ++i) os << " ";
  return os;
}


} // namespace mirtk

#endif // MIRTK_Indent_H

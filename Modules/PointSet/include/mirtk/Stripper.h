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

#ifndef MIRTK_Stripper_H
#define MIRTK_Stripper_H

#include "mirtk/MeshFilter.h"


namespace mirtk {


/**
 * Creates line strips from contiguous lines and adjacent triangles
 *
 * This filter is the equivalent of the vtkStripper, but with a more efficient
 * line stripping algorithm. Moreover, VTK <7 does not have the
 * vtkStripper::SetJoinContigousLines option.
 */
class Stripper : public MeshFilter
{
  mirtkObjectMacro(Stripper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Convert contiguous lines to poly-lines
  mirtkPublicAttributeMacro(bool, StripLines);

  /// Convert adjacent triangles to triangle strips
  mirtkPublicAttributeMacro(bool, StripTriangles);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const Stripper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  Stripper();

  /// Copy constructor
  Stripper(const Stripper &);

  /// Assignment operator
  Stripper &operator =(const Stripper &);

  /// Destructor
  virtual ~Stripper();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Execute filter
  virtual void Execute();

  // ---------------------------------------------------------------------------
  // Alternative VTK-like API

public:

  /// Enable/disable joining of contiguous lines
  mirtkOnOffMacro(StripLines);

  /// Enable/disable generation of triangle strips
  mirtkOnOffMacro(StripTriangles);

};


} // namespace mirtk

#endif // MIRTK_Stripper_H

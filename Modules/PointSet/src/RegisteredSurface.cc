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

#include "mirtk/RegisteredSurface.h"
#include "mirtk/Transformation.h"

#include "vtkPolyData.h"


namespace mirtk {


// -----------------------------------------------------------------------------
RegisteredSurface::RegisteredSurface(vtkPolyData                *data,
                                     const class Transformation *transform)
:
  RegisteredPointSet(data, transform)
{
}

// -----------------------------------------------------------------------------
RegisteredSurface::RegisteredSurface(const RegisteredSurface &other)
:
  RegisteredPointSet(other)
{
}

// -----------------------------------------------------------------------------
RegisteredSurface &RegisteredSurface::operator =(const RegisteredSurface &other)
{
  RegisteredPointSet::operator =(other);
  return *this;
}

// -----------------------------------------------------------------------------
RegisteredSurface::~RegisteredSurface()
{
}

// -----------------------------------------------------------------------------
void RegisteredSurface::Initialize()
{
  // Ensure that input dataset is of valid type
  _InputSurface = vtkPolyData::SafeDownCast(_InputPointSet);
  if (_InputSurface == NULL) {
    cerr << "RegisteredSurface::Initialize: Input dataset must be a vtkPolyData" << endl;
    exit(1);
  }

  // Build cells and links
  _InputSurface->BuildLinks();

  // Initialize base class -- makes shallow copy of input surface, incl. links
  RegisteredPointSet::Initialize();
}


} // namespace mirtk

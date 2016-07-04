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

#include "mirtk/SurfaceFilter.h"

#include "mirtk/PointSetUtils.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void SurfaceFilter::CopyAttributes(const SurfaceFilter &other)
{
}

// -----------------------------------------------------------------------------
SurfaceFilter::SurfaceFilter()
{
}

// -----------------------------------------------------------------------------
SurfaceFilter::SurfaceFilter(const SurfaceFilter &other)
:
  MeshFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
SurfaceFilter &SurfaceFilter::operator =(const SurfaceFilter &other)
{
  if (this != &other) {
    MeshFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
SurfaceFilter::~SurfaceFilter()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void SurfaceFilter::Initialize()
{
  // Initialize base class
  MeshFilter::Initialize();

  // Check input mesh
  if (!IsSurfaceMesh(_Input)) {
    cerr << this->NameOfType() << "::Initialize: Input mesh must be a surface mesh with 2D faces" << endl;
    exit(1);
  }
}


} // namespace mirtk

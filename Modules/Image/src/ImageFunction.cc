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

#include "mirtk/ImageFunction.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
ImageFunction::ImageFunction()
:
  _Input(NULL),
  _DefaultValue(.0)
{
}

// -----------------------------------------------------------------------------
ImageFunction::ImageFunction(const ImageFunction &other)
:
  Object(other),
  _Input(other._Input),
  _DefaultValue(other._DefaultValue)
{
}

// -----------------------------------------------------------------------------
ImageFunction::~ImageFunction()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void ImageFunction::Initialize()
{
  if (_Input == NULL) {
    cerr << this->NameOfClass() << "::Initialize: Function has no input image" << endl;
    exit(1);
  }
}


} // namespace mirtk

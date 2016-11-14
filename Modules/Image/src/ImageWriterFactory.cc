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

#include "mirtk/ImageWriterFactory.h"

#include "mirtk/Assert.h"
#include "mirtk/Path.h"
#include "mirtk/Array.h"


namespace mirtk {


// -----------------------------------------------------------------------------
ImageWriterFactory &ImageWriterFactory::Instance()
{
  static ImageWriterFactory instance;
  return instance;
}

// -----------------------------------------------------------------------------
ImageWriterFactory::ImageWriterFactory()
{
}

// -----------------------------------------------------------------------------
ImageWriterFactory::~ImageWriterFactory()
{
}

// -----------------------------------------------------------------------------
bool ImageWriterFactory::Register(const Array<string> &exts, ImageWriterCreator creator)
{
  #ifndef NDEBUG
    UniquePtr<ImageWriter> writer(creator());
    mirtkAssert(writer != nullptr, "ImageWriterCreator produces object");
  #endif
  for (auto ext = exts.begin(); ext != exts.end(); ++ext) {
    _Associations[*ext] = creator;
  }
  return true;
}

// -----------------------------------------------------------------------------
ImageWriter *ImageWriterFactory::New(const char *fname) const
{
  auto it = _Associations.find(Extension(fname));
  return (it == _Associations.end() ? nullptr : (it->second)());
}


} // namespace mirtk

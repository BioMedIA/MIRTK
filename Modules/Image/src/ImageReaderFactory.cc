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

#include <mirtkImageReaderFactory.h>

#include <mirtkAssert.h>


namespace mirtk {


// -----------------------------------------------------------------------------
ImageReaderFactory &ImageReaderFactory::Instance()
{
  static ImageReaderFactory instance;
  return instance;
}

// -----------------------------------------------------------------------------
ImageReaderFactory::ImageReaderFactory()
{
}

// -----------------------------------------------------------------------------
ImageReaderFactory::~ImageReaderFactory()
{
}

// -----------------------------------------------------------------------------
bool ImageReaderFactory::Register(ImageReaderCreator creator)
{
  mirtkAssert(creator() != nullptr, "ImageReaderCreator produces object");
  _Creators.push_back(creator);
  return true;
}

// -----------------------------------------------------------------------------
ImageReader *ImageReaderFactory::New(const char *fname) const
{
  ImageReader *reader = nullptr;
  for (auto it = _Creators.begin(); it != _Creators.end(); ++it) {
    reader = (*it)();
    if (reader->CanRead(fname)) break;
    delete reader;
    reader = nullptr;
  }
  return reader;
}


} // namespace mirtk

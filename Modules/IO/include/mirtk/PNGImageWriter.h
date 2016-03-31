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

#ifndef MIRTK_PNGImageWriter_H
#define MIRTK_PNGImageWriter_H

#include "mirtk/ImageWriter.h"

#include "mirtk/Array.h"


namespace mirtk {


/**
 * Class for image to PNG file filter.
 *
 * This is a class which takes an image as input and produces an image file
 * in PNG file format. Note that PNG file formats support only 2D images!!!
 */
class PNGImageWriter : public ImageWriter
{
  mirtkObjectMacro(PNGImageWriter);

protected:

  /// Initialize filter
  virtual void Initialize();

public:

  /// List of file name extensions
  static Array<string> Extensions();

  /// Write entire image
  virtual void Run();

};


} // namespace mirtk

#endif // MIRTK_PNGImageWriter_H

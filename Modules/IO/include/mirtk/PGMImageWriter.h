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

#ifndef MIRTK_PGMImageWriter_H
#define MIRTK_PGMImageWriter_H

#include "mirtk/ImageWriter.h"

#include "mirtk/Array.h"


namespace mirtk {


/**
 * Class for image to PGM file filter.
 *
 * This is a class which takes an image as input and produces an image file
 * in PGM file format. Note that PGM file formats support only 2D images!!!
 */
class PGMImageWriter : public ImageWriter
{
  mirtkObjectMacro(PGMImageWriter);

public:

  /// List of file name extensions
  static Array<string> Extensions();

  /// Write entire image
  virtual void Run();

protected:

  /// Initialize filter
  virtual void Initialize();

};


} // namespace mirtk

#endif // MIRTK_PGMImageWriter_H

/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2019 Imperial College London
 * Copyright 2019 Andreas Schuh
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

#ifndef MIRTK_MetaImageWriter_H
#define MIRTK_MetaImageWriter_H

#include "mirtk/ImageWriter.h"

#include "mirtk/Memory.h"
#include "mirtk/Array.h"


namespace mirtk {


/**
 * Writes images in MetaImage format.
 *
 * \sa http://www.itk.org/Wiki/MetaIO
 */
class MetaImageWriter : public ImageWriter
{
  mirtkObjectMacro(MetaImageWriter);

public:

  /// Constructor
  MetaImageWriter();

  /// Destructor
  virtual ~MetaImageWriter();

  /// List of file name extensions
  static Array<string> Extensions();

  /// Write image
  virtual void Run();

};


} // namespace mirtk

#endif // MIRTK_MetaImageWriter_H

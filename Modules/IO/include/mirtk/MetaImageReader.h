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

#ifndef MIRTK_MetaImageReader_H
#define MIRTK_MetaImageReader_H

#include "mirtk/ImageReader.h"
#include "mirtk/Memory.h"


class MetaImage;


namespace mirtk {


/**
 * Read images saved in MetaImage format.
 *
 * \sa http://www.itk.org/Wiki/MetaIO
 */
class MetaImageReader : public ImageReader
{
  mirtkObjectMacro(MetaImageReader);

  /// Path/name of image data file
  mirtkReadOnlyAttributeMacro(string, ImageName);

  /// MetaImage instance
  UniquePtr<MetaImage> _MetaImage;

public:

  /// Returns whether file has correct header
  static bool CheckHeader(const char *);

  /// Check if this reader can read a given image file
  virtual bool CanRead(const char *) const;

  /// Constructor
  MetaImageReader();

  /// Destructor
  virtual ~MetaImageReader();

  /// Open file and read image header
  virtual void Initialize();

  /// Read image from file
  ///
  /// \returns Newly read image. Must be deleted by caller.
  virtual BaseImage *Run();

protected:

  /// Copy header information from MetaImage instance
  void CopyHeader(const MetaImage &);

  /// Read header of MetaImage file
  virtual void ReadHeader();

};


} // namespace mirtk

#endif // MIRTK_MetaImageReader_H

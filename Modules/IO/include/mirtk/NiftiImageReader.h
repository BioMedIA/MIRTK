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

#ifndef MIRTK_NiftiImageReader_H
#define MIRTK_NiftiImageReader_H

#include "mirtk/ImageReader.h"
#include "mirtk/Memory.h"


namespace mirtk {


// Forward declaration of internal NIfTI-1 image wrapper
class NiftiImage;


/**
 * Class for reading images in NIFTI file format.
 *
 * This is a class which reads images in NIFTI file format and converts them
 * into images. The NIFTI file format is a file format for 3D and 4D images.
 *
 * \sa http://nifti.nimh.nih.gov/nifti-1/
 */
class NiftiImageReader : public ImageReader
{
  mirtkObjectMacro(NiftiImageReader);

  /// Path/name of image data file
  mirtkReadOnlyAttributeMacro(string, ImageName);

  /// NIfTI image
  UniquePtr<NiftiImage> _Nifti;

public:

  /// Returns whether file has correct header
  static bool CheckHeader(const char *);

  /// Check if this reader can read a given image file
  virtual bool CanRead(const char *) const;

  /// Constructor
  NiftiImageReader();

  /// Destructor
  virtual ~NiftiImageReader();

  /// Open file and read image header
  virtual void Initialize();

  /// Print image file information
  virtual void Print() const;

protected:

  /// Read header of NIFTI file
  virtual void ReadHeader();

};


} // namespace mirtk

#endif // MIRTK_NiftiImageReader_H

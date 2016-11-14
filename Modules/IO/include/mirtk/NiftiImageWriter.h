/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_NiftiImageWriter_H
#define MIRTK_NiftiImageWriter_H

#include "mirtk/ImageWriter.h"

#include "mirtk/Memory.h"
#include "mirtk/Array.h"


namespace mirtk {


// Forward declaration of internal NIfTI-1 image wrapper
class NiftiImage;


/**
 * Class for image to NIFTI file writer.
 *
 * This is a class which takes an image as input and produces an image file
 * in NIFTI file format.
 */
class NiftiImageWriter : public ImageWriter
{
  mirtkObjectMacro(NiftiImageWriter);

  /// NIfTI image
  UniquePtr<NiftiImage> _Nifti;

public:

  /// Constructor
  NiftiImageWriter();

  /// Destructor
  virtual ~NiftiImageWriter();

  /// List of file name extensions
  static Array<string> Extensions();

  /// Write image
  virtual void Run();

protected:

  /// Open file and write header
  virtual void Initialize();

  /// Close file
  virtual void Finalize();

};


} // namespace mirtk

#endif // MIRTK_NiftiImageWriter_H

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

#ifndef MIRTK_ImageWriter_H
#define MIRTK_ImageWriter_H

#include "mirtk/Cofstream.h"
#include "mirtk/BaseImage.h"


namespace mirtk {


/**
 * Abstract base class for any general image to file filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take an image as input and produce an image file as output.
 * Each derived class has to implement all of the abstract member functions.
 */
class ImageWriter : protected Cofstream
{
  mirtkAbstractMacro(ImageWriter);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Input image
  mirtkPublicAggregateMacro(const BaseImage, Input);

  /// Output file name
  mirtkPublicAttributeMacro(string, FileName);

protected:

  /// Start address of the data in the image file.
  /// Should be initialized by overridden Initialize() function.
  int _Start;

  /// Flag whether to reflect X axis
  int _ReflectX;

  /// Flag whether to reflect Y axis
  int _ReflectY;

  /// Flag whether to reflect Z axis
  int _ReflectZ;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Constructor
  ImageWriter();

public:

  /// Destructor
  virtual ~ImageWriter();

  /// Static constructor. This constructor allocates a derived class which
  /// can be used to write the image file in a format corresponding to the
  /// given file name extension.
  static ImageWriter *New(const char *);

  // ---------------------------------------------------------------------------
  // Execution

  /// Write image
  virtual void Run();

protected:

  /// Initialize filter
  virtual void Initialize();

  /// Finalize filter
  virtual void Finalize();

};


} // namespace mirtk

#endif // MIRTK_ImageWriter_H

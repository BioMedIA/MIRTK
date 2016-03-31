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

#ifndef MIRTK_ImageReader_H
#define MIRTK_ImageReader_H

#include "mirtk/Cifstream.h"
#include "mirtk/ImageAttributes.h"
#include "mirtk/BaseImage.h"


namespace mirtk {


/**
 * Abstract base class for any general file to image filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take a filename (referrencing an image file) as input and
 * produce an image as output. For each file format a derived class should
 * be created. Each derived class has to implement abstract member functions
 * ReadHeader() and CheckHeader().
 */
class ImageReader : protected Cifstream
{
  mirtkAbstractMacro(ImageReader);

  /// Path/name of image file
  mirtkPublicAttributeMacro(string, FileName);

  /// Image attributes
  mirtkReadOnlyAttributeMacro(ImageAttributes, Attributes);

  /// Type of voxels
  mirtkReadOnlyAttributeMacro(int, DataType);

  /// No. of bytes per voxel
  mirtkReadOnlyAttributeMacro(int, Bytes);
 
  /// Intensity scaling parameter - slope (default: 1)
  mirtkReadOnlyAttributeMacro(double, Slope);

  /// Intensity scaling parameter -  intercept (default: 0)
  mirtkReadOnlyAttributeMacro(double, Intercept);

protected:

  /// Flag whether to reflect X axis
  int _ReflectX;

  /// Flag whether to reflect Y axis
  int _ReflectY;

  /// Flag whether to reflect Z axis
  int _ReflectZ;

  /// Start of image data
  int _Start;

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Contructor
  ImageReader();

  /// Destructor
  virtual ~ImageReader();

  /// Static constructor. This constructor allocates a derived class which
  /// is then used to read the image file. This involves checking the file
  /// format of the file and dynamically allocating the corresonding derived
  /// class. The reader instance must be deleted by the caller.
  ///
  /// \returns Image reader instance or \c nullptr.
  static ImageReader *TryNew(const char *);

  /// Static constructor. This constructor allocates a derived class which
  /// is then used to read the image file. This involves checking the file
  /// format of the file and dynamically allocating the corresonding derived
  /// class. The reader instance must be deleted by the caller.
  static ImageReader *New(const char *);

  // ---------------------------------------------------------------------------
  // Execution

  /// Check if this reader can read a given image file
  virtual bool CanRead(const char *) const = 0;

  /// Open image file and read header information
  virtual void Initialize();

  /// Print image header information
  virtual void Print() const;

  /// Read image from file
  ///
  /// \returns Newly read image. Must be deleted by caller.
  virtual BaseImage *Run();

protected:

  /// Read header. This is an abstract function. Each derived class has to
  /// implement this function in order to initialize the read-only attributes
  /// of this class.
  virtual void ReadHeader() = 0;

};


} // namespace mirtk

#endif // MIRTK_ImageReader_H

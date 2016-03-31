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

#ifndef MIRTK_ImageFunction_H
#define MIRTK_ImageFunction_H

#include "mirtk/Object.h"
#include "mirtk/BaseImage.h"


namespace mirtk {


/**
 * Abstract base class for any general image interpolation function filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take an image as input and sample that image at arbitrary
 * location. Each derived class has to implement all abstract member functions.
 */
class ImageFunction : public Object
{
  mirtkAbstractMacro(ImageFunction);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Input image
  mirtkPublicAggregateMacro(const BaseImage, Input);

  /// Default value to return
  mirtkPublicAttributeMacro(double, DefaultValue);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Default constructor
  ImageFunction();

  /// Copy constructor
  ImageFunction(const ImageFunction &);

public:

  /// Destructor
  virtual ~ImageFunction();

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize the filter. This function must be called by any derived
  /// filter class to perform some initialize tasks.
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Evaluate the filter at an arbitrary image location (in pixels)
  virtual double Evaluate(double, double, double, double = 0) const = 0;

};


} // namespace mirtk

#endif // MIRTK_ImageFunction_H

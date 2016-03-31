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

#ifndef MIRTK_MirrorExtrapolateImageFunction_H
#define MIRTK_MirrorExtrapolateImageFunction_H

#include "mirtk/ExtrapolateImageFunction.h"
#include "mirtk/BaseImage.h"


namespace mirtk {


/**
 * Extrapolation of generic image by mirroring along the boundaries
 */
template <class TImage>
class GenericMirrorExtrapolateImageFunction
: public IndexExtrapolateImageFunction<TImage>
{
  mirtkExtrapolatorMacro(
    GenericMirrorExtrapolateImageFunction,
    Extrapolation_Mirror
  );

public:

  /// Constructor
  GenericMirrorExtrapolateImageFunction() {}

  /// Destructor
  virtual ~GenericMirrorExtrapolateImageFunction() {}

  /// Mirror index at boundary such that it is inside the range [0, max]
  /// \note Use static function as MirrorExtrapolateImageFunction::Apply.
  static void Apply(int &index, int max)
  {
    if (max == 0) {
      index = 0;
    } else if (index < 0) {
      index = -index;
      int n = index / max;
      int m = index - n * max;
      if (n & 1) index = max - m;
      else       index = m;
    } else if (index > max) {
      index -= max;
      int n = index / max;
      int m = index - n * max;
      if (n & 1) index = m;
      else       index = max - m;
    }
  }

  /// Mirror index at boundary such that it is inside the range [0, max]
  virtual void TransformIndex(int &index, int max) const
  {
    Apply(index, max);
  }

};


/**
 * Extrapolation of any image by mirroring along the boundaries
 */
class MirrorExtrapolateImageFunction
: public GenericMirrorExtrapolateImageFunction<BaseImage>
{
  mirtkObjectMacro(MirrorExtrapolateImageFunction);

public:

  /// Constructor
  MirrorExtrapolateImageFunction() {}

  /// Destructor
  virtual ~MirrorExtrapolateImageFunction() {}

};


} // namespace mirtk

#endif // MIRTK_MirrorExtrapolateImageFunction_H

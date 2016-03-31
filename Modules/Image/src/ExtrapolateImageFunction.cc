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

#include "mirtk/ExtrapolateImageFunction.h"

#include "mirtk/ConstExtrapolateImageFunction.h"
#include "mirtk/ConstExtrapolateImageFunctionWithPeriodicTime.h"
#include "mirtk/NearestNeighborExtrapolateImageFunction.h"
#include "mirtk/RepeatExtrapolateImageFunction.h"
#include "mirtk/MirrorExtrapolateImageFunction.h"


namespace mirtk {


// -----------------------------------------------------------------------------
ExtrapolateImageFunction::ExtrapolateImageFunction()
{
}

// -----------------------------------------------------------------------------
ExtrapolateImageFunction::~ExtrapolateImageFunction()
{
}

// -----------------------------------------------------------------------------
ExtrapolateImageFunction *
ExtrapolateImageFunction::New(enum ExtrapolationMode mode, const BaseImage *image)
{
  ExtrapolateImageFunction *p = NULL;
  switch (mode) {
    case Extrapolation_None:   { p = NULL;                                          break; }
    case Extrapolation_Const:  { p = new ConstExtrapolateImageFunction();           break; }
    case Extrapolation_NN:     { p = new NearestNeighborExtrapolateImageFunction(); break; }
    case Extrapolation_Repeat: { p = new RepeatExtrapolateImageFunction();          break; }
    case Extrapolation_Mirror: { p = new MirrorExtrapolateImageFunction();          break; }
    case Extrapolation_ConstWithPeriodicTime:
      p = new ConstExtrapolateImageFunctionWithPeriodicTime();
      break;
    default:
   	  cerr << "ExtrapolateImageFunction::New: Unknwon extrapolation mode: " << mode << endl;
      exit(1);
  }
  if (p) p->Input(image);
  return p;
}


} // namespace mirtk

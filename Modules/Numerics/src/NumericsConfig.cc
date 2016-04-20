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

#include "mirtk/NumericsConfig.h"
#include "mirtk/ObjectFactory.h"

#ifndef MIRTK_AUTO_REGISTER
  #include "mirtk/GradientDescent.h"
  #include "mirtk/ConjugateGradientDescent.h"
  #if MIRTK_Numerics_WITH_LBFGS
    #include "mirtk/LimitedMemoryBFGSDescent.h"
  #endif
#endif


namespace mirtk {


// -----------------------------------------------------------------------------
static void RegisterOptimizers()
{
  #ifndef MIRTK_AUTO_REGISTER
    mirtkRegisterOptimizerMacro(GradientDescent);
    mirtkRegisterOptimizerMacro(ConjugateGradientDescent);
    #if MIRTK_Numerics_WITH_LBFGS
      mirtkRegisterOptimizerMacro(LimitedMemoryBFGSDescent);
    #endif
  #endif
}

// -----------------------------------------------------------------------------
void InitializeNumericsLibrary()
{
  static bool initialized = false;
  if (!initialized) {
    RegisterOptimizers();
    initialized = true;
  }
}


} // namespace mirtk

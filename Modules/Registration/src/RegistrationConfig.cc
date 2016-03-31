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

#include "mirtk/RegistrationConfig.h"
#include "mirtk/ObjectFactory.h"

#ifndef MIRTK_AUTO_REGISTER
  #include "mirtk/SumOfSquaredIntensityDifferences.h"
  #include "mirtk/MutualImageInformation.h"
  #include "mirtk/NormalizedMutualImageInformation.h"
  #include "mirtk/IntensityCrossCorrelation.h"
  #include "mirtk/NormalizedIntensityCrossCorrelation.h"
  #include "mirtk/CosineOfNormalizedGradientField.h"
  #if MIRTK_Registration_WITH_PointSet
    #include "mirtk/CurrentsDistance.h"
    #include "mirtk/PointCorrespondenceDistance.h"
    #include "mirtk/FiducialRegistrationError.h"
  #endif
#endif


namespace mirtk {


// -----------------------------------------------------------------------------
static void RegisterImageSimilarities()
{
  #ifndef MIRTK_AUTO_REGISTER
    mirtkRegisterEnergyTermMacro(SumOfSquaredIntensityDifferences);
    mirtkRegisterEnergyTermMacro(MutualImageInformation);
    mirtkRegisterEnergyTermMacro(NormalizedMutualImageInformation);
    mirtkRegisterEnergyTermMacro(IntensityCrossCorrelation);
    mirtkRegisterEnergyTermMacro(NormalizedIntensityCrossCorrelation);
    mirtkRegisterEnergyTermMacro(CosineOfNormalizedGradientField);
    #if MIRTK_Registration_WITH_PointSet
      mirtkRegisterEnergyTermMacro(CurrentsDistance);
      mirtkRegisterEnergyTermMacro(PointCorrespondenceDistance);
      mirtkRegisterEnergyTermMacro(FiducialRegistrationError);
    #endif
  #endif
}

// -----------------------------------------------------------------------------
static void RegisterPointSetDistances()
{
  #if !defined(MIRTK_AUTO_REGISTER) && MIRTK_Registration_WITH_PointSet
    mirtkRegisterEnergyTermMacro(CurrentsDistance);
    mirtkRegisterEnergyTermMacro(PointCorrespondenceDistance);
    mirtkRegisterEnergyTermMacro(FiducialRegistrationError);
  #endif
}


// -----------------------------------------------------------------------------
void InitializeRegistrationLibrary()
{
  static bool initialized = false;
  if (!initialized) {
    RegisterImageSimilarities();
    RegisterPointSetDistances();
    initialized = true;
  }
}


} // namespace mirtk

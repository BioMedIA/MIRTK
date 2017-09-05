/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2017 Andreas Schuh
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

#include "mirtk/TransformationConfig.h"
#include "mirtk/ObjectFactory.h"

#ifndef MIRTK_AUTO_REGISTER
  #include "mirtk/SmoothnessConstraint.h"
  #include "mirtk/LinearElasticityConstraint.h"
  #include "mirtk/LogJacobianConstraint.h"
  #include "mirtk/NegJacobianConstraint.h"
  #include "mirtk/TopologyPreservationConstraint.h"
  #include "mirtk/VolumePreservationConstraint.h"
  #include "mirtk/SparsityConstraint.h"
#endif


namespace mirtk {


// -----------------------------------------------------------------------------
static void RegisterTransformationConstraints()
{
  #ifndef MIRTK_AUTO_REGISTER
    mirtkRegisterEnergyTermMacro(SmoothnessConstraint);
    mirtkRegisterEnergyTermMacro(LinearElasticityConstraint);
    mirtkRegisterEnergyTermMacro(LogJacobianConstraint);
    mirtkRegisterEnergyTermMacro(NegJacobianConstraint);
    mirtkRegisterEnergyTermMacro(TopologyPreservationConstraint);
    mirtkRegisterEnergyTermMacro(VolumePreservationConstraint);
    mirtkRegisterEnergyTermMacro(SparsityConstraint);
  #endif
}

// -----------------------------------------------------------------------------
void InitializeTransformationLibrary()
{
  static bool initialized = false;
  if (!initialized) {
    RegisterTransformationConstraints();
    initialized = true;
  }
}


} // namespace mirtk

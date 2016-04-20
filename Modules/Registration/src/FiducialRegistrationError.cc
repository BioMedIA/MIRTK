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

#include "mirtk/FiducialRegistrationError.h"

#include "mirtk/FiducialMatch.h"
#include "mirtk/ObjectFactory.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(FiducialRegistrationError);


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
FiducialRegistrationError::FiducialRegistrationError(const char *name, double weight)
:
  PointCorrespondenceDistance(name, weight, new FiducialMatch())
{
  _ParameterPrefix.push_back("FRE ");
  _ParameterPrefix.push_back("Fiducial registration error ");
  _ParameterPrefix.push_back("Fiducial error ");
  _ParameterPrefix.push_back("Landmark registration error ");
  _ParameterPrefix.push_back("Landmark error ");
}

// -----------------------------------------------------------------------------
FiducialRegistrationError::FiducialRegistrationError(const FiducialRegistrationError &other)
:
  PointCorrespondenceDistance(other)
{
}

// -----------------------------------------------------------------------------
FiducialRegistrationError &FiducialRegistrationError::operator =(const FiducialRegistrationError &other)
{
  PointCorrespondenceDistance::operator =(other);
  return *this;
}

// -----------------------------------------------------------------------------
FiducialRegistrationError::~FiducialRegistrationError()
{
}


} // namespace mirtk

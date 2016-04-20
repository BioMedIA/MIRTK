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

#ifndef MIRTK_FiducialRegistrationError_H
#define MIRTK_FiducialRegistrationError_H

#include "mirtk/PointCorrespondenceDistance.h"


namespace mirtk {


/**
 * Fiducial registration error (FRE) measure
 */
class FiducialRegistrationError : public PointCorrespondenceDistance
{
  mirtkEnergyTermMacro(FiducialRegistrationError, EM_FRE);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  FiducialRegistrationError(const char * = "", double = 1.0);

  /// Copy constructor
  FiducialRegistrationError(const FiducialRegistrationError &);

  /// Assignment operator
  FiducialRegistrationError &operator =(const FiducialRegistrationError &);

  /// Destructor
  virtual ~FiducialRegistrationError();

};


} // namespace mirtk

#endif // MIRTK_FiducialRegistrationError_H

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

#ifndef MIRTK_VolumePreservationConstraint_H
#define MIRTK_VolumePreservationConstraint_H

#include "mirtk/LogJacobianConstraint.h"


namespace mirtk {


/**
 * Volume preservation constraint for deformable image registration
 *
 * Unlike the base class, NegJacobianConstraint, this constraint is
 * by default always applied to the Jacobian of the displacement field.
 *
 *
 */
class VolumePreservationConstraint : public LogJacobianConstraint
{
  mirtkEnergyTermMacro(VolumePreservationConstraint, EM_VolumePreservation);

public:

  /// Constructor
  VolumePreservationConstraint(const char * = "");

  /// Destructor
  virtual ~VolumePreservationConstraint();

};


} // namespace mirtk

#endif // MIRTK_VolumePreservationConstraint_H

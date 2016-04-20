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

#ifndef MIRTK_DisplacementFieldLog_H
#define MIRTK_DisplacementFieldLog_H

#include "mirtk/DisplacementToVelocityFieldBCH.h"


namespace mirtk {


// ----------------------------------------------------------------------------
/// Compute logarithmic map of displacement field
template <class VoxelType>
void Log(GenericImage<VoxelType> *d)
{
  DisplacementToVelocityFieldBCH<VoxelType> dtov;
  dtov.SetInput (d);
  dtov.SetOutput(d);
  dtov.Run();
}


} // namespace mirtk

#endif // MIRTK_DisplacementFieldLog_H

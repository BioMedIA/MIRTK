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

#ifndef MIRTK_FreeFormTransformationIntegration_H
#define MIRTK_FreeFormTransformationIntegration_H


// ============================================================================
// Integration methods
// ============================================================================

// Enumeration of integration methods
#include "mirtk/FFDIntegrationMethod.h"

// Runge-Kutta integration methods
#include "FreeFormTransformationRungeKutta.h"

// ============================================================================
// Auxiliary macros
// ============================================================================

// ---------------------------------------------------------------------------
/// Define integration method for particular FFD
#define MIRTK_FFDIM2(NAME, FFD) \
  typedef mirtk::FreeFormTransformationIntegration##NAME<FFD> NAME;

// ---------------------------------------------------------------------------
/// Define integration method for particular FFD
#define MIRTK_FFDIM3(NAME, FFD, METHOD) \
  typedef mirtk::FreeFormTransformationIntegration##METHOD<FFD> NAME;


#endif // MIRTK_FreeFormTransformationIntegration_H

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

#ifndef MIRTK_Transformations_H
#define MIRTK_Transformations_H


// Homogeneous transformations
#include <mirtkHomogeneousTransformation.h>
#include <mirtkRigidTransformation.h>
#include <mirtkSimilarityTransformation.h>
#include <mirtkAffineTransformation.h>

// Free-form transformations
#include <mirtkFreeFormTransformation.h>
#include <mirtkFreeFormTransformation3D.h>
#include <mirtkFreeFormTransformation4D.h>

#include <mirtkBSplineFreeFormTransformation3D.h>
#include <mirtkBSplineFreeFormTransformation4D.h>
#include <mirtkBSplineFreeFormTransformationSV.h>
#include <mirtkBSplineFreeFormTransformationTD.h>
#include <mirtkBSplineFreeFormTransformationStatistical.h>
#include <mirtkLinearFreeFormTransformation3D.h>
#include <mirtkLinearFreeFormTransformation4D.h>
#include <mirtkLinearFreeFormTransformationTD.h>

// Composite transformations
#include <mirtkMultiLevelTransformation.h>
#include <mirtkMultiLevelFreeFormTransformation.h>
#include <mirtkMultiLevelStationaryVelocityTransformation.h>
#include <mirtkFluidFreeFormTransformation.h>

// Decorators (i.e., wrappers)
#include <mirtkInverseAffineTransformation.h>
#include <mirtkPartialAffineTransformation.h>
#include <mirtkPartialBSplineFreeFormTransformationSV.h>
#include <mirtkPartialMultiLevelStationaryVelocityTransformation.h>

// Image transformation filters
#include <mirtkImageTransformation.h>


#endif // MIRTK_Transformations_H

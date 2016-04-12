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


// Transformation model enumerations and auxiliary functions
#include "mirtk/TransformationModel.h"

// Homogeneous transformations
#include "mirtk/HomogeneousTransformation.h"
#include "mirtk/RigidTransformation.h"
#include "mirtk/SimilarityTransformation.h"
#include "mirtk/AffineTransformation.h"

// Free-form transformations
#include "mirtk/FreeFormTransformation.h"
#include "mirtk/FreeFormTransformation3D.h"
#include "mirtk/FreeFormTransformation4D.h"

#include "mirtk/BSplineFreeFormTransformation3D.h"
#include "mirtk/BSplineFreeFormTransformation4D.h"
#include "mirtk/BSplineFreeFormTransformationSV.h"
#include "mirtk/BSplineFreeFormTransformationTD.h"
#include "mirtk/BSplineFreeFormTransformationStatistical.h"
#include "mirtk/LinearFreeFormTransformation3D.h"
#include "mirtk/LinearFreeFormTransformation4D.h"
#include "mirtk/LinearFreeFormTransformationTD.h"

// Composite transformations
#include "mirtk/MultiLevelTransformation.h"
#include "mirtk/MultiLevelFreeFormTransformation.h"
#include "mirtk/MultiLevelStationaryVelocityTransformation.h"
#include "mirtk/FluidFreeFormTransformation.h"

// Decorators (i.e., wrappers)
#include "mirtk/InverseAffineTransformation.h"
#include "mirtk/PartialAffineTransformation.h"
#include "mirtk/PartialBSplineFreeFormTransformationSV.h"
#include "mirtk/PartialMultiLevelStationaryVelocityTransformation.h"

// Image transformation filters
#include "mirtk/ImageTransformation.h"


#endif // MIRTK_Transformations_H

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

#ifndef MIRTK_Registration_H
#define MIRTK_Registration_H

#include "mirtk/RegistrationConfig.h"

// Base classes
#include "mirtk/DataFidelity.h"
#include "mirtk/ImageSimilarity.h"
#include "mirtk/HistogramImageSimilarity.h"
#include "mirtk/GradientFieldSimilarity.h"

#if MIRTK_Registration_WITH_PointSet
#  include "mirtk/PointSetDistance.h"
#  include "mirtk/SurfaceDistance.h"
#endif // MIRTK_Registration_WITH_PointSet

// Image (dis-)similarities
#include "mirtk/SimilarityMeasure.h"
#include "mirtk/CosineOfNormalizedGradientField.h"
#include "mirtk/IntensityCrossCorrelation.h"
#include "mirtk/MutualImageInformation.h"
#include "mirtk/NormalizedGradientFieldSimilarity.h"
#include "mirtk/NormalizedIntensityCrossCorrelation.h"
#include "mirtk/NormalizedMutualImageInformation.h"
#include "mirtk/SumOfSquaredIntensityDifferences.h"

// Point set distances
#if MIRTK_Registration_WITH_PointSet
#  include "mirtk/PointSetDistanceMeasure.h"
#  include "mirtk/PointCorrespondenceDistance.h"
#  include "mirtk/FiducialRegistrationError.h"
#  include "mirtk/CurrentsDistance.h"
#endif // MIRTK_Registration_WITH_PointSet

// Generic registration filter
#include "mirtk/RegistrationEnergy.h"
#include "mirtk/RegistrationFilter.h"
#include "mirtk/GenericRegistrationFilter.h"
#include "mirtk/GenericRegistrationLogger.h"
#include "mirtk/GenericRegistrationDebugger.h"


#endif // MIRTK_Registration_H

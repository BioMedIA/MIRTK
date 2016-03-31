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

/**
 * \file  mirtkInitialization.h
 * \brief Auxiliary functions used by applications to initialize object factories.
 *
 * The object factories have to be explicitly initialized by an application
 * when the MIRTK modules are build as static libraries. Otherwise, each
 * shared library registers its object types automatically during static
 * initialization. Because this file is common to all applications, it is
 * part of the Common module. As this core module cannot depend on all other
 * modules, it must be a header only to be included by the translation unit
 * that defines the main function of the application.
 */

#ifndef MIRTK_Initialization_H
#define MIRTK_Initialization_H

#include "mirtk/ObjectFactory.h"

#ifndef NDEBUG
#  include "mirtk/Stream.h"
#endif


// =============================================================================
// Iterative optimizers
// =============================================================================

#if defined(HAVE_MIRTK_Numerics) && !defined(MIRTK_AUTO_REGISTER)
  #include "mirtk/GradientDescent.h"
  #include "mirtk/ConjugateGradientDescent.h"
  #ifdef HAVE_MIRTK_ThirdPartyLBFGS
    #include "mirtk/LimitedMemoryBFGSDescent.h"
  #endif
#endif


namespace mirtk {


void RegisterOptimizers()
{
  #ifdef HAVE_MIRTK_Numerics
    #ifndef MIRTK_AUTO_REGISTER
      mirtkRegisterOptimizerMacro(GradientDescent);
      mirtkRegisterOptimizerMacro(ConjugateGradientDescent);
      #ifdef HAVE_MIRTK_ThirdPartyLBFGS
        mirtkRegisterOptimizerMacro(LimitedMemoryBFGSDescent);
      #endif
    #endif
  #elif !defined(NDEBUG)
    cerr << "Warning: Iterative optimizers not available. Rebuild with Numerics module enabled." << endl;
  #endif
}


} // namespace mirtk

// =============================================================================
// Image readers/writers
// =============================================================================

#if defined(HAVE_MIRTK_Image) && !defined(MIRTK_AUTO_REGISTER)
  #include "mirtk/ImageReaderFactory.h"
  #include "mirtk/ImageWriterFactory.h"
  #include "mirtk/GIPLImageReader.h"
  #include "mirtk/GIPLImageWriter.h"
  #include "mirtk/PGMImageReader.h"
  #include "mirtk/PGMImageWriter.h"
  #ifdef HAVE_PNG
    #include "mirtk/PNGImageWriter.h"
  #endif
  #ifdef HAVE_MIRTK_ThirdPartyNIfTI
    #include "mirtk/NiftiImageReader.h"
    #include "mirtk/NiftiImageWriter.h"
  #endif
#endif


namespace mirtk {


void RegisterImageReaders()
{
  #ifdef HAVE_MIRTK_Image
    #ifndef MIRTK_AUTO_REGISTER
      mirtkRegisterImageReaderMacro(GIPLImageReader);
      mirtkRegisterImageReaderMacro(PGMImageReader);
      #ifdef HAVE_MIRTK_ThirdPartyNIfTI
        mirtkRegisterImageReaderMacro(NiftiImageReader);
      #endif
    #endif
  #elif !defined(NDEBUG)
    cerr << "Warning: Image readers not available. Rebuild with Image module enabled." << endl;
  #endif
}

void RegisterImageWriters()
{
  #ifdef HAVE_MIRTK_Image
    #ifndef MIRTK_AUTO_REGISTER
      mirtkRegisterImageWriterMacro(GIPLImageWriter);
      mirtkRegisterImageWriterMacro(PGMImageWriter);
      #ifdef HAVE_PNG
        mirtkRegisterImageWriterMacro(PNGImageWriter);
      #endif
      #ifdef HAVE_MIRTK_ThirdPartyNIfTI
        mirtkRegisterImageWriterMacro(NiftiImageWriter);
      #endif
    #endif
  #elif !defined(NDEBUG)
    cerr << "Warning: Image writers not available. Rebuild with Image module enabled." << endl;
  #endif
}


} // namespace mirtk

// =============================================================================
// Image (dis-)similarities
// =============================================================================

#if defined(HAVE_MIRTK_Image) && defined(HAVE_MIRTK_Registration) && !defined(MIRTK_AUTO_REGISTER)
  #include "mirtk/SumOfSquaredIntensityDifferences.h"
  #include "mirtk/MutualImageInformation.h"
  #include "mirtk/NormalizedMutualImageInformation.h"
  #include "mirtk/IntensityCrossCorrelation.h"
  #include "mirtk/NormalizedIntensityCrossCorrelation.h"
  #include "mirtk/CosineOfNormalizedGradientField.h"
#endif


namespace mirtk {


void RegisterImageSimilarities()
{
  #if defined(HAVE_MIRTK_Image) && defined(HAVE_MIRTK_Registration)
    #ifndef MIRTK_AUTO_REGISTER
      mirtkRegisterEnergyTermMacro(SumOfSquaredIntensityDifferences);
      mirtkRegisterEnergyTermMacro(MutualImageInformation);
      mirtkRegisterEnergyTermMacro(NormalizedMutualImageInformation);
      mirtkRegisterEnergyTermMacro(IntensityCrossCorrelation);
      mirtkRegisterEnergyTermMacro(NormalizedIntensityCrossCorrelation);
      mirtkRegisterEnergyTermMacro(CosineOfNormalizedGradientField);
    #endif
  #elif !defined(NDEBUG)
    cerr << "Warning: Image (dis-)similarities not available. Rebuild with Image and Registration modules enabled." << endl;
  #endif
}


} // namespace mirtk

// =============================================================================
// Point set distances
// =============================================================================

#if defined(HAVE_MIRTK_PointSet) && defined(HAVE_MIRTK_Registration) && !defined(MIRTK_AUTO_REGISTER)
  #include "mirtk/CurrentsDistance.h"
  #include "mirtk/PointCorrespondenceDistance.h"
  #include "mirtk/FiducialRegistrationError.h"
#endif

namespace mirtk {

void RegisterPointSetDistances()
{
  #if defined(HAVE_MIRTK_PointSet) && defined(HAVE_MIRTK_Registration)
    #ifndef MIRTK_AUTO_REGISTER
      mirtkRegisterEnergyTermMacro(CurrentsDistance);
      mirtkRegisterEnergyTermMacro(PointCorrespondenceDistance);
      mirtkRegisterEnergyTermMacro(FiducialRegistrationError);
    #endif
  #elif !defined(NDEBUG)
    cerr << "Warning: Point set distances not available. Rebuild with PointSet and Registration modules enabled." << endl;
  #endif
}


} // namespace mirtk

// =============================================================================
// Transformation constraints
// =============================================================================

#if defined(HAVE_MIRTK_Transformation) && !defined(MIRTK_AUTO_REGISTER)
  #include "mirtk/SmoothnessConstraint.h"
  #include "mirtk/LogJacobianConstraint.h"
  #include "mirtk/MinJacobianConstraint.h"
  #include "mirtk/TopologyPreservationConstraint.h"
  #include "mirtk/VolumePreservationConstraint.h"
  #include "mirtk/SparsityConstraint.h"
#endif


namespace mirtk {


void RegisterTransformationConstraints()
{
  #ifdef HAVE_MIRTK_Transformation
    #ifndef MIRTK_AUTO_REGISTER
      mirtkRegisterEnergyTermMacro(SmoothnessConstraint);
      mirtkRegisterEnergyTermMacro(LogJacobianConstraint);
      mirtkRegisterEnergyTermMacro(MinJacobianConstraint);
      mirtkRegisterEnergyTermMacro(TopologyPreservationConstraint);
      mirtkRegisterEnergyTermMacro(VolumePreservationConstraint);
      mirtkRegisterEnergyTermMacro(SparsityConstraint);
    #endif
  #elif !defined(NDEBUG)
    cerr << "Warning: Transformation constraints not available. Rebuild with Tranformation module enabled." << endl;
  #endif
}


} // namespace mirtk

// =============================================================================
// External forces
// =============================================================================

#if defined(HAVE_MIRTK_Image) && defined(HAVE_MIRTK_PointSet) && defined(HAVE_MIRTK_Deformable) && !defined(MIRTK_AUTO_REGISTER)
  #include "mirtk/BalloonForce.h"
  #include "mirtk/ImageEdgeForce.h"
  #include "mirtk/ImplicitSurfaceDistance.h"
#endif


namespace mirtk {


void RegisterExternalForces()
{
  #if defined(HAVE_MIRTK_Image) && defined(HAVE_MIRTK_PointSet) && defined(HAVE_MIRTK_Deformable)
    #ifndef MIRTK_AUTO_REGISTER
      mirtkRegisterEnergyTermMacro(BalloonForce);
      mirtkRegisterEnergyTermMacro(ImageEdgeForce);
      mirtkRegisterEnergyTermMacro(ImplicitSurfaceDistance);
    #endif
  #elif !defined(NDEBUG)
    cerr << "Warning: External point set forces not available. Rebuild with Image and Deformable modules enabled." << endl;
  #endif
}


} // namespace mirtk

// =============================================================================
// Internal forces
// =============================================================================

#if defined(HAVE_MIRTK_PointSet) && defined(HAVE_MIRTK_Deformable) && !defined(MIRTK_AUTO_REGISTER)
  #include "mirtk/CurvatureConstraint.h"
  #include "mirtk/InflationForce.h"
  #include "mirtk/MetricDistortion.h"
  #include "mirtk/NonSelfIntersectionConstraint.h"
  #include "mirtk/RepulsiveForce.h"
  #include "mirtk/StretchingForce.h"
#endif


namespace mirtk {


void RegisterInternalForces()
{
  #if defined(HAVE_MIRTK_PointSet) && defined(HAVE_MIRTK_Deformable)
    #ifndef MIRTK_AUTO_REGISTER
      mirtkRegisterEnergyTermMacro(CurvatureConstraint);
      mirtkRegisterEnergyTermMacro(InflationForce);
      mirtkRegisterEnergyTermMacro(MetricDistortion);
      mirtkRegisterEnergyTermMacro(NonSelfIntersectionConstraint);
      mirtkRegisterEnergyTermMacro(RepulsiveForce);
      mirtkRegisterEnergyTermMacro(StretchingForce);
    #endif
  #elif !defined(NDEBUG)
    cerr << "Warning: Internal point set forces not available. Rebuild with Deformable module enabled." << endl;
  #endif
}


} // namespace mirtk

// =============================================================================
// Deformable point set forces
// =============================================================================

namespace mirtk {

void RegisterPointSetForces()
{
  RegisterExternalForces();
  RegisterInternalForces();
}

} // namespace mirtk

// =============================================================================
// All energy terms
// =============================================================================

namespace mirtk {

void RegisterEnergyTerms()
{
  RegisterImageSimilarities();
  RegisterPointSetDistances();
  RegisterPointSetForces();
  RegisterTransformationConstraints();
}

} // namespace mirtk

#endif // MIRTK_Initialization_H

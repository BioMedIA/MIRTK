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

#ifndef MIRTK_VoxelFunction_H
#define MIRTK_VoxelFunction_H

#include "mirtk/ImageAttributes.h"
#include "mirtk/VoxelDomain.h"
#include "mirtk/Parallel.h"
#include "mirtk/Stream.h"


namespace mirtk {


////////////////////////////////////////////////////////////////////////////////
// Base classes
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Base class for voxel functions
// =============================================================================

/**
 * Base class for voxel functions
 *
 * The basic voxel functions which are part of this library and included by the
 * mirtkVoxelFunction.h module demonstrate the implementation and use of the
 * ForEachVoxel function templates. These templates can be well optimized by the
 * compiler and provide the fastest way of iterating over all voxels in a given
 * image region (or the entire image) of one or more images. When processing
 * multiple images at the same time, these must have the same image size.
 *
 * Example usage:
 * \code
 * // Determine intensity range of grayscale image
 * GreyImage image(attr);
 * UnaryVoxelFunction::GetMinMax minmax;
 * ForEachVoxel(image, minmax);
 * double min = minmax.GetMinAsDouble();
 * double max = minmax.GetMaxAsDouble();
 *
 * // Add three images within a given region and store result in another
 * NaryVoxelFunction::Sum sum;
 * GreyImage input1(attr);
 * GreyImage input2(attr);
 * GreyImage input3(attr);
 * GreyImage output(attr);
 * blocked_range3d<int> region(1, 2, 0, attr._y, 0, attr._x);
 * ForEachVoxel(region, input1, input2, input3, output, sum);
 *
 * // Subtract one image from another
 * // Note: Voxel function passed in this case by copy instead of reference!
 * //       This is only possible if it is not a voxel reduction.
 * GreyImage minuend   (attr);
 * GreyImage subtrahend(attr);
 * ParallelForEachVoxel(BinaryVoxelFunction::Sub(), subtrahend, minuend);
 * \endcode
 *
 * If a voxel function needs to consider also the intensities of neighboring
 * voxels, it can access these by manipulating the given pointer to the current
 * voxel to be processed. Useful helpers for this are the NeighborhoodOffsets
 * and ImageRegion classes.
 *
 * \sa NeighborhoodOffsets, ImageRegion, ForEachVoxelDomain
 */
struct VoxelFunction
{
  /// Finite image domain of images processed by this voxel function
  ///
  /// Attributes of first image which should be identical to other images,
  /// at least the number of voxels must be identical. This const pointer,
  /// when not set before by the user of the voxel function or its constructor,
  /// is set by the ForEachVoxelBody constructor called by the ForEachVoxel
  /// and ParallelForEachVoxel template functions.
  ///
  /// These image attributes give access to commonly needed functions within
  /// the voxel function operator() implementation such as:
  /// - ImageAttributes::LatticeToWorld
  /// - ImageAttributes::WorldToLattice
  /// - ImageAttributes::IsBoundary
  /// and of course the actual data members of the ImageAttributes structure.
  const ImageAttributes *_Domain;

  /// Default constructor
  VoxelFunction() : _Domain(nullptr) {}

  /// Used by ParallelForEachVoxel to determine if voxel function has to
  /// be executed using parallel_reduce
  static bool IsReduction() { return false; }

  /// Split "constructor"
  ///
  /// \note A method is used instead of an actual split constructor such that
  ///       subclasses which are not used for reduction do not need to
  ///       implement such constructor.
  void split(VoxelFunction &) {}

  /// Join results
  void join(VoxelFunction &) {}
};


/**
 * Base class for voxel functions which implement a reduction of voxel values
 *
 * Voxel functions derived from this base class are run by ParallelForEachVoxel
 * using parallel_reduce instead of parallel_for. Typical voxel reduction
 * functions are unary functions computing image statistics such as the min,
 * max, mean, or variance of voxel values within an image region.
 */
struct VoxelReduction : public VoxelFunction
{
  /// Used by ParallelForEachVoxel to determine if voxel function has to
  /// be executed using parallel_reduce
  static bool IsReduction() { return true; }

  /// Split "constructor"
  void split(VoxelFunction &)
  {
    cerr << "VoxelReduction::split must be overriden by each subclass!" << endl;
    cerr << "Otherwise you should use VoxelFunction with parallel_for instead." << endl;
    exit(1);
  }

  /// Join results
  void join(VoxelFunction &)
  {
    cerr << "VoxelReduction::join must be overriden by each subclass!" << endl;
    cerr << "Otherwise you should use VoxelFunction with parallel_for instead." << endl;
    exit(1);
  }
};

// =============================================================================
// Base class for ForEachVoxel body with single voxel function
// =============================================================================

/**
 * Base class for ForEachVoxel template function body with single voxel
 * function for each voxel
 */
template <class VoxelFunc>
struct ForEachVoxelBody
{
  // ---------------------------------------------------------------------------
  // Members

  VoxelFunc _VoxelFunc; ///< Functor executed for each voxel
  int       _k, _l;     ///< Indices for fixed dimensions

  // ---------------------------------------------------------------------------
  // Construction

  /// Constructor
  ForEachVoxelBody(const VoxelFunc &vf, const ImageAttributes &attr)
  :
    _VoxelFunc(vf), _k(0), _l(0)
  {
    if (!_VoxelFunc.VoxelFunction::_Domain) {
      _VoxelFunc.VoxelFunction::_Domain = &attr;
    }
  }

  /// Copy constructor
  ForEachVoxelBody(const ForEachVoxelBody &o)
  :
    _VoxelFunc(o._VoxelFunc), _k(o._k), _l(o._l)
  {}

  // ---------------------------------------------------------------------------
  // Parallel reduction

  /// Split constructor
  ForEachVoxelBody(ForEachVoxelBody &o, split s)
  :
    _VoxelFunc(o._VoxelFunc), _k(o._k), _l(o._l)
  {
    _VoxelFunc.split(o._VoxelFunc);
  }

  /// Join results
  void join(ForEachVoxelBody &rhs)
  {
    _VoxelFunc.join(rhs._VoxelFunc);
  }
};

// =============================================================================
// Base class for ForEachVoxelIf body with separate inside/outside voxel functions
// =============================================================================

/**
 * Base class for ForEachVoxelIf template function body with separate voxel
 * function for inside and outside voxels
 */
template <class VoxelFunc, class OutsideFunc>
struct ForEachVoxelIfBody : public ForEachVoxelBody<VoxelFunc>
{
  // ---------------------------------------------------------------------------
  // Members

  OutsideFunc _OutsideFunc; ///< Functor executed for each background voxel

  // ---------------------------------------------------------------------------
  // Construction

  /// Constructor
  ForEachVoxelIfBody(const VoxelFunc &vf, const OutsideFunc &of, const ImageAttributes &attr)
  :
    ForEachVoxelBody<VoxelFunc>(vf, attr), _OutsideFunc(of)
  {}

  /// Copy constructor
  ForEachVoxelIfBody(const ForEachVoxelIfBody &o)
  :
    ForEachVoxelBody<VoxelFunc>(o), _OutsideFunc(o._OutsideFunc)
  {}

  // ---------------------------------------------------------------------------
  // Parallel reduction

  /// Split constructor
  ForEachVoxelIfBody(ForEachVoxelIfBody &o, split s)
  :
    ForEachVoxelBody<VoxelFunc>(o, s),
    _OutsideFunc(o._OutsideFunc)
  {
    _OutsideFunc.split(o._OutsideFunc);
  }

  /// Join results
  void join(ForEachVoxelIfBody &rhs)
  {
    ForEachVoxelBody<VoxelFunc>::join(rhs);
    _OutsideFunc.join(rhs._OutsideFunc);
  }
};

// =============================================================================
// NOP voxel function
// =============================================================================

namespace NaryVoxelFunction {


/**
 * NOP voxel function, e.g., used as default outside function by ForEachVoxelIf
 */
struct NOP : public VoxelFunction
{
  /// Unary voxel function operator
  template <class TImage>
  void operator()(const TImage&, int, const void *) const {}
  /// Unary voxel function operator
  void operator()(int, int, int, int, const void *) const {}

  /// Binary voxel function operator
  template <class TImage>
  void operator()(const TImage&, int, const void *, const void *) const {}
  /// Binary voxel function operator
  void operator()(int, int, int, int, const void *, const void *) const {}

  /// Ternary voxel function operator
  template <class TImage>
  void operator()(const TImage&, int, const void *, const void *, const void *) const {}
  /// Ternary voxel function operator
  void operator()(int, int, int, int, const void *, const void *, const void *) const {}

  /// Quaternary voxel function operator
  template <class TImage>
  void operator()(const TImage&, int, const void *, const void *, const void *, const void *) const {}
  /// Quaternary voxel function operator
  void operator()(int, int, int, int, const void *, const void *, const void *, const void *) const {}

  /// Quinery voxel function operator
  template <class TImage>
  void operator()(const TImage&, int, const void *, const void *, const void *, const void *, const void *) const {}
  /// Quinery voxel function operator
  void operator()(int, int, int, int, const void *, const void *, const void *, const void *, const void *) const {}

  /// Senary voxel function operator
  template <class TImage>
  void operator()(const TImage&, int, const void *, const void *, const void *, const void *, const void *, const void *) const {}
  /// Senary voxel function operator
  void operator()(int, int, int, int, const void *, const void *, const void *, const void *, const void *, const void *) const {}

  /// Septenary voxel function operator
  template <class TImage>
  void operator()(const TImage&, int, const void *, const void *, const void *, const void *, const void *, const void *, const void *) const {}
  /// Septenary voxel function operator
  void operator()(int, int, int, int, const void *, const void *, const void *, const void *, const void *, const void *, const void *) const {}

  /// Octary voxel function operator
  template <class TImage>
  void operator()(const TImage&, int, const void *, const void *, const void *, const void *, const void *, const void *, const void *, const void *) const {}
  /// Octary voxel function operator
  void operator()(int, int, int, int, const void *, const void *, const void *, const void *, const void *, const void *, const void *, const void *) const {}

  /// Nonary voxel function operator
  template <class TImage>
  void operator()(const TImage&, int, const void *, const void *, const void *, const void *, const void *, const void *, const void *, const void *, const void *) const {}
  /// Nonary voxel function operator
  void operator()(int, int, int, int, const void *, const void *, const void *, const void *, const void *, const void *, const void *, const void *, const void *) const {}
};


} // namespace NaryVoxelFunction
} // namespace mirtk

// =============================================================================
// ForEachVoxel template functions for various number of image arguments
// =============================================================================

#include "mirtk/ForEachUnaryVoxelFunction.h"
#include "mirtk/ForEachBinaryVoxelFunction.h"
#include "mirtk/ForEachTernaryVoxelFunction.h"
#include "mirtk/ForEachQuaternaryVoxelFunction.h"
#include "mirtk/ForEachQuinaryVoxelFunction.h"
#include "mirtk/ForEachSenaryVoxelFunction.h"
#include "mirtk/ForEachSeptenaryVoxelFunction.h"
#include "mirtk/ForEachOctaryVoxelFunction.h"
#include "mirtk/ForEachNonaryVoxelFunction.h"


#endif // MIRTK_VoxelFunction_H

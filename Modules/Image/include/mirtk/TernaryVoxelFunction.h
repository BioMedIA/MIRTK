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

#ifndef MIRTK_TernaryVoxelFunction_H
#define MIRTK_TernaryVoxelFunction_H

#include "mirtk/VoxelFunction.h"


namespace mirtk {


/**
 * These basic ternary voxel functions can be used as VoxelFunc template parameter
 * of the ternary ForEachVoxel function templates with three image arguments, e.g.:
 *
 * \code
 * GreyImage input1(attr);
 * GreyImage input2(attr);
 * GreyImage ouptut(attr);
 * // Compute voxel-wise sum of input1 and input2
 * TernaryVoxelFunction::Sum sum;
 * ForEachVoxel(input1, input2, output, sum);
 * // Compute voxel-wise input1/input2 for each foreground voxel
 * TernaryVoxelFunction::Div div;
 * ParallelForEachVoxelInside(input1, input2, output, div);
 * \endcode
 */
namespace TernaryVoxelFunction {


// -----------------------------------------------------------------------------
/**
 * Sums up the voxel values of two images
 */
struct Sum : public VoxelFunction
{
  template <class T1, class T2, class T3>
  void operator ()(const T1 *in1, const T2 *in2, T3 *out)
  {
    *out = static_cast<T3>(static_cast<double>(*in1) + static_cast<double>(*in2));
  }

  template <class TImage, class T1, class T2, class T3>
  void operator ()(const TImage&, int, const T1 *in1, const T2 *in2, T3 *out)
  {
    this->operator ()(in1, in2, out);
  }

  template <class T1, class T2, class T3>
  void operator ()(int, int, int, int, const T1 *in1, const T2 *in2, T3 *out)
  {
    this->operator ()(in1, in2, out);
  }
};

// -----------------------------------------------------------------------------
/**
 * Computes the voxel intensity difference of two images
 */
struct Diff : public VoxelFunction
{
  template <class T1, class T2, class T3>
  void operator ()(const T1 *in1, const T2 *in2, T3 *out)
  {
    *out = static_cast<T3>(static_cast<double>(*in1) - static_cast<double>(*in2));
  }

  template <class TImage, class T1, class T2, class T3>
  void operator ()(const TImage&, int, const T1 *in1, const T2 *in2, T3 *out)
  {
    this->operator ()(in1, in2, out);
  }

  template <class T1, class T2, class T3>
  void operator ()(int, int, int, int, const T1 *in1, const T2 *in2, T3 *out)
  {
    this->operator ()(in1, in2, out);
  }
};

// -----------------------------------------------------------------------------
/**
 * Calculates the product of the voxel values of two images
 */
struct Mul : public VoxelFunction
{
  template <class T1, class T2, class T3>
  void operator ()(const T1 *in1, const T2 *in2, T3 *out)
  {
    *out = static_cast<T3>(static_cast<double>(*in1) * static_cast<double>(*in2));
  }

  template <class TImage, class T1, class T2, class T3>
  void operator ()(const TImage&, int, const T1 *in1, const T2 *in2, T3 *out)
  {
    this->operator ()(in1, in2, out);
  }

  template <class T1, class T2, class T3>
  void operator ()(int, int, int, int, const T1 *in1, const T2 *in2, T3 *out)
  {
    this->operator ()(in1, in2, out);
  }
};

/// Alternative name for Mul voxel function
typedef Mul Prod;

// -----------------------------------------------------------------------------
/**
 * Calculates the divison of the voxel values of two images
 */
struct Div : public VoxelFunction
{
  template <class T1, class T2, class T3>
  void operator ()(const T1 *in1, const T2 *in2, T3 *out)
  {
    double divisor = static_cast<double>(*in2);
    if (divisor == .0) {
      *out = static_cast<T3>(0);
    } else {
      *out = static_cast<T3>(static_cast<double>(*in1) / divisor);
    }
  }

  template <class TImage, class T1, class T2, class T3>
  void operator ()(const TImage&, int, const T1 *in1, const T2 *in2, T3 *out)
  {
    this->operator ()(in1, in2, out);
  }

  template <class T1, class T2, class T3>
  void operator ()(int, int, int, int, const T1 *in1, const T2 *in2, T3 *out)
  {
    this->operator ()(in1, in2, out);
  }
};


} } // namespace mirtk::TernaryVoxelFunction

#endif // MIRTK_TernaryVoxelFunction_H

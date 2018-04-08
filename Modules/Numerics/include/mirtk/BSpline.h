/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2018 Andreas Schuh
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

#ifndef MIRTK_BSpline_H
#define MIRTK_BSpline_H

#include "mirtk/Config.h"
#include "mirtk/Math.h"
#include "mirtk/Stream.h"

#include "mirtk/NumericsExport.h"


namespace mirtk {


/**
 * Cubic B-spline, its basis functions, and their derivatives.
 *
 * This static class provides methods to evaluate the cubic B-spline function,
 * each of its four basis functions, and the corresponding first and second
 * derivatives in the interval [0, 1]. It is in particular intended for use
 * by classes implementing the irtkFreeFormTransformation interface using
 * a cubic B-spline basis to represent the displacement or velocity field,
 * respectively.
 *
 * The cubic B-spline pieces are indexed here in reverse order!
 * This is especially important when using the derivatives, because otherwise
 * the sign may be incorrect when the intervals are exchanged.
 * Spline function B(t) consists of pieces:
 * - B3(t + 2) for t in [-2, -1)
 * - B2(t + 1) for t in [-1,  0)
 * - B1(t - 0) for t in [ 0,  1)
 * - B0(t - 1) for t in [ 1,  2]
 *
 * \note Though not explicitly enforced by protecting these members,
 *       do <b>not</b> modify the lookup tables of precomputed function
 *       values as it will affect all code which makes use of them!
 *       Moreover, make sure to call the static Initialize() method at least
 *       once in your program before accessing any of these precomputed values.
 *
 * Usage:
 * \code
 * BSplineFunction kernel;
 * // evaluate function explicitly
 * double exact = kernel.B1(0.2);
 * // or use lookup table; note that the Initialize() method must only be
 * // called once during the life-time of the program
 * BSplineFunction::Initialize();
 * double approx = kernel.LookupTable[kernel.VariableToIndex(0.2)][1];
 * \endcode
 */

template <class TReal>
class BSpline
{
public:

  // ---------------------------------------------------------------------------
  // Generel 1D interpolation kernel interface

  /// Floating point precision type used
  typedef TReal RealType;

  /// Lookup table of B-spline basis function values for t=0
  /// (i.e., kernel centered at lattice coordinate)
  static const TReal LatticeWeights[4];

  /// Lookup table of B-spline basis function 1st derivative values for t=0
  /// (i.e., kernel centered at lattice coordinate)
  static const TReal LatticeWeights_I[4];

  /// Lookup table of B-spline basis function 2nd derivative values for t=0
  /// (i.e., kernel centered at lattice coordinate)
  static const TReal LatticeWeights_II[4];

  /// Returns the value of the B-spline function
  MIRTKCU_API static TReal Weight(TReal);

  /// Returns the values of the B-spline function for
  /// - [0]: t in ]-2, -1[
  /// - [1]: t in [-1, 0]
  /// - [2]: t in ]0, 1]
  /// - [3]: t in ]1, 2[
  MIRTKCU_API static void Weights(TReal, TReal[4]);

  /// Size of lookup tables of pre-computed B-spline values
  static const unsigned int LookupTableSize = 1000;

  /// Initialize lookup tables
  static void Initialize(bool = false);

  /// Returns the lookup table index for a given value in [0,1]
  static int VariableToIndex(TReal);

  /// Returns the lookup table indices for any t value
  ///
  /// \param[in]  t Cubic B-spline function parameter.
  /// \param[out] i First lookup table index in [0, LookupTableSize - 1].
  /// \param[out] j Second lookup table index in [0, 3].
  ///
  /// \returns Whether t is within local support of cubic B-spline function.
  ///          When the return value is false, lookup table indices corresponding
  ///          to a boundary of the local support region is returned.
  static bool VariableToIndex(TReal, int &i, int &j);

  /// Lookup table of B-spline function values
  MIRTK_Numerics_EXPORT static TReal WeightLookupTable[LookupTableSize];

  /// Lookup table of B-spline basis function values
  MIRTK_Numerics_EXPORT static TReal LookupTable[LookupTableSize][4];

  /// Lookup table of B-spline basis function 1st derivative values
  MIRTK_Numerics_EXPORT static TReal LookupTable_I[LookupTableSize][4];

  /// Lookup table of B-spline basis function 2nd derivative values
  MIRTK_Numerics_EXPORT static TReal LookupTable_II[LookupTableSize][4];

  // ---------------------------------------------------------------------------
  // B-spline basis functions

  /// Returns the value of the B-spline function
  MIRTKCU_API static TReal B(TReal);

  /// Returns the value of the i-th B-spline basis function
  MIRTKCU_API static TReal B(int, TReal);

  /// Returns the value of the first B-spline basis function
  MIRTKCU_API static TReal B0(TReal);

  /// Returns the value of the second B-spline basis function
  MIRTKCU_API static TReal B1(TReal);

  /// Returns the value of the third B-spline basis function
  MIRTKCU_API static TReal B2(TReal);

  /// Returns the value of the fourth B-spline basis function
  MIRTKCU_API static TReal B3(TReal);

  /// Returns the 1st derivative value of the B-spline function
  MIRTKCU_API static TReal B_I(TReal);

  /// Returns the 1st derivative value of the i-th B-spline basis function
  MIRTKCU_API static TReal B_I(int, TReal);
  
  /// Returns the 1st derivative value of the first B-spline basis function
  MIRTKCU_API static TReal B0_I(TReal);

  /// Returns the 1st derivative value of the second B-spline basis function
  MIRTKCU_API static TReal B1_I(TReal);

  /// Returns the 1st derivative value of the third B-spline basis function
  MIRTKCU_API static TReal B2_I(TReal);

  /// Returns the 1st derivative value of the fourth B-spline basis function
  MIRTKCU_API static TReal B3_I(TReal);

  /// Returns the 2nd derivative value of the B-spline function
  MIRTKCU_API static TReal B_II(TReal);

  /// Returns the 2nd derivative value of the i-th B-spline basis function
  MIRTKCU_API static TReal B_II(int, TReal);

  /// Returns the 2nd derivative value of the first B-spline basis function
  MIRTKCU_API static TReal B0_II(TReal);

  /// Returns the 2nd derivative value of the second B-spline basis function
  MIRTKCU_API static TReal B1_II(TReal);

  /// Returns the 2nd derivative value of the third B-spline basis function
  MIRTKCU_API static TReal B2_II(TReal);

  /// Returns the 2nd derivative value of the fourth B-spline basis function
  MIRTKCU_API static TReal B3_II(TReal);

  /// Returns the 3rd derivative value of the B-spline function
  MIRTKCU_API static TReal B_III(TReal);

  /// Returns the 3rd derivative value of the i-th B-spline basis function
  MIRTKCU_API static TReal B_III(int, TReal);

  /// Returns the 3rd derivative value of the first B-spline basis function
  MIRTKCU_API static TReal B0_III(TReal);

  /// Returns the 3rd derivative value of the second B-spline basis function
  MIRTKCU_API static TReal B1_III(TReal);

  /// Returns the 3rd derivative value of the third B-spline basis function
  MIRTKCU_API static TReal B2_III(TReal);

  /// Returns the 3rd derivative value of the fourth B-spline basis function
  MIRTKCU_API static TReal B3_III(TReal);

  /// Returns the n-th derivative value of the B-spline function
  MIRTKCU_API static TReal B_nI(int, TReal);

  /// Returns the n-th derivative value of the i-th B-spline basis function
  MIRTKCU_API static TReal B_nI(int, int, TReal);

  /// Returns the n-th derivative value of the first B-spline basis function
  MIRTKCU_API static TReal B0_nI(int, TReal);

  /// Returns the n-th derivative value of the second B-spline basis function
  MIRTKCU_API static TReal B1_nI(int, TReal);

  /// Returns the n-th derivative value of the third B-spline basis function
  MIRTKCU_API static TReal B2_nI(int, TReal);

  /// Returns the n-th derivative value of the fourth B-spline basis function
  MIRTKCU_API static TReal B3_nI(int, TReal);

protected:

  /// Flag which indicates whether the lookup tables are initialized
  MIRTK_Numerics_EXPORT static bool _initialized;
};

// -----------------------------------------------------------------------------
// Fix Clang warning regarding missing definition of static template members
#ifdef __clang__

// Lookup table of B-spline function values
template <class TReal> TReal BSpline<TReal>::WeightLookupTable[BSpline<TReal>::LookupTableSize];

// Lookup table of B-spline basis function values
template <class TReal> TReal BSpline<TReal>::LookupTable   [BSpline<TReal>::LookupTableSize][4];
template <class TReal> TReal BSpline<TReal>::LookupTable_I [BSpline<TReal>::LookupTableSize][4];
template <class TReal> TReal BSpline<TReal>::LookupTable_II[BSpline<TReal>::LookupTableSize][4];

// Wether lookup tables of B-spline kernel were initialized
template <class TReal> bool BSpline<TReal>::_initialized = false;

#endif // defined(__clang__)

// -----------------------------------------------------------------------------
// Weights often used by B-spline transformations for evaluation of derivatives
// at lattice points only where B-spline parameter t is zero
template <class TReal>
const TReal BSpline<TReal>::LatticeWeights[4] = {
  TReal(1.0/6.0),  TReal(2.0/3.0), TReal(1.0/6.0), TReal(0.0)
};

template <class TReal>
const TReal BSpline<TReal>::LatticeWeights_I[4] = {
  TReal(-0.5), TReal(0.0), TReal(0.5), TReal(0.0)
};

template <class TReal>
const TReal BSpline<TReal>::LatticeWeights_II[4] = {
  TReal(1.0), TReal(-2.0), TReal(1.0), TReal(0.0)
};

// -----------------------------------------------------------------------------
typedef BSpline<double> BSplineFunction;

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <class TReal>
inline void BSpline<TReal>::Weights(TReal t, TReal value[4])
{
  TReal c = static_cast<TReal>(1.0 - t);
  TReal b = t * t;
  TReal a = b * t;
  value[0] = static_cast<TReal>(c * c * c / 6.0);
  value[1] = static_cast<TReal>(( 3.0 * a - 6.0 * b + 4.0) / 6.0);
  value[2] = static_cast<TReal>((-3.0 * a + 3.0 * b + 3.0 * t + 1.0) / 6.0);
  value[3] = static_cast<TReal>(a / 6.0);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::Weight(TReal t)
{
  double value = 0.0;
  t = fabs(t);
  if (t < 2.0) {
    if (t < 1.0) {
      value = 2.0/3.0 + (0.5*t-1.0)*t*t;  // B1
    } else {
      t -= 2.0;
      value = -t*t*t/6.0;  // B3
    }
  }
  return static_cast<TReal>(value);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B(TReal t)
{
  return Weight(t);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B(int i, TReal t)
{
  switch (i) {
    case 0:  return B0(t);
    case 1:  return B1(t);
    case 2:  return B2(t);
    case 3:  return B3(t);
    default: return TReal(0);
  }
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B0(TReal t)
{
  return static_cast<TReal>((1.0-t)*(1.0-t)*(1.0-t)/6.0);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B1(TReal t)
{
  return static_cast<TReal>((3.0*t*t*t - 6.0*t*t + 4.0)/6.0);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B2(TReal t)
{
  return static_cast<TReal>((-3.0*t*t*t + 3.0*t*t + 3.0*t + 1.0)/6.0);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B3(TReal t)
{
  return static_cast<TReal>((t*t*t)/6.0);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B_I(TReal t)
{
  double y     = abs(t);
  double value = 0.0;
  if (y < 2.0) {
    if (y < 1.0) {
      value = (1.5 * y - 2.0) * t;
    } else {
      y -= 2.0;
      value = -0.5 * y * y;
      if (t < 0.0) value = -value;
    }
  }
  return static_cast<TReal>(value);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B_I(int i, TReal t)
{
  switch (i) {
    case 0:  return B0_I(t);
    case 1:  return B1_I(t);
    case 2:  return B2_I(t);
    case 3:  return B3_I(t);
    default: return TReal(0);
  }
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B0_I(TReal t)
{
  return static_cast<TReal>(-(1.0-t)*(1.0-t)/2.0);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B1_I(TReal t)
{
  return static_cast<TReal>((9.0*t*t - 12.0*t)/6.0);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B2_I(TReal t)
{
  return static_cast<TReal>((-9.0*t*t + 6.0*t + 3.0)/6.0);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B3_I(TReal t)
{
  return static_cast<TReal>((t*t)/2.0);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B_II(TReal t)
{
  double value = 0.0;
  t = fabs(t);
  if (t < 2.0) {
    if (t < 1.0) {
      value = 3.0 * t - 2.0;
    } else {
      value = -t + 2.0;
    }
  }
  return static_cast<TReal>(value);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B_II(int i, TReal t)
{
  switch (i) {
    case 0:  return B0_II(t);
    case 1:  return B1_II(t);
    case 2:  return B2_II(t);
    case 3:  return B3_II(t);
    default: return TReal(0);
  }
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B0_II(TReal t)
{
  return static_cast<TReal>(1.0 - t);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B1_II(TReal t)
{
  return static_cast<TReal>(3.0*t - 2.0);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B2_II(TReal t)
{
  return static_cast<TReal>(-3.0*t + 1.0);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B3_II(TReal t)
{
  return t;
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B_III(TReal t)
{
  TReal y = fabs(t), value;
  if (y < 2.0) {
    if (y < 1.0) {
      value = copysign(TReal(3.0),  t);
    } else {
      value = copysign(TReal(1.0), -t);
    }
  } else {
    value = TReal(0);
  }
  return value;
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B_III(int i, TReal t)
{
  switch (i) {
    case 0:  return B0_III(t);
    case 1:  return B1_III(t);
    case 2:  return B2_III(t);
    case 3:  return B3_III(t);
    default: return TReal(0);
  }
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B0_III(TReal)
{
  return TReal(-1.0);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B1_III(TReal)
{
  return TReal(3.0);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B2_III(TReal)
{
  return TReal(-3.0);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B3_III(TReal)
{
  return TReal(1.0);
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B_nI(int n, TReal t)
{
  switch(n) {
    case 0:  return B    (t);
    case 1:  return B_I  (t);
    case 2:  return B_II (t);
    case 3:  return B_III(t);
    default: return TReal(0);
  }
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B_nI(int i, int n, TReal t)
{
  switch(n) {
    case 0:  return B    (i, t);
    case 1:  return B_I  (i, t);
    case 2:  return B_II (i, t);
    case 3:  return B_III(i, t);
    default: return TReal(0);
  }
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B0_nI(int n, TReal t)
{
  switch(n) {
    case 0:  return B0    (t);
    case 1:  return B0_I  (t);
    case 2:  return B0_II (t);
    case 3:  return B0_III(t);
    default: return TReal(0);
  }
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B1_nI(int n, TReal t)
{
  switch(n) {
    case 0:  return B1    (t);
    case 1:  return B1_I  (t);
    case 2:  return B1_II (t);
    case 3:  return B1_III(t);
    default: return TReal(0);
  }
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B2_nI(int n, TReal t)
{
  switch(n) {
    case 0:  return B2    (t);
    case 1:  return B2_I  (t);
    case 2:  return B2_II (t);
    case 3:  return B2_III(t);
    default: return TReal(0);
  }
}

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal BSpline<TReal>::B3_nI(int n, TReal t)
{
  switch(n) {
    case 0:  return B3    (t);
    case 1:  return B3_I  (t);
    case 2:  return B3_II (t);
    case 3:  return B3_III(t);
    default: return TReal(0);
  }
}

// -----------------------------------------------------------------------------
template <class TReal>
inline void BSpline<TReal>::Initialize(bool force)
{
  if (!_initialized || force) {
    for (unsigned int i = 0; i < LookupTableSize; i++) {
      WeightLookupTable[i] = Weight(static_cast<TReal>(i) / static_cast<TReal>(LookupTableSize-1));
      for (int l = 0; l < 4; l++) {
        LookupTable   [i][l] = B   (l, static_cast<TReal>(i) / static_cast<TReal>(LookupTableSize-1));
        LookupTable_I [i][l] = B_I (l, static_cast<TReal>(i) / static_cast<TReal>(LookupTableSize-1));
        LookupTable_II[i][l] = B_II(l, static_cast<TReal>(i) / static_cast<TReal>(LookupTableSize-1));
      }
    }
    _initialized = true;
  }
}

// -----------------------------------------------------------------------------
template <class TReal>
inline int BSpline<TReal>::VariableToIndex(TReal t)
{
  if      (t <  .0) return 0;
  else if (t > 1.0) return LookupTableSize-1;
  else              return iround(t * static_cast<TReal>(LookupTableSize-1));
}

// -----------------------------------------------------------------------------
template <class TReal>
inline bool BSpline<TReal>::VariableToIndex(TReal t, int &i, int &j)
{
  if (t < 0.) {
    if (t <= -2.) {
      j = 3;
      i = 0;
      return false;
    }
    if (t < -1.) {
      j = 3;
      i = iround((t + 2.) * static_cast<TReal>(LookupTableSize - 1));
    } else {
      j = 2;
      i = iround((t + 1.) * static_cast<TReal>(LookupTableSize - 1));
    }
  } else {
    if (t >= 2.) {
      j = 0;
      i = LookupTableSize - 1;
      return false;
    }
    if (t < 1.) {
      j = 1;
      i = iround(t * static_cast<TReal>(LookupTableSize - 1));
    } else {
      j = 0;
      i = iround((t - 1.) * static_cast<TReal>(LookupTableSize - 1));
    }
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// B-spline interpolation weights
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
/// Compute indices and weights required for 2D B-spline interpolation
template <class TReal>
inline void ComputeBSplineIndicesAndWeights(double x, double y, int spline_degree,
                                            int   xIndex [6], int   yIndex [6],
                                            TReal xWeight[6], TReal yWeight[6])
{
  // Compute the interpolation indices
  int i, j;
  if (spline_degree & 1) {
    i = static_cast<int>(floor(x      )) - spline_degree / 2;
    j = static_cast<int>(floor(y      )) - spline_degree / 2;
  } else {
    i = static_cast<int>(floor(x + 0.5)) - spline_degree / 2;
    j = static_cast<int>(floor(y + 0.5)) - spline_degree / 2;
  }
  for (int m = 0; m <= spline_degree; m++) {
    xIndex[m] = i++;
    yIndex[m] = j++;
  }

  // Compute the interpolation weights
  double w, w2, w4, t, t0, t1;

  switch (spline_degree) {
    case 2:
      /* x */
      w = x - static_cast<double>(xIndex[1]);
      xWeight[1] = TReal(3.0 / 4.0 - w * w);
      xWeight[2] = TReal((1.0 / 2.0) * (w - xWeight[1] + 1.0));
      xWeight[0] = TReal(1.0 - xWeight[1] - xWeight[2]);
      /* y */
      w = y - static_cast<double>(yIndex[1]);
      yWeight[1] = TReal(3.0 / 4.0 - w * w);
      yWeight[2] = TReal((1.0 / 2.0) * (w - yWeight[1] + 1.0));
      yWeight[0] = TReal(1.0 - yWeight[1] - yWeight[2]);
      break;
    case 3:
      /* x */
      w = x - static_cast<double>(xIndex[1]);
      xWeight[3] = TReal((1.0 / 6.0) * w * w * w);
      xWeight[0] = TReal((1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3]);
      xWeight[2] = TReal(w + xWeight[0] - 2.0 * xWeight[3]);
      xWeight[1] = TReal(1.0 - xWeight[0] - xWeight[2] - xWeight[3]);
      /* y */
      w = y - static_cast<double>(yIndex[1]);
      yWeight[3] = TReal((1.0 / 6.0) * w * w * w);
      yWeight[0] = TReal((1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - yWeight[3]);
      yWeight[2] = TReal(w + yWeight[0] - 2.0 * yWeight[3]);
      yWeight[1] = TReal(1.0 - yWeight[0] - yWeight[2] - yWeight[3]);
      break;
    case 4:
      /* x */
      w = x - static_cast<double>(xIndex[2]);
      w2 = w * w;
      t = (1.0 / 6.0) * w2;
      xWeight[0] = TReal(1.0 / 2.0 - w);
      xWeight[0] *= xWeight[0];
      xWeight[0] *= TReal(1.0 / 24.0) * xWeight[0];
      t0 = w * (t - 11.0 / 24.0);
      t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
      xWeight[1] = TReal(t1 + t0);
      xWeight[3] = TReal(t1 - t0);
      xWeight[4] = TReal(xWeight[0] + t0 + (1.0 / 2.0) * w);
      xWeight[2] = TReal(1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4]);
      /* y */
      w = y - static_cast<double>(yIndex[2]);
      w2 = w * w;
      t = (1.0 / 6.0) * w2;
      yWeight[0] = TReal(1.0 / 2.0 - w);
      yWeight[0] *= yWeight[0];
      yWeight[0] *= TReal(1.0 / 24.0) * yWeight[0];
      t0 = w * (t - 11.0 / 24.0);
      t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
      yWeight[1] = TReal(t1 + t0);
      yWeight[3] = TReal(t1 - t0);
      yWeight[4] = TReal(yWeight[0] + t0 + (1.0 / 2.0) * w);
      yWeight[2] = TReal(1.0 - yWeight[0] - yWeight[1] - yWeight[3] - yWeight[4]);
      break;
    case 5:
      /* x */
      w = x - static_cast<double>(xIndex[2]);
      w2 = w * w;
      xWeight[5] = TReal((1.0 / 120.0) * w * w2 * w2);
      w2 -= w;
      w4 = w2 * w2;
      w -= 1.0 / 2.0;
      t = w2 * (w2 - 3.0);
      xWeight[0] = TReal((1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5]);
      t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
      t1 = (-1.0 / 12.0) * w * (t + 4.0);
      xWeight[2] = TReal(t0 + t1);
      xWeight[3] = TReal(t0 - t1);
      t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
      t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
      xWeight[1] = TReal(t0 + t1);
      xWeight[4] = TReal(t0 - t1);
      /* y */
      w = y - static_cast<double>(yIndex[2]);
      w2 = w * w;
      yWeight[5] = TReal((1.0 / 120.0) * w * w2 * w2);
      w2 -= w;
      w4 = w2 * w2;
      w -= 1.0 / 2.0;
      t = w2 * (w2 - 3.0);
      yWeight[0] = TReal((1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - yWeight[5]);
      t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
      t1 = (-1.0 / 12.0) * w * (t + 4.0);
      yWeight[2] = TReal(t0 + t1);
      yWeight[3] = TReal(t0 - t1);
      t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
      t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
      yWeight[1] = TReal(t0 + t1);
      yWeight[4] = TReal(t0 - t1);
      break;
    default:
      cerr << "ComputeBSplineIndicesAndWeights: Unsupported B-spline degree: " << spline_degree << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
/// Compute indices and weights required for 3D B-spline interpolation
template <class TReal>
inline void ComputeBSplineIndicesAndWeights(double x, double y, double z, int spline_degree,
                                            int   xIndex [6], int   yIndex [6], int   zIndex [6],
                                            TReal xWeight[6], TReal yWeight[6], TReal zWeight[6])
{
  // Compute the interpolation indices
  int i, j, k, m;

  if (spline_degree & 1) {
    i = static_cast<int>(floor(x)) - spline_degree / 2;
    j = static_cast<int>(floor(y)) - spline_degree / 2;
    k = static_cast<int>(floor(z)) - spline_degree / 2;
    for (m = 0; m <= spline_degree; m++) {
      xIndex[m] = i++;
      yIndex[m] = j++;
      zIndex[m] = k++;
    }
  } else {
    i = static_cast<int>(floor(x + 0.5)) - spline_degree / 2;
    j = static_cast<int>(floor(y + 0.5)) - spline_degree / 2;
    k = static_cast<int>(floor(z + 0.5)) - spline_degree / 2;
    for (m = 0; m <= spline_degree; m++) {
      xIndex[m] = i++;
      yIndex[m] = j++;
      zIndex[m] = k++;
    }
  }

  // Compute the interpolation weights
  double w, w2, w4, t, t0, t1;

  switch (spline_degree) {
    case 2:
      /* x */
      w = x - static_cast<double>(xIndex[1]);
      xWeight[1] = TReal(3.0 / 4.0 - w * w);
      xWeight[2] = TReal((1.0 / 2.0) * (w - xWeight[1] + 1.0));
      xWeight[0] = TReal(1.0 - xWeight[1] - xWeight[2]);
      /* y */
      w = y - static_cast<double>(yIndex[1]);
      yWeight[1] = TReal(3.0 / 4.0 - w * w);
      yWeight[2] = TReal((1.0 / 2.0) * (w - yWeight[1] + 1.0));
      yWeight[0] = TReal(1.0 - yWeight[1] - yWeight[2]);
      /* z */
      w = z - static_cast<double>(zIndex[1]);
      zWeight[1] = TReal(3.0 / 4.0 - w * w);
      zWeight[2] = TReal((1.0 / 2.0) * (w - zWeight[1] + 1.0));
      zWeight[0] = TReal(1.0 - zWeight[1] - zWeight[2]);
      break;
    case 3:
      /* x */
      w = x - static_cast<double>(xIndex[1]);
      xWeight[3] = TReal((1.0 / 6.0) * w * w * w);
      xWeight[0] = TReal((1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3]);
      xWeight[2] = TReal(w + xWeight[0] - 2.0 * xWeight[3]);
      xWeight[1] = TReal(1.0 - xWeight[0] - xWeight[2] - xWeight[3]);
      /* y */
      w = y - static_cast<double>(yIndex[1]);
      yWeight[3] = TReal((1.0 / 6.0) * w * w * w);
      yWeight[0] = TReal((1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - yWeight[3]);
      yWeight[2] = TReal(w + yWeight[0] - 2.0 * yWeight[3]);
      yWeight[1] = TReal(1.0 - yWeight[0] - yWeight[2] - yWeight[3]);
      /* z */
      w = z - static_cast<double>(zIndex[1]);
      zWeight[3] = TReal((1.0 / 6.0) * w * w * w);
      zWeight[0] = TReal((1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - zWeight[3]);
      zWeight[2] = TReal(w + zWeight[0] - 2.0 * zWeight[3]);
      zWeight[1] = TReal(1.0 - zWeight[0] - zWeight[2] - zWeight[3]);
      break;
    case 4:
      /* x */
      w = x - static_cast<double>(xIndex[2]);
      w2 = w * w;
      t = (1.0 / 6.0) * w2;
      xWeight[0] = TReal(1.0 / 2.0 - w);
      xWeight[0] *= xWeight[0];
      xWeight[0] *= TReal(1.0 / 24.0) * xWeight[0];
      t0 = w * (t - 11.0 / 24.0);
      t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
      xWeight[1] = TReal(t1 + t0);
      xWeight[3] = TReal(t1 - t0);
      xWeight[4] = TReal(xWeight[0] + t0 + (1.0 / 2.0) * w);
      xWeight[2] = TReal(1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4]);
      /* y */
      w = y - static_cast<double>(yIndex[2]);
      w2 = w * w;
      t = (1.0 / 6.0) * w2;
      yWeight[0] = TReal(1.0 / 2.0 - w);
      yWeight[0] *= yWeight[0];
      yWeight[0] *= TReal(1.0 / 24.0) * yWeight[0];
      t0 = w * (t - 11.0 / 24.0);
      t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
      yWeight[1] = TReal(t1 + t0);
      yWeight[3] = TReal(t1 - t0);
      yWeight[4] = TReal(yWeight[0] + t0 + (1.0 / 2.0) * w);
      yWeight[2] = TReal(1.0 - yWeight[0] - yWeight[1] - yWeight[3] - yWeight[4]);
      /* z */
      w = z - static_cast<double>(zIndex[2]);
      w2 = w * w;
      t = (1.0 / 6.0) * w2;
      zWeight[0] = TReal(1.0 / 2.0 - w);
      zWeight[0] *= zWeight[0];
      zWeight[0] *= TReal(1.0 / 24.0) * zWeight[0];
      t0 = w * (t - 11.0 / 24.0);
      t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
      zWeight[1] = TReal(t1 + t0);
      zWeight[3] = TReal(t1 - t0);
      zWeight[4] = TReal(zWeight[0] + t0 + (1.0 / 2.0) * w);
      zWeight[2] = TReal(1.0 - zWeight[0] - zWeight[1] - zWeight[3] - zWeight[4]);
      break;
    case 5:
      /* x */
      w = x - static_cast<double>(xIndex[2]);
      w2 = w * w;
      xWeight[5] = TReal((1.0 / 120.0) * w * w2 * w2);
      w2 -= w;
      w4 = w2 * w2;
      w -= 1.0 / 2.0;
      t = w2 * (w2 - 3.0);
      xWeight[0] = TReal((1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5]);
      t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
      t1 = (-1.0 / 12.0) * w * (t + 4.0);
      xWeight[2] = TReal(t0 + t1);
      xWeight[3] = TReal(t0 - t1);
      t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
      t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
      xWeight[1] = TReal(t0 + t1);
      xWeight[4] = TReal(t0 - t1);
      /* y */
      w = y - static_cast<double>(yIndex[2]);
      w2 = w * w;
      yWeight[5] = TReal((1.0 / 120.0) * w * w2 * w2);
      w2 -= w;
      w4 = w2 * w2;
      w -= 1.0 / 2.0;
      t = w2 * (w2 - 3.0);
      yWeight[0] = TReal((1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - yWeight[5]);
      t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
      t1 = (-1.0 / 12.0) * w * (t + 4.0);
      yWeight[2] = TReal(t0 + t1);
      yWeight[3] = TReal(t0 - t1);
      t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
      t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
      yWeight[1] = TReal(t0 + t1);
      yWeight[4] = TReal(t0 - t1);
      /* z */
      w = z - static_cast<double>(zIndex[2]);
      w2 = w * w;
      zWeight[5] = TReal((1.0 / 120.0) * w * w2 * w2);
      w2 -= w;
      w4 = w2 * w2;
      w -= 1.0 / 2.0;
      t = w2 * (w2 - 3.0);
      zWeight[0] = TReal((1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - zWeight[5]);
      t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
      t1 = (-1.0 / 12.0) * w * (t + 4.0);
      zWeight[2] = TReal(t0 + t1);
      zWeight[3] = TReal(t0 - t1);
      t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
      t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
      zWeight[1] = TReal(t0 + t1);
      zWeight[4] = TReal(t0 - t1);
      break;
    default:
      cerr << "ComputeBSplineIndicesAndWeights: Unsupported B-spline degree: " << spline_degree << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
/// Compute indices and weights required for 4D B-spline interpolation
template <class TReal>
inline void ComputeBSplineIndicesAndWeights(double x, double y, double z, double t, int spline_degree,
                                            int   xIndex [6], int   yIndex [6], int   zIndex [6], int   tIndex[6],
                                            TReal xWeight[6], TReal yWeight[6], TReal zWeight[6], TReal tWeight[6])
{
  // Compute the interpolation indices
  int i, j, k, l, m;

  if (spline_degree & 1) {
    i = ifloor(x) - spline_degree / 2;
    j = ifloor(y) - spline_degree / 2;
    k = ifloor(z) - spline_degree / 2;
    l = ifloor(t) - spline_degree / 2;
    for (m = 0; m <= spline_degree; m++) {
      xIndex[m] = i++;
      yIndex[m] = j++;
      zIndex[m] = k++;
      tIndex[m] = l++;
    }
  } else {
    i = ifloor(x + .5) - spline_degree / 2;
    j = ifloor(y + .5) - spline_degree / 2;
    k = ifloor(z + .5) - spline_degree / 2;
    l = ifloor(t + .5) - spline_degree / 2;
    for (m = 0; m <= spline_degree; m++) {
      xIndex[m] = i++;
      yIndex[m] = j++;
      zIndex[m] = k++;
      tIndex[m] = l++;
    }
  }

  // Compute the interpolation weights
  double w, w2, w4, t0, t1;

  switch (spline_degree) {
    case 2:
      /* x */
      w = x - static_cast<double>(xIndex[1]);
      xWeight[1] = TReal(3.0 / 4.0 - w * w);
      xWeight[2] = TReal((1.0 / 2.0) * (w - xWeight[1] + 1.0));
      xWeight[0] = TReal(1.0 - xWeight[1] - xWeight[2]);
      /* y */
      w = y - static_cast<double>(yIndex[1]);
      yWeight[1] = TReal(3.0 / 4.0 - w * w);
      yWeight[2] = TReal((1.0 / 2.0) * (w - yWeight[1] + 1.0));
      yWeight[0] = TReal(1.0 - yWeight[1] - yWeight[2]);
      /* z */
      w = z - static_cast<double>(zIndex[1]);
      zWeight[1] = TReal(3.0 / 4.0 - w * w);
      zWeight[2] = TReal((1.0 / 2.0) * (w - zWeight[1] + 1.0));
      zWeight[0] = TReal(1.0 - zWeight[1] - zWeight[2]);
      /* t */
      w = t - static_cast<double>(tIndex[1]);
      tWeight[1] = TReal(3.0 / 4.0 - w * w);
      tWeight[2] = TReal((1.0 / 2.0) * (w - tWeight[1] + 1.0));
      tWeight[0] = TReal(1.0 - tWeight[1] - tWeight[2]);
      break;
    case 3:
      /* x */
      w = x - static_cast<double>(xIndex[1]);
      xWeight[3] = TReal((1.0 / 6.0) * w * w * w);
      xWeight[0] = TReal((1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3]);
      xWeight[2] = TReal(w + xWeight[0] - 2.0 * xWeight[3]);
      xWeight[1] = TReal(1.0 - xWeight[0] - xWeight[2] - xWeight[3]);
      /* y */
      w = y - static_cast<double>(yIndex[1]);
      yWeight[3] = TReal((1.0 / 6.0) * w * w * w);
      yWeight[0] = TReal((1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - yWeight[3]);
      yWeight[2] = TReal(w + yWeight[0] - 2.0 * yWeight[3]);
      yWeight[1] = TReal(1.0 - yWeight[0] - yWeight[2] - yWeight[3]);
      /* z */
      w = z - static_cast<double>(zIndex[1]);
      zWeight[3] = TReal((1.0 / 6.0) * w * w * w);
      zWeight[0] = TReal((1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - zWeight[3]);
      zWeight[2] = TReal(w + zWeight[0] - 2.0 * zWeight[3]);
      zWeight[1] = TReal(1.0 - zWeight[0] - zWeight[2] - zWeight[3]);
      /* t */
      w = t - static_cast<double>(tIndex[1]);
      tWeight[3] = TReal((1.0 / 6.0) * w * w * w);
      tWeight[0] = TReal((1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - tWeight[3]);
      tWeight[2] = TReal(w + tWeight[0] - 2.0 * tWeight[3]);
      tWeight[1] = TReal(1.0 - tWeight[0] - tWeight[2] - tWeight[3]);
      break;
    case 4:
      /* x */
      w = x - static_cast<double>(xIndex[2]);
      w2 = w * w;
      t = (1.0 / 6.0) * w2;
      xWeight[0] = TReal(1.0 / 2.0 - w);
      xWeight[0] *= xWeight[0];
      xWeight[0] *= TReal(1.0 / 24.0) * xWeight[0];
      t0 = w * (t - 11.0 / 24.0);
      t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
      xWeight[1] = TReal(t1 + t0);
      xWeight[3] = TReal(t1 - t0);
      xWeight[4] = TReal(xWeight[0] + t0 + (1.0 / 2.0) * w);
      xWeight[2] = TReal(1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4]);
      /* y */
      w = y - static_cast<double>(yIndex[2]);
      w2 = w * w;
      t = (1.0 / 6.0) * w2;
      yWeight[0] = TReal(1.0 / 2.0 - w);
      yWeight[0] *= yWeight[0];
      yWeight[0] *= TReal(1.0 / 24.0) * yWeight[0];
      t0 = w * (t - 11.0 / 24.0);
      t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
      yWeight[1] = TReal(t1 + t0);
      yWeight[3] = TReal(t1 - t0);
      yWeight[4] = TReal(yWeight[0] + t0 + (1.0 / 2.0) * w);
      yWeight[2] = TReal(1.0 - yWeight[0] - yWeight[1] - yWeight[3] - yWeight[4]);
      /* z */
      w = z - static_cast<double>(zIndex[2]);
      w2 = w * w;
      t = (1.0 / 6.0) * w2;
      zWeight[0] = TReal(1.0 / 2.0 - w);
      zWeight[0] *= zWeight[0];
      zWeight[0] *= TReal(1.0 / 24.0) * zWeight[0];
      t0 = w * (t - 11.0 / 24.0);
      t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
      zWeight[1] = TReal(t1 + t0);
      zWeight[3] = TReal(t1 - t0);
      zWeight[4] = TReal(zWeight[0] + t0 + (1.0 / 2.0) * w);
      zWeight[2] = TReal(1.0 - zWeight[0] - zWeight[1] - zWeight[3] - zWeight[4]);
      /* t */
      w = t - static_cast<double>(tIndex[2]);
      w2 = w * w;
      t = (1.0 / 6.0) * w2;
      tWeight[0] = TReal(1.0 / 2.0 - w);
      tWeight[0] *= tWeight[0];
      tWeight[0] *= TReal(1.0 / 24.0) * tWeight[0];
      t0 = w * (t - 11.0 / 24.0);
      t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
      tWeight[1] = TReal(t1 + t0);
      tWeight[3] = TReal(t1 - t0);
      tWeight[4] = TReal(tWeight[0] + t0 + (1.0 / 2.0) * w);
      tWeight[2] = TReal(1.0 - tWeight[0] - tWeight[1] - tWeight[3] - tWeight[4]);
      break;
    case 5:
      /* x */
      w = x - static_cast<double>(xIndex[2]);
      w2 = w * w;
      xWeight[5] = TReal((1.0 / 120.0) * w * w2 * w2);
      w2 -= w;
      w4 = w2 * w2;
      w -= 1.0 / 2.0;
      t = w2 * (w2 - 3.0);
      xWeight[0] = TReal((1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5]);
      t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
      t1 = (-1.0 / 12.0) * w * (t + 4.0);
      xWeight[2] = TReal(t0 + t1);
      xWeight[3] = TReal(t0 - t1);
      t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
      t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
      xWeight[1] = TReal(t0 + t1);
      xWeight[4] = TReal(t0 - t1);
      /* y */
      w = y - static_cast<double>(yIndex[2]);
      w2 = w * w;
      yWeight[5] = TReal((1.0 / 120.0) * w * w2 * w2);
      w2 -= w;
      w4 = w2 * w2;
      w -= 1.0 / 2.0;
      t = w2 * (w2 - 3.0);
      yWeight[0] = TReal((1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - yWeight[5]);
      t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
      t1 = (-1.0 / 12.0) * w * (t + 4.0);
      yWeight[2] = TReal(t0 + t1);
      yWeight[3] = TReal(t0 - t1);
      t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
      t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
      yWeight[1] = TReal(t0 + t1);
      yWeight[4] = TReal(t0 - t1);
      /* z */
      w = z - static_cast<double>(zIndex[2]);
      w2 = w * w;
      zWeight[5] = TReal((1.0 / 120.0) * w * w2 * w2);
      w2 -= w;
      w4 = w2 * w2;
      w -= 1.0 / 2.0;
      t = w2 * (w2 - 3.0);
      zWeight[0] = TReal((1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - zWeight[5]);
      t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
      t1 = (-1.0 / 12.0) * w * (t + 4.0);
      zWeight[2] = TReal(t0 + t1);
      zWeight[3] = TReal(t0 - t1);
      t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
      t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
      zWeight[1] = TReal(t0 + t1);
      zWeight[4] = TReal(t0 - t1);
      /* t */
      w = t - static_cast<double>(tIndex[2]);
      w2 = w * w;
      tWeight[5] = TReal((1.0 / 120.0) * w * w2 * w2);
      w2 -= w;
      w4 = w2 * w2;
      w -= 1.0 / 2.0;
      t = w2 * (w2 - 3.0);
      tWeight[0] = TReal((1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - tWeight[5]);
      t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
      t1 = (-1.0 / 12.0) * w * (t + 4.0);
      tWeight[2] = TReal(t0 + t1);
      tWeight[3] = TReal(t0 - t1);
      t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
      t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
      tWeight[1] = TReal(t0 + t1);
      tWeight[4] = TReal(t0 - t1);
      break;
    default:
      cerr << "ComputeBSplineIndicesAndWeights: Unsupported B-spline degree: " << spline_degree << endl;
      exit(1);
  }
}


} // namespace mirtk

#endif // MIRTK_BSpline_H

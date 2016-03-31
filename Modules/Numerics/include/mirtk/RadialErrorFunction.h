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

#ifndef MIRTK_RadialErrorFunction_H
#define MIRTK_RadialErrorFunction_H

#include "mirtk/Object.h"


namespace mirtk {


/**
 * Abstract radial fiducial registration error (FRE) function
 */
class RadialErrorFunction : public Object
{
  mirtkAbstractMacro(RadialErrorFunction);

protected:

  /// Constructor
  RadialErrorFunction();

public:

  /// Enumeration of available error functions
  enum TypeId {
    Unknown,     ///< Unknown/invalid/default function
    Distance,    ///< (Euclidean) distance
    Squared,     ///< Squared (Euclidean) distance
    Gaussian,    ///< Gaussian error function
    Charbonnier, ///< Charbonnier fiducial registration error
    PeronaMalik  ///< Perona-Malik fiducial registration error
  };

  /// Construct a new instance of specified type
  static RadialErrorFunction *New(TypeId);

  /// Construct a new instance of specified type
  static RadialErrorFunction *New(const char *);

  /// Copy construct a new instance
  virtual RadialErrorFunction *NewInstance() const = 0;

  /// Destructor
  virtual ~RadialErrorFunction() = 0;

  /// Type enumeration value
  virtual TypeId Type() const = 0;

  /// Evaluate radial registration error
  ///
  /// \param[in] d Squared (Euclidean) distance.
  virtual double Value(double d) const = 0;

  /// Evaluate derivative of radial registration error
  ///
  /// \param[in] d Squared (Euclidean) distance.
  virtual double Derivative(double d) const = 0;

};


// -----------------------------------------------------------------------------
template <> bool FromString(const char *str, RadialErrorFunction::TypeId &);
template <> string ToString(const RadialErrorFunction::TypeId &, int, char, bool);


} // namespace mirtk

#endif // MIRTK_RadialErrorFunction_H

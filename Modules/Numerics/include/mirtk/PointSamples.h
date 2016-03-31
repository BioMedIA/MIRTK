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

#ifndef MIRTK_PointSamples_H
#define MIRTK_PointSamples_H

#include "mirtk/PointSet.h"


namespace mirtk {


/**
 * Auxiliary class for generation of point samples
 */
class PointSamples : public PointSet
{
  mirtkObjectMacro(PointSamples);

  /// Random number generator
  void *_RandomNumberGenerator;

public:

  // ---------------------------------------------------------------------------
  // Construction/destruction

  /// Constructor
  ///
  /// \param[in] n    Number of point samples.
  /// \param[in] seed Seed of random number generator. Pass a negative value to
  ///                 use the current timestamp as seed value.
  PointSamples(int n = 0, int seed = 0);

  /// Destructor
  virtual ~PointSamples();

  // ---------------------------------------------------------------------------
  // Regular grid sampling

  /// Sample axes-aligned uniform grid
  void SampleGrid(const Point &p1, const Point &p2,
                  int nx, int ny, int nz);

  /// Sample axes-aligned uniform grid
  void SampleGrid(const Point &p1, const Point &p2,
                  double dx, double dy, double dz);

  /// Sample axes-aligned uniform grid
  void SampleGrid(double x1, double y1, double z1,
                  double x2, double y2, double z2,
                  int    nx, int    ny, int    nz);

  /// Sample axes-aligned uniform grid
  void SampleGrid(double x1, double y1, double z1,
                  double x2, double y2, double z2,
                  double dx, double dy, double dz);

  // ---------------------------------------------------------------------------
  // Uniform spherical distribution

  /// Add uniform spherical point samples
  void SampleSphere(double r = 1.0);

  /// Add uniform spherical point samples
  ///
  /// \param[in] c Center in each dimension.
  /// \param[in] r Radius in each dimension.
  void SampleSphere(double c, double r);

  /// Add uniform spherical point samples
  ///
  /// \param[in] c Center point.
  /// \param[in] r Radius in each dimension.
  void SampleSphere(const Point & c, double r = 1.0);

  /// Add uniform spherical point samples
  ///
  /// \param[in] c  Center point.
  /// \param[in] rx Radius in x direction.
  /// \param[in] ry Radius in x direction.
  /// \param[in] rz Radius in x direction.
  void SampleSphere(const Point & c, double rx, double ry, double rz);

  /// Add uniform spherical point samples
  ///
  /// \param[in] cx Center in x direction.
  /// \param[in] cy Center in y direction.
  /// \param[in] cz Center in z direction.
  /// \param[in] r  Radius in each dimension.
  void SampleSphere(double cx, double cy, double cz, double r);

  /// Add uniform spherical point samples
  ///
  /// \param[in] cx Center in x direction.
  /// \param[in] cy Center in y direction.
  /// \param[in] cz Center in z direction.
  /// \param[in] rx Radius in x direction.
  /// \param[in] ry Radius in x direction.
  /// \param[in] rz Radius in x direction.
  void SampleSphere(double cx, double cy, double cz,
                    double rx, double ry, double rz);

  // ---------------------------------------------------------------------------
  // Regular spherical samples

  /// Add regular spherical point samples
  void SampleRegularSphere(double r = 1.0);

  /// Add regular spherical point samples
  ///
  /// \param[in] c Center in each dimension.
  /// \param[in] r Radius in each dimension.
  void SampleRegularSphere(double c, double r);

  /// Add regular spherical point samples
  ///
  /// \param[in] c Center point.
  /// \param[in] r Radius in each dimension.
  void SampleRegularSphere(const Point & c, double r = 1.0);

  /// Add regular spherical point samples
  ///
  /// \param[in] c  Center point.
  /// \param[in] rx Radius in x direction.
  /// \param[in] ry Radius in x direction.
  /// \param[in] rz Radius in x direction.
  void SampleRegularSphere(const Point & c, double rx, double ry, double rz);

  /// Add regular spherical point samples
  ///
  /// \param[in] cx Center in x direction.
  /// \param[in] cy Center in y direction.
  /// \param[in] cz Center in z direction.
  /// \param[in] r  Radius in each dimension.
  void SampleRegularSphere(double cx, double cy, double cz, double r);

  /// Add regular spherical point samples
  ///
  /// \param[in] cx Center in x direction.
  /// \param[in] cy Center in y direction.
  /// \param[in] cz Center in z direction.
  /// \param[in] rx Radius in x direction.
  /// \param[in] ry Radius in x direction.
  /// \param[in] rz Radius in x direction.
  void SampleRegularSphere(double cx, double cy, double cz,
                           double rx, double ry, double rz);

  // ---------------------------------------------------------------------------
  // Regular spherical sampling of half sphere

  /// Add regular spherical point samples
  void SampleRegularHalfSphere(double r = 1.0);

  /// Add regular spherical point samples
  ///
  /// \param[in] c Center in each dimension.
  /// \param[in] r Radius in each dimension.
  void SampleRegularHalfSphere(double c, double r);

  /// Add regular spherical point samples
  ///
  /// \param[in] c Center point.
  /// \param[in] r Radius in each dimension.
  void SampleRegularHalfSphere(const Point & c, double r = 1.0);

  /// Add regular spherical point samples
  ///
  /// \param[in] c  Center point.
  /// \param[in] rx Radius in x direction.
  /// \param[in] ry Radius in x direction.
  /// \param[in] rz Radius in x direction.
  void SampleRegularHalfSphere(const Point & c, double rx, double ry, double rz);

  /// Add regular spherical point samples
  ///
  /// \param[in] cx Center in x direction.
  /// \param[in] cy Center in y direction.
  /// \param[in] cz Center in z direction.
  /// \param[in] r  Radius in each dimension.
  void SampleRegularHalfSphere(double cx, double cy, double cz, double r);

  /// Add regular spherical point samples
  ///
  /// \param[in] cx Center in x direction.
  /// \param[in] cy Center in y direction.
  /// \param[in] cz Center in z direction.
  /// \param[in] rx Radius in x direction.
  /// \param[in] ry Radius in x direction.
  /// \param[in] rz Radius in x direction.
  void SampleRegularHalfSphere(double cx, double cy, double cz,
                               double rx, double ry, double rz);

  // ---------------------------------------------------------------------------
  // Normal distribution

  /// Add normally distributed point samples
  ///
  /// \param[in] s Standard deviation in each dimension.
  void SampleGaussian(double s = 1.0);

  /// Add normally distributed point samples
  ///
  /// \param[in] m Mean in each dimension.
  /// \param[in] s Standard deviation in each dimension.
  void SampleGaussian(double m, double s);

  /// Add normally distributed point samples
  ///
  /// \param[in] mx Mean in x direction.
  /// \param[in] my Mean in y direction.
  /// \param[in] mz Mean in z direction.
  /// \param[in] s  Standard deviation in each dimension.
  void SampleGaussian(double mx, double my, double mz, double s);

  /// Add normally distributed point samples
  ///
  /// \param[in] m Mean of normal distribution.
  /// \param[in] s Standard deviation in each dimension.
  void SampleGaussian(const Point & m, double s);

  /// Add normally distributed point samples
  ///
  /// \param[in] m  Mean of normal distribution.
  /// \param[in] sx Standard deviation in x direction.
  /// \param[in] sy Standard deviation in y direction.
  /// \param[in] sz Standard deviation in z direction.
  void SampleGaussian(const Point & m, double sx, double sy, double sz);

  /// Add normally distributed point samples
  ///
  /// \param[in] mx Mean in x direction.
  /// \param[in] my Mean in y direction.
  /// \param[in] mz Mean in z direction.
  /// \param[in] sx Standard deviation in x direction.
  /// \param[in] sy Standard deviation in y direction.
  /// \param[in] sz Standard deviation in z direction.
  void SampleGaussian(double mx, double my, double mz,
                      double sx, double sy, double sz);

};


} // namespace mirtk

#endif // MIRTK_PointSamples_H

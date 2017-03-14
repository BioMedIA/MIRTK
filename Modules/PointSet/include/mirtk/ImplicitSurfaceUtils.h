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

#ifndef MIRTK_ImplicitSurfaceUtils_H
#define MIRTK_ImplicitSurfaceUtils_H

#include "mirtk/Assert.h"
#include "mirtk/Math.h"
#include "mirtk/Array.h"
#include "mirtk/Algorithm.h"
#include "mirtk/GenericImage.h"
#include "mirtk/RegisteredImage.h"
#include "mirtk/InterpolateImageFunction.h"
#include "mirtk/ImageGradientFunction.h"
#include "mirtk/PointSamples.h"
#include "mirtk/PointSetUtils.h"

// Include inline definitions of image interpolation functions
#include "mirtk/LinearInterpolateImageFunction.hxx"
#include "mirtk/FastLinearImageGradientFunction.hxx"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"


namespace mirtk {


// TODO: Wrap these functions in a class named ImplicitSurface, which is
//       a subclass of GenericImage and maintains its own interpolators.
namespace ImplicitSurfaceUtils {


// =============================================================================
// Typedefs
// =============================================================================

typedef GenericImage<RegisteredImage::VoxelType>              DistanceImage;
typedef GenericLinearInterpolateImageFunction<DistanceImage>  DistanceFunction;
typedef GenericFastLinearImageGradientFunction<DistanceImage> DistanceGradient;

// =============================================================================
// Contouring
// =============================================================================

// -----------------------------------------------------------------------------
/// Extract isosurface from distance image
///
/// @param[in] dmap      Distance image.
/// @param[in] offset    Isovalue.
/// @param[in] blurring  Standard deviation of Gaussian kernel with which to
///                      blur the distance image before resampling/contouring.
/// @param[in] isotropic Whether to resample distance image to isotropic voxel size.
/// @param[in] close     Whether to close surface at boundary of image domain.
/// @param[in] normals   Whether to compute normal vectors at mesh points.
/// @param[in] gradients Whether to compute gradient vectors at mesh points.
///
/// @returns Discrete isosurface mesh.
vtkSmartPointer<vtkPolyData> Isosurface(const DistanceImage &dmap, double offset = .0,
                                        double blurring = .0, bool isotropic = true, bool close = true,
                                        bool normals = true, bool gradients = false);

// =============================================================================
// Implicit surface distance function
// =============================================================================

// -----------------------------------------------------------------------------
/// Get distance value at given world position
inline double Evaluate(const DistanceFunction &dist, const double p[3], double offset = .0)
{
  double d, x = p[0], y = p[1], z = p[2];
  dist.WorldToImage(x, y, z);
  if (dist.IsInside(x, y, z)) {
    d = dist.GetInside(x, y, z);
  } else {
    double xmin, ymin, zmin, xmax, ymax, zmax;
    dist.Inside(xmin, ymin, zmin, xmax, ymax, zmax);
    x = clamp(x, xmin, xmax);
    y = clamp(y, ymin, ymax);
    z = clamp(z, zmin, zmax);
    d = max(0., static_cast<double>(dist.GetInside(x, y, z)));
    dist.ImageToWorld(x, y, z);
    x -= p[0], y -= p[1], z -= p[2];
    d += sqrt(x*x + y*y + z*z);
  }
  return d - offset;
}

// -----------------------------------------------------------------------------
/// Get (normalized) distance gradient at given world position
inline void Evaluate(const DistanceGradient &gradient, const double p[3], double g[3], bool normalize = true)
{
  double x = p[0], y = p[1], z = p[2];
  gradient.WorldToImage(x, y, z);
  if (gradient.IsInside(x, y, z)) {
    DistanceGradient::GradientType v = gradient.GetInside(x, y, z);
    if (normalize) v.Normalize();
    g[0] = static_cast<double>(v._x);
    g[1] = static_cast<double>(v._y);
    g[2] = static_cast<double>(v._z);
  } else {
    g[0] = g[1] = g[2] = .0;
  }
}

// =============================================================================
// Auxiliary functions for finding implicit surface intersections
// =============================================================================

// -----------------------------------------------------------------------------
/// Bracket zero crossing along specified ray
inline bool BracketZeroCrossing(const double p[3], const double e[3],
                                double minh, double maxh,
                                double &a, double &da, double &b, double &db,
                                const DistanceFunction &distance, double offset,
                                double tol = 1e-3)
{
  if (fequal(da, .0, tol)) {
    b = a, db = da;
    return true;
  }
  if (abs(da) < maxh) {
    double x[3], step = max(minh, abs(da));
    for (b = step; b <= maxh; b += step) {
      x[0] = p[0] + b * e[0];
      x[1] = p[1] + b * e[1];
      x[2] = p[2] + b * e[2];
      db = Evaluate(distance, x, offset);
      if (fequal(db, .0, tol) || da * db <= .0) return true; // zero crossing
      a = b, da = db;
      step = max(minh, abs(da));
    }
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Find zero crossing in [a, b] using binary search
inline double BinarySearch(const double p[3], const double e[3],
                           double a, double da, double b, double db,
                           const DistanceFunction &distance, double offset = .0,
                           double tol = 1e-3)
{
  if (fequal(da, .0, tol)) return a;
  if (fequal(db, .0, tol)) return b;
  mirtkAssert(da * db < .0, "[a, b] brackets zero crossing");
  double x[3], d, h = a;
  while (a < b) {
    h = .5 * (a + b);
    x[0] = p[0] + h * e[0];
    x[1] = p[1] + h * e[1];
    x[2] = p[2] + h * e[2];
    d = Evaluate(distance, x, offset);
    if (fequal(d, .0, tol)) break;
    if (da * d < .0) b = h, db = d;
    else             a = h, da = d;
  }
  return h;
}

// -----------------------------------------------------------------------------
/// Find zero crossing in [a, b] using minimum surface distance as adaptive step length
inline double AdaptiveSearch(const double p[3], const double e[3],
                             double a, double da, double b, double db,
                             const DistanceFunction &distance, double offset = .0,
                             double tol = 1e-3)
{
  if (fequal(da, .0, tol)) return a;
  if (fequal(db, .0, tol)) return b;
  mirtkAssert(da * db < .0, "[a, b] brackets zero crossing");
  double x[3], h, d = da, step = abs(da), prev;
  for (h = a + step; h <= b; h += step) {
    prev = abs(d);
    x[0] = p[0] + h * e[0];
    x[1] = p[1] + h * e[1];
    x[2] = p[2] + h * e[2];
    d = Evaluate(distance, x, offset);
    if (prev - abs(d) < tol) return h;
    step = copysign(abs(d), da * d);
  }
  return b;
}

// =============================================================================
// Implicit surface intersections
// =============================================================================

// -----------------------------------------------------------------------------
/// Find closest intersection of line segment with implicit surface
///
/// \returns Distance of intersection point from point \p p. If the line segment
///          does not intersect the surface, \p maxh is returned.
inline double IntersectWithLine(const double p[3], const double e[3],
                                double d0, double minh, double maxh,
                                const DistanceFunction &distance, double offset = .0,
                                double tol = 1e-3)
{
  if (fequal(d0, .0, tol)) return .0;
  double a = .0, da = d0, b, db;
  if (!BracketZeroCrossing(p, e, minh, maxh, a, da, b, db, distance, offset, tol)) {
    return maxh;
  }
#if 1
  return BinarySearch(p, e, a, da, b, db, distance, offset, tol);
#else
  return AdaptiveSearch(p, e, a, da, b, db, distance, offset, tol);
#endif
}

// -----------------------------------------------------------------------------
/// Determine distances to all isosurface intersections with a given line segment
inline int Intersections(Array<double> &distances, const double p[3], const double e1[3],
                         double d0, double minh, double l,
                         const DistanceFunction &distance, double offset = .0,
                         double tol = 1e-3)
{
  double x[3] = {p[0], p[1], p[2]}, d = .0, h, maxh = l;
  distances.clear();
  if (fequal(d0, .0, tol)) {
    distances.push_back(.0);
    while (fequal(d0, .0, tol)) {
      d    += minh;
      x[0] += minh * e1[0];
      x[1] += minh * e1[1];
      x[2] += minh * e1[2];
      maxh -= minh;
      if (maxh <= minh) break;
      d0 = Evaluate(distance, x, offset);
    }
  }
  while (maxh > minh) {
    d += (h = IntersectWithLine(x, e1, d0, minh, maxh, distance, offset, tol));
    if (h == maxh) break;
    distances.push_back(d);
    if (h != .0) {
      x[0] += h * e1[0];
      x[1] += h * e1[1];
      x[2] += h * e1[2];
      maxh -= h;
      d0 = Evaluate(distance, x, offset);
    }
    while (fequal(d0, .0, tol)) {
      d    += minh;
      x[0] += minh * e1[0];
      x[1] += minh * e1[1];
      x[2] += minh * e1[2];
      maxh -= minh;
      if (maxh <= minh) break;
      d0 = Evaluate(distance, x, offset);
    }
  }
  return static_cast<int>(distances.size());
}

// -----------------------------------------------------------------------------
/// Determine distances to all isosurface intersections with a given line segment
inline int Intersections(Array<double> &d, const double p[3], const double e[3],
                         double minh, double l,
                         const DistanceFunction &distance, double offset = .0,
                         double tol = 1e-3)
{
  const double d0 = Evaluate(distance, p, offset);
  return Intersections(d, p, e, d0, minh, l, distance, offset, tol);
}

// =============================================================================
// Directional surface distance (vs. minimum distance)
// =============================================================================

// -----------------------------------------------------------------------------
/// Determine surface distance at a point in a specified direction
///
/// \returns Surface distance in specified direction, clamped to \p maxd.
inline double ForwardDistance(const double p[3], const double e1[3],
                              double mind, double minh, double maxd,
                              const DistanceFunction &distance, double offset = .0,
                              double tol = 1e-3)
{
  return IntersectWithLine(p, e1, mind, minh, maxd, distance, offset, tol);
}

// -----------------------------------------------------------------------------
/// Determine surface distance at a point in a specified direction
///
/// \returns Surface distance in specified direction, clamped to \p maxd.
inline double ForwardDistance(const double p[3], const double e1[3],
                              double minh, double maxd,
                              const DistanceFunction &distance, double offset = .0,
                              double tol = 1e-3)
{
  const double mind = Evaluate(distance, p, offset);
  return ForwardDistance(p, e1, mind, minh, maxd, distance, offset, tol);
}

// -----------------------------------------------------------------------------
/// Determine surface distance at a point in a specified direction
///
/// \returns Surface distance in specified direction, clamped to \p maxd.
inline double BackwardDistance(const double p[3], const double e1[3],
                               double mind, double minh, double maxd,
                               const DistanceFunction &distance, double offset = .0,
                               double tol = 1e-3)
{
  const double e2[3] = {-e1[0], -e1[1], -e1[2]};
  return IntersectWithLine(p, e2, mind, minh, maxd, distance, offset, tol);
}

// -----------------------------------------------------------------------------
/// Determine surface distance at a point in a specified direction
///
/// \returns Surface distance in specified direction, clamped to \p maxd.
inline double BackwardDistance(const double p[3], const double e1[3],
                               double minh, double maxd,
                               const DistanceFunction &distance, double offset = .0,
                               double tol = 1e-3)
{
  const double mind = Evaluate(distance, p, offset);
  return BackwardDistance(p, e1, mind, minh, maxd, distance, offset, tol);
}

// -----------------------------------------------------------------------------
/// Determine surface distance at a point in a specified direction
///
/// \returns Surface distance in specified direction, clamped to \p maxd.
inline double Distance(const double p[3], const double e1[3],
                       double mind, double minh, double maxd,
                       const DistanceFunction &distance, double offset = .0,
                       double tol = 1e-3)
{
  double d1 = ForwardDistance (p, e1, mind, minh, maxd, distance, offset, tol);
  double d2 = BackwardDistance(p, e1, mind, minh, maxd, distance, offset, tol);
  return min(d1, d2);
}

// -----------------------------------------------------------------------------
/// Determine surface distance at a point in a specified direction
///
/// \returns Surface distance in specified direction, clamped to \p maxd.
inline double Distance(const double p[3], const double e1[3],
                       double minh, double maxd,
                       const DistanceFunction &distance, double offset = .0,
                       double tol = 1e-3)
{
  const double mind = Evaluate(distance, p, offset);
  return Distance(p, e1, mind, minh, maxd, distance, offset, tol);
}

// -----------------------------------------------------------------------------
/// Determine signed surface distance at a point in normal direction
///
/// \returns Signed surface distance in normal direction, clamped to \p maxd.
inline double SignedDistance(const double p[3], const double n[3],
                             double mind, double minh, double maxd,
                             const DistanceFunction &distance, double offset = .0,
                             double tol = 1e-3)
{
  if (mind < .0) return -ForwardDistance(p, n, mind, minh, maxd, distance, offset, tol);
  else           return BackwardDistance(p, n, mind, minh, maxd, distance, offset, tol);
}

// -----------------------------------------------------------------------------
/// Determine signed surface distance at a point in normal direction
///
/// \returns Signed surface distance in normal direction, clamped to \p maxw.
inline double SignedDistance(const double p[3], const double n[3],
                             double minh, double maxd,
                             const DistanceFunction &distance, double offset = .0,
                             double tol = 1e-3)
{
  const double mind = Evaluate(distance, p, offset);
  return SignedDistance(p, n, mind, minh, maxd, distance, offset, tol);
}

// -----------------------------------------------------------------------------
/// Determine width of implicit surface at a point in a specified direction
///
/// The width of the implicit surface is the sum of the distances to the nearest
/// surface points along a ray cast in opposing directions. The starting point
/// can be either in a convex region outside the implicit surface or a concave
/// region inside of it.
///
/// \returns Width of implicit surface in given direction, clamped to \p maxw.
///          If either ray does not intersect the implicit surface, \p maxw is returned.
inline double Width(const double p[3], const double e1[3],
                    double mind, double minh, double maxw,
                    const DistanceFunction &distance, double offset = .0,
                    double tol = 1e-3)
{
  double d1 = ForwardDistance(p, e1, mind, minh, maxw, distance, offset, tol);
  if (d1 >= maxw) return maxw;
  double d2 = BackwardDistance(p, e1, mind, minh, maxw, distance, offset, tol);
  return min(d1 + d2, maxw);
}

// -----------------------------------------------------------------------------
/// Determine width of implicit surface at a point in a specified direction
///
/// The width of the implicit surface is the sum of the distances to the nearest
/// surface points along a ray cast in opposing directions. The starting point
/// can be either in a convex region outside the implicit surface or a concave
/// region inside of it.
///
/// \returns Width of implicit surface in given direction, clamped to \p maxw.
///          If either ray does not intersect the implicit surface, \p maxw is returned.
inline double Width(const double p[3], const double e1[3],
                    double minh, double maxw,
                    const DistanceFunction &distance, double offset = .0,
                    double tol = 1e-3)
{
  const double mind = Evaluate(distance, p, offset);
  return Width(p, e1, mind, minh, maxw, distance, offset, tol);
}

// =============================================================================
// Miscellaneous distance measurements (e.g., width of cortical sulcus)
// =============================================================================

// -----------------------------------------------------------------------------
/// Base class of distance measurement functors
class DistanceMeasurement
{
  /// Distance measurements
  mirtkAttributeMacro(Array<double>, Values);

protected:

  /// Construct distance measurement for given number of values
  DistanceMeasurement(int nvalues = 1)
  :
    _Values(nvalues, numeric_limits<double>::quiet_NaN())
  {}

  /// Destructor
  virtual ~DistanceMeasurement() {}

public:

  /// Get number of distance measurements
  int NumberOfValues() const
  {
    return static_cast<int>(_Values.size());
  }

  /// Set distance value(s)
  void Set(double v)
  {
    for (size_t i = 0; i < _Values.size(); ++i) {
      _Values[i] = v;
    }
  }

  /// Set distance value(s)
  void Set(int i, double v)
  {
    _Values[i] = v;
  }

  /// Get i-th distance measurement
  double Get(int i = 0) const
  {
    return _Values[i];
  }

  /// Determine distance value(s) in given directions
  ///
  /// \param[in] p        Center point at which to evaluate width of gap.
  /// \param[in] dirs     Point samples on unit sphere centered at the origin.
  /// \param[in] mind     Minimum implicit surface distance at point \p p.
  /// \param[in] minh     Minimum step length for bracketing of intersections.
  /// \param[in] maxd     Maximum distance considered for ray casting.
  /// \param[in] distance Implicit surface distance function.
  /// \param[in] offset   Isovalue of implicit surface.
  /// \param[in] tol      Tolerance by which distance value may differ from \p offset.
  virtual void Evaluate(const double p[3], const PointSamples &dirs,
                        double mind, double minh, double maxd,
                        const DistanceFunction &distance, double offset = .0,
                        double tol = 1e-3) = 0;

  /// Determine distance value(s) in given directions
  ///
  /// \param[in] p        Center point at which to evaluate width of gap.
  /// \param[in] dirs     Point samples on unit sphere centered at the origin.
  /// \param[in] minh     Minimum step length for bracketing of intersections.
  /// \param[in] maxd     Maximum distance considered for ray casting.
  /// \param[in] distance Implicit surface distance function.
  /// \param[in] offset   Isovalue of implicit surface.
  /// \param[in] tol      Tolerance by which distance value may differ from \p offset.
  virtual void Evaluate(const double p[3], const PointSamples &dirs,
                        double minh, double maxd,
                        const DistanceFunction &distance, double offset = .0,
                        double tol = 1e-3)
  {
    const double mind = ImplicitSurfaceUtils::Evaluate(distance, p, offset);
    if (abs(mind) < tol) {
      Set(.0);   // point lies on isosurface
    } else if (abs(mind) >= maxd) {
      Set(maxd); // point is too far from isosurface
    } else {
      this->Evaluate(p, dirs, mind, minh, maxd, distance, offset, tol);
    }
  }
};

// -----------------------------------------------------------------------------
/// Determine minimum of the widths measured in multiple directions
struct MinWidth : public DistanceMeasurement
{
  using DistanceMeasurement::Evaluate;
  void Evaluate(const double p[3], const PointSamples &dirs,
                double mind, double minh, double maxw,
                const DistanceFunction &distance, double offset = .0,
                double tol = 1e-3)
  {
    double e[3], w = maxw;
    for (int i = 0; i < dirs.Size(); ++i) {
      dirs.GetPoint(i, e);
      w = min(w, Width(p, e, mind, minh, maxw, distance, offset, tol));
    }
    Set(0, w);
  }
};

// -----------------------------------------------------------------------------
/// Determine maximum of the widths measured in multiple directions
struct MaxWidth : public DistanceMeasurement
{
  using DistanceMeasurement::Evaluate;
  void Evaluate(const double p[3], const PointSamples &dirs,
                double mind, double minh, double maxw,
                const DistanceFunction &distance, double offset = .0,
                double tol = 1e-3)
  {
    double e[3], w = .0;
    for (int i = 0; i < dirs.Size(); ++i) {
      dirs.GetPoint(i, e);
      w = max(w, Width(p, e, mind, minh, maxw, distance, offset, tol));
    }
    Set(0, w);
  }
};

// -----------------------------------------------------------------------------
/// Determine maximum of the widths measured in multiple directions
struct WidthExtrema : public DistanceMeasurement
{
  WidthExtrema() : DistanceMeasurement(2) {}

  using DistanceMeasurement::Evaluate;
  void Evaluate(const double p[3], const PointSamples &dirs,
                double mind, double minh, double maxw,
                const DistanceFunction &distance, double offset = .0,
                double tol = 1e-3)
  {
    double e[3], w, wmin = maxw, wmax = .0;
    for (int i = 0; i < dirs.Size(); ++i) {
      dirs.GetPoint(i, e);
      w = Width(p, e, mind, minh, maxw, distance, offset, tol);
      if (w < wmin) wmin = w;
      if (w > wmax) wmax = w;
    }
    Set(0, wmin);
    Set(1, wmax);
  }

  double Min() const
  {
    return Get(0);
  }

  double Max() const
  {
    return Get(1);
  }
};

/// Alternative name for WidthExtrema
typedef WidthExtrema MinMaxWidth;

// -----------------------------------------------------------------------------
/// Determine average of the widths measured in multiple directions
struct MeanWidth : public DistanceMeasurement
{
  using DistanceMeasurement::Evaluate;
  void Evaluate(const double p[3], const PointSamples &dirs,
                double mind, double minh, double maxw,
                const DistanceFunction &distance, double offset = .0,
                double tol = 1e-3)
  {
    double e[3], w = .0;
    for (int i = 0; i < dirs.Size(); ++i) {
      dirs.GetPoint(i, e);
      w += Width(p, e, mind, minh, maxw, distance, offset, tol);
    }
    Set(0, w / dirs.Size());
  }
};

// -----------------------------------------------------------------------------
/// Determine average of the widths measured in multiple directions
struct MedianWidth : public DistanceMeasurement
{
  using DistanceMeasurement::Evaluate;
  void Evaluate(const double p[3], const PointSamples &dirs,
                double mind, double minh, double maxw,
                const DistanceFunction &distance, double offset = .0,
                double tol = 1e-3)
  {
    double e[3];
    Array<double> w(dirs.Size());
    for (int i = 0; i < dirs.Size(); ++i) {
      dirs.GetPoint(i, e);
      w[i] = Width(p, e, mind, minh, maxw, distance, offset, tol);
    }
    sort(w.begin(), w.end());
    Set(0, w[w.size() / 2]);
  }
};

// =============================================================================
// Auxiliary functions to evaluate distance measurements
//
// TODO: Consider to make these member functions of the base class DistanceMeasurement.
// =============================================================================

// -----------------------------------------------------------------------------
/// Obtain value of distance measurement evaluated in multiple directions
///
/// \param[in] p        Center point at which to evaluate width of gap.
/// \param[in] dirs     Point samples on unit sphere centered at the origin.
/// \param[in] mind     Minimum implicit surface distance at point \p p.
/// \param[in] minh     Minimum step length for bracketing of intersections.
/// \param[in] maxh     Maximum distance considered for ray casting.
/// \param[in] distance Implicit surface distance function.
/// \param[in] offset   Isovalue of implicit surface.
/// \param[in] tol      Tolerance by which distance value may differ from \p offset.
///
/// \returns Measured distance value.
///
/// \see MinWidth, MaxWidth, MeanWidth, MedianWidth
template <class DistanceMeasurement>
inline double Evaluate(const double p[3], const PointSamples &dirs,
                       double mind, double minh, double maxh,
                       const DistanceFunction &distance, double offset = .0,
                       double tol = 1e-3)
{
  DistanceMeasurement d;
  d.Evaluate(p, dirs, mind, minh, maxh, distance, offset, tol);
  return d.Get(0);
}

// -----------------------------------------------------------------------------
/// Obtain value of distance measurement evaluated in multiple directions
///
/// \param[in] p        Center point at which to evaluate width of gap.
/// \param[in] dirs     Point samples on unit sphere centered at the origin.
/// \param[in] minh     Minimum step length for bracketing of intersections.
/// \param[in] maxh     Maximum distance considered for ray casting.
/// \param[in] distance Implicit surface distance function.
/// \param[in] offset   Isovalue of implicit surface.
/// \param[in] tol      Tolerance by which distance value may differ from \p offset.
///
/// \returns Measured distance value clamped to \p maxw.
///          When the point \p p lies on the isosurface, zero is returned.
///
/// \see MinWidth, MaxWidth, MeanWidth, MedianWidth
template <class DistanceMeasurement>
inline double Evaluate(const double p[3], const PointSamples &dirs,
                       double minh, double maxh,
                       const DistanceFunction &distance, double offset = .0,
                       double tol = 1e-3)
{
  DistanceMeasurement d;
  d.Evaluate(p, dirs, minh, maxh, distance, offset, tol);
  return d.Get(0);
}

// -----------------------------------------------------------------------------
/// Obtain value(s) of distance measurement evaluated in spherical directions
///
/// \param[out] d        Distance measurement(s).
/// \param[in]  p        Center point at which to evaluate width of gap.
/// \param[in]  mind     Minimum implicit surface distance at point \p p.
/// \param[in]  minh     Minimum step length for bracketing of intersections.
/// \param[in]  maxh     Maximum distance considered for ray casting.
/// \param[in]  distance Implicit surface distance function.
/// \param[in]  offset   Isovalue of implicit surface.
/// \param[in]  tol      Tolerance by which distance value may differ from \p offset.
/// \param[in]  ndirs    Desired number of spherical direction samples.
///
/// \see MinWidth, MaxWidth, MeanWidth, MedianWidth
///
/// \note When evaluating the distance measure at multiple points, it is
///       more efficient to precompute the spherical direction samples using
///       the PointSamples::SampleRegularHalfSphere function and to call
///       DistanceMeasurement::Evaluate directly with these as argument.
template <class DistanceMeasurement>
inline void EvaluateSpherical(DistanceMeasurement &d, const double p[3],
                              double mind, double minh, double maxh,
                              const DistanceFunction &distance, double offset = .0,
                              double tol = 1e-3, int ndirs = 20)
{
  PointSamples dirs(ndirs);
  dirs.SampleRegularHalfSphere();
  d.Evaluate(p, dirs, mind, minh, maxh, distance, offset, tol);
}

// -----------------------------------------------------------------------------
/// Obtain value of distance measurement evaluated in spherical directions
///
/// \param[in] p        Center point at which to evaluate width of gap.
/// \param[in] n        Face normal vector.
/// \param[in] mind     Minimum implicit surface distance at point \p p.
/// \param[in] minh     Minimum step length for bracketing of intersections.
/// \param[in] maxh     Maximum distance considered for ray casting.
/// \param[in] distance Implicit surface distance function.
/// \param[in] offset   Isovalue of implicit surface.
/// \param[in] tol      Tolerance by which distance value may differ from \p offset.
/// \param[in] ndirs    Desired number of spherical direction samples.
///
/// \returns Measured distance value in tangent plane.
///
/// \see MinWidth, MaxWidth, MeanWidth, MedianWidth
///
/// \note When evaluating the distance measure at multiple points, it is
///       more efficient to precompute the spherical direction samples using
///       the PointSamples::SampleRegularHalfSphere function and to call
///       DistanceMeasurement::Evaluate directly with these as argument.
template <class DistanceMeasurement>
inline double EvaluateSpherical(const double p[3], const double n[3],
                                double mind, double minh, double maxh,
                                const DistanceFunction &distance, double offset = .0,
                                double tol = 1e-3, int ndirs = 20)
{
  DistanceMeasurement d;
  EvaluateSpherical(d, p, n, mind, minh, maxh, distance, offset, tol, ndirs);
  return d.Get(0);
}

// -----------------------------------------------------------------------------
/// Obtain value(s) of distance measurement evaluated in tangent directions
///
/// \param[out] d        Distance measurement(s).
/// \param[in]  p        Center point at which to evaluate width of gap.
/// \param[in]  n        Face normal vector.
/// \param[in]  mind     Minimum implicit surface distance at point \p p.
/// \param[in]  minh     Minimum step length for bracketing of intersections.
/// \param[in]  maxh     Maximum distance considered for ray casting.
/// \param[in]  distance Implicit surface distance function.
/// \param[in]  offset   Isovalue of implicit surface.
/// \param[in]  tol      Tolerance by which distance value may differ from \p offset.
/// \param[in]  ndirs    Number of direction samples in tangent plane.
/// \param[in]  beta     Maximum angular deviation of normal vector.
///
/// \see MinWidth, MaxWidth, MeanWidth, MedianWidth
template <class DistanceMeasurement>
inline void EvaluateTangential(DistanceMeasurement &d,
                               const double p[3], const double n[3],
                               double mind, double minh, double maxh,
                               const DistanceFunction &distance, double offset = .0,
                               double tol = 1e-3, int ndirs = 4, double beta = .0)
{
  double e1[3], e2[3];
  if (!ComputeTangents(n, e1, e2)) {
    d.Set(numeric_limits<double>::quiet_NaN());
    return;
  }
  const int nbetas = (beta == .0 ? 1 : 3);
  PointSamples dirs(ndirs * nbetas);
  switch (ndirs) {
    case 1: {
      dirs(0) = Point(e1);
    } break;
    case 2: {
      dirs(0) = Point(e1);
      dirs(1) = Point(e2);
    } break;
    case 4: {
      const double sqrt2 = sqrt(2.0);
      dirs(0) = Point(e1);
      dirs(1) = Point(e2);
      dirs(2) = (Point(e1) + Point(e2)) / sqrt2;
      dirs(3) = (Point(e1) - Point(e2)) / sqrt2;
    } break;
    default:
      double alpha = .0, delta = pi / ndirs;
      for (int i = 0; i < ndirs; ++i, alpha += delta) {
        dirs(i) = Point(e1) * cos(alpha) + Point(e2) * sin(alpha);
      }
      break;
  }
  if (nbetas > 1) {
    beta *= rad_per_deg;
    Point sin_beta_n(n);
    sin_beta_n *= sin(beta);
    const double cos_beta = cos(beta);
    for (int i = 0, j = ndirs; i < ndirs; ++i) {
      dirs(j++) = dirs(i) * cos_beta + sin_beta_n;
      dirs(j++) = dirs(i) * cos_beta - sin_beta_n;
    }
  }
  d.Evaluate(p, dirs, mind, minh, maxh, distance, offset, tol);
}

// -----------------------------------------------------------------------------
/// Obtain value of distance measurement evaluated in tangent directions
///
/// \param[in] p        Center point at which to evaluate width of gap.
/// \param[in] n        Face normal vector.
/// \param[in] mind     Minimum implicit surface distance at point \p p.
/// \param[in] minh     Minimum step length for bracketing of intersections.
/// \param[in] maxh     Maximum distance considered for ray casting.
/// \param[in] distance Implicit surface distance function.
/// \param[in] offset   Isovalue of implicit surface.
/// \param[in] tol      Tolerance by which distance value may differ from \p offset.
/// \param[in] ndirs    Number of direction samples in tangent plane.
/// \param[in] beta     Maximum angular deviation of normal vector.
///
/// \returns Measured distance value in tangent plane.
///
/// \see MinWidth, MaxWidth, MeanWidth, MedianWidth
template <class DistanceMeasurement>
inline double EvaluateTangential(const double p[3], const double n[3],
                                 double mind, double minh, double maxh,
                                 const DistanceFunction &distance, double offset = .0,
                                 double tol = 1e-3, int ndirs = 4, double beta = .0)
{
  DistanceMeasurement d;
  EvaluateTangential(d, p, n, mind, minh, maxh, distance, offset, tol, ndirs, beta);
  return d.Get(0);
}

// -----------------------------------------------------------------------------
/// Obtain value of distance measurement evaluated in tangent directions
///
/// \param[out] d        Distance measurement(s).
/// \param[in]  p        Center point at which to evaluate width of gap.
/// \param[in]  n        Face normal vector.
/// \param[in]  minh     Minimum step length for bracketing of intersections.
/// \param[in]  maxh     Maximum distance considered for ray casting.
/// \param[in]  distance Implicit surface distance function.
/// \param[in]  offset   Isovalue of implicit surface.
/// \param[in]  tol      Tolerance by which distance value may differ from \p offset.
/// \param[in]  ndirs    Number of direction samples in tangent plane.
/// \param[in]  beta     Maximum angular deviation of normal vector.
///
/// \see MinWidth, MaxWidth, MeanWidth, MedianWidth
template <class DistanceMeasurement>
inline void EvaluateTangential(DistanceMeasurement &d,
                               const double p[3], const double n[3],
                               double minh, double maxh,
                               const DistanceFunction &distance, double offset = .0,
                               double tol = 1e-3, int ndirs = 4, double beta = .0)
{
  const double mind = Evaluate(distance, p, offset);
  if (abs(mind) >= maxh) {
    d.Set(maxh); // point is too far from isosurface
  // Note: abs(mind) < tol is not used in this case because the tangent vectors
  //       can be orthogonal to the tangent plane of the isosurface. In this case,
  //       we still might want to measure some distance in the direction of the
  //       outwards (or inwards) pointing isosurface normal.
  } else if (Distance(p, n, mind, minh, 2.0 * tol + tol, distance, offset, .5 * tol) <= tol) {
    d.Set(.0); // point is either very close to or lies on the isosurface
  } else {
    EvaluateTangential(d, p, n, mind, minh, maxh, distance, offset, tol, ndirs, beta);
  }
}

// -----------------------------------------------------------------------------
/// Obtain value of distance measurement evaluated in tangent directions
///
/// \param[in] p        Center point at which to evaluate width of gap.
/// \param[in] n        Face normal vector.
/// \param[in] minh     Minimum step length for bracketing of intersections.
/// \param[in] maxh     Maximum distance considered for ray casting.
/// \param[in] distance Implicit surface distance function.
/// \param[in] offset   Isovalue of implicit surface.
/// \param[in] tol      Tolerance by which distance value may differ from \p offset.
/// \param[in] ndirs    Number of direction samples in tangent plane.
/// \param[in] beta     Maximum angular deviation of normal vector.
///
/// \returns Measured distance value clamped to \p maxh.
///          When the point \p p lies on the isosurface, zero is returned.
///
/// \see MinWidth, MaxWidth, MeanWidth, MedianWidth
template <class DistanceMeasurement>
inline double EvaluateTangential(const double p[3], const double n[3],
                                 double minh, double maxh,
                                 const DistanceFunction &distance, double offset = .0,
                                 double tol = 1e-3, int ndirs = 4, double beta = .0)
{
  DistanceMeasurement d;
  EvaluateTangential(d, p, n, minh, maxh, distance, offset, tol, ndirs, beta);
  return d.Get(0);
}

// -----------------------------------------------------------------------------
/// Obtain for all voxels a distance measurement in spherical directions
///
/// What distance measurement is obtained depends on the type of the template
/// function argument, such as for example the MinWidth.
///
/// \returns Image with measured distance values for each voxel in the distance image.
///
/// \see MinWidth, MaxWidth, MeanWidth, MedianWidth
template <class DistanceMeasurement>
DistanceImage Evaluate(double maxh, const DistanceImage &dmap,
                       double offset = .0, double tol = 1e-3, int ndirs = 20)
{
  const double ds   = min(min(dmap.GetXSize(), dmap.GetYSize()), dmap.GetZSize());
  const double minh = .1 * ds;

  PointSamples dirs(ndirs);
  dirs.SampleRegularHalfSphere();

  DistanceMeasurement d;
  DistanceImage output(dmap.Attributes(), d.NumberOfValues());

  DistanceFunction distance;
  distance.Input(&dmap);
  distance.Initialize();

  double p[3];
  for (int k = 0; k < output.Z(); ++k)
  for (int j = 0; j < output.Y(); ++j)
  for (int i = 0; i < output.X(); ++i) {
    p[0] = i, p[1] = j, p[2] = k;
    output.ImageToWorld(p[0], p[1], p[2]);
    d.Evaluate(p, dirs, minh, maxh, distance, offset, tol);
    for (int l = 0; l < output.T(); ++l) {
      output(i, j, k, l) = static_cast<DistanceImage::VoxelType>(d.Get(l));
    }
  }

  return output;
}


} } // namespace mirtk::ImplicitSurfaceUtils

#endif // MIRTK_ImplicitSurfaceUtils_H

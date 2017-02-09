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

#ifndef MIRTK_LinearFreeFormTransformationTD_H
#define MIRTK_LinearFreeFormTransformationTD_H

#include "mirtk/LinearFreeFormTransformation4D.h"

#include <cmath>


namespace mirtk {


// Forward declaration of other TD FFD
class BSplineFreeFormTransformationTD;


/**
 * Temporal diffeomorphic free-form transformation.
 *
 * This class implements a free-form transformation which is represented
 * by a time-varying velocity field (3D+t). The 3D displacement field at a
 * specific time is obtained by integrating the velocity field starting at the
 * reference time point. The integration steps are adjusted if necessary in
 * order to ensure that the resulting spatial transformation is diffeomorphic.
 *
 * For more details about the implementation see De Craene et al. (2012).
 * Temporal diffeomorphic free-form deformation: application to motion and
 * strain estimation from 3D echocardiography.
 * Medical image analysis, 16(2), 427, 2012. doi:10.1016/j.media.2011.10.006
 */
class LinearFreeFormTransformationTD : public LinearFreeFormTransformation4D
{
  mirtkTransformationMacro(LinearFreeFormTransformationTD);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Minimum length of temporal steps (in ms)
  mirtkPublicAttributeMacro(double, MinTimeStep);

  /// Maximum length of temporal steps (in ms)
  mirtkPublicAttributeMacro(double, MaxTimeStep);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  LinearFreeFormTransformationTD();

  /// Construct free-form transformation for given domain and lattice spacing
  explicit LinearFreeFormTransformationTD(const ImageAttributes &,
                                          double, double, double, double);

  /// Construct free-form transformation for given target image and lattice spacing
  explicit LinearFreeFormTransformationTD(const BaseImage &,
                                          double, double, double, double);

  /// Constructor
  explicit LinearFreeFormTransformationTD(const BSplineFreeFormTransformationTD &);

  /// Copy constructor
  LinearFreeFormTransformationTD(const LinearFreeFormTransformationTD &);

  /// Destructor
  virtual ~LinearFreeFormTransformationTD();

  // ---------------------------------------------------------------------------
  // Approximation/Interpolation

  using LinearFreeFormTransformation4D::ApproximateAsNew;

  /// Approximate displacements: This function takes a set of 3D displacement fields
  /// and corresponding time points and temporal intervals. Given these inputs,
  /// it finds a time-varying velocity field which approximates these displacements.
  virtual void ApproximateDOFs(const GenericImage<double> * const *,
                               const double *, const double *, int,
                               bool = false, int = 3, int = 8);

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds !new! parameters such that the resulting
  /// transformation approximates the displacements as good as possible.
  virtual void ApproximateDOFs(const double *, const double *, const double *, const double *,
                               const double *, const double *, const double *, int);

  /// Finds gradient of approximation error: This function takes a set of points
  /// and a set of errors. It finds a gradient w.r.t. the transformation parameters
  /// which minimizes the L2 norm of the approximation error and adds it to the
  /// input gradient with the given weight.
  virtual void ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                                       const double *, const double *, const double *, int,
                                       double *, double = 1.0) const;

  /// Approximate displacements: This function takes a set of 3D displacement fields
  /// and corresponding time points and temporal intervals. Given these inputs,
  /// it finds a time-varying velocity field which approximates these displacements.
  /// The displacements are replaced by the residual displacements of the newly
  /// approximated transformation. Returns the approximation error of the resulting FFD.
  virtual double ApproximateAsNew(GenericImage<double> **,
                                  const double *, const double *, int,
                                  bool = false, int = 3, int = 8);

  /// Interpolates displacements: This function takes a set of displacements defined at
  /// the control points and finds a time-varying velocity field such that the
  /// resulting transformation interpolates these displacements.
  virtual void Interpolate(const double *, const double *, const double *);

  // ---------------------------------------------------------------------------
  // Point transformation

  /// Transforms a single point
  virtual void LocalTransform(double &, double &, double &, double, double) const;

  /// Transforms a single point using the inverse transformation
  virtual bool LocalInverse(double &, double &, double &, double, double) const;

  // ---------------------------------------------------------------------------
  // I/O

  // Do not hide methods of base class
  using LinearFreeFormTransformation4D::Print;

  /// Prints the parameters of the transformation
  virtual void Print(ostream &, Indent = 0) const;

  /// Whether this transformation can read a file of specified type (i.e. format)
  virtual bool CanRead(TransformationType) const;

protected:

  /// Reads transformation parameters from a file stream
  virtual Cifstream &ReadDOFs(Cifstream &, TransformationType);

  /// Writes transformation parameters to a file stream
  virtual Cofstream &WriteDOFs(Cofstream &) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
inline void LinearFreeFormTransformationTD::LocalTransform(double &x, double &y, double &z, double t1, double t2) const
{
  double u, v, w, s;
  Matrix sjac(3, 3);

  // Number of integration steps
  int minsteps = abs(iround((t2 - t1) / _MaxTimeStep));
  int maxsteps = abs(iround((t2 - t1) / _MinTimeStep));
  if (minsteps < 2) minsteps = 2;
  if (maxsteps < 2) maxsteps = 2;

  // Initial step size
  int    N  = minsteps;
  double dt = (t2 - t1) / N;
  // Transform point by integrating the velocity over time
  double t = t1;
  for (int n = 0; n < N; n++) {
    u = x;
    v = y;
    w = z;
    this->WorldToLattice(u, v, w);
    s = this->TimeToLattice(t);
    // Check if still inside FFD domain
    if (-2 <= u && u <= _x + 1 &&
        -2 <= v && v <= _y + 1 &&
        -2 <= w && w <= _z + 1 &&
        -2 <= s && s <= _t + 1) {
      // Reduce step size if necessary to ensure invertibility
      while (N < maxsteps) {
        // Compute current factor of Jacobian of deformation w.r.t the spatial variable.
        //
        // See equation (5) of De Craene et al., Medical Image Analysis 16 (2012).
        // The partial derivatives are calculated w.r.t the spatial variables.
        // The goal is to avoid foldings and singularities in the spatial domain
        // while computing the temporal trajectory of a point.
        EvaluateJacobian(sjac, u, v, w, s);
        sjac = sjac * _matW2L(0, 0, 3, 3) * dt;
        sjac(0, 0) += 1;
        sjac(1, 1) += 1;
        sjac(2, 2) += 1;
        // In case of negative Jacobian determinant, reduce step size
        if (sjac.Det3x3() < 0) {
          N *= 2;
          if (N > maxsteps) {
            N  = maxsteps;
            dt = (t2 - t1) / N;
          } else {
            dt /= 2;
          }
          continue;
        }
        break;
      }
      // Integration step
      EvaluateInside(u, v, w, s); // velocity
      x += u * dt;
      y += v * dt;
      z += w * dt;
      t += dt;
    // Outside FFD domain, velocity is zero
    } else {
      break;
    }
  }
}

// -----------------------------------------------------------------------------
inline bool LinearFreeFormTransformationTD::LocalInverse(double &x, double &y, double &z, double t1, double t2) const
{
  this->LocalTransform(x, y, z, t2, t1);
  return true;
}


} // namespace mirtk

#endif // MIRTK_LinearFreeFormTransformationTD_H

/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
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

// TODO: Change velocities to be given in lattice units rather than world units.
//       This reduces the number of conversions from world to lattice coordinates
//       during the numerical integration.
//       -as12312

#ifndef MIRTK_FreeFormTransformationRungeKutta_H
#define MIRTK_FreeFormTransformationRungeKutta_H

#include "mirtk/Math.h"
#include "mirtk/Matrix.h"
#include "mirtk/Vector3D.h"
#include "mirtk/TransformationJacobian.h"


namespace mirtk {


// ============================================================================
// Base class of Runge-Kutta integration methods
// ============================================================================

/**
 * Base class of integration methods used by FFDs which are parameterized by a
 * velocity field to compute the displacements.
 */
class FreeFormTransformationRungeKutta
{
protected:

  // -------------------------------------------------------------------------
  /// Helper for computation of Jacobian w.r.t spatial coordinates
  static void dkdx(Matrix &dk, const Matrix &Dv, const Matrix &dx, double h)
  {
    dk(0, 0) = (Dv(0, 0) * dx(0, 0) + Dv(0, 1) * dx(1, 0) + Dv(0, 2) * dx(2, 0)) * h;
    dk(0, 1) = (Dv(0, 0) * dx(0, 1) + Dv(0, 1) * dx(1, 1) + Dv(0, 2) * dx(2, 1)) * h;
    dk(0, 2) = (Dv(0, 0) * dx(0, 2) + Dv(0, 1) * dx(1, 2) + Dv(0, 2) * dx(2, 2)) * h;
    dk(1, 0) = (Dv(1, 0) * dx(0, 0) + Dv(1, 1) * dx(1, 0) + Dv(1, 2) * dx(2, 0)) * h;
    dk(1, 1) = (Dv(1, 0) * dx(0, 1) + Dv(1, 1) * dx(1, 1) + Dv(1, 2) * dx(2, 1)) * h;
    dk(1, 2) = (Dv(1, 0) * dx(0, 2) + Dv(1, 1) * dx(1, 2) + Dv(1, 2) * dx(2, 2)) * h;
    dk(2, 0) = (Dv(2, 0) * dx(0, 0) + Dv(2, 1) * dx(1, 0) + Dv(2, 2) * dx(2, 0)) * h;
    dk(2, 1) = (Dv(2, 0) * dx(0, 1) + Dv(2, 1) * dx(1, 1) + Dv(2, 2) * dx(2, 1)) * h;
    dk(2, 2) = (Dv(2, 0) * dx(0, 2) + Dv(2, 1) * dx(1, 2) + Dv(2, 2) * dx(2, 2)) * h;
  }

  // -------------------------------------------------------------------------
  /// Helper for computation of Jacobian w.r.t control point
  static void dkdp(Matrix &dk, const Matrix &Dv, const Matrix &dx, const double dv[3], double h)
  {
    dk(0, 0) = (Dv(0, 0) * dx(0, 0) + Dv(0, 1) * dx(1, 0) + Dv(0, 2) * dx(2, 0) + dv[0]) * h;
    dk(0, 1) = (Dv(0, 0) * dx(0, 1) + Dv(0, 1) * dx(1, 1) + Dv(0, 2) * dx(2, 1)        ) * h;
    dk(0, 2) = (Dv(0, 0) * dx(0, 2) + Dv(0, 1) * dx(1, 2) + Dv(0, 2) * dx(2, 2)        ) * h;
    dk(1, 0) = (Dv(1, 0) * dx(0, 0) + Dv(1, 1) * dx(1, 0) + Dv(1, 2) * dx(2, 0)        ) * h;
    dk(1, 1) = (Dv(1, 0) * dx(0, 1) + Dv(1, 1) * dx(1, 1) + Dv(1, 2) * dx(2, 1) + dv[1]) * h;
    dk(1, 2) = (Dv(1, 0) * dx(0, 2) + Dv(1, 1) * dx(1, 2) + Dv(1, 2) * dx(2, 2)        ) * h;
    dk(2, 0) = (Dv(2, 0) * dx(0, 0) + Dv(2, 1) * dx(1, 0) + Dv(2, 2) * dx(2, 0)        ) * h;
    dk(2, 1) = (Dv(2, 0) * dx(0, 1) + Dv(2, 1) * dx(1, 1) + Dv(2, 2) * dx(2, 1)        ) * h;
    dk(2, 2) = (Dv(2, 0) * dx(0, 2) + Dv(2, 1) * dx(1, 2) + Dv(2, 2) * dx(2, 2) + dv[2]) * h;
  }

}; // FreeFormTransformationRungeKutta

// ============================================================================
// Generic implementation of explicit Runge-Kutta
// ============================================================================

/**
 * Explicit Runge-Kutta integration method for FFD parameterized by velocity field.
 * 
 * The first template argument is the type of the FFD which has to implement
 * the following methods:
 *
 * - Evaluate
 * - EvaluateJacobianWorld
 * - EvaluateJacobianDOFs
 *
 * The second template argument is a struct of the respective Butcher tableau
 * of the Runge-Kutta method which has been defined using the macro
 * MIRTK_DEFINE_FFDRK_EXPLICIT.
 *
 * \sa BSplineFreeFormTransformationTD
 */
template <class TFreeFormTransformation, class TButcherTableau>
class FreeFormTransformationExplicitRungeKutta
:
  public FreeFormTransformationRungeKutta
{
public:

  typedef TButcherTableau BT; ///< Short-hand for Butcher tableau template argument

  // -------------------------------------------------------------------------
  static void Transform(const TFreeFormTransformation *v,
                        double &x, double &y, double &z,
                        double t1, double t2, double dt)
  {
    if (t1 == t2) return;

    Vector3D<double> k[BT::s];                        // Intermediate evaluations
    const double     d = copysign(1.0, t2 - t1); // Direction of integration
    double           h = d * abs(dt);            // Initial step size
    int              i, j;                            // Butcher tableau indices
    double           l;                               // Temporal lattice coordinate

    // Integrate from t=t1 to t=t2
    double t = t1;
    while (d * t < d * t2) {
      // Ensure that last step ends at t2
      if (d * (t + h) > d * t2) h = t2 - t;
      // Evaluate velocities at intermediate steps
      if (BT::fsal && t != t1) k[0] = k[BT::s - 1], i = 1;
      else                                          i = 0;
      for (/*i = 0|1*/; i < BT::s; i++) {
        // Lattice coordinates of current intermediate step
        k[i]._x = x, k[i]._y = y, k[i]._z = z;
        for (j = 0; j < i; j++) k[i] += k[j] * BT::a[i][j];
        v->WorldToLattice(k[i]._x, k[i]._y, k[i]._z);
        l = v->TimeToLattice(t + BT::c[i] * h);
        // Evaluate velocity at intermediate point
        v->Evaluate(k[i]._x, k[i]._y, k[i]._z, l);
        k[i] *= h;
      }
      // Perform step
      for (i = 0; i < BT::s; i++) {
        x += k[i]._x * BT::b[i];
        y += k[i]._y * BT::b[i];
        z += k[i]._z * BT::b[i];
      }
      t += h;
    }
  }

  // -------------------------------------------------------------------------
  static void Jacobian(const TFreeFormTransformation *v,
                       Matrix &jac,
                       double &x, double &y, double &z,
                       double t1, double t2, double dt)
  {
    if (t1 == t2) return;

    Vector3D<double> k [BT::s];                       // Intermediate evaluations
    Matrix           dk[BT::s];                       // Derivative of k_i
    Matrix           dx(3, 3);                        // Derivative of intermediate location
    Matrix           Dv(3, 3);                        // Partial derivative of velocity field
    const double     d = copysign(1.0, t2 - t1); // Direction of integration
    double           h = d * abs(dt);            // Initial step size
    int              i, j;                            // Butcher tableau indices
    double           l;                               // Temporal lattice coordinate

    for (i = 0; i < BT::s; ++i) dk[i].Initialize(3, 3);

    // Integrate from t=t1 to t=t2
    double t = t1;
    while (d * t < d * t2) {
      // Ensure that last step ends at t2
      if (d * (t + h) > d * t2) h = t2 - t;
      // Evaluate at intermediate steps
      if (BT::fsal && t != t1) {
        k [0] = k [BT::s-1];
        dk[0] = dk[BT::s-1];
        i     = 1;
      } else {
        i     = 0;
      }
      for (/*i = 0|1*/; i < BT::s; i++) {
        // Lattice coordinates
        k[i]._x = x, k[i]._y = y, k[i]._z = z;
        for (j = 0; j < i; j++) k[i] += k[j] * BT::a[i][j];
        v->WorldToLattice(k[i]._x, k[i]._y, k[i]._z);
        l = v->TimeToLattice(t + BT::c[i] * h);
        // Evaluate partial derivatives of velocity
        v->EvaluateJacobianWorld(Dv, k[i]._x, k[i]._y, k[i]._z, l);
        // Evaluate velocity
        v->Evaluate(k[i]._x, k[i]._y, k[i]._z, l);
        k[i] *= h;
        // Calculate derivatives of k_i
        dx.Ident();
        for (j = 0; j < i; j++) dx += dk[j] * BT::a[i][j];
        dkdx(dk[i], Dv, dx, h); // dk_i = Dv dx h
      }
      // Perform step with local extrapolation
      dx.Ident();
      for (i = 0; i < BT::s; i++) {
        x  +=  k[i]._x * BT::b[i];
        y  +=  k[i]._y * BT::b[i];
        z  +=  k[i]._z * BT::b[i];
        dx += dk[i]    * BT::b[i];
      }
      jac = dx * jac;
      t  += h;
    }
  }

  // -------------------------------------------------------------------------
  static void JacobianDOFs(const TFreeFormTransformation *v,
                           Matrix &jac, int ci, int cj, int ck, int cl,
                           double &x, double &y, double &z,
                           double t1, double t2, double dt)
  {
    if (t1 == t2) return;

    Vector3D<double> k [BT::s];                       // Intermediate evaluations
    Matrix           dk[BT::s];                       // Derivative of k_i w.r.t control point
    Matrix           dx;                              // Derivative of intermediate location w.r.t control point
    Matrix           Dv;                              // Partial derivative of velocity field w.r.t spatial coordinates
    double           dv[3];                           // Partial derivative of velocity field w.r.t control point
    const double     d = copysign(1.0, t2 - t1); // Direction of integration
    double           h = d * abs(dt);            // Initial step size
    int              i, j;                            // Butcher tableau indices
    double           l;                               // Temporal lattice coordinate

    for (i = 0; i < BT::s; i++) dk[i].Initialize(3, 3);

    // Integrate from t=t1 to t=t2
    double t = t1;
    while (d * t < d * t2) {
      // Ensure that last step ends at t2
      if (d * (t + h) > d * t2) h = t2 - t;
      // Evaluate at intermediate steps
      if (BT::fsal && t != t1) {
        k [0] = k [BT::s-1];
        dk[0] = dk[BT::s-1];
        i     = 1;
      } else {
        i     = 0;
      }
      for (/*i = 0|1*/; i < BT::s; i++) {
        // Lattice coordinates
        k[i]._x = x, k[i]._y = y, k[i]._z = z;
        for (j = 0; j < i; j++) k[i] += k[j] * BT::a[i][j];
        v->WorldToLattice(k[i]._x, k[i]._y, k[i]._z);
        l = v->TimeToLattice(t + BT::c[i] * h);
        // Evaluate partial derivatives of velocity
        v->EvaluateJacobianWorld(Dv,                 k[i]._x, k[i]._y, k[i]._z, l);
        v->EvaluateJacobianDOFs (dv, ci, cj, ck, cl, k[i]._x, k[i]._y, k[i]._z, l);
        // Evaluate velocity
        v->Evaluate(k[i]._x, k[i]._y, k[i]._z, l);
        k[i] *= h;
        // Jacobian at intermediate step
        dx = jac;
        for (j = 0; j < i; j++) dx += dk[j] * BT::a[i][j];
        // Calculate derivatives of k_i
        dkdp(dk[i], Dv, dx, dv, h); // dk = (Dv * dx + dv) * h
      }
      // Perform step with local extrapolation
      for (i = 0; i < BT::s; i++) {
        x   +=  k[i]._x * BT::b[i];
        y   +=  k[i]._y * BT::b[i];
        z   +=  k[i]._z * BT::b[i];
        jac += dk[i]    * BT::b[i];
      }
      t += h;
    }
  }

  // -------------------------------------------------------------------------
  static void JacobianDOFs(const TFreeFormTransformation *v,
                           TransformationJacobian &jac,
                           double &x, double &y, double &z,
                           double t1, double t2, double dt)
  {
    if (t1 == t2) return;

    Vector3D<double>       k [BT::s];                       // Intermediate evaluations
    TransformationJacobian dk[BT::s];                       // Derivative of k_i w.r.t transformation parameters
                                                            // and temporarily used to store partial derivative of
                                                            // velocity field w.r.t control points
    TransformationJacobian dx;                              // Current derivative w.r.t transformation parameters
    Matrix                 Dv;                              // Partial derivative of velocity field w.r.t spatial coordinates
    const double           d = copysign(1.0, t2 - t1); // Direction of integration
    double                 h = d * abs(dt);            // Initial step size
    int                    i, j;                            // Butcher tableau indices
    double                 l;                               // Temporal lattice coordinate

    // Integrate from t=t1 to t=t2
    double t = t1;
    while (d * t < d * t2) {
      // Ensure that last step ends at t2
      if (d * (t + h) > d * t2) h = t2 - t;
      // Evaluate at intermediate steps
      if (BT::fsal && t != t1) {
        k [0] = k [BT::s - 1];
        dk[0] = dk[BT::s - 1];
        i     = 1;
      } else {
        i     = 0;
      }
      for (/*i = 0|1*/; i < BT::s; i++) {
        // Lattice coordinates
        k[i]._x = x, k[i]._y = y, k[i]._z = z;
        for (j = 0; j < i; j++) k[i] += k[j] * BT::a[i][j];
        v->WorldToLattice(k[i]._x, k[i]._y, k[i]._z);
        l = v->TimeToLattice(t + BT::c[i] * h);
        // Evaluate partial derivatives of velocity
        v->EvaluateJacobianWorld(Dv,           k[i]._x, k[i]._y, k[i]._z, l);
        v->EvaluateJacobianDOFs (dk[i]/*=dv*/, k[i]._x, k[i]._y, k[i]._z, l);
        // Evaluate velocity
        v->Evaluate(k[i]._x, k[i]._y, k[i]._z, l);
        k[i] *= h;
        // Jacobian at intermediate step
        dx = jac;
        for (j = 0; j < i; j++) dx.add(dk[j], BT::a[i][j]);
        // Calculate derivatives of k_i
        dx *= (Dv *= h);  // += Dv * dx * h
        dx.add(dk[i], h); // += dv * h
        dk[i] = dx;       // dk = (Dv * dx + dv) * h
      }
      // Perform step with local extrapolation
      for (i = 0; i < BT::s; i++) {
        x += k[i]._x * BT::b[i];
        y += k[i]._y * BT::b[i];
        z += k[i]._z * BT::b[i];
        jac.add(dk[i], BT::b[i]);
      }
      t += h;
    }
  }

}; // FreeFormTransformationExplicitRungeKutta

// ===========================================================================
// Generic implementation of embedded Runge-Kutta
// ===========================================================================

/**
 * Embedded Runge-Kutta integration method for FFD parameterized by velocity field.
 * 
 * The first template argument is the type of the FFD which has to implement
 * the following methods:
 *
 * - Evaluate
 * - EvaluateJacobianWorld
 * - EvaluateJacobianDOFs
 *
 * The second template argument is a struct of the respective Butcher tableau
 * of the Runge-Kutta method which has been defined using the macro
 * MIRTK_DEFINE_FFDRK_EMBEDDED.
 *
 * \sa BSplineFreeFormTransformationTD
 */
template <class TFreeFormTransformation, class TButcherTableau>
class FreeFormTransformationEmbeddedRungeKutta
:
  public FreeFormTransformationRungeKutta
{
public:

  typedef TButcherTableau BT; ///< Short-hand for Butcher tableau template argument

  // -------------------------------------------------------------------------
  static void Transform(const TFreeFormTransformation *v,
                        double &x, double &y, double &z,
                        double t1, double t2, double mindt, double maxdt, double tol)
  {
    if (t1 == t2) return;

    Vector3D<double> k[BT::s];                             // k_i = h * v(t + c_i * h, x + sum_{j=0}^{i-1} a_j * k_j)
    Vector3D<double> tmp;                                  // Solution of order p-1
    Vector3D<double> next;                                 // Solution of order p
    double           error;                                // Local error estimate
    double           h = t2 - t1;                          // Initial step size
    double           hnext;                                // Next step size
    const double     d = copysign(1.0, h);            // Direction of integration
    const double     e = 1.0 / static_cast<double>(BT::p); // Exponent for step size scaling factor
    int              i, j;                                 // Butcher tableau indices
    double           l;                                    // Temporal lattice coordinate

    // Decrease initial step size if necessary
    if (abs(h) > maxdt) h = copysign(maxdt, d);

    // Integrate from t=t1 to t=t2
    double t = t1;
    while (d * t < d * t2) {
      // Ensure that last step ends at t2
      if (d * (t + h) > d * t2) h = t2 - t;
      // Evaluate velocities at intermediate steps
      if (BT::fsal && t != t1) k[0] = k[BT::s - 1], i = 1;
      else                                          i = 0;
      for (/*i = 0|1*/; i < BT::s; i++) {
        k[i]._x = x, k[i]._y = y, k[i]._z = z;
        for (j = 0; j < i; j++) k[i] += k[j] * BT::a[i][j];
        v->WorldToLattice(k[i]._x, k[i]._y, k[i]._z);
        l = v->TimeToLattice(t + BT::c[i] * h);
        v->Evaluate(k[i]._x, k[i]._y, k[i]._z, l);
        k[i] *= h;
      }
      // Calculate solution of order p
      next._x = x, next._y = y, next._z = z;
      for (i = 0; i < BT::s; i++) next += k[i] * BT::b[1][i];
      // Adapt step size
      if (mindt < maxdt) {
        // Calculate solution of order p-1
        tmp._x = x, tmp._y = y, tmp._z = z;
        for (i = 0; i < BT::s; i++) tmp += k[i] * BT::b[0][i];
        // Estimate local error
        error = max(max(abs(next._x - tmp._x),
                                  abs(next._y - tmp._y)),
                                  abs(next._z - tmp._z));
        // If local error exceeds tolerance...
        if (abs(h) > mindt && error > tol) {
          // ...decrease step size
          h *= 0.8 * pow(tol / error, e);
          if (abs(h) < mindt) h = copysign(mindt, d);
          // ...and redo current step
          continue;
        // Otherwise, increase step size
        } else {
          hnext = 0.8 * pow(tol / error, e) * h;
          if      (abs(hnext) < mindt) hnext = copysign(mindt, d);
          else if (abs(hnext) > maxdt) hnext = copysign(maxdt, d);
        }
      } else {
        hnext = h;
      }
      // Perform step with local extrapolation
      x  = next._x;
      y  = next._y;
      z  = next._z;
      t += h;
      // Update step size
      h = hnext;
    }
  }

  // -------------------------------------------------------------------------
  static void Jacobian(const TFreeFormTransformation *v,
                       Matrix &jac, double &x, double &y, double &z,
                       double t1, double t2, double mindt, double maxdt, double tol)
  {
    if (t1 == t2) return;

    Vector3D<double> k [BT::s];                            // Intermediate evaluations
    Matrix           dk[BT::s];                            // Derivative of k_i
    Matrix           dx(3, 3);                             // Derivative of intermediate location
    Matrix           Dv(3, 3);                             // Partial derivative of velocity field
    double           h = t2 - t1;                          // Initial step size
    double           hnext;                                // Next step size
    const double     d = copysign(1.0, h);            // Direction of integration
    Vector3D<double> tmp;                                  // Solution of order p-1
    Vector3D<double> next;                                 // Solution of order p
    double           error;                                // Local error estimate
    const double     e = 1.0 / static_cast<double>(BT::p); // Exponent for step size scaling factor
    int              i, j;                                 // Butcher tableau indices
    double           l;                                    // Temporal lattice coordinate

    for (i = 0; i < BT::s; i++) dk[i].Initialize(3, 3);

    // Decrease initial step size if necessary
    if (abs(h) > maxdt) h = copysign(maxdt, d);

    // Integrate from t=t1 to t=t2
    double t = t1;
    while (d * t < d * t2) {
      // Ensure that last step ends at t2
      if (d * (t + h) > d * t2) h = t2 - t;
      // Evaluate at intermediate steps
      if (BT::fsal && t != t1) {
        k [0] = k [BT::s-1];
        dk[0] = dk[BT::s-1];
        i     = 1;
      } else {
        i     = 0;
      }
      for (/*i = 0|1*/; i < BT::s; i++) {
        // Lattice coordinates
        k[i]._x = x, k[i]._y = y, k[i]._z = z;
        for (j = 0; j < i; j++) k[i] += k[j] * BT::a[i][j];
        v->WorldToLattice(k[i]._x, k[i]._y, k[i]._z);
        l = v->TimeToLattice(t + BT::c[i] * h);
        // Evaluate partial derivatives of velocity
        v->EvaluateJacobianWorld(Dv, k[i]._x, k[i]._y, k[i]._z, l);
        // Evaluate velocity
        v->Evaluate(k[i]._x, k[i]._y, k[i]._z, l);
        k[i] *= h;
        // Jacobian at current intermediate step
        dx.Ident();
        for (j = 0; j < i; j++) dx += dk[j] * BT::a[i][j];
        // Calculate derivatives of k_i
        dkdx(dk[i], Dv, dx, h); // dk = Dv * dx * h
      }
      // Calculate solution of order p
      next._x = x, next._y = y, next._z = z;
      for (i = 0; i < BT::s; i++) next += k[i] * BT::b[1][i];
      // Adapt step size
      if (mindt < maxdt) {
        // Calculate solution of order p-1
        tmp._x = x, tmp._y = y, tmp._z = z;
        for (i = 0; i < BT::s; i++) tmp += k[i] * BT::b[0][i];
        // Estimate local error
        error = max(max(abs(next._x - tmp._x),
                                  abs(next._y - tmp._y)),
                                  abs(next._z - tmp._z));
        // If local error exceeds tolerance...
        if (abs(h) > mindt && error > tol) {
          // ...decrease step size
          h *= 0.8 * pow(tol / error, e);
          if (abs(h) < mindt) h = copysign(mindt, d);
          // ...and redo current step
          continue;
          // Otherwise, increase step size
        } else {
          hnext = 0.8 * pow(tol / error, e) * h;
          if      (abs(hnext) < mindt) hnext = copysign(mindt, d);
          else if (abs(hnext) > maxdt) hnext = copysign(maxdt, d);
        }
      } else {
        hnext = h;
      }
      // Update Jacobian
      dx.Ident();
      for (i = 0; i < BT::s; i++) dx += dk[i] * BT::b[1][i];
      jac = dx * jac;
      // Perform step with local extrapolation
      x  = next._x;
      y  = next._y;
      z  = next._z;
      t += h;
      // Update step size
      h = hnext;
    }
  }

  // -------------------------------------------------------------------------
  static void JacobianDOFs(const TFreeFormTransformation *v,
                           Matrix &jac, int ci, int cj, int ck, int cl,
                           double &x, double &y, double &z, double t1, double t2,
                           double mindt, double maxdt, double tol)
  {
    if (t1 == t2) return;

    Vector3D<double> k [BT::s];                            // Intermediate evaluations
    Matrix           dk[BT::s];                            // Derivative of k_i w.r.t control point
    Matrix           dx;                                   // Derivative of intermediate location w.r.t control point
    Matrix           Dv;                                   // Partial derivative of velocity field w.r.t spatial coordinates
    double           dv[3];                                // Partial derivative of velocity field w.r.t control point
    double           h = t2 - t1;                          // Initial step size
    double           hnext;                                // Next step size
    const double     d = copysign(1.0, h);            // Direction of integration
    Vector3D<double> tmp;                                  // Solution of order p-1
    Vector3D<double> next;                                 // Solution of order p
    double           error;                                // Local error estimate
    const double     e = 1.0 / static_cast<double>(BT::p); // Exponent for step size scaling factor
    int              i, j;                                 // Butcher tableau indices
    double           l;                                    // Temporal lattice coordinate

    for (i = 0; i < BT::s; i++) dk[i].Initialize(3, 3);

    // Decrease initial step size if necessary
    if (abs(h) > maxdt) h = copysign(maxdt, d);

    // Integrate from t=t1 to t=t2
    double t = t1;
    while (d * t < d * t2) {
      // Ensure that last step ends at t2
      if (d * (t + h) > d * t2) h = t2 - t;
      // Evaluate at intermediate steps
      if (BT::fsal && t != t1) {
        k [0] = k [BT::s-1];
        dk[0] = dk[BT::s-1];
        i     = 1;
      } else {
        i     = 0;
      }
      for (/*i = 0|1*/; i < BT::s; i++) {
        // Lattice coordinates
        k[i]._x = x, k[i]._y = y, k[i]._z = z;
        for (j = 0; j < i; j++) k[i] += k[j] * BT::a[i][j];
        v->WorldToLattice(k[i]._x, k[i]._y, k[i]._z);
        l = v->TimeToLattice(t + BT::c[i] * h);
        // Evaluate partial derivatives of velocity
        v->EvaluateJacobianWorld(Dv,                 k[i]._x, k[i]._y, k[i]._z, l);
        v->EvaluateJacobianDOFs (dv, ci, cj, ck, cl, k[i]._x, k[i]._y, k[i]._z, l);
        // Evaluate velocity
        v->Evaluate(k[i]._x, k[i]._y, k[i]._z, l);
        k[i] *= h;
        // Jacobian at current intermediate step
        dx = jac;
        for (j = 0; j < i; j++) dx += dk[j] * BT::a[i][j];
        // Calculate derivatives of k_i
        dkdp(dk[i], Dv, dx, dv, h); // dk = (Dv * dx + dv) * h
      }
      // Calculate solution of order p
      next._x = x, next._y = y, next._z = z;
      for (i = 0; i < BT::s; i++) next += k[i] * BT::b[1][i];
      // Adapt step size
      if (mindt < maxdt) {
        // Calculate solution of order p-1
        tmp._x = x, tmp._y = y, tmp._z = z;
        for (i = 0; i < BT::s; i++) tmp += k[i] * BT::b[0][i];
        // Estimate local error
        error = max(max(abs(next._x - tmp._x),
                                  abs(next._y - tmp._y)),
                                  abs(next._z - tmp._z));
        // If local error exceeds tolerance...
        if (abs(h) > mindt && error > tol) {
          // ...decrease step size
          h *= 0.8 * pow(tol / error, e);
          if (abs(h) < mindt) h = copysign(mindt, d);
          // ...and redo current step
          continue;
        // Otherwise, increase step size
        } else {
          hnext = 0.8 * pow(tol / error, e) * h;
          if      (abs(hnext) < mindt) hnext = copysign(mindt, d);
          else if (abs(hnext) > maxdt) hnext = copysign(maxdt, d);
        }
      } else {
        hnext = h;
      }
      // Update Jacobian
      for (i = 0; i < BT::s; i++) jac += dk[i] * BT::b[1][i];
      // Perform step with local extrapolation
      x  = next._x;
      y  = next._y;
      z  = next._z;
      t += h;
      // Update step size
      h = hnext;
    }
  }

  // -------------------------------------------------------------------------
  static void JacobianDOFs(const TFreeFormTransformation *v,
                           TransformationJacobian &jac,
                           double &x, double &y, double &z,
                           double t1, double t2, double mindt, double maxdt, double tol)
  {
    if (t1 == t2) return;

    Vector3D<double>       k [BT::s];                  // Intermediate evaluations
    TransformationJacobian dk[BT::s];                  // Derivative of k_i w.r.t transformation parameters
                                                       // and temporarily used to store partial derivative of
                                                       // velocity field w.r.t control points
    TransformationJacobian dx;                         // Current derivative w.r.t transformation parameters
    Matrix                 Dv;                         // Partial derivative of velocity field w.r.t spatial coordinates
    double                 h = t2 - t1;                // Initial step size
    const double           d = copysign(1.0, h);       // Direction of integration
    double                 hnext;                      // Next step size
    Vector3D<double>       tmp;                        // Solution of order p-1
    Vector3D<double>       next;                       // Solution of order p
    double                 error;                      // Local error estimate
    const double           e = 1.0 / double(BT::p);    // Exponent for step size scaling factor
    int                    i, j;                       // Butcher tableau indices
    double                 l;                          // Temporal lattice coordinate

    // Decrease initial step size if necessary
    if (abs(h) > maxdt) h = copysign(maxdt, d);

    // Integrate from t=t1 to t=t2
    double t = t1;
    while (d * t < d * t2) {
      // Ensure that last step ends at t2
      if (d * (t + h) > d * t2) h = t2 - t;
      // Evaluate at intermediate steps
      if (BT::fsal && t != t1) {
        k [0] = k [BT::s-1];
        dk[0] = dk[BT::s-1];
        i     = 1;
      } else {
        i     = 0;
      }
      for (/*i = 0|1*/; i < BT::s; i++) {
        // Lattice coordinates
        k[i]._x = x, k[i]._y = y, k[i]._z = z;
        for (j = 0; j < i; j++) k[i] += k[j] * BT::a[i][j];
        v->WorldToLattice(k[i]._x, k[i]._y, k[i]._z);
        l = v->TimeToLattice(t + BT::c[i] * h);
        // Evaluate partial derivatives of velocity
        v->EvaluateJacobianWorld(Dv,           k[i]._x, k[i]._y, k[i]._z, l);
        v->EvaluateJacobianDOFs (dk[i]/*=dv*/, k[i]._x, k[i]._y, k[i]._z, l);
        // Evaluate velocity
        v->Evaluate(k[i]._x, k[i]._y, k[i]._z, l);
        k[i] *= h;
        // Jacobian at intermediate step
        dx = jac;
        for (j = 0; j < i; j++) dx.add(dk[j], BT::a[i][j]);
        // Calculate derivatives of k_i
        dx *= (Dv *= h);  // += Dv * dx * h
        dx.add(dk[i], h); // += dv * h
        dk[i] = dx;       // dk = (Dv * dx + dv) * h
      }
      // Calculate solution of order p
      next._x = x, next._y = y, next._z = z;
      for (i = 0; i < BT::s; i++) next += k[i] * BT::b[1][i];
      // Adapt step size
      if (mindt < maxdt) {
        // Calculate solution of order p-1
        tmp._x = x, tmp._y = y, tmp._z = z;
        for (i = 0; i < BT::s; i++) tmp += k[i] * BT::b[0][i];
        // Estimate local error
        error = max(max(abs(next._x - tmp._x),
                                  abs(next._y - tmp._y)),
                                  abs(next._z - tmp._z));
        // If local error exceeds tolerance...
        if (abs(h) > mindt && error > tol) {
          // ...decrease step size
          h *= 0.8 * pow(tol / error, e);
          if (abs(h) < mindt) h = copysign(mindt, d);
          // ...and redo current step
          continue;
        // Otherwise, increase step size
        } else {
          hnext = 0.8 * pow(tol / error, e) * h;
          if      (abs(hnext) < mindt) hnext = copysign(mindt, d);
          else if (abs(hnext) > maxdt) hnext = copysign(maxdt, d);
        }
      } else {
        hnext = h;
      }
      // Update Jacobian
      for (i = 0; i < BT::s; i++) jac.add(dk[i], BT::b[1][i]);
      // Perform step with local extrapolation
      x  = next._x;
      y  = next._y;
      z  = next._z;
      t += h;
      // Update step size
      h = hnext;
    }
  }

}; // FreeFormTransformationEmbeddedRungeKutta

// =============================================================================
// Auxiliary Macros
// =============================================================================

// -----------------------------------------------------------------------------
/// Declare explicit Runge-Kutta method
#define MIRTK_DECLARE_FFDRK_EXPLICIT(NAME, S)                                  \
  struct FreeFormTransformationButcherTableau##NAME                            \
  {                                                                            \
    static const int    s = S;                                                 \
    static const int    p;                                                     \
    static const bool   fsal;                                                  \
    static const double a[S][S];                                               \
    static const double b[S];                                                  \
    static const double c[S];                                                  \
  };                                                                           \
                                                                               \
  template <class TFreeFormTransformation>                                     \
  class FreeFormTransformationIntegration##NAME                                \
    : public mirtk::FreeFormTransformationExplicitRungeKutta                   \
             <TFreeFormTransformation, FreeFormTransformationButcherTableau##NAME> \
  { }

// -----------------------------------------------------------------------------
/// Declare embedded Runge-Kutta method
#define MIRTK_DECLARE_FFDRK_EMBEDDED(NAME, S)                                  \
  struct FreeFormTransformationButcherTableau##NAME                            \
  {                                                                            \
    static const int    s = S;                                                 \
    static const int    p;                                                     \
    static const bool   fsal;                                                  \
    static const double a[S][S];                                               \
    static const double b[2][S];                                               \
    static const double c[S];                                                  \
  };                                                                           \
                                                                               \
  template <class TFreeFormTransformation>                                     \
  class FreeFormTransformationIntegration##NAME                                \
    : public mirtk::FreeFormTransformationEmbeddedRungeKutta                   \
             <TFreeFormTransformation, FreeFormTransformationButcherTableau##NAME> \
  { }

// -----------------------------------------------------------------------------
/// Define explicit Runge-Kutta method
#define MIRTK_DEFINE_FFDRK_EXPLICIT(NAME, S, P, FSAL, C, A, B)                 \
  /*const int    FreeFormTransformationButcherTableau##NAME::s       = S;*/    \
  const int    FreeFormTransformationButcherTableau##NAME::p       = P;        \
  const bool   FreeFormTransformationButcherTableau##NAME::fsal    = FSAL;     \
  const double FreeFormTransformationButcherTableau##NAME::a[S][S] = A;        \
  const double FreeFormTransformationButcherTableau##NAME::b[S]    = B;        \
  const double FreeFormTransformationButcherTableau##NAME::c[S]    = C

// -----------------------------------------------------------------------------
/// Define embedded Runge-Kutta method
#define MIRTK_DEFINE_FFDRK_EMBEDDED(NAME, S, P, FSAL, C, A, BY, BZ)            \
  /*const int    FreeFormTransformationButcherTableau##NAME::s       = S;*/    \
  const int    FreeFormTransformationButcherTableau##NAME::p       = P;        \
  const bool   FreeFormTransformationButcherTableau##NAME::fsal    = FSAL;     \
  const double FreeFormTransformationButcherTableau##NAME::a[S][S] = A;        \
  const double FreeFormTransformationButcherTableau##NAME::b[2][S] = {BY, BZ}; \
  const double FreeFormTransformationButcherTableau##NAME::c[S]    = C

// ============================================================================
// Runge-Kutta integration methods
// ============================================================================

//                          name    steps    description
// ----------------------------------------------------------------------------
MIRTK_DECLARE_FFDRK_EXPLICIT(RKE1,   1); ///< Forward Euler method
MIRTK_DECLARE_FFDRK_EMBEDDED(RKEH12, 2); ///< Euler-Heun method of order 2(1)
MIRTK_DECLARE_FFDRK_EXPLICIT(RKE2,   2); ///< Modified Euler method (Heun's method, explicit midpoint rule)
MIRTK_DECLARE_FFDRK_EXPLICIT(RKH2,   2); ///< Improved Euler method (Heun's method, explicit trapezoidal rule)
MIRTK_DECLARE_FFDRK_EMBEDDED(RKBS23, 4); ///< Bogacki-Shampine method of order 3(2)
MIRTK_DECLARE_FFDRK_EXPLICIT(RK4,    4); ///< Classical Runge-Kutta method
MIRTK_DECLARE_FFDRK_EMBEDDED(RKF45,  6); ///< Fehlberg method of order 5(4)
MIRTK_DECLARE_FFDRK_EMBEDDED(RKCK45, 6); ///< Cash-Karp method of order 5(4)
MIRTK_DECLARE_FFDRK_EMBEDDED(RKDP45, 7); ///< Dormand–Prince method of order 5(4)


} // namespace mirtk

#endif // MIRTK_FreeFormTransformationRungeKutta_H

/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2017 Andreas Schuh
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

#include "mirtk/ConjugateGradientDescent.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/ObjectFactory.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterOptimizerMacro(ConjugateGradientDescent);


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
ConjugateGradientDescent::ConjugateGradientDescent(ObjectiveFunction *f)
:
  GradientDescent(f),
  _UseConjugateGradient(true),
  _ConjugateTotalGradient(true),
  _g(nullptr), _h(nullptr)
{
}

// -----------------------------------------------------------------------------
void ConjugateGradientDescent::CopyAttributes(const ConjugateGradientDescent &other)
{
  _UseConjugateGradient   = other._UseConjugateGradient;
  _ConjugateTotalGradient = other._ConjugateTotalGradient;

  Deallocate(_g);
  if (other._g && _Function) {
    Allocate(_g, _Function->NumberOfDOFs());
    memcpy(_g, other._g, _Function->NumberOfDOFs() * sizeof(double));
  }

  Deallocate(_h);
  if (other._h && _Function) {
    Allocate(_h, _Function->NumberOfDOFs());
    memcpy(_h, other._h, _Function->NumberOfDOFs() * sizeof(double));
  }
}

// -----------------------------------------------------------------------------
ConjugateGradientDescent::ConjugateGradientDescent(const ConjugateGradientDescent &other)
:
  GradientDescent(other),
  _g(nullptr), _h(nullptr)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
ConjugateGradientDescent &ConjugateGradientDescent::operator =(const ConjugateGradientDescent &other)
{
  if (this != &other) {
    GradientDescent::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
ConjugateGradientDescent::~ConjugateGradientDescent()
{
  Deallocate(_g);
  Deallocate(_h);
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool ConjugateGradientDescent::Set(const char *name, const char *value)
{
  if (strcmp(name, "Conjugate total energy gradient") == 0 ||
      strcmp(name, "Conjugate total gradient") == 0) {
    return FromString(value, _ConjugateTotalGradient);
  }
  return GradientDescent::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList ConjugateGradientDescent::Parameter() const
{
  ParameterList params = GradientDescent::Parameter();
  Insert(params, "Conjugate total energy gradient", _ConjugateTotalGradient);
  return params;
}

// =============================================================================
// Optimization
// =============================================================================

// -----------------------------------------------------------------------------
void ConjugateGradientDescent::Initialize()
{
  GradientDescent::Initialize();
  Deallocate(_g), Allocate(_g, _Function->NumberOfDOFs());
  Deallocate(_h), Allocate(_h, _Function->NumberOfDOFs());
  ResetConjugateGradient();
}

// -----------------------------------------------------------------------------
void ConjugateGradientDescent::Finalize()
{
  GradientDescent::Finalize();
  Deallocate(_g);
  Deallocate(_h);
}

// -----------------------------------------------------------------------------
void ConjugateGradientDescent::Gradient(double *gradient, double step, bool *sgn_chg)
{
  // Compute gradient of objective function
  Function()->DataFidelityGradient(gradient, step, sgn_chg);
  if (!_UseConjugateGradient || _ConjugateTotalGradient) {
    Function()->AddConstraintGradient(gradient, step, sgn_chg);
  }

  // Conjugate gradient
  if (_UseConjugateGradient) {
    ConjugateGradient(gradient);
  } else {
    ResetConjugateGradient();
  }

  // Add non-conjugated constraint gradient
  if (_UseConjugateGradient && !_ConjugateTotalGradient) {
    Function()->AddConstraintGradient(gradient, step, sgn_chg);
  }
}

// -----------------------------------------------------------------------------
void ConjugateGradientDescent::ConjugateGradient(double *gradient)
{
  const int ndofs = _Function->NumberOfDOFs();
  if (IsNaN(_g[0])) {
    for (int i = 0; i < ndofs; ++i) _g[i] = -gradient[i];
    memcpy(_h, _g, ndofs * sizeof(double));
  } else {
    double gg  = .0;
    double dgg = .0;
    for (int i = 0; i < ndofs; ++i) {
      gg  += _g[i] * _g[i];
      dgg += (gradient[i] + _g[i]) * gradient[i];
    }
    double gamma = max(dgg / gg, .0);
    for (int i = 0; i < ndofs; ++i) {
      _g[i] = -gradient[i];
      _h[i] = _g[i] + gamma * _h[i];
      gradient[i] = -_h[i];
    }
  }
}


} // namespace mirtk

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

#include "mirtk/EnergyTerm.h"

#include "mirtk/Math.h"
#include "mirtk/Algorithm.h" // transform


namespace mirtk {


// =============================================================================
// Factory method
// =============================================================================

// Define singleton factory class for instantiation of energy term
mirtkDefineObjectFactory(EnergyMeasure, EnergyTerm);

// -----------------------------------------------------------------------------
EnergyTerm::FactoryType &EnergyTerm::Factory()
{
  return EnergyTermFactory::Instance();
}

// -----------------------------------------------------------------------------
EnergyTerm *EnergyTerm::TryNew(enum EnergyMeasure em, const char *name, double w)
{
  EnergyTerm *term = Factory().New(em);
  if (term) {
    term->Name(name);
    term->Weight(w);
  }
  return term;
}

// -----------------------------------------------------------------------------
EnergyTerm *EnergyTerm::New(enum EnergyMeasure em, const char *name, double w)
{
  EnergyTerm *term = TryNew(em, name, w);
  if (term) return term;
  cerr << NameOfType() << "::New: Energy term unknown or not available: " << em << endl;
  exit(1);
  return nullptr;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
EnergyTerm::EnergyTerm(const char *name, double weight)
:
  Configurable(name),
  _Weight(weight),
  _Transformation(NULL),
  _DivideByInitialValue(false)
{
  ResetValue();
  ResetInitialValue();
}

// -----------------------------------------------------------------------------
void EnergyTerm::CopyAttributes(const EnergyTerm &other)
{
  _Weight               = other._Weight;
  _Transformation       = other._Transformation;
  _DivideByInitialValue = other._DivideByInitialValue;
  _InitialValue         = other._InitialValue;
  _Value                = other._Value;
}

// -----------------------------------------------------------------------------
EnergyTerm::EnergyTerm(const EnergyTerm &other)
:
  Configurable(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
EnergyTerm &EnergyTerm::operator =(const EnergyTerm &other)
{
  if (this != &other) {
    Configurable::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
EnergyTerm::~EnergyTerm()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool EnergyTerm::SetWithPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Divide energy terms by initial value") == 0) {
    return FromString(value, _DivideByInitialValue);
  }
  return Configurable::SetWithPrefix(param, value);
}

// -----------------------------------------------------------------------------
bool EnergyTerm::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Relative to initial value") == 0) {
    return FromString(value, _DivideByInitialValue);
  }
  if (strcmp(param, "Weight") == 0) {
    double weight;
    if (!FromString(value, weight)) return false;
    // Keep sign of weight by default, as it depends on whether the energy term
    // is being minimized or maximized; only change the relative weighting of
    // the term; a negative value can be used to flip the sign
    _Weight = copysign(1.0, _Weight) * weight;
    return true;
  }
  if (strcmp(param, "Weight (signed)") == 0) {
    return FromString(value, _Weight);
  }
  return Configurable::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList EnergyTerm::Parameter() const
{
  ParameterList params;
  InsertWithPrefix(params, "Weight (signed)",           _Weight);
  InsertWithPrefix(params, "Relative to initial value", _DivideByInitialValue);
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void EnergyTerm::Initialize()
{
  ResetValue();
  ResetInitialValue();
}

// -----------------------------------------------------------------------------
void EnergyTerm::Update(bool)
{
  ResetValue();
}

// -----------------------------------------------------------------------------
bool EnergyTerm::Upgrade()
{
  return false;
}

// -----------------------------------------------------------------------------
void EnergyTerm::ResetValue()
{
  _Value = numeric_limits<double>::quiet_NaN();
}

// -----------------------------------------------------------------------------
void EnergyTerm::ResetInitialValue()
{
  _InitialValue = numeric_limits<double>::quiet_NaN();
}

// -----------------------------------------------------------------------------
double EnergyTerm::InitialValue()
{
  if (IsNaN(_InitialValue)) {
    _InitialValue = _Value = (_Weight != .0 ? this->Evaluate() : .0);
  }
  return _InitialValue;
}

// -----------------------------------------------------------------------------
double EnergyTerm::Value()
{
  if (IsNaN(_Value))        _Value        = (_Weight != .0 ? this->Evaluate() : .0);
  if (IsNaN(_InitialValue)) _InitialValue = _Value;

  double value = _Value;
  if (_DivideByInitialValue && _InitialValue != .0) {
    value /= abs(_InitialValue);
  }

  return _Weight * value;
}

// -----------------------------------------------------------------------------
double EnergyTerm::RawValue(double value) const
{
  if (_Weight != .0) value /= _Weight;
  if (_DivideByInitialValue) {
    double abs_init_value = abs(_InitialValue);
    if (abs_init_value > 0) value *= abs_init_value;
  }
  return value;
}

// -----------------------------------------------------------------------------
double EnergyTerm::RawValue()
{
  return this->RawValue(this->Value());
}

// -----------------------------------------------------------------------------
void EnergyTerm::Gradient(double *gradient, double step)
{
  double weight = _Weight;
  if (_DivideByInitialValue) {
    if (IsNaN(_InitialValue)) _InitialValue = this->InitialValue();
    if (_InitialValue != .0) weight /= abs(_InitialValue);
  }
  if (weight != .0) this->EvaluateGradient(gradient, step, weight);
}

// -----------------------------------------------------------------------------
void EnergyTerm::NormalizedGradient(double *gradient, double step)
{
  if (_Weight != .0) {
    const int ndofs = _Transformation->NumberOfDOFs();
    double *grad = CAllocate<double>(ndofs);
    this->EvaluateGradient(grad, step, 1.0);
    double norm = _Transformation->DOFGradientNorm(grad);
    if (norm > .0) norm = _Weight / norm;
    for (int dof = 0; dof < ndofs; ++dof) gradient[dof] += norm * grad[dof];
    Deallocate(grad);
  }
}

// -----------------------------------------------------------------------------
void EnergyTerm::GradientStep(const double *, double &, double &) const
{
  // By default, step length range chosen by user/optimizer
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void EnergyTerm::Print(Indent indent) const
{
  cout << indent << "Name:    " << _Name   << "\n";
  cout << indent << "Weight:  " << _Weight << "\n";
  cout << indent << "Initial: " << _InitialValue << "\n";
  cout << indent << "Value:   " << _Value << "\n";
}

// -----------------------------------------------------------------------------
string EnergyTerm::Prefix(const char *prefix) const
{
  if (_Name.empty()) return string(prefix);
  string name = ToLower(_Name);
  replace(name.begin(), name.end(), '.',  '_');
  replace(name.begin(), name.end(), ' ',  '_');
  replace(name.begin(), name.end(), '\t', '_');
  if (prefix) name.insert(0, prefix);
  name += "_";
  return name;
}

// -----------------------------------------------------------------------------
void EnergyTerm::WriteDataSets(const char *, const char *, bool) const
{
}

// -----------------------------------------------------------------------------
void EnergyTerm::WriteGradient(const char *, const char *) const
{
}


} // namespace mirtk

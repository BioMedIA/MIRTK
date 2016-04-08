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

#ifndef MIRTK_EnergyTerm_H
#define MIRTK_EnergyTerm_H

#include "mirtk/Configurable.h"

#include "mirtk/Indent.h"
#include "mirtk/EnergyMeasure.h"
#include "mirtk/Transformation.h"
#include "mirtk/ObjectFactory.h"


namespace mirtk {


/**
 * Base class for one term of an objective function
 *
 * In particular, this is the base class for both the data similarity term
 * and transformation regularization term commonly seen in objective functions
 * used for image/surface registration.
 */
class EnergyTerm : public Configurable
{
  mirtkAbstractMacro(EnergyTerm);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Weight of energy term
  mirtkPublicAttributeMacro(double, Weight);

  /// Transformation with free parameters of energy function
  mirtkPublicAggregateMacro(class Transformation, Transformation);

  /// Whether to divide energy term by its initial value
  mirtkPublicAttributeMacro(bool, DivideByInitialValue);

  /// Initial unweighted value of energy term
  double _InitialValue;

  /// Cached unweighted value of energy term
  double _Value;

  /// Copy attributes of this class from another instance
  void CopyAttributes(const EnergyTerm &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  EnergyTerm(const char * = "", double = 1.0);

  /// Copy constructor
  EnergyTerm(const EnergyTerm &);

  /// Assignment operator
  EnergyTerm &operator =(const EnergyTerm &);

public:

  /// Type of energy term factory
  typedef ObjectFactory<enum EnergyMeasure, EnergyTerm> FactoryType;

  /// Get global energy term factory instance
  static FactoryType &Factory();

  /// Construct new energy term or return nullptr if term not available
  static EnergyTerm *TryNew(EnergyMeasure, const char * = "", double = 1.0);

  /// Construct new energy term
  static EnergyTerm *New(EnergyMeasure, const char * = "", double = 1.0);

  /// Destructor
  virtual ~EnergyTerm();

  /// Energy measure implemented by this term
  virtual enum EnergyMeasure EnergyMeasure() const = 0;

  // ---------------------------------------------------------------------------
  // Parameters

protected:

  /// Set parameter value from string
  virtual bool SetWithPrefix(const char *, const char *);

  /// Set parameter value from string
  virtual bool SetWithoutPrefix(const char *, const char *);

public:

  // Import other overloads
  using Configurable::Parameter;

  /// Get parameter key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Initialize energy term once input and parameters have been set
  virtual void Initialize();

  /// Update internal state after change of DoFs
  ///
  /// \param[in] gradient Whether to also update internal state for evaluation
  ///                     of energy gradient. If \c false, only the internal state
  ///                     required for the energy evaluation need to be updated.
  virtual void Update(bool gradient = true);

  /// Update energy term after convergence
  virtual bool Upgrade();

  /// Reset initial value of energy term
  void ResetInitialValue();

  /// Reset cached value of energy term
  void ResetValue();

  /// Returns initial value of energy term
  double InitialValue();

  /// Evaluate energy term
  double Value();

  /// Evaluate gradient of energy term
  ///
  /// \param[in,out] gradient Gradient to which the evaluated gradient of this
  ///                         energy term is added to with its resp. weight.
  /// \param[in]     step     Step length for finite differences.
  void Gradient(double *gradient, double step);

  /// Evaluate and normalize gradient of energy term
  ///
  /// \param[in,out] gradient Gradient to which the evaluated normalized gradient
  ///                         of this energy term is added to with its resp. weight.
  /// \param[in]     step     Step length for finite differences.
  void NormalizedGradient(double *gradient, double step);

  /// Adjust step length range
  ///
  /// \param[in]     gradient Gradient of objective function.
  /// \param[in,out] min      Minimum step length.
  /// \param[in,out] max      Maximum step length.
  virtual void GradientStep(const double *gradient, double &min, double &max) const;

protected:

  /// Evaluate unweighted energy term
  virtual double Evaluate() = 0;

  /// Evaluate and add gradient of energy term
  ///
  /// \param[in,out] gradient Gradient to which the computed gradient of the
  ///                         energy term should be added to.
  /// \param[in]     step     Step length for finite differences.
  /// \param[in]     weight   Weight to use when adding the gradient.
  virtual void EvaluateGradient(double *gradient, double step, double weight) = 0;

  // ---------------------------------------------------------------------------
  // Debugging

public:

  /// Return unweighted and unnormalized raw energy term value
  /// \remarks Use for progress reporting only.
  virtual double RawValue(double) const;

  /// Return unweighted and unnormalized raw energy term value
  /// \remarks Use for progress reporting only.
  double RawValue();

  /// Print debug information
  virtual void Print(Indent = 0) const;

  /// Prefix to be used for debug output files
  string Prefix(const char * = NULL) const;

  /// Write input of data fidelity term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

  /// Write gradient of data fidelity term w.r.t each transformed input
  virtual void WriteGradient(const char *, const char *) const;

};

////////////////////////////////////////////////////////////////////////////////
// Auxiliary macros for optimizer implementation
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
#define mirtkEnergyTermMacro(name, id)                                         \
  mirtkObjectMacro(name);                                                      \
public:                                                                        \
  /** Energy measure implemented by this term */                               \
  static mirtk::EnergyMeasure ID() { return id; }                              \
  /** Energy measure implemented by this term */                               \
  virtual mirtk::EnergyMeasure EnergyMeasure() const { return id; }            \
private:

// -----------------------------------------------------------------------------
/// Register object type with factory singleton
#define mirtkRegisterEnergyTermMacro(type)                                     \
  mirtk::EnergyTerm::Factory().Register(type::ID(), type::NameOfType(),        \
                                        mirtk::New<mirtk::EnergyTerm, type>)

// -----------------------------------------------------------------------------
/// Register object type with factory singleton at static initialization time
#define mirtkAutoRegisterEnergyTermMacro(type)                                 \
  mirtkAutoRegisterObjectTypeMacro(mirtk::EnergyTerm::Factory(),               \
                                   mirtk::EnergyMeasure, type::ID(),           \
                                   mirtk::EnergyTerm,    type)


} // namespace mirtk

#endif // MIRTK_EnergyTerm_H

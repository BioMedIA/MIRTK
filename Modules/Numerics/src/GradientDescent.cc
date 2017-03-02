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

#include "mirtk/GradientDescent.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/ObjectFactory.h"

#include <algorithm>


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterOptimizerMacro(GradientDescent);


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
GradientDescent::GradientDescent(ObjectiveFunction *f)
:
  LocalOptimizer(f),
  _NumberOfRestarts      (-1),
  _NumberOfFailedRestarts(-1),
  _LineSearchStrategy    (LS_Adaptive),
  _LineSearch            (NULL),
  _LineSearchOwner       (false),
  _Gradient              (NULL),
  _AllowSignChange       (NULL)
{
  _EventDelegate.Bind(MakeDelegate(this, &Observable::Broadcast));
  this->Epsilon(- this->Epsilon()); // relative to current best value
}

// -----------------------------------------------------------------------------
void GradientDescent::CopyAttributes(const GradientDescent &other)
{
  _NumberOfRestarts       = other._NumberOfRestarts;
  _NumberOfFailedRestarts = other._NumberOfFailedRestarts;
  _LineSearchStrategy     = other._LineSearchStrategy;
  _LineSearchParameter    = other._LineSearchParameter;
  Deallocate(_Gradient);
  Deallocate(_AllowSignChange);
  if (_LineSearchOwner) Delete(_LineSearch);
  if (!other._LineSearchOwner) {
    _LineSearch      = other._LineSearch;
    _LineSearchOwner = false;
  }
}

// -----------------------------------------------------------------------------
GradientDescent::GradientDescent(const GradientDescent &other)
:
  LocalOptimizer(other),
  _Gradient       (NULL),
  _AllowSignChange(NULL)
{
  _EventDelegate.Bind(MakeDelegate(this, &Observable::Broadcast));
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
GradientDescent &GradientDescent::operator =(const GradientDescent &other)
{
  if (this != &other) {
    LocalOptimizer::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
GradientDescent::~GradientDescent()
{
  Deallocate(_Gradient);
  Deallocate(_AllowSignChange);
  if (_LineSearchOwner) Delete(_LineSearch);
}

// -----------------------------------------------------------------------------
void GradientDescent::Function(ObjectiveFunction *f)
{
  // Reallocation of possibly previously allocated gradient vector required
  // if new function has differing number of DoFs
  if (!f || !this->Function() || this->Function()->NumberOfDOFs() != f->NumberOfDOFs()) {
    Deallocate(_Gradient);
    Deallocate(_AllowSignChange);
  }
  LocalOptimizer::Function(f);
}

// -----------------------------------------------------------------------------
void GradientDescent::LineSearch(class LineSearch *s, bool transfer_ownership)
{
  if (_LineSearchOwner) {
    Delete(_LineSearch);
    _LineSearchOwner = false;
  }
  if (s) {
    _LineSearchStrategy = s->Strategy();
    _LineSearch         = s;
    _LineSearchOwner    = transfer_ownership;
  }
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool GradientDescent::Set(const char *name, const char *value)
{
  if (strcmp(name, "Maximum no. of restarts")    == 0 ||
      strcmp(name, "Maximum number of restarts") == 0 ||
      strcmp(name, "No. of restarts")            == 0 ||
      strcmp(name, "Number of restarts")         == 0) {
    return FromString(value, _NumberOfRestarts);
  }
  if (strcmp(name, "Maximum no. of failed restarts")          == 0 ||
      strcmp(name, "Maximum number of failed restarts")       == 0 ||
      strcmp(name, "No. of failed restarts")                  == 0 ||
      strcmp(name, "Number of failed restarts")               == 0 ||
      strcmp(name, "Maximum no. of unsuccessful restarts")    == 0 ||
      strcmp(name, "Maximum number of unsuccessful restarts") == 0 ||
      strcmp(name, "No. of unsuccessful restarts")            == 0 ||
      strcmp(name, "Number of unsuccessful restarts")         == 0 ||
      strcmp(name, "Maximum streak of failed restarts")       == 0 ||
      strcmp(name, "Maximum streak of unsuccessful restarts") == 0 ||
      strcmp(name, "Maximum streak of rejected restarts")     == 0) {
    return FromString(value, _NumberOfFailedRestarts);
  }
  if (strcmp(name, "Line search strategy") == 0) {
    return FromString(value, _LineSearchStrategy);
  }
  if (strstr(name, "line search") != NULL ||
      strstr(name, "line iterations") != NULL ||
      strcmp(name, "Maximum streak of rejected steps") == 0    ||
      strcmp(name, "Length of steps")  == 0 ||
      strcmp(name, "Minimum length of steps") == 0 ||
      strcmp(name, "Maximum length of steps") == 0 ||
      strcmp(name, "Strict step length range") == 0 ||
      strcmp(name, "Strict total step length range") == 0 ||
      strcmp(name, "Strict incremental step length range") == 0 ||
      strcmp(name, "Strict accumulated step length range") == 0 ||
      strcmp(name, "Step length rise") == 0 ||
      strcmp(name, "Step length drop") == 0 ||
      strcmp(name, "Reuse previous step length") == 0) {
    Insert(_LineSearchParameter, name, value);
    return true;
  }
  return LocalOptimizer::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList GradientDescent::Parameter() const
{
  ParameterList params = LocalOptimizer::Parameter();
  if (_LineSearch && _LineSearch->Strategy() == _LineSearchStrategy) {
    Insert(params, _LineSearch->Parameter());
  }
  Insert(params, _LineSearchParameter);
  Insert(params, "Maximum no. of restarts",        _NumberOfRestarts);
  Insert(params, "Maximum no. of failed restarts", _NumberOfFailedRestarts);
  Insert(params, "Line search strategy",           _LineSearchStrategy);
  return params;
}

// =============================================================================
// Optimization
// =============================================================================

// -----------------------------------------------------------------------------
void GradientDescent::Initialize()
{
  // Initialize base class -- checks that _Function is valid
  LocalOptimizer::Initialize();
  // Default values
  if (_NumberOfRestarts       < 0) _NumberOfRestarts       = 100;
  if (_NumberOfFailedRestarts < 0) _NumberOfFailedRestarts =   5;
  // Use existing line search object which was created when user called
  // Initialize before Run in order to be able to set the line search
  // parameters directly through the setters instead of the generic Set.
  if (!_LineSearch || _LineSearch->Strategy() != _LineSearchStrategy) {
    if (_LineSearchOwner) Delete(_LineSearch);
    _LineSearch      = LineSearch::New(_LineSearchStrategy);
    _LineSearchOwner = true;
  }
  // Allocate memory for gradient vector if not done before
  if (!_Gradient)        Allocate(_Gradient,        Function()->NumberOfDOFs());
  if (!_AllowSignChange) Allocate(_AllowSignChange, Function()->NumberOfDOFs());
  // Initialize line search object
  _LineSearch->Function   (Function());
  _LineSearch->Parameter  (_LineSearchParameter);
  _LineSearch->Delta      (_Delta);
  _LineSearch->Epsilon    (max(.0, _Epsilon));
  _LineSearch->Direction  (_Gradient);
  _LineSearch->Revert     (true);
  _LineSearch->AddObserver(_EventDelegate);
  _LineSearch->Initialize();
  // Check line search parameters
  if (_LineSearch->MaxStepLength() == .0) {
    cerr << this->NameOfClass() << "::Initialize: Line search interval length is zero!" << endl;
    cerr << "  Check the \"Minimum/Maximum length of steps\" line search parameter." << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
double GradientDescent::Run()
{
  double value     = numeric_limits<double>::max();
  int    nrestarts = 0; // Number of restarts after convergence
  int    nfailed   = 0; // Number of consecutive restarts without improvement

  // Initialize
  this->Initialize();

  // Initial line search parameters
  const double initial_delta    = _LineSearch->Delta();
  const double initial_epsilon  = _LineSearch->Epsilon();
  const double initial_step     = _LineSearch->StepLength();
  const double initial_min_step = _LineSearch->MinStepLength();
  const double initial_max_step = _LineSearch->MaxStepLength();

  // Perform initial update of energy function before StartEvent because
  // it may trigger some further delayed initialization with LogEvent's
  Function()->Update(true);

  // Notify observers about start of optimization
  Broadcast(StartEvent);

  // Total number of performed gradient steps
  Iteration step(0, _NumberOfSteps);

  // Repeat gradient descent optimization for each restart with modified
  // energy function. If energy remains fixed, only one iteration is done.
  while (true) {

    // Get initial value of (modified) energy function
    value = _Function->Value();
    _LastValues.clear();
    _LastValues.push_back(value);

    // Current number of iterations
    const int iter = step.Iter();

    // Descent along computed gradient
    _Converged = false;
    while (!_Converged && step.Next()) {

      // Notify observers about start of gradient descent iteration
      Broadcast(IterationStartEvent, &step);

      // Update current best value
      _LineSearch->CurrentValue(value);

      // Compute gradient of objective function
      if (step.Iter() > 1) Function()->Update(true);
      this->Gradient(_Gradient, _LineSearch->StepLength(), _AllowSignChange);

      // Adjust step length range if necessary
      //
      // This is required for the optimization of the registration cost function
      // with L1-norm sparsity constraint on the multi-level free-form deformation
      // parameters (Wenzhe et al.'s Sparse FFD). Furthermore, if the gradient of
      // an energy term is approximated using finite differences, this energy term
      // will set min_step = max_step such that the step length corresponds to
      // the one with which the gradient was approximated.
      double max_norm = Function()->GradientNorm(_Gradient);
      double min_step = _LineSearch->MinStepLength();
      double max_step = _LineSearch->MaxStepLength();
      Function()->GradientStep(_Gradient, min_step, max_step);

      // Set line search range
      _LineSearch->MinStepLength (min_step);
      _LineSearch->MaxStepLength (max_step);
      _LineSearch->StepLengthUnit(max_norm);

      // Perform line search along computed gradient direction
      value = _LineSearch->Run();

      // Adjust epsilon if relative to current best value, i.e.,
      // epsilon parameter is set to a negative value
      if (_Epsilon < .0) _LineSearch->Epsilon(abs(_Epsilon * value));

      // Check convergence
      if (_LineSearch->StepLength() > 0.) {
        _Converged = this->Converged(step.Iter(), value, _Gradient);
      } else {
        _Converged = true;
      }

      // Notify observers about end of gradient descent iteration
      Broadcast(IterationEndEvent, &step);
    }

    // Stop if previous restart did not bring any improvement
    if (step.Iter() == (iter + 1) && value == _LineSearch->CurrentValue()) {
      ++nfailed;
      if (nfailed >= _NumberOfFailedRestarts) break;
    } else {
      nfailed = 0;
    }

    // Update current best value
    _LineSearch->CurrentValue(value);

    // Stop if maximum number of iterations exceeded
    if (step.End()) break;

    // Stop if maximum number of allowed restarts exceeded
    if (nrestarts >= _NumberOfRestarts) break;
    ++nrestarts;

    // If there was no improvement compared to the previous converged
    // iterative gradient descent, or the objective function remains
    // unmodified, stop here. Otherwise, start another optimization of
    // the amended objective (cf. FiducialRegistrationError).
    // Such restart realizes an alternating optimization, where the
    // ObjectiveFunction::Upgrade is performing the optimization w.r.t.
    // some function parameters different from the "public" DoFs.
    if (!_Function->Upgrade()) break;

    // Update energy function for initial value evaluation as well as for
    // the first gradient computation
    Function()->Update(true);

    // Reset line search parameters
    _LineSearch->Delta        (initial_delta);
    _LineSearch->Epsilon      (initial_epsilon);
    _LineSearch->StepLength   (initial_step);
    _LineSearch->MinStepLength(initial_min_step);
    _LineSearch->MaxStepLength(initial_max_step);

    // Notify observers about restart
    Broadcast(RestartEvent);
  }

  // Notify observers about end of optimization
  Broadcast(EndEvent, &value);

  // Finalize
  this->Finalize();

  return value;
}

// -----------------------------------------------------------------------------
void GradientDescent::Gradient(double *gradient, double step, bool *sgn_chg)
{
  Function()->Gradient(gradient, step, sgn_chg);
}

// -----------------------------------------------------------------------------
void GradientDescent::Finalize()
{
  Deallocate(_Gradient);
  Deallocate(_AllowSignChange);
  _LineSearch->DeleteObserver(_EventDelegate);
  if (_LineSearchOwner) Delete(_LineSearch);
}


} // namespace mirtk

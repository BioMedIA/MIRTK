/*
 * Medical Image Registration ToolKit (MIRTK) LBFGS Library
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include "mirtk/LimitedMemoryBFGSDescent.h"

#include "mirtk/ObjectFactory.h"

#include "lbfgs.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterOptimizerMacro(LimitedMemoryBFGSDescent);


// =============================================================================
// liblbfgs error messages
// =============================================================================

// -----------------------------------------------------------------------------
const char *lbfgs_error_message(int err)
{
  switch (err) {
    case LBFGS_ALREADY_MINIMIZED:            return "The initial variables already minimize the objective function.";
    case LBFGSERR_LOGICERROR:                return "Logic error.";
    case LBFGSERR_OUTOFMEMORY:               return "Insufficient memory.";
    case LBFGSERR_CANCELED:                  return "The minimization process has been canceled.";
    case LBFGSERR_INVALID_N:                 return "Invalid number of variables specified.";
    case LBFGSERR_INVALID_N_SSE:             return "Invalid number of variables (for SSE) specified.";
    case LBFGSERR_INVALID_X_SSE:             return "The array x must be aligned to 16 (for SSE).";
    case LBFGSERR_INVALID_EPSILON:           return "Invalid parameter lbfgs_parameter_t::epsilon specified.";
    case LBFGSERR_INVALID_TESTPERIOD:        return "Invalid parameter lbfgs_parameter_t::past specified.";
    case LBFGSERR_INVALID_DELTA:             return "Invalid parameter lbfgs_parameter_t::delta specified.";
    case LBFGSERR_INVALID_LINESEARCH:        return "Invalid parameter lbfgs_parameter_t::linesearch specified.";
    case LBFGSERR_INVALID_MINSTEP:           return "Invalid parameter lbfgs_parameter_t::max_step specified.";
    case LBFGSERR_INVALID_MAXSTEP:           return "Invalid parameter lbfgs_parameter_t::max_step specified.";
    case LBFGSERR_INVALID_FTOL:              return "Invalid parameter lbfgs_parameter_t::ftol specified.";
    case LBFGSERR_INVALID_WOLFE:             return "Invalid parameter lbfgs_parameter_t::wolfe specified.";
    case LBFGSERR_INVALID_GTOL:              return "Invalid parameter lbfgs_parameter_t::gtol specified.";
    case LBFGSERR_INVALID_XTOL:              return "Invalid parameter lbfgs_parameter_t::xtol specified.";
    case LBFGSERR_INVALID_MAXLINESEARCH:     return "Invalid parameter lbfgs_parameter_t::max_linesearch specified.";
    case LBFGSERR_INVALID_ORTHANTWISE:       return "Invalid parameter lbfgs_parameter_t::orthantwise_c specified.";
    case LBFGSERR_INVALID_ORTHANTWISE_START: return "Invalid parameter lbfgs_parameter_t::orthantwise_start specified.";
    case LBFGSERR_INVALID_ORTHANTWISE_END:   return "Invalid parameter lbfgs_parameter_t::orthantwise_end specified.";
    case LBFGSERR_OUTOFINTERVAL:             return "The line-search step went out of the interval of uncertainty.";
    case LBFGSERR_INCORRECT_TMINMAX:         return "A logic error occurred; alternatively, the interval of uncertainty became too small.";
    case LBFGSERR_ROUNDING_ERROR:            return "A rounding error occurred; alternatively, no line-search step satisfies the sufficient decrease and curvature conditions.";
    case LBFGSERR_MINIMUMSTEP:               return "The line-search step became smaller than lbfgs_parameter_t::min_step.";
    case LBFGSERR_MAXIMUMSTEP:               return "The line-search step became larger than lbfgs_parameter_t::max_step.";
    case LBFGSERR_MAXIMUMLINESEARCH:         return "The line-search routine reaches the maximum number of evaluations.";
    case LBFGSERR_MAXIMUMITERATION:          return "The algorithm routine reaches the maximum number of iterations.";
    case LBFGSERR_WIDTHTOOSMALL:             return "Relative width of the interval of uncertainty is at most lbfgs_parameter_t::xtol.";
    case LBFGSERR_INVALIDPARAMETERS:         return "A logic error (negative line-search step) occurred.";
    case LBFGSERR_INCREASEGRADIENT:          return "The current search direction increases the objective function value.";
    default:                                 return "Unknown error.";
  }
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
LimitedMemoryBFGSDescent::LimitedMemoryBFGSDescent(ObjectiveFunction *f)
:
  LocalOptimizer(f),
  _NumberOfIterations(100),
  _MinStepLength     (0.01),
  _MaxStepLength     (1.0)
{
}

// -----------------------------------------------------------------------------
void LimitedMemoryBFGSDescent::CopyAttributes(const LimitedMemoryBFGSDescent &other)
{
  _NumberOfIterations = other._NumberOfIterations;
  _MinStepLength      = other._MinStepLength;
  _MaxStepLength      = other._MaxStepLength;
  _CurrentStep        = other._CurrentStep;
}

// -----------------------------------------------------------------------------
LimitedMemoryBFGSDescent::LimitedMemoryBFGSDescent(const LimitedMemoryBFGSDescent &other)
:
  LocalOptimizer(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
LimitedMemoryBFGSDescent &LimitedMemoryBFGSDescent::operator =(const LimitedMemoryBFGSDescent &other)
{
  if (this != &other) {
    LocalOptimizer::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
LimitedMemoryBFGSDescent::~LimitedMemoryBFGSDescent()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool LimitedMemoryBFGSDescent::Set(const char *name, const char *value)
{
  if (strcmp(name, "No. of iterations") == 0) {
    return FromString(value, _NumberOfIterations) && _NumberOfIterations > 0;
  }
  if (strcmp(name, "Minimum length of steps") == 0) {
    return FromString(value, _MinStepLength) && _MinStepLength > .0;
  }
  if (strcmp(name, "Maximum length of steps") == 0) {
    return FromString(value, _MaxStepLength) && _MaxStepLength > .0;
  }
  if (strcmp(name, "Length of steps") == 0) {
    if (!FromString(value, _MaxStepLength) || _MaxStepLength <= .0) return false;
    _MinStepLength = _MaxStepLength;
    return true;
  }
  return LocalOptimizer::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList LimitedMemoryBFGSDescent::Parameter() const
{
  ParameterList params = LocalOptimizer::Parameter();
  Insert(params, "No. of iterations", ToString(_NumberOfIterations));
  if (_MinStepLength == _MaxStepLength) {
    Insert(params, "Length of steps", ToString(_MinStepLength));
  } else {
    Insert(params, "Minimum length of steps", ToString(_MinStepLength));
    Insert(params, "Maximum length of steps", ToString(_MaxStepLength));
  }
  return params;
}

// =============================================================================
// Optimization
// =============================================================================

// L-BFGS helper
struct LBFGSCallback
{
  static lbfgsfloatval_t Evaluate(
      void                  *o,    // Instance of LimitedMemoryBFGSDescent
      const lbfgsfloatval_t *x,    // Current parameters
      lbfgsfloatval_t       *g,    // Gradient of objective function
      const int              n,    // Number of parameters
      const lbfgsfloatval_t  step) // Line search step length
  {
    LimitedMemoryBFGSDescent *_this = reinterpret_cast<LimitedMemoryBFGSDescent *>(o);
    _this->Function()->Put(x);
    lbfgsfloatval_t m = _this->Function()->Evaluate(g, step);
    _this->_CurrentStep._Value = m;
    return m;
  }

  static int Progress(
      void                  *o,        // Instance of LimitedMemoryBFGSDescent
      const lbfgsfloatval_t *x,        // Current parameters
      const lbfgsfloatval_t * /*g*/,   // Current gradient
      const lbfgsfloatval_t fx,        // Current value
      const lbfgsfloatval_t xnorm,     // Norm of parameter vector
      const lbfgsfloatval_t /*gnorm*/, // Norm of gradient vector
      const lbfgsfloatval_t step,      // Current step length
      int n,                           // Number of parameters
      int k,                           // Current iteration
      int ls                           // Current line search iteration
      )
  {
    LimitedMemoryBFGSDescent *_this = reinterpret_cast<LimitedMemoryBFGSDescent *>(o);
    string msg = "Current best metric value is ";
    msg += ToString(fx);
    msg += ", step ";
    msg += ToString(step);
    msg += ", after ";
    msg += ToString(ls);
    msg += " line search steps\n";
    _this->Broadcast(LogEvent, msg.c_str());
    return 0;
  }
};

// -----------------------------------------------------------------------------
double LimitedMemoryBFGSDescent::Run()
{
  lbfgs_parameter_t param;
  lbfgs_parameter_init(&param);

  param.past           = 1;
  param.delta          = _Epsilon;
  param.epsilon        = _Delta;
  param.max_iterations = _NumberOfIterations;
  param.min_step       = _MinStepLength;
  param.max_step       = _MaxStepLength;

  const int        n = Function()->NumberOfDOFs();
  lbfgsfloatval_t *x = lbfgs_malloc(n);
  Function()->Get(x);

  int ret = lbfgs(n, x, NULL, LBFGSCallback::Evaluate, LBFGSCallback::Progress, this, &param);
  if (ret < 0 && ret != LBFGSERR_MAXIMUMITERATION && ret != LBFGSERR_MAXIMUMLINESEARCH) {
    cerr << "L-BFGS optimization failed: " << lbfgs_error_message(ret) << " (error code: " << ret << ")" << endl;
    exit(1);
  }

  lbfgs_free(x);
  return _CurrentStep._Value;
}


} // namespace mirtk

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

#ifndef MIRTK_Event_H
#define MIRTK_Event_H

#include "mirtk/Stream.h"

#include <functional>


namespace mirtk {


/// Events that can be observed
enum Event
{
  // Generic
  AnyEvent,                      ///< Any event
  ModifiedEvent,                 ///< Observable has modified state/parameters
  StatusEvent,                   ///< Status message given as event data (const char *)
  LogEvent,                      ///< Log message given as event data (const char *)
  InitEvent,                     ///< Before initialization of filter execution
  StartEvent,                    ///< Start of filter execution (after initialization)
  RestartEvent,                  ///< Restart of filter execution
  EndEvent,                      ///< End of filter execution (before finalization)
  FinishEvent,                   ///< After finalization of filter execution
  IterationEvent,                ///< Iteration with line search about to start
  IterationStartEvent,           ///< Iteration started
  IterationEndEvent,             ///< Iteration finished
  // Observable
  RegisteredEvent,               ///< Observer has been registered
  UnregisteredEvent,             ///< Observer has been unregistered
  // LineSearch
  LineSearchStartEvent,          ///< Line search started
  LineSearchIterationStartEvent, ///< Line search iteration
  LineSearchIterationEndEvent,   ///< Line search iteration
  LineSearchEndEvent,            ///< Line search finished
  AcceptedStepEvent,             ///< Accepted objective function value
                                 ///< Event data is the new value as double
  RejectedStepEvent              ///< Rejected objective function value
                                 ///< Event data is the rejected value as double
};


/// Data of IterationStartEvent and IterationEndEvent
///
/// This event data structure can further be used as iteration counter
/// directly to avoid the need for unnecessary duplicate counter.
///
/// \code
/// irtkIteration iter(0, 10);
/// while (iter.Next()) {
///   Broadcast(IterationEvent, &iter);
///   // ...
/// }
/// \endcode
struct Iteration
{
  int _Iter;  ///< Iteration counter
  int _Begin; ///< Start value of iteration counter
  int _End;   ///< End value of iteration counter
  int _Inc;   ///< Increment after each iteration

  Iteration() : _Iter(0), _Begin(0), _End(1), _Inc(1) {}

  Iteration(int b, int e, int s = 1)
    : _Iter(b), _Begin(b), _End(e), _Inc((_End >= _Begin ? 1 : -1) * (s < 0 ? -s : s))
  {
    if (_Inc == 0) {
      cerr << "Iteration::Iteration: Invalid increment value" << endl;
      exit(1);
    }
  }

  bool Start() const { return _Iter == _Begin; }
  bool End()   const { return _Iter == _End; }
  int  Total() const { return (_End  - _Begin) / _Inc; }
  int  Count() const { return (_Iter - _Begin) / _Inc; }
  int  Iter()  const { return _Iter; }
  void Reset()       { _Iter = _Begin; }
  void Break()       { _Iter = _End; }

  int operator ++()
  {
    _Iter += _Inc;
    if (_Inc > 0) {
      if (_Iter > _End) _Iter = _End;
    } else {
      if (_Iter < _End) _Iter = _End;
    }
    return _Iter;
  }

  int operator --()
  {
    _Iter -= _Inc;
    if (_Inc > 0) {
      if (_Iter < _Begin) _Iter = _Begin;
    } else {
      if (_Iter > _Begin) _Iter = _Begin;
    }
    return _Iter;
  }

  bool Next()
  {
    if (_Iter != _End) {
      _Iter += _Inc;
      if (_Inc > 0) {
        if (_Iter > _End) {
          _Iter = _End;
          return false;
        }
      } else {
        if (_Iter < _End) {
          _Iter = _End;
          return false;
        }
      }
      return true;
    }
    return false;
  }
};

/// Data of AcceptedStepEvent and RejectedStepEvent
struct LineSearchStep
{
  const char *_Info;        ///< Adjective detailing step for log output
  double     *_Direction;   ///< Step direction (i.e., function gradient)
  double      _Current;     ///< Current objective function value
  double      _Value;       ///< Accepted/rejected function value
  double      _Length;      ///< Accepted/rejected step length
  double      _TotalLength; ///< Accumulated length of accepted steps
  double      _MinLength;   ///< Minimum allowed step length
  double      _MaxLength;   ///< Maximum allowed step length
  double      _Unit;        ///< Step length unit
  double      _Delta;       ///< Maximum change of parameters (DoFs)
  double      _TotalDelta;  ///< Accumulated change of parameters (DoFs)

  LineSearchStep()
  :
    _Info(NULL), _Direction(NULL),
    _Current(.0), _Value(.0), _Length(.0), _TotalLength(.0),
    _MinLength(.0), _MaxLength(.0), _Unit(.0),
    _Delta(.0), _TotalDelta(.0)
  {}
};


} // namespace mirtk


namespace std {

/// Compute hash value from Event enumeration value
template<>
struct hash<mirtk::Event> {
    size_t operator()(const mirtk::Event &enum_value) const {
        return std::hash<int>()(enum_value);
    }
};


} // namespace std

#endif // MIRTK_Event_H

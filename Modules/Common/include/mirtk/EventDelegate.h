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

#ifndef MIRTK_EventDelegate_H
#define MIRTK_EventDelegate_H

#include "mirtk/Observer.h"

#include "mirtk/FastDelegate.h"
#include "mirtk/UnorderedMap.h"


namespace mirtk {


/**
 * Invokes a callback function for observed events
 *
 * This class can be registered as observer of an observable object derived
 * from Observable. It receives the event notifications from the observed
 * object(s) and dispatches the corresponding registered callback (member)
 * function (if any). It therefore uses instances of FastDelegate.
 *
 * \note Do not derive from this class. Instead add instances of it as attributes
 *       to classes which observe other objects and handle the events.
 *
 * \code
 * EventDelegate delegate;
 * // Bind non-member function to any event
 * delegate.Bind(&function);
 * // Bind static member function to EndEvent
 * delegate.Bind(EndEvent, &function);
 * // Bind non-static member function to IterationEvent
 * delegate.Bind(IterationEvent, MakeDelegate(&obj, &method));
 * \endcode
 */
class EventDelegate : public Observer
{
  mirtkObjectMacro(EventDelegate);

  // ---------------------------------------------------------------------------
  // Types

  typedef FastDelegate0<>                                  Delegate0;
  typedef FastDelegate1<Event>                             Delegate1;
  typedef FastDelegate2<Event, const void *>               Delegate2;
  typedef FastDelegate3<Observable *, Event, const void *> Delegate3;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Copy constructor
  /// \note Intentionally not implemented
  EventDelegate(const EventDelegate &);

  /// Assignment operator
  /// \note Intentionally not implemented
  EventDelegate &operator =(const EventDelegate &);

public:

  /// Default constructor
  EventDelegate();

  /// Destructor
  ~EventDelegate();

  // ---------------------------------------------------------------------------
  // Bind callback functions

  /// Bind callback (member) function for any event with no arguments
  void Bind(Delegate0);

  /// Bind callback (member) function for any event with one argument
  void Bind(Delegate1);

  /// Bind callback (member) function for any event with two arguments
  void Bind(Delegate2);

  /// Bind callback (member) function for any event with four arguments
  void Bind(Delegate3);

  /// Bind callback (member) function to specified event
  void Bind(Event, Delegate0);

  /// Bind callback (member) function to specified event
  void Bind(Event, Delegate1);

  /// Bind callback (member) function to specified event
  void Bind(Event, Delegate2);

  /// Bind callback (member) function to specified event
  void Bind(Event, Delegate3);

  // ---------------------------------------------------------------------------
  // Event forwarding

  /// Forward event notification to registered callback (member) function
  void HandleEvent(Observable *, Event, const void * = NULL);

  // ---------------------------------------------------------------------------
  // Attributes
private:

  /// Auxiliary structure to store disperse callback delegates
  struct Callback
  {
    int             _N;
    DelegateMemento _Memento;
    Callback() : _N(0), _Memento() {}
    Callback(DelegateMemento memento, int n) : _N(n), _Memento(memento) {}
    Callback(const Callback &o) : _N(o._N), _Memento(o._Memento) {}
  };

  void Forward(Callback &func, Observable *obj, Event event, const void *data);

  UnorderedMap<Event, Callback> _SingleEventDelegate; ///< Event delegates
  Callback                      _AnyEventDelegate;    ///< AnyEvent delegate
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline EventDelegate::EventDelegate()
{
}

// -----------------------------------------------------------------------------
inline EventDelegate::~EventDelegate()
{
}

// -----------------------------------------------------------------------------
inline void EventDelegate::Bind(Delegate0 delegate)
{
  _AnyEventDelegate = Callback(delegate.GetMemento(), 0);
}

// -----------------------------------------------------------------------------
inline void EventDelegate::Bind(Delegate1 delegate)
{
  _AnyEventDelegate = Callback(delegate.GetMemento(), 1);
}

// -----------------------------------------------------------------------------
inline void EventDelegate::Bind(Delegate2 delegate)
{
  _AnyEventDelegate = Callback(delegate.GetMemento(), 2);
}

// -----------------------------------------------------------------------------
inline void EventDelegate::Bind(Delegate3 delegate)
{
  _AnyEventDelegate = Callback(delegate.GetMemento(), 3);
}

// -----------------------------------------------------------------------------
inline void EventDelegate::Bind(Event event, Delegate0 delegate)
{
  _SingleEventDelegate[event] = Callback(delegate.GetMemento(), 0);
}

// -----------------------------------------------------------------------------
inline void EventDelegate::Bind(Event event, Delegate1 delegate)
{
  _SingleEventDelegate[event] = Callback(delegate.GetMemento(), 1);
}

// -----------------------------------------------------------------------------
inline void EventDelegate::Bind(Event event, Delegate2 delegate)
{
  _SingleEventDelegate[event] = Callback(delegate.GetMemento(), 2);
}

// -----------------------------------------------------------------------------
inline void EventDelegate::Bind(Event event, Delegate3 delegate)
{
  _SingleEventDelegate[event] = Callback(delegate.GetMemento(), 3);
}

// -----------------------------------------------------------------------------
inline void EventDelegate::Forward(Callback &func, Observable *obj, Event event, const void *data)
{
  switch (func._N) {
    case 0: {
      Delegate0 delegate;
      delegate.SetMemento(func._Memento);
      if (delegate) delegate();
      break;
    }
    case 1: {
      Delegate1 delegate;
      delegate.SetMemento(func._Memento);
      if (delegate) delegate(event);
      break;
    }
    case 2: {
      Delegate2 delegate;
      delegate.SetMemento(func._Memento);
      if (delegate) delegate(event, data);
      break;
    }
    case 3: {
      Delegate3 delegate;
      delegate.SetMemento(func._Memento);
      if (delegate) delegate(obj, event, data);
      break;
    }
  }
}

// -----------------------------------------------------------------------------
inline void EventDelegate::HandleEvent(Observable *obj, Event event, const void *data)
{
  Forward(_AnyEventDelegate, obj, event, data);
  UnorderedMap<Event, Callback>::iterator it = _SingleEventDelegate.find(event);
  if (it != _SingleEventDelegate.end()) Forward(it->second, obj, event, data);
}


} // namespace mirtk

#endif // MIRTK_EventDelegate_H

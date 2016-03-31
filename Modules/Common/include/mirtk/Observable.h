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

#ifndef MIRTK_Observable_H
#define MIRTK_Observable_H

#include "mirtk/Object.h"

#include "mirtk/Observer.h"
#include "mirtk/OrderedSet.h"
#include "mirtk/Event.h"


namespace mirtk {


/**
 * Base class for an observable object
 *
 * Any object which should be observable must derive from irtkObservable
 * instead of irtkObject. A client which wants to observe events emitted
 * by an observable must use an instance of irtkObservable and register
 * callback functions at this instance for the events of interest.
 *
 * \attention Adding and removing observers is not thread-safe!
 */
class Observable : public Object
{
  mirtkAbstractMacro(Observable);

  /// Container storing pointers to observers
  typedef OrderedSet<Observer *> ObserverSet;

  /// Iterator for registered observers
  typedef ObserverSet::iterator ObserverIterator;

  /// Whether this object has changed and should notify observers upon request
  mirtkPublicAttributeMacro(bool, Changed);

  /// Registered observers
  mirtkAttributeMacro(ObserverSet, Observers);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const Observable &);

protected:

  /// Default constructor
  Observable();

  /// Copy constructor
  Observable(const Observable &);

  /// Assignment operator
  Observable &operator =(const Observable &);

public:

  /// Destructor
  virtual ~Observable();

  /// Number of current observers
  int NumberOfObservers() const;

  /// Add observer
  void AddObserver(Observer &);

  /// Delete observer
  void DeleteObserver(Observer &);

  /// Delete all observers
  void ClearObservers();

  /// Broadcast event to observers
  void Broadcast(Event, const void * = NULL);

  /// Notify all observers about given event if this object has changed
  void NotifyObservers(Event, const void * = NULL);

};

////////////////////////////////////////////////////////////////////////////////
// Auxiliary macro for broadcasting log messages
////////////////////////////////////////////////////////////////////////////////

/// Broadcast a LogEvent message
#define MIRTK_LOG_EVENT(msg) \
  do { \
    ostringstream ss; \
    ss << msg; \
    Broadcast(mirtk::LogEvent, ss.str().c_str()); \
  } while(false)

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
inline Observable::Observable()
:
  _Changed(false)
{
}

// -----------------------------------------------------------------------------
inline void Observable::CopyAttributes(const Observable &other)
{
  _Changed = other._Changed;
  // Do not copy observers!
}

// -----------------------------------------------------------------------------
inline Observable::Observable(const Observable &other)
:
  Object(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
inline Observable &Observable::operator =(const Observable &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
    Changed(true);
  }
  return *this;
}

// -----------------------------------------------------------------------------
inline void Observable::ClearObservers()
{
  for (ObserverIterator it = _Observers.begin(); it != _Observers.end(); ++it) {
    (*it)->StopObserving(this);
  }
  _Observers.clear();
}

// -----------------------------------------------------------------------------
inline Observable::~Observable()
{
  ClearObservers();
}

// =============================================================================
// Add/Remove observer
// =============================================================================

// -----------------------------------------------------------------------------
inline int Observable::NumberOfObservers() const
{
  return static_cast<int>(_Observers.size());
}

// -----------------------------------------------------------------------------
inline void Observable::AddObserver(Observer &observer)
{
  _Observers.insert(&observer);
  observer.StartObserving(this);
}

// -----------------------------------------------------------------------------
inline void Observable::DeleteObserver(Observer &observer)
{
  observer.StopObserving(this);
  _Observers.erase(&observer);
}

// =============================================================================
// Events
// =============================================================================

// -----------------------------------------------------------------------------
inline void Observable::Broadcast(Event event, const void *data)
{
  for (ObserverIterator it = _Observers.begin(); it != _Observers.end(); ++it) {
    (*it)->HandleEvent(this, event, data);
  }
}

// -----------------------------------------------------------------------------
inline void Observable::NotifyObservers(Event event, const void *data)
{
  if (_Changed) {
    Broadcast(event, data);
    Changed(false);
  }
}


} // namespace mirtk

#endif // MIRTK_Observable_H

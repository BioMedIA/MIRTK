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

#ifndef MIRTK_Observer_H
#define MIRTK_Observer_H

#include "mirtk/Object.h"

#include "mirtk/OrderedSet.h"
#include "mirtk/Event.h"


namespace mirtk {


// Type of objects that can be observed
class Observable;


/**
 * Observer of an observable object
 */
class Observer : public Object
{
  mirtkAbstractMacro(Observer);

  // ---------------------------------------------------------------------------
  // Observable/Observer interaction

  friend class Observable;

  void StartObserving(Observable *);
  void StopObserving (Observable *);

  // ---------------------------------------------------------------------------
  // Types

  /// Set of objects observed by this instance
  typedef OrderedSet<Observable *> ObservableSet;

  /// Iterator of observable objects
  typedef ObservableSet::iterator ObservableIterator;
 
  // ---------------------------------------------------------------------------
  // Attributes

  /// Observed objects
  mirtkAttributeMacro(ObservableSet, Observables);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const Observer &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Default constructor
  Observer();

  /// Copy constructor
  Observer(const Observer &);

  /// Assignment operator
  Observer &operator =(const Observer &);

  /// Destructor
  virtual ~Observer();

  // ---------------------------------------------------------------------------
  // Observation
public:

  /// Stop observing any of the currently monitored observables
  void ClearObservables();

  /// Receives event messages from observed objects
  virtual void HandleEvent(Observable *, Event, const void * = NULL) = 0;

};


} // namespace mirtk

#endif // MIRTK_Observer_H

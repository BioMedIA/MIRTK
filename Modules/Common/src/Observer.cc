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

#include "mirtk/Observer.h"
#include "mirtk/Observable.h"


namespace mirtk {


// -----------------------------------------------------------------------------
Observer::Observer()
{
}

// -----------------------------------------------------------------------------
void Observer::CopyAttributes(const Observer &other)
{
  ObservableIterator it = other._Observables.begin();
  while (it != other._Observables.end()) {
    (*it)->AddObserver(*this);
    ++it;
  }
}

// -----------------------------------------------------------------------------
Observer::Observer(const Observer &other)
:
  Object(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
Observer &Observer::operator =(const Observer &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
void Observer::ClearObservables()
{
  while (!_Observables.empty()) {
    Observable *observable = *_Observables.begin();
    // Calls this->StopObserving(observable) which removes the observable
    // and thus would invalidate any iterators to _Observable
    observable->DeleteObserver(*this);
  }
}

// -----------------------------------------------------------------------------
Observer::~Observer()
{
  ClearObservables();
}

// -----------------------------------------------------------------------------
void Observer::StartObserving(Observable *obj)
{
  if (_Observables.find(obj) == _Observables.end()) {
    _Observables.insert(obj);
    this->HandleEvent(obj, RegisteredEvent);
  }
}

// -----------------------------------------------------------------------------
void Observer::StopObserving(Observable *obj)
{
  if (_Observables.find(obj) != _Observables.end()) {
    _Observables.erase(obj);
    this->HandleEvent(obj, UnregisteredEvent);
  }
}


} // namespace mirtk

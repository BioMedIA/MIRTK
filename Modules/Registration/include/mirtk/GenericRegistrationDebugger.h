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

#ifndef MIRTK_GenericRegistrationDebugger_H
#define MIRTK_GenericRegistrationDebugger_H

#include "mirtk/Observer.h"

#include "mirtk/Array.h"


namespace mirtk {


class ImageSimilarity;
class GenericRegistrationFilter;


/**
 * Writes intermediate registration data to the current working directory
 *
 * Usage:
 * \code
 * GenericRegistrationFilter   registration;
 * GenericRegistrationDebugger debugger;
 * registration.AddObserver(debugger);
 * \endcode
 */
class GenericRegistrationDebugger : public Observer
{
  mirtkObjectMacro(GenericRegistrationDebugger);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Prefix for output file names
  mirtkPublicAttributeMacro(string, Prefix);

  /// Whether to use level specific file name prefix
  mirtkPublicAttributeMacro(bool, LevelPrefix);

  /// Current level
  mirtkAttributeMacro(int, Level);

  /// Current iteration
  mirtkAttributeMacro(int, Iteration);

  /// Current line iteration
  mirtkAttributeMacro(int, LineIteration);

  /// Similarity terms
  mirtkAttributeMacro(Array<ImageSimilarity *>, Similarity);

  /// Reference to the registration filter object
  mirtkAggregateMacro(GenericRegistrationFilter, Registration);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
private:

  /// Copy construction
  /// \note Intentionally not implemented.
  GenericRegistrationDebugger(const GenericRegistrationDebugger &);

  /// Assignment operator
  /// \note Intentionally not implemented.
  GenericRegistrationDebugger &operator =(const GenericRegistrationDebugger &);

public:

  /// Constructor
  GenericRegistrationDebugger(const char * = "");

  /// Destructor
  ~GenericRegistrationDebugger();

  /// Handle event and print message to output stream
  void HandleEvent(Observable *, Event, const void *);

};


} // namespace mirtk

#endif // MIRTK_GenericRegistrationDebugger_H

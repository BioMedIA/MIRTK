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

#ifndef MIRTK_Configurable_H
#define MIRTK_Configurable_H

#include "mirtk/Observable.h"

#include "mirtk/Array.h"


namespace mirtk {


/**
 * Base class of observable objects which have a primary name and optionally
 * one or more alternative names
 */
class Configurable : public Observable
{
  mirtkAbstractMacro(Configurable);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Primary name of object
  mirtkPublicAttributeMacro(string, Name);

  /// Parameter name prefixes recognized by Set
  mirtkPublicAttributeMacro(Array<string>, ParameterPrefix);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const Configurable &other);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  Configurable(const char * = "");

  /// Copy constructor
  Configurable(const Configurable &);

  /// Assignment operator
  Configurable &operator =(const Configurable &);

public:

  /// Destructor
  virtual ~Configurable();

  /// Default name derived from class name
  ///
  /// Example usage:
  /// \code
  /// irtkConfigurable *obj = new MyConfigurableObject();
  /// obj->Name(obj->DefaultName());
  /// \endcode
  string DefaultName() const;

  // ---------------------------------------------------------------------------
  // Parameters

protected:

  /// Whether this object has an explicit name
  bool HasName() const;

  /// Whether this object has either an explicit name or default prefix
  bool HasPrefix() const;

  /// Get default object name prefix (if any)
  string DefaultPrefix() const;

  /// Get name of parameter with default object name prefix
  string ParameterNameWithPrefix(const string &) const;

  /// Get name of parameter with default object name prefix
  string ParameterNameWithPrefix(const char *) const;

  /// Get name of parameter without object name prefix
  string ParameterNameWithoutPrefix(const char *) const;

  /// Insert parameter into name/value list with object name prefix
  template <class T> bool InsertWithPrefix(ParameterList &, string, T) const;

  /// Insert parameters into name/value list with object name prefix
  bool InsertWithPrefix(ParameterList &, const ParameterList &) const;

  /// Set parameter value from string
  virtual bool SetWithPrefix(const char *, const char *);

  /// Set parameter value from string
  virtual bool SetWithoutPrefix(const char *, const char *);

public:

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <class T>
bool Configurable::InsertWithPrefix(ParameterList &params, string name, T value) const
{
  if (name.empty() || !HasPrefix()) return false;
  Insert(params, ParameterNameWithPrefix(name), value);
  return true;
}

// -----------------------------------------------------------------------------
inline bool Configurable::InsertWithPrefix(ParameterList &params, const ParameterList &other) const
{
  string prefix = DefaultPrefix();
  if (prefix.empty()) return false;
  if (prefix.back() == ' ') prefix.back() = '\0';
  Insert(params, other, prefix.c_str());
  return true;
}


} // namespace mirtk

#endif // MIRTK_Configurable_H

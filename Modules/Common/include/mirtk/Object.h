/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#ifndef MIRTK_Object_H
#define MIRTK_Object_H

#include "mirtk/Array.h"
#include "mirtk/Exception.h"
#include "mirtk/Pair.h"
#include "mirtk/String.h"
#include "mirtk/Stream.h"


namespace mirtk {


// =============================================================================
// Basic object interface
// =============================================================================

/// Ordered list of parameter name/value pairs
typedef Array<Pair<string, string> >    ParameterList;
typedef ParameterList::iterator         ParameterIterator;
typedef ParameterList::const_iterator   ParameterConstIterator;


/**
 * Base class for all MIRTK object classes
 *
 * \note This base class must be a virtual interface without any own data
 *       members to not change the type size of subclasses! Derive another
 *       intermediate abstract base class from it to add data members shared
 *       by a more specific class of objects.
 */
class Object
{
public:

  /// Get name of this class type
  static const char *NameOfType();

  /// Get name of class, which this object is an instance of
  virtual const char *NameOfClass() const = 0;

  /// Destructor
  virtual ~Object();

  /// Set parameter value from string
  ///
  /// \param[in] name  Parameter name.
  /// \param[in] value Parameter value.
  ///
  /// \returns Whether the named parameter is valid and a valid value was
  ///          set from the provided string representation of the value.
  virtual bool Set(const char *name, const char *value);

  /// Get parameter name/value pairs
  virtual ParameterList Parameter() const;

  /// Set parameters from name/value pairs
  bool Parameter(const ParameterList &);

  // ---------------------------------------------------------------------------
  // Error handling
protected:

  /// Raise error in member function
  ///
  /// The current implementation prints the error message to STDERR and terminates
  /// the program with exit code 1. In future releases, when all library code has
  /// been rewritten to use this function, a suitable runtime exception may be
  /// thrown instead.
  ///
  /// \param[in] err  Error type. Unused at the moment, but may be used in
  ///                 future releases to throw the appropriate exception type.
  /// \param[in] cls  Name of class that defines the member function.
  /// \param[in] func Name of member function which is throwing the error (i.e., __func__).
  /// \param[in] args Error message. The given arguments are converted to strings
  ///                 using the ToString template function. These strings are then
  ///                 concatenated to produce the complete error message.
  template <typename... Args>
  static void ThrowStatic(ErrorType err, const char *cls, const char *func, Args... args);

  /// Raise error in member function
  ///
  /// The current implementation prints the error message to STDERR and terminates
  /// the program with exit code 1. In future releases, when all library code has
  /// been rewritten to use this function, a suitable runtime exception may be
  /// thrown instead.
  ///
  /// \param[in] err  Error type. Unused at the moment, but may be used in
  ///                 future releases to throw the appropriate exception type.
  /// \param[in] func Name of member function which is throwing the error (i.e., __func__).
  /// \param[in] args Error message. The given arguments are converted to strings
  ///                 using the ToString template function. These strings are then
  ///                 concatenated to produce the complete error message.
  template <typename... Args>
  void Throw(ErrorType err, const char *func, Args... args) const;
};

// =============================================================================
// Inline definitions
// =============================================================================

// -----------------------------------------------------------------------------
inline Object::~Object()
{
}

// -----------------------------------------------------------------------------
inline const char *Object::NameOfType()
{
  return "mirtk::Object";
}

// -----------------------------------------------------------------------------
inline bool Object::Set(const char *, const char *)
{
  return false;
}

// -----------------------------------------------------------------------------
inline ParameterList Object::Parameter() const
{
  return ParameterList();
}

// -----------------------------------------------------------------------------
inline bool Object::Parameter(const ParameterList &param)
{
  bool ok = true;
  for (ParameterConstIterator it = param.begin(); it != param.end(); ++it) {
    ok = this->Set(it->first.c_str(), it->second.c_str()) && ok;
  }
  return ok;
}

// -----------------------------------------------------------------------------
template <typename... Args>
void Object::ThrowStatic(ErrorType err, const char *cls, const char *func, Args... args)
{
  string member(cls);
  member += "::";
  member += func;
  mirtk::Throw(err, member.c_str(), args...);
}

// -----------------------------------------------------------------------------
template <typename... Args>
void Object::Throw(ErrorType err, const char *func, Args... args) const
{
  ThrowStatic(err, this->NameOfClass(), func, args...);
}

// =============================================================================
// Auxiliary functions for subclass implementation
// =============================================================================

// -----------------------------------------------------------------------------
/// Find parameter in parameters list
inline ParameterConstIterator Find(const ParameterList &params, string name)
{
  ParameterConstIterator it = params.begin();
  while (it != params.end() && it->first != name) ++it;
  return it;
}

// -----------------------------------------------------------------------------
/// Find parameter in parameters list
inline ParameterIterator Find(ParameterList &params, string name)
{
  ParameterIterator it = params.begin();
  while (it != params.end() && it->first != name) ++it;
  return it;
}

// -----------------------------------------------------------------------------
/// Whether parameter is in parameters list
inline bool Contains(const ParameterList &params, string name)
{
  return Find(params, name) != params.end();
}

// -----------------------------------------------------------------------------
/// Get parameter value from parameters list
inline string Get(const ParameterList &params, string name)
{
  ParameterConstIterator pos = Find(params, name);
  if (pos == params.end()) return string("");
  return pos->second;
}

// -----------------------------------------------------------------------------
/// Insert/replace value into/in parameters list
template <class T>
inline ParameterList &Insert(ParameterList &params, string name, T value)
{
  ParameterIterator pos = Find(params, name);
  if (pos == params.end()) {
    params.push_back(make_pair(name, ToString(value)));
  } else {
    pos->second = ToString(value);
  }
  return params;
}

// -----------------------------------------------------------------------------
/// Insert/replace string value into/in parameters list
template <>
inline ParameterList &Insert(ParameterList &params, string name, const char *value)
{
  ParameterIterator pos = Find(params, name);
  if (pos == params.end()) {
    params.push_back(make_pair(name, string(value)));
  } else {
    pos->second = value;
  }
  return params;
}

// -----------------------------------------------------------------------------
/// Insert/replace string value into/in parameters list
template <>
inline ParameterList &Insert(ParameterList &params, string name, string value)
{
  ParameterIterator pos = Find(params, name);
  if (pos == params.end()) {
    params.push_back(make_pair(name, value));
  } else {
    pos->second = value;
  }
  return params;
}

// -----------------------------------------------------------------------------
/// Insert/replace values into/in parameters list
inline ParameterList &Insert(ParameterList       &params,
                             const ParameterList &other,
                             const char          *prefix = NULL)
{
  if (prefix) {
    string name;
    for (ParameterConstIterator it = other.begin(); it != other.end(); ++it) {
      name    = it->first;
      name[0] = ::tolower(name[0]);
      Insert(params, string(prefix) + " " + name, it->second);
    }
  } else {
    for (ParameterConstIterator it = other.begin(); it != other.end(); ++it) {
      Insert(params, it->first, it->second);
    }
  }
  return params;
}

// -----------------------------------------------------------------------------
/// Remove parameter from parameters list
inline ParameterList &Remove(ParameterList &params, string name)
{
  ParameterIterator pos = Find(params, name);
  if (pos != params.end()) params.erase(pos);
  return params;
}

// =============================================================================
// Auxiliary macros for subclass implementation
// =============================================================================

// -----------------------------------------------------------------------------
/// Declare abstract base class derived from Object
#define mirtkAbstractMacro(name)                                               \
  public:                                                                      \
    /** Get name of this class type */                                         \
    inline static const char *NameOfType() { return #name; }                   \
    /** Get name of class, which this object is an instance of */              \
    virtual const char *NameOfClass() const = 0;                               \
  private:

// -----------------------------------------------------------------------------
/// Declare class derived from Object
#define mirtkObjectMacro(name)                                                 \
  public:                                                                      \
    /** Get name of this class type */                                         \
    inline static  const char *NameOfType() { return #name; }                  \
    /** Get name of class, which this object is an instance of */              \
    inline virtual const char *NameOfClass() const { return #name; }           \
  private:

// -----------------------------------------------------------------------------
/// Declare class of mutable objects, i.e., ones which define their own
/// NameOfClass implementation that returns a different type identifier
/// depending on the state of the object.
#define mirtkMutableObjectMacro(name)                                          \
  public:                                                                      \
    /** Get name of this class type */                                         \
    inline static const char *NameOfType() { return #name; }                   \
    /** Get name of class, which this object is an instance of */              \
    virtual const char *NameOfClass() const;                                   \
  private:

// -----------------------------------------------------------------------------
/// Define setter for class member variable
/// \sa mirtkPublicAttributeMacro
#define mirtkSetMacro(name, type)                                              \
    virtual void Set##name(type arg) { this->_##name = arg; }

// -----------------------------------------------------------------------------
/// Define getter for class member variable
/// \sa mirtkPublicAttributeMacro, mirtkReadOnlyAttributeMacro
#define mirtkGetMacro(name, type)                                              \
    type Get##name() const { return this->_##name; }

// -----------------------------------------------------------------------------
/// Define VTK-like On/Off setter for boolean class member variable
/// \sa mirtkPublicAttributeMacro
#define mirtkOnOffMacro(name)                                                  \
    virtual void name##On()  { this->_##name = true;  }                        \
    virtual void name##Off() { this->_##name = false; }

// -----------------------------------------------------------------------------
/// Define read-only class attribute and corresponding accessors
#define mirtkDefineReadOnlyAttributeMacro(typedecl, type, name)                \
  protected:                                                                   \
    typedecl _##name;                                                          \
  public:                                                                      \
    /** Get value of _##name attribute */                                      \
    inline type &name() { return _##name; }                                    \
    /** Get value of _##name attribute */                                      \
    inline const type &name() const { return _##name; }                        \
  private:

// -----------------------------------------------------------------------------
/// Define class attribute and corresponding accessors
#define mirtkDefineAttributeMacro(access, typedecl, type, name)                \
  protected:                                                                   \
    typedecl _##name;                                                          \
  access:                                                                      \
    /** Set value of _##name attribute */                                      \
    inline virtual void name(type arg) { _##name = arg; }                      \
    /** Get value of _##name attribute */                                      \
    inline type &name() { return _##name; }                                    \
    /** Get value of _##name attribute */                                      \
    inline const type &name() const { return _##name; }                        \
  private:

// -----------------------------------------------------------------------------
/// Define public attribute
#define mirtkPublicAttributeMacro(type, name)                                  \
  mirtkDefineAttributeMacro(public, type, type, name)
/// Define public read-only attribute
#define mirtkReadOnlyAttributeMacro(type, name)                                \
  mirtkDefineReadOnlyAttributeMacro(type, type, name)
/// Define public mutable attribute
#define mirtkPublicMutableAttributeMacro(type, name)                           \
  mirtkDefineAttributeMacro(public, mutable type, type, name)

// -----------------------------------------------------------------------------
/// Define protected attribute
#define mirtkAttributeMacro(type, name)                                        \
  mirtkDefineAttributeMacro(protected, type, type, name)
/// Define protected mutable attribute
#define mirtkMutableAttributeMacro(type, name)                                 \
  mirtkDefineAttributeMacro(protected, mutable type, type, name)

// -----------------------------------------------------------------------------
/// Define pointer to aggregate (cf. UML aggregation) and corresponding accessors
#define mirtkDefineAggregateMacro(access, type, name)                          \
  protected:                                                                   \
    type *_##name;                                                             \
  access:                                                                      \
    /** Set value of _##name attribute */                                      \
    inline virtual void name(type *arg) { _##name = arg; }                     \
    /** Get value of _##name attribute */                                      \
    inline type *name() const { return _##name; }                              \
  private:

// -----------------------------------------------------------------------------
/// Define pointer to aggregate (cf. UML aggregation) and corresponding accessors
#define mirtkDefineReadOnlyAggregateMacro(access, type, name)                  \
  protected:                                                                   \
    type *_##name;                                                             \
  access:                                                                      \
    /** Get value of _##name attribute */                                      \
    inline type *name() const { return _##name; }                              \
  private:

// -----------------------------------------------------------------------------
/// Define public pointer to aggregate (cf. UML aggregation)
#define mirtkPublicAggregateMacro(type, name)                                  \
  mirtkDefineAggregateMacro(public, type, name)
/// Define public read-only pointer to aggregate (cf. UML aggregation)
#define mirtkReadOnlyAggregateMacro(type, name)                                \
  mirtkDefineReadOnlyAggregateMacro(public, type, name)

// -----------------------------------------------------------------------------
/// Define protected pointer to aggregate (cf. UML aggregation)
#define mirtkAggregateMacro(type, name)                                        \
  mirtkDefineAggregateMacro(protected, type, name)

// -----------------------------------------------------------------------------
/// Define pointer to component (cf. UML composition) and corresponding accessors
#define mirtkDefineComponentMacro(access, type, name)                          \
  protected:                                                                   \
    type *_##name;                                                             \
  access:                                                                      \
    /** Set pointer to _##name attribute */                                    \
    inline virtual void name(type *arg) { delete _##name; _##name = arg; }     \
    /** Get pointer to _##name attribute */                                    \
    inline type *name() { return _##name; }                                    \
    /** Get const pointer to _##name attribute */                              \
    inline const type *name() const { return _##name; }                        \
  private:

// -----------------------------------------------------------------------------
/// Define public pointer to component (cf. UML composition)
#define mirtkPublicComponentMacro(type, name)                                  \
  mirtkDefineComponentMacro(public, type, name)
/// Define public read-only pointer to component (cf. UML composition)
#define mirtkReadOnlyComponentMacro(type, name)                                \
  mirtkDefineReadOnlyAggregateMacro(public, type, name)

// -----------------------------------------------------------------------------
/// Define protected pointer to component (cf. UML composition)
#define mirtkComponentMacro(type, name)                                        \
  mirtkDefineComponentMacro(protected, type, name)

// -----------------------------------------------------------------------------
/// Define pointer to component (cf. UML composition) which can also be set to
/// an externally managed instance instead using the corresponding accessors
///
/// A loose component is a pointer that must never by a nullptr. It must be
/// initialized during construction of an instance of the class which this
/// component attribute belongs.
///
/// Unlike a regular component, a loose component can alternatively be set to an
/// externally instantiated object of the component type. In this case, the class
/// revokes ownership of the object which must be deleted by the client code.
/// When a nullptr is given as argument to the setter of the loose component,
/// the current pointer is replaced by a new copy that is owned by this class.
/// It can thus be used to release the external object again before destruction
/// in cases where the class that the loose component belongs to remains in use.
#define mirtkLooseComponentMacro(type, name)                                   \
  protected:                                                                   \
    type *_##name;                                                             \
    bool  _##name##Owner;                                                      \
  public:                                                                      \
    /** Set pointer to _##name attribute */                                    \
    inline virtual void name(type *arg) {                                      \
      if (_##name##Owner) delete _##name;                                      \
      if (arg) {                                                               \
        _##name        = arg;                                                  \
        _##name##Owner = false;                                                \
      } else {                                                                 \
        _##name        = new type(*_##name);                                   \
        _##name##Owner = true;                                                 \
      }                                                                        \
    }                                                                          \
    /** Get pointer to _##name attribute */                                    \
    inline type *name() { return _##name; }                                    \
    /** Get const pointer to _##name attribute */                              \
    inline const type *name() const { return _##name; }                        \
  private:


} // namespace mirtk

#endif // MIRTK_Object_H

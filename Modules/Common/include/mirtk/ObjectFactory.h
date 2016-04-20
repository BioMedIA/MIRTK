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

#ifndef MIRTK_ObjectFactory_H
#define MIRTK_ObjectFactory_H

#include "mirtk/UnorderedMap.h"
#include "mirtk/Stream.h"

#include "mirtk/CommonExport.h"


namespace mirtk {


// Global "debug" flag (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int debug;


/**
 * Constructs a new instance of a class derived from the template type
 * given a type identification such as an enumeration value or string
 */
template <typename TId, typename TInterface>
class ObjectFactory
{
  // ---------------------------------------------------------------------------
  // Types
public:

  /// Type of object type identifier
  typedef TId Id;

  /// Type of object base class
  typedef TInterface Interface;

  /// Type of object creator
  typedef Interface *(*Creator)();

  /// Type of object factory instance
  typedef ObjectFactory<Id, Interface> InstanceType;

  // ---------------------------------------------------------------------------
  // Abstract base class for singleton
protected:

  /// Constructor
  ObjectFactory() {}

  /// Destructor
  virtual ~ObjectFactory() {}

  /// Copy constructor. Intentionally not implemented.
  ObjectFactory(const ObjectFactory &);

  /// Assignment operator. Intentionally not implemented.
  void operator =(const ObjectFactory &);

public:

  // ---------------------------------------------------------------------------
  // Object creation
private:

  /// Type of associative map
  typedef UnorderedMap<Id, Creator> TypeIdAssociations;

  /// Type of associative map
  typedef UnorderedMap<const char *, Creator> TypeNameAssociations;

  /// Registered object type creators
  TypeIdAssociations _TypeIdAssociations;

  /// Registered object type creators
  TypeNameAssociations _TypeNameAssociations;

public:

  /// Register new object creator
  ///
  /// \param[in] type_id    Object type identifier.
  /// \param[in] type_name  Object type name.
  /// \param[in] creator    Object creator function.
  bool Register(Id type_id, const char *type_name, Creator creator)
  {
    _TypeIdAssociations  [type_id  ] = creator;
    _TypeNameAssociations[type_name] = creator;
    if (debug) {
      cout << "Registered object type with name " << type_name << " and ID " << type_id << endl;
    }
    return true;
  }

  /// Create new object
  ///
  /// \param[in] type_id Object type identifier.
  ///
  /// \returns Object of given type which must be deleted by caller or nullptr.
  Interface *New(Id type_id) const
  {
    auto it = _TypeIdAssociations.find(type_id);
    return (it == _TypeIdAssociations.end() ? nullptr : (it->second)());
  }

  /// Create new object
  ///
  /// \param[in] type_name Object type name.
  ///
  /// \returns Object of given type which must be deleted by caller or nullptr.
  Interface *New(const char *type_name) const
  {
    auto it = _TypeNameAssociations.find(type_name);
    return (it == _TypeNameAssociations.end() ? nullptr : (it->second)());
  }

};

// -----------------------------------------------------------------------------
/// Default object creation function
template <class BaseType, class ObjectType> BaseType *New()
{
  return new ObjectType();
}

// -----------------------------------------------------------------------------
/// Define singleton object factory class in .cc file of object base class
#define mirtkDefineObjectFactory(id_type, base)                                \
    class base##Factory : public mirtk::ObjectFactory<id_type, base>           \
    {                                                                          \
    private:                                                                   \
      /** Constructor */                                                       \
      base##Factory() {}                                                       \
      /** Destructor */                                                        \
      virtual ~base##Factory() {}                                              \
      /** Copy constructor. Intentionally not implemented. */                  \
      base##Factory(const base##Factory &);                                    \
      /** Assignment operator. Intentionally not implemented. */               \
      void operator =(const base##Factory &);                                  \
    public:                                                                    \
      /** Get singleton instance */                                            \
      static base##Factory &Instance()                                         \
      {                                                                        \
        static base##Factory instance;                                         \
        return instance;                                                       \
      }                                                                        \
    }

// -----------------------------------------------------------------------------
/// Register object type with factory at static initialization time
#ifdef MIRTK_AUTO_REGISTER
  #define mirtkAutoRegisterObjectTypeMacro(factory, id_type, id, base, type)   \
    namespace {                                                                \
      static auto _##type##Registered =                                        \
        (factory).Register(id, type::NameOfType(), mirtk::New<base, type>);    \
    }
#else // MIRTK_AUTO_REGISTER
  #define mirtkAutoRegisterObjectTypeMacro(factory, id_type, id, base, type)
#endif // MIRTK_AUTO_REGISTER


} // namespace mirtk

#endif // MIRTK_ObjectFactory_H

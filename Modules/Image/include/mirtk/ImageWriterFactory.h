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

#ifndef MIRTK_ImageWriterFactory_H
#define MIRTK_ImageWriterFactory_H

#include "mirtk/ObjectFactory.h" // New<BaseType, ObjectType>()
#include "mirtk/ImageWriter.h"

#include "mirtk/Array.h"
#include "mirtk/UnorderedMap.h"



namespace mirtk {


/**
 * Factory for instantiation of image writers
 */
class ImageWriterFactory
{
  // ---------------------------------------------------------------------------
  // Types
public:

  /// Type of object creator
  typedef ImageWriter *(*ImageWriterCreator)();

  // ---------------------------------------------------------------------------
  // Singleton
private:

  /// Constructor
  ImageWriterFactory();

  /// Destructor
  ~ImageWriterFactory();

  /// Copy constructor. Intentionally not implemented.
  ImageWriterFactory(const ImageWriterFactory &);

  /// Assignment operator. Intentionally not implemented.
  void operator =(const ImageWriterFactory &);

public:

  /// Singleton instance
  /// \attention This function is not thread-safe!
  static ImageWriterFactory &Instance();


private:

  /// Type of associative map
  typedef UnorderedMap<string, ImageWriterCreator> Associations;

  /// Associates file name extensions with image writer creators
  Associations _Associations;

public:

  /// Register new image writer
  ///
  /// \param[in] exts    File name extensions.
  /// \param[in] creator Image writer instantiation function.
  bool Register(const Array<string> &exts, ImageWriterCreator creator);

  /// Construct new image writer for given output image file name
  ///
  /// \returns Image writer which is able to write an image in the format
  ///          corresponding to the given file name extension or nullptr
  ///          when no such writer is registered.
  ImageWriter *New(const char *fname) const;

};

// -----------------------------------------------------------------------------
/// Register image reader with factory singleton
#define mirtkRegisterImageWriterMacro(type)                                    \
  mirtk::ImageWriterFactory::Instance()                                        \
      .Register(type::Extensions(), mirtk::New<mirtk::ImageWriter, type>)

// -----------------------------------------------------------------------------
/// Register image reader with factory singleton at static initialization time
#ifdef MIRTK_AUTO_REGISTER
  #define mirtkAutoRegisterImageWriterMacro(type)                              \
    namespace {                                                                \
      static auto _##type##Registered =                                        \
        mirtk::ImageWriterFactory::Instance()                                  \
          .Register(type::Extensions(), mirtk::New<mirtk::ImageWriter, type>); \
    }
#else // MIRTK_AUTO_REGISTER
  #define mirtkAutoRegisterImageWriterMacro(type)
#endif // MIRTK_AUTO_REGISTER



} // namespace mirtk

#endif // MIRTK_ImageWriterFactory_H

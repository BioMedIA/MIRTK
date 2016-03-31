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

#ifndef MIRTK_ImageReaderFactory_H
#define MIRTK_ImageReaderFactory_H

#include "mirtk/ObjectFactory.h" // New<BaseType, ObjectType>()
#include "mirtk/ImageReader.h"
#include "mirtk/List.h"


namespace mirtk {


/**
 * Factory for instantiation of image readers
 */
class ImageReaderFactory
{
  // ---------------------------------------------------------------------------
  // Types
public:

  /// Type of object creator
  typedef ImageReader *(*ImageReaderCreator)();

  // ---------------------------------------------------------------------------
  // Singleton
private:

  /// Constructor
  ImageReaderFactory();

  /// Destructor
  ~ImageReaderFactory();

  /// Copy constructor. Intentionally not implemented.
  ImageReaderFactory(const ImageReaderFactory &);

  /// Assignment operator. Intentionally not implemented.
  void operator =(const ImageReaderFactory &);

public:

  /// Singleton instance
  /// \attention This function is not thread-safe!
  static ImageReaderFactory &Instance();

  // ---------------------------------------------------------------------------
  // Object creation
private:

  /// Type of associative map
  typedef List<ImageReaderCreator> Creators;

  /// Registered object type creators
  Creators _Creators;

public:

  /// Register new object creator
  ///
  /// \param[in] creator Object creator function.
  bool Register(ImageReaderCreator creator);

  /// Construct new image reader for given image file
  ///
  /// \param[in] fname Input image file path/name incl. file name extension.
  ///
  /// \returns First found image reader which is able to read the given image file
  ///          or nullptr when file cannot be read by any registered reader.
  ImageReader *New(const char *fname) const;

};

// -----------------------------------------------------------------------------
/// Register image reader with factory singleton
#define mirtkRegisterImageReaderMacro(type)                                    \
  mirtk::ImageReaderFactory::Instance()                                        \
      .Register(mirtk::New<mirtk::ImageReader, type>)

// -----------------------------------------------------------------------------
/// Register image reader with factory singleton at static initialization time
#ifdef MIRTK_AUTO_REGISTER
  #define mirtkAutoRegisterImageReaderMacro(type)                              \
    namespace {                                                                \
      static auto _##type##Registered =                                        \
        mirtk::ImageReaderFactory::Instance()                                  \
          .Register(mirtk::New<mirtk::ImageReader, type>);                     \
    }
#else // MIRTK_AUTO_REGISTER
  #define mirtkAutoRegisterImageReaderMacro(type)
#endif // MIRTK_AUTO_REGISTER


} // namespace mirtk

#endif // MIRTK_ImageReaderFactory_H

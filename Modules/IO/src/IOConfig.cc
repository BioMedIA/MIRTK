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

#include "mirtk/IOConfig.h"

#if !defined(MIRTK_AUTO_REGISTER)
  #include "mirtk/ImageReaderFactory.h"
  #include "mirtk/ImageWriterFactory.h"
  #include "mirtk/GIPLImageReader.h"
  #include "mirtk/GIPLImageWriter.h"
  #include "mirtk/PGMImageReader.h"
  #include "mirtk/PGMImageWriter.h"
  #include "mirtk/PNGImageWriter.h"
  #if MIRTK_IO_WITH_NIfTI
    #include "mirtk/NiftiImageReader.h"
    #include "mirtk/NiftiImageWriter.h"
  #endif
#endif


namespace mirtk {



// -----------------------------------------------------------------------------
static void RegisterImageReaders()
{
  #ifndef MIRTK_AUTO_REGISTER
    mirtkRegisterImageReaderMacro(GIPLImageReader);
    mirtkRegisterImageReaderMacro(PGMImageReader);
    #if MIRTK_IO_WITH_NIfTI
      mirtkRegisterImageReaderMacro(NiftiImageReader);
    #endif
  #endif
}

// -----------------------------------------------------------------------------
static void RegisterImageWriters()
{
  #ifndef MIRTK_AUTO_REGISTER
    mirtkRegisterImageWriterMacro(GIPLImageWriter);
    mirtkRegisterImageWriterMacro(PGMImageWriter);
    #if MIRTK_IO_WITH_PNG
      mirtkRegisterImageWriterMacro(PNGImageWriter);
    #endif
    #if MIRTK_IO_WITH_NIfTI
      mirtkRegisterImageWriterMacro(NiftiImageWriter);
    #endif
  #endif
}

// -----------------------------------------------------------------------------
void InitializeIOLibrary()
{
  static bool initialized = false;
  if (!initialized) {
    RegisterImageReaders();
    RegisterImageWriters();
    initialized = true;
  }
}


} // namespace mirtk

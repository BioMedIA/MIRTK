/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2017 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
 * Copyright 2015-2017 Andreas Schuh
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

#include "mirtk/Config.h"

#include "mirtk/PNGImageWriter.h"
#include "png.h"

#include "mirtk/ImageWriterFactory.h"


namespace mirtk {


// Register image reader with object factory during static initialization
mirtkAutoRegisterImageWriterMacro(PNGImageWriter);


// -----------------------------------------------------------------------------
Array<string> PNGImageWriter::Extensions()
{
  Array<string> exts(1);
  exts[0] = ".png";
  return exts;
}

// -----------------------------------------------------------------------------
void PNGImageWriter::Initialize()
{
  // Initialize base class
  ImageWriter::Initialize();

  // Check input
  if (_Input->Z() != 1 && _Input->Z() != 3) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Supports only images with z=1 (grayscale) or z=3 (RGB)");
  }
}

// -----------------------------------------------------------------------------
void PNGImageWriter::Run()
{
  // Initialize filter
  this->Initialize();

  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, (png_voidp)nullptr, nullptr, nullptr);
  if (!png_ptr) {
    Throw(ERR_RuntimeError, __FUNCTION__, "Unable to write PNG file");
  }

  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    png_destroy_write_struct(&png_ptr, (png_infopp)nullptr);
    Throw(ERR_RuntimeError, __FUNCTION__, "Unable to write PNG file");
  }

  if (setjmp(png_jmpbuf(png_ptr)) != 0) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    Throw(ERR_RuntimeError, __FUNCTION__, "Unable to write PNG file");
  }

  // Open file
#ifdef WINDOWS
  FILE *fp;
  errno_t err = fopen_s(&fp, _FileName.c_str(), "wb");
  if (err != 0) fp = nullptr;
#else
  FILE *fp = fopen(_FileName.c_str(), "wb");
#endif
  if (fp == nullptr) {
    Throw(ERR_RuntimeError, __FUNCTION__, "Failed to open PNG file for writing: ", _FileName);
  }

  // Initialize PNG I/O
  png_init_io(png_ptr, fp);

  // Initialize header, and convert and write image data
  UniquePtr<png_byte[]> data;
  UniquePtr<png_bytep[]> image;
  data.reset(new png_byte[_Input->X() * _Input->Y() * _Input->Z()]);
  image.reset(new png_bytep[_Input->Y()]);

  double min_val, max_val;
  _Input->GetMinMaxAsDouble(min_val, max_val);
  double slope = 1., intercept = 0.;
  if (min_val < 0. || max_val > 255.) {
    slope = 255. / (max_val - min_val);
    intercept = - slope * min_val;
  }
  if (_Input->Z() == 1) {
    png_set_IHDR(png_ptr, info_ptr, _Input->X(), _Input->Y(),
                 8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);
    for (int j = 0; j < _Input->Y(); ++j) {
      png_bytep ptr = &data[j * _Input->X()];
      image[_Input->Y() - j - 1] = ptr;
      for (int i = 0; i < _Input->X(); ++i, ++ptr) {
        *ptr = voxel_cast<png_byte>(slope * _Input->GetAsDouble(i, j) + intercept);
      }
    }
  } else if (_Input->Z() == 3) {
    png_set_IHDR(png_ptr, info_ptr, _Input->X(), _Input->Y(),
                 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);
    for (int j = 0; j < _Input->Y(); ++j) {
      png_bytep ptr = &data[j * _Input->X()];
      image[_Input->Y() - j - 1] = ptr;
      for (int i = 0; i < _Input->X(); ++i) {
        *ptr++ = voxel_cast<png_byte>(slope * _Input->GetAsDouble(i, j, 0) + intercept);
        *ptr++ = voxel_cast<png_byte>(slope * _Input->GetAsDouble(i, j, 1) + intercept);
        *ptr++ = voxel_cast<png_byte>(slope * _Input->GetAsDouble(i, j, 2) + intercept);
      }
    }
  } else {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Unsupported input image type/size");
  }
  png_write_info(png_ptr, info_ptr);
  png_write_image(png_ptr, image.get());
  png_write_end(png_ptr, info_ptr);

  // Destroy PNG data
  png_destroy_write_struct(&png_ptr, &info_ptr);

  // Close file
  fclose(fp);

  // Finalize filter
  this->Finalize();
}


} // namespace mirtk

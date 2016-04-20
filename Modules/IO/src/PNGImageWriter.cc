/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
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
  if (_Input->Z() != 3) {
    cerr << this->NameOfClass() << " supports only images with (z = 3, e.g. R-G-B)" << endl;
    exit(1);
  }
  if (dynamic_cast<const GenericImage<unsigned char> *>(_Input) == NULL) {
    cerr << this->NameOfClass() << " supports only images of voxel type unsigned char" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void PNGImageWriter::Run()
{
  // Initialize filter
  this->Initialize();

  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, (png_voidp)NULL, NULL, NULL);
  if (!png_ptr) {
    cerr << this->NameOfClass() << "::Run: Unable to write PNG file!" << endl;
    exit(1);
  }

  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    cerr << this->NameOfClass() << "::Run: Unable to write PNG file!" << endl;;
    png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
    exit(1);
  }

  if (setjmp(png_jmpbuf(png_ptr)) != 0) {
    cerr << this->NameOfClass() << "::Run: Unable to write PNG file!" << endl;;
    png_destroy_write_struct(&png_ptr, &info_ptr);
    exit(1);
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
    cerr << this->NameOfClass() << "::Run: Failed to open PNG file for writing: " << _FileName << endl;
    exit(1);
  }

  // Initialize PNG I/O
  png_init_io(png_ptr, fp);

  // Initialize header
  png_set_IHDR(png_ptr, info_ptr, _Input->X(), _Input->Y(),
               8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_DEFAULT,
               PNG_FILTER_TYPE_DEFAULT);

  // Write header
  png_write_info(png_ptr, info_ptr);

  // Copy image
  png_byte *data = new png_byte  [3*_Input->X()*_Input->Y()];
  png_byte **ptr = new png_byte *[  _Input->Y()];
  png_byte *ptr2data = data;

  const GenericImage<unsigned char> *image;
  image = dynamic_cast<const GenericImage<unsigned char> *>(_Input);
  for (int j = 0; j < image->Y(); ++j) {
    for (int i = 0; i < image->X(); ++i) {
      *ptr2data++ = image->Get(i, j, 0, 0);
      *ptr2data++ = image->Get(i, j, 1, 0);
      *ptr2data++ = image->Get(i, j, 2, 0);
    }
    ptr[_Input->Y() - j - 1] = &(data[_Input->X() * j * 3]);
  }
  png_write_image(png_ptr, ptr);
  png_write_end(png_ptr, info_ptr);

  // Delete pointers
  delete[] data;
  delete[] ptr;

  // Destroy PNG data
  png_destroy_write_struct(&png_ptr, &info_ptr);

  // Close file
  fclose(fp);

  // Finalize filter
  this->Finalize();
}


} // namespace mirtk

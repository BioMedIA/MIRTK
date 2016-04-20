/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2016 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_NiftiImage_H
#define MIRTK_NiftiImage_H

#include "nifti/nifti2_io.h"

#ifndef LSB_FIRST
#  define LSB_FIRST 1
#endif
#ifndef MSB_FIRST
#  define MSB_FIRST 2
#endif

#include "mirtk/Matrix.h"


namespace mirtk {


/**
 * The NIfTI-1 image class.
 *
 * This is a wrapper around the nifti_image struct.
 */
class NiftiImage
{
public:

  enum DataOrder { RADIOLOGICAL = -1, INCONSISTENT = 0, NEUROLOGICAL = 1 };

  /// The "NIFTI-1" image storage struct.
  nifti_image *nim;

  /// Constructor
  NiftiImage(const char *fname = NULL);

  /// Destructor
  ~NiftiImage();

  /// Read header
  void Read(const char *);

  // Initialize header with minimal set of fields and orientation info.
  void Initialize(int, int, int, int,
                  double, double, double, double,
                  int, Matrix &qmat, Matrix *smat = 0,
                  double torigin = 0, const void *data = NULL);

  // Initialize header with minimal set of fields and orientation info.
  void Initialize(int, int, int, int, int,
                  double, double, double, double,
                  int, Matrix &qmat, Matrix *smat = 0,
                  double torigin = 0, const void *data = NULL);

  /// Print header (for debugging purposes)
  void Print();

};


} // namespace mirtk

#endif // MIRTK_NiftiImage_H

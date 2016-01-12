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

#ifndef MIRTK_GIPL_H
#define MIRTK_GIPL_H


// GIPL magic number
#define GIPL_MAGIC1         719555000u
#define GIPL_MAGIC2        4026526128u
#define GIPL_MAGIC_EXT            815u

// GIPL header size
#define GIPL_HEADERSIZE   256

// GIPL filter types
#define GIPL_BINARY       1
#define GIPL_CHAR         7
#define GIPL_U_CHAR       8
#define GIPL_SHORT        15
#define GIPL_U_SHORT      16
#define GIPL_U_INT        31
#define GIPL_INT          32
#define GIPL_FLOAT        64
#define GIPL_DOUBLE       65
#define GIPL_VECTOR_3D_CHAR   66
#define GIPL_VECTOR_3D_SHORT  67
#define GIPL_VECTOR_3D_FLOAT  68
#define GIPL_VECTOR_3D_DOUBLE 69
#define GIPL_C_SHORT      144
#define GIPL_C_INT        160
#define GIPL_C_FLOAT      192
#define GIPL_C_DOUBLE     193


#endif // MIRTK_GIPL_H

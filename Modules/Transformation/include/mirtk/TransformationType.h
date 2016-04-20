/*
 * Medical Image Registration ToolKit ()
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_TransformationType_H
#define MIRTK_TransformationType_H


namespace mirtk {


/// Enumeration of transformation types
///
/// Each transformation class has its own enumeration value which is written
/// to the transformation file right after the magic number. Different versions
/// are distinguished by different type IDs.
///
/// Enumeration values 50-60 were assigned to the refactored FFD types
/// which store the control point data in an instance of GenericImage.
///
/// Enumeration values 70-80 were assigned to the refactored FFD types
/// which store also the method used to extrapolate control point coefficients
/// outside the finite discrete lattice on which the FFD is defined.
///
/// \attention Do not change the enumeration value of existing entries as then
///            already saved transformation files are not identified correctly.
///            To yet allow a better ordering of the entries, the enumeration
///            values are thus assigned explicitly, also to remind of this.
enum TransformationType
{
  TRANSFORMATION_MAGIC                           = 815007,
  TRANSFORMATION_UNKNOWN                         =      0,
  // linear transformations
  TRANSFORMATION_HOMOGENEOUS                     =      1,
  TRANSFORMATION_RIGID                           =      2,
  TRANSFORMATION_SIMILARITY                      =     22,
  TRANSFORMATION_AFFINE                          =      3,
  TRANSFORMATION_HOMO_TEMPORAL                   =     30,
  TRANSFORMATION_RIGID_TEMPORAL                  =     31,
  TRANSFORMATION_AFFINE_TEMPORAL                 =     32,
  // linear FFD
  TRANSFORMATION_LINEAR_FFD_2D_v1                =     70,
  TRANSFORMATION_LINEAR_FFD_2D = TRANSFORMATION_LINEAR_FFD_2D_v1,
  TRANSFORMATION_LINEAR_FFD_3D_v1                =      5,
  TRANSFORMATION_LINEAR_FFD_3D_v2                =     13,
  TRANSFORMATION_LINEAR_FFD_3D_v3                =     51,
  TRANSFORMATION_LINEAR_FFD_3D_v4                =     71,
  TRANSFORMATION_LINEAR_FFD_3D = TRANSFORMATION_LINEAR_FFD_3D_v4,
  TRANSFORMATION_LINEAR_FFD_4D_v1                =     17,
  TRANSFORMATION_LINEAR_FFD_4D_v2                =     52,
  TRANSFORMATION_LINEAR_FFD_4D_v3                =     72,
  TRANSFORMATION_LINEAR_FFD_4D = TRANSFORMATION_LINEAR_FFD_4D_v3,
  TRANSFORMATION_LINEAR_FFD_SV_v1                =     73,
  TRANSFORMATION_LINEAR_FFD_SV = TRANSFORMATION_LINEAR_FFD_SV_v1,
  TRANSFORMATION_LINEAR_FFD_TD_v1                =     18,
  TRANSFORMATION_LINEAR_FFD_TD_v2                =     54,
  TRANSFORMATION_LINEAR_FFD_TD_v3                =     74,
  TRANSFORMATION_LINEAR_FFD_TD = TRANSFORMATION_LINEAR_FFD_TD_v3,
  // B-spline FFD
  TRANSFORMATION_BSPLINE_FFD_2D_v1               =     75,
  TRANSFORMATION_BSPLINE_FFD_3D_v1               =      4,
  TRANSFORMATION_BSPLINE_FFD_3D_v2               =     12,
  TRANSFORMATION_BSPLINE_FFD_3D_v3               =     56,
  TRANSFORMATION_BSPLINE_FFD_3D_v4               =     76,
  TRANSFORMATION_BSPLINE_FFD_3D = TRANSFORMATION_BSPLINE_FFD_3D_v4,
  TRANSFORMATION_BSPLINE_FFD_4D_v1               =     14,
  TRANSFORMATION_BSPLINE_FFD_4D_v2               =     57,
  TRANSFORMATION_BSPLINE_FFD_4D_v3               =     77,
  TRANSFORMATION_BSPLINE_FFD_4D = TRANSFORMATION_BSPLINE_FFD_4D_v3,
  TRANSFORMATION_BSPLINE_FFD_SV_v1               =     16,
  TRANSFORMATION_BSPLINE_FFD_SV_v2               =     23,
  TRANSFORMATION_BSPLINE_FFD_SV_v3               =     24,
  TRANSFORMATION_BSPLINE_FFD_SV_v4               =     25,
  TRANSFORMATION_BSPLINE_FFD_SV_v5               =     27,
  TRANSFORMATION_BSPLINE_FFD_SV_v6               =     58,
  TRANSFORMATION_BSPLINE_FFD_SV_v7               =     65,
  TRANSFORMATION_BSPLINE_FFD_SV_v8               =     78,
  TRANSFORMATION_BSPLINE_FFD_SV = TRANSFORMATION_BSPLINE_FFD_SV_v8,
  TRANSFORMATION_BSPLINE_FFD_TD_v1               =     15,
  TRANSFORMATION_BSPLINE_FFD_TD_v2               =     21,
  TRANSFORMATION_BSPLINE_FFD_TD_v3               =     59,
  TRANSFORMATION_BSPLINE_FFD_TD_v4               =     79,
  TRANSFORMATION_BSPLINE_FFD_TD = TRANSFORMATION_BSPLINE_FFD_TD_v4,
  // "decorating" transformations
  TRANSFORMATION_BSPLINE_FFD_STATISTICAL         =     61,
  // composite transformations
  TRANSFORMATION_MFFD                            =      7,
  TRANSFORMATION_FLUID_v1                        =      8,
  TRANSFORMATION_FLUID_v2                        =     81,
  TRANSFORMATION_FLUID = TRANSFORMATION_FLUID_v2,
  TRANSFORMATION_MFFD_SV                         =     26,
  // others not included in MIRTK (cf. BioMedIA/IRTK)
  TRANSFORMATION_LATTICE_FFD                     =      9,
  TRANSFORMATION_MULTI_FRAME_LATTICE_FFD         =     10,
  TRANSFORMATION_QUATERNION                      =     11,
  TRANSFORMATION_EIGEN_FFD_3D_v1                 =      6,
  TRANSFORMATION_EIGEN_FFD_3D_v2                 =     60,
  TRANSFORMATION_EIGEN_FFD_3D_v3                 =     80,
  TRANSFORMATION_EIGEN_FFD_3D = TRANSFORMATION_EIGEN_FFD_3D_v3,
  TRANSFORMATION_PERIODIC_v1                     =     20, // obsolete
  TRANSFORMATION_PERIODIC = TRANSFORMATION_PERIODIC_v1
};


} // namespace mirtk

#endif // MIRTK_TransformationType_H

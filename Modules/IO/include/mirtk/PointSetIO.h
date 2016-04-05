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

#ifndef MIRTK_PointSetIO_H
#define MIRTK_PointSetIO_H

#include "mirtk/PointSet.h"

#include "vtkSmartPointer.h"
#include "vtkDataSet.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"


namespace mirtk {


// =============================================================================
// File name extension
// =============================================================================

/// Default extension for given data set
const char *DefaultExtension(vtkDataSet *);

// =============================================================================
// Generic I/O functions
// =============================================================================

/// Read point set from file
///
/// @param[in]  fname           File name.
/// @param[out] ftype           File type (VTK_ASCII or VTK_BINARY) of legacy VTK input file.
/// @param[in]  exit_on_failure Call exit when point set could not be read.
///
/// @return Point set. Dataset is empty if file could not be read and @p exit_on_failure is @c false.
vtkSmartPointer<vtkPointSet> ReadPointSet(const char *fname, int *ftype = nullptr, bool exit_on_failure = true);

/// Read polygonal dataset from file
///
/// @param[in]  fname           File name.
/// @param[out] ftype           File type (VTK_ASCII or VTK_BINARY) of legacy VTK input file.
/// @param[in]  exit_on_failure Call exit when point set could not be read.
///
/// @return Polygonal dataset. Dataset is empty if file could not be read and @p exit_on_failure is @c false.
vtkSmartPointer<vtkPolyData> ReadPolyData(const char *fname, int *ftype = nullptr, bool exit_on_failure = true);

/// Write point set to file
///
/// @param fname    File name. The extension determines the output format.
/// @param pointset Point set to write.
/// @param compress Whether to use compression when writing to VTK XML format (.vtp).
/// @param ascii    Whether to use ASCII format when writing to legacy VTK format (.vtk).
///
/// @return Whether point set was written successfully to the specified file.
bool WritePointSet(const char *fname, vtkPointSet *pointset, bool compress = true, bool ascii = false);

/// Write polygonal dataset to file
///
/// @param fname    File name. The extension determines the output format.
/// @param polydata Polydata to write.
/// @param compress Whether to use compression when writing to VTK XML format (.vtp).
/// @param ascii    Whether to use ASCII format when writing to legacy VTK format (.vtk).
///
/// @return Whether point set was written successfully to the specified file.
bool WritePolyData(const char *fname, vtkPolyData *polydata, bool compress = true, bool ascii = false);

// =============================================================================
// TetGen I/O functions
// =============================================================================

/// Write point set to TetGen .node file
///
/// @param fname    File name.
/// @param pointset Point set.
///
/// @return Whether dataset was written successfully to the specified file.
bool WriteTetGenNode(const char *fname, vtkPointSet *pointset);

/// Write polygonal dataset to TetGen .poly file
///
/// @param fname    File name.
/// @param polydata Polygonal dataset.
/// @param holes    Hole list.
///
/// @return Whether dataset was written successfully to the specified file.
bool WriteTetGenPoly(const char *fname, vtkPolyData *polydata, const PointSet *holes = nullptr);

/// Write polygonal dataset to TetGen .smesh file
///
/// @param fname    File name.
/// @param polydata Polygonal dataset.
/// @param holes    Hole list.
///
/// @return Whether dataset was written successfully to the specified file.
bool WriteTetGenSMesh(const char *fname, vtkPolyData *polydata, const PointSet *holes = nullptr);

// =============================================================================
// BrainSuite I/O functions
// =============================================================================

/// Read BrainSuite surface from .dfs file
///
/// @param[in] fname File name.
///
/// @return Polygonal dataset. Dataset is empty if file could not be read.
vtkSmartPointer<vtkPolyData> ReadDFS(const char *fname);

/// Write surface in BrainSuite .dfs format
///
/// @param[in] fname    File name.
/// @param[in] polydata Polygonal dataset.
///
/// @return Whether surface was written successfully to the specified file.
bool WriteDFS(const char *fname, vtkPolyData *polydata);

// =============================================================================
// Object File Format I/O functions
// =============================================================================

/// Read polygonal dataset from Object File Format (.off) file
///
/// @param[in] fname File name.
///
/// @return Polygonal dataset. Dataset is empty if file could not be read.
vtkSmartPointer<vtkPolyData> ReadOFF(const char *fname);

/// Write polygonal dataset to Object File Format (.off) file
///
/// @param[in] fname    File name.
/// @param[in] polydata Polygonal dataset.
///
/// @return Whether dataset was written successfully to the specified file.
bool WriteOFF(const char *fname, vtkPolyData *polydata);


} // namespace mirtk

#endif // MIRTK_VtkPointSetIO_H

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

#include "mirtk/IOConfig.h"
#include "mirtk/PointSet.h"

#include "vtkSmartPointer.h"
#include "vtkDataSet.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"

#if MIRTK_IO_WITH_GIFTI
  #include "vtkPoints.h"
  #include "vtkCellArray.h"
  #include "vtkDataArray.h"
  #include "vtkPointData.h"
  #include "vtkInformationStringKey.h"
  #include "vtkInformationIntegerKey.h"
  #include "vtkInformationDoubleKey.h"
#endif


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

// =============================================================================
// GIFTI I/O functions -- https://www.nitrc.org/projects/gifti/
// =============================================================================
#if MIRTK_IO_WITH_GIFTI

/// GIFTI meta data keys
///
/// @sa Section 3.0 Standard MetaData of GIFTI file format specification at
///     https://www.nitrc.org/frs/download.php/2871/GIFTI_Surface_Format.pdf
class GiftiMetaData
{
  // ---------------------------------------------------------------------------
  // vtkInformation key instances for standard GIFTI meta data entries
public:

  /// Date and possibly time when the GIFTI file was written
  static vtkInformationStringKey *DATE();

  /// Name of user that wrote the GIFTI file
  static vtkInformationStringKey *USER_NAME();

  /// ID of subject whose anatomical structure the point set models
  static vtkInformationStringKey *SUBJECT_ID();

  /// A unique string that identifies a surface
  static vtkInformationStringKey *SURFACE_ID();

  /// A unique string that identifies a data array (UUID)
  static vtkInformationStringKey *UNIQUE_ID();

  /// Name of GIFTI data array
  static vtkInformationStringKey *NAME();

  /// Description of GIFTI file or data array
  static vtkInformationStringKey *DESCRIPTION();

  /// Included in a GIFTI time series file, specifies repetition time (TR)
  /// and is equivalent to the slice_duration of a NIfTI volume
  static vtkInformationDoubleKey *TIME_STEP();

  /// Data space of point set before any coordinate system transformation
  static vtkInformationStringKey *DATA_SPACE();

  /// Anatomical structure that the point set models
  static vtkInformationStringKey *ANATOMICAL_STRUCTURE_PRIMARY();

  /// Further describe anatomical structure that the point set models
  static vtkInformationStringKey *ANATOMICAL_STRUCTURE_SECONDARY();

  /// Describes geometry of the point set
  static vtkInformationStringKey *GEOMETRIC_TYPE();

  /// Topology of surface model
  static vtkInformationStringKey *TOPOLOGICAL_TYPE();

  /// Intent code of functional data array
  static vtkInformationIntegerKey *INTENT_CODE();

  /// First parameter of statistical test
  static vtkInformationDoubleKey *INTENT_P1();

  /// Second parameter of statistical test
  static vtkInformationDoubleKey *INTENT_P2();

  /// Third parameter of statistical test
  static vtkInformationDoubleKey *INTENT_P3();

};

/// Read point set coordinates from GIFTI ([.coord].gii) file
///
/// @param[in]     fname  File name.
/// @param[in,out] info   vtkInformation to which to add geometric meta data.
///                       If nullptr, the meta data is discarded.
/// @param[in]     errmsg Whether to print error messages if any.
///
/// @return Point set coordinates or empty set if file could not be read
///         or has no valid GIFTI data array with intent @c NIFTI_INTENT_POINTSET.
vtkSmartPointer<vtkPoints> ReadGIFTICoordinates(const char     *fname,
                                                vtkInformation *info   = nullptr,
                                                bool            errmsg = false);

/// Read surface topology from GIFTI ([.surf|.topo].gii) file
///
/// @param[in]     fname  File name.
/// @param[in,out] info   vtkInformation to which to add topological meta data.
///                       If nullptr, the meta data is discarded.
/// @param[in]     errmsg Whether to print error messages if any.
///
/// @return Triangle list or nullptr if file could not be read or has no valid
///         GIFTI data array with intent @c NIFTI_INTENT_TRIANGLE.
vtkSmartPointer<vtkCellArray> ReadGIFTITopology(const char     *fname,
                                                vtkInformation *info   = nullptr,
                                                bool            errmsg = false);

/// Read point data arrays from GIFTI (.gii) file
///
/// @param[in] fname  File name.
/// @param[in] errmsg Whether to print error messages if any.
///
/// @return Point data arrays or nullptr if file could not be read.
vtkSmartPointer<vtkPointData> ReadGIFTIPointData(const char *fname,
                                                 bool errmsg = false);

/// Read polygonal dataset from GIFTI (.gii) file
///
/// Standard GIFTI meta data is stored in the vtkInformation of the returned
/// vtkPolyData instance and the respective information of the corresponding
/// vtkDataArray instances of its point data.
///
/// @param[in]     fname   File name.
/// @param[in,out] surface Polygonal dataset to which to add GIFTI data arrays.
///                        If nullptr, a new polygonal dataset is allocated
///                        and a GIFTI data array with @c NIFTI_INTENT_POINTSET
///                        must be present in the GIFTI file. If input surface
///                        has a previously read point set, the GIFTI file from
///                        which to read additional point data can use sparse
///                        storage using a data array with @c NIFTI_INTENT_NODE_INDEX
///                        that specifies the points to which the read data array
///                        values are added. Missing data values default to zero.
/// @param[in]     errmsg  Whether to print error messages if any.
///
/// @return Polygonal dataset. Dataset is empty if file could not be read.
///
/// Example:
/// @code
/// // Read geometry and topology of cortical surface model
/// vtkSmartPointer<vtkPolyData> surface = ReadGIFTI("cortex.surf.gii");
/// // Add shape measurements and functional statistics as surface point data
/// surface = ReadGIFTI("cortex.shape.gii", surface);
/// surface = ReadGIFTI("cortex.func.gii",  surface);
/// @endcode
///
/// @sa https://www.nitrc.org/projects/gifti/
vtkSmartPointer<vtkPolyData> ReadGIFTI(const char  *fname,
                                       vtkPolyData *surface = nullptr,
                                       bool         errmsg  = false);

/// Write polygonal dataset to GIFTI ([.surf].gii) file
///
/// @param[in] fname    File name. Based on the extension either all or only
///                     certain data arrays of the polygonal dataset are saved.
///                     - .surf.gii:   Point coordinates and topology.
///                     - .coords.gii: Point coordinates.
///                     - .topo.gii:   Topology.
///                     - otherwise:   Point coordinates, topology, and point data.
/// @param[in] polydata Polygonal dataset.
/// @param[in] compress Whether to compress binary data.
/// @param[in] ascii    Whether to write data in ASCII format.
///
/// @return Whether dataset was written successfully to the specified file.
///
/// @sa https://www.nitrc.org/projects/gifti/
bool WriteGIFTI(const char *fname, vtkPolyData *polydata,
                bool compress = true, bool ascii = false);

#endif // MIRTK_IO_WITH_GIFTI

} // namespace mirtk

#endif // MIRTK_PointSetIO_H

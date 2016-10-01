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

#ifndef MIRTK_Vtk_H
#define MIRTK_Vtk_H

#include "mirtk/String.h"

#include "vtkConfigure.h"
#include "vtkNew.h"
#include "vtkSmartPointer.h"
#include "vtkDataArray.h"
#include "vtkDataSetAttributes.h"


namespace mirtk {


// =============================================================================
// Defines used by MIRTK image I/O filters
// =============================================================================

// VTK Legacy file format magic header
#define VTK_MAGIC1 "# vtk DataFile Version 1.0"
#define VTK_MAGIC2 "# vtk DataFile Version 2.0"
#define VTK_MAGIC3 "# vtk DataFile Version 3.0"

// VTK data types
#define VTK_DATA_CHAR    "char"
#define VTK_DATA_U_CHAR  "unsigned_char"
#define VTK_DATA_SHORT   "short"
#define VTK_DATA_U_SHORT "unsigned_short"
#define VTK_DATA_FLOAT   "float"

// =============================================================================
// VTK 5/6 transition
// =============================================================================

// Auxiliary macros to set/add VTK filter input (connection)
#if VTK_MAJOR_VERSION >= 6
#  define SetVTKInput(filter, dataset) (filter)->SetInputData(dataset);
#  define AddVTKInput(filter, dataset) (filter)->AddInputData(dataset);
#  define SetVTKConnection(filter2, filter1) (filter2)->SetInputConnection((filter1)->GetOutputPort());
#  define SetNthVTKInput(filter, n, dataset) (filter)->SetInputData(n, dataset);
#  define AddNthVTKInput(filter, n, dataset) (filter)->AddInputData(n, dataset);
#  define SetNthVTKConnection(filter2, n2, filter1, n1) (filter2)->SetInputConnection(n2, (filter1)->GetOutputPort(n1));
#else
#  define SetVTKInput(filter, dataset) (filter)->SetInput(dataset);
#  define AddVTKInput(filter, dataset) (filter)->AddInput(dataset);
#  define SetVTKConnection(filter2, filter1) (filter2)->SetInput((filter1)->GetOutput());
#  define SetNthVTKInput(filter, n, dataset) (filter)->SetInput(n, dataset);
#  define AddNthVTKInput(filter, n, dataset) (filter)->AddInput(n, dataset);
#  define SetNthVTKConnection(filter2, n2, filter1, n1) (filter2)->SetInput(n2, (filter1)->GetOutput(n1));
#endif

// =============================================================================
// Data set attributes
// =============================================================================

// -----------------------------------------------------------------------------
/// Instantiate new VTK data array of given type
///
/// \param[in] type   VTK data type ID, e.g., VTK_FLOAT. When VTK_VOID, a floating
///                   point data array with default precision, i.e., either single
///                   or double is returned.
/// \param[in] tuples Number of tuples. The array is uninitialized when non-positive.
/// \param[in] comps  Number of components per tuple.
/// \param[in] name   Data array name.
///
/// \returns New VTK data array instance.
vtkSmartPointer<vtkDataArray>
NewVtkDataArray(int type = VTK_VOID, int tuples = 0, int comps = 1, const char *name = nullptr);

/// \deprecated Use NewVtkDataArray instead
vtkSmartPointer<vtkDataArray> NewVTKDataArray(int type = VTK_VOID);

// -----------------------------------------------------------------------------
/// Convert string to vtkDataSetAttributes::AttributeType
template <>
inline bool FromString(const char *str, vtkDataSetAttributes::AttributeTypes &type)
{
  string name = ToLower(str);
  if      (name == "scalars")     type = vtkDataSetAttributes::SCALARS;
  else if (name == "vectors")     type = vtkDataSetAttributes::VECTORS;
  else if (name == "normals")     type = vtkDataSetAttributes::NORMALS;
  else if (name == "tcoords")     type = vtkDataSetAttributes::TCOORDS;
  else if (name == "tensors")     type = vtkDataSetAttributes::TENSORS;
  else if (name == "globalids")   type = vtkDataSetAttributes::GLOBALIDS;
  else if (name == "pedigreeids") type = vtkDataSetAttributes::PEDIGREEIDS;
  else if (name == "edgeflag")    type = vtkDataSetAttributes::EDGEFLAG;
  else if (name == "other")       type = vtkDataSetAttributes::NUM_ATTRIBUTES;
  else return false;
  return true;
}

// -----------------------------------------------------------------------------
/// Convert vtkDataSetAttributes::AttributeType to string
template <>
inline string ToString(const vtkDataSetAttributes::AttributeTypes &type,
                       int w, char c, bool left)
{
  const char *str;
  switch (type) {
    case vtkDataSetAttributes::SCALARS:     str = "scalars";     break;
    case vtkDataSetAttributes::VECTORS:     str = "vectors";     break;
    case vtkDataSetAttributes::NORMALS:     str = "normals";     break;
    case vtkDataSetAttributes::TCOORDS:     str = "tcoords";     break;
    case vtkDataSetAttributes::TENSORS:     str = "tensors";     break;
    case vtkDataSetAttributes::GLOBALIDS:   str = "globalids";   break;
    case vtkDataSetAttributes::PEDIGREEIDS: str = "pedigreeids"; break;
    case vtkDataSetAttributes::EDGEFLAG:    str = "edgeflag";    break;
    default:                                str = "other";       break;
  }
  return ToString(str, w, c, left);
}

// -----------------------------------------------------------------------------
/// Convert vtkDataSetAttributes::AttributeType to string
inline string VtkAttributeTypeString(int type)
{
  if (type < 0 || type >= vtkDataSetAttributes::NUM_ATTRIBUTES) return "other";
  return ToString(static_cast<vtkDataSetAttributes::AttributeTypes>(type));
}

// -----------------------------------------------------------------------------
/// Convert vtkDataSetAttributes::AttributeType to string
inline string VtkDataTypeString(int type)
{
  switch (type)
  {
    case VTK_VOID:               return "void";
    case VTK_CHAR:               return "char";
    case VTK_SHORT:              return "short";
    case VTK_INT:                return "int";
    case VTK_LONG:               return "long";
    case VTK_LONG_LONG:          return "int64";
    case VTK_UNSIGNED_CHAR:      return "uchar";
    case VTK_UNSIGNED_SHORT:     return "ushort";
    case VTK_UNSIGNED_INT:       return "uint";
    case VTK_UNSIGNED_LONG:      return "ulong";
    case VTK_UNSIGNED_LONG_LONG: return "uint64";
    case VTK_FLOAT:              return "float";
    case VTK_DOUBLE:             return "double";
    default:                     return "unknown";
  }
}


} // namespace mirtk

#endif // MIRTK_Vtk_H

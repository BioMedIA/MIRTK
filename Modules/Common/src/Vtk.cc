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

#include "mirtk/Config.h" // MIRTK_USE_FLOAT_BY_DEFAULT

#include "mirtk/Vtk.h"
#include "mirtk/Stream.h"
#include "mirtk/Exception.h"

#include "vtkCharArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkShortArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkIntArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"


namespace mirtk {


// -----------------------------------------------------------------------------
vtkSmartPointer<vtkDataArray> NewVtkDataArray(int type, int tuples, int comps, const char *name)
{
  vtkSmartPointer<vtkDataArray> arr;
  switch (type) {
    case VTK_VOID: {
      #if MIRTK_USE_FLOAT_BY_DEFAULT
        arr = vtkSmartPointer<vtkFloatArray>::New();
      #else
        arr = vtkSmartPointer<vtkDoubleArray>::New();
      #endif
    } break;
    case VTK_CHAR:           { arr = vtkSmartPointer<vtkCharArray>::New(); } break;
    case VTK_UNSIGNED_CHAR:  { arr = vtkSmartPointer<vtkUnsignedCharArray>::New(); } break;
    case VTK_SHORT:          { arr = vtkSmartPointer<vtkShortArray>::New(); } break;
    case VTK_UNSIGNED_SHORT: { arr = vtkSmartPointer<vtkUnsignedShortArray>::New(); } break;
    case VTK_INT:            { arr = vtkSmartPointer<vtkIntArray>::New(); } break;
    case VTK_UNSIGNED_INT:   { arr = vtkSmartPointer<vtkUnsignedIntArray>::New(); } break;
    case VTK_FLOAT:          { arr = vtkSmartPointer<vtkFloatArray>::New(); } break;
    case VTK_DOUBLE:         { arr = vtkSmartPointer<vtkDoubleArray>::New(); } break;
    case VTK_ID_TYPE:        { arr = vtkSmartPointer<vtkIdTypeArray>::New(); } break;
    default:
      Throw(ERR_LogicError, __FUNCTION__, "Invalid VTK data type: ", type);
  }
  if (comps  > 0) arr->SetNumberOfComponents(comps);
  if (tuples > 0) arr->SetNumberOfTuples(tuples);
  if (name != nullptr) arr->SetName(name);
  return arr;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkDataArray> NewVTKDataArray(int vtkType)
{
  return NewVtkDataArray(vtkType);
}


} // namespace mirtk

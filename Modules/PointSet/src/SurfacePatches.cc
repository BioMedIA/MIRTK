/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#include "mirtk/SurfacePatches.h"

#include "mirtk/Queue.h"
#include "mirtk/Algorithm.h"

#include "vtkNew.h"
#include "vtkSmartPointer.h"
#include "vtkDataArray.h"
#include "vtkIdList.h"
#include "vtkCellData.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void SurfacePatches::CopyAttributes(const SurfacePatches &other)
{
  _Ordering        = other._Ordering;
  _NumberOfPatches = other._NumberOfPatches;
  _PatchSize       = other._PatchSize;
}

// -----------------------------------------------------------------------------
SurfacePatches::SurfacePatches()
:
  _Ordering(LargestFirst),
  _NumberOfPatches(0)
{
}

// -----------------------------------------------------------------------------
SurfacePatches::SurfacePatches(const SurfacePatches &other)
:
  SurfaceFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
SurfacePatches &SurfacePatches::operator =(const SurfacePatches &other)
{
  if (this != &other) {
    SurfaceFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
SurfacePatches::~SurfacePatches()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
vtkDataArray *SurfacePatches::GetLabelsArray() const
{
  return _Output->GetCellData()->GetArray("PatchLabel");
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void SurfacePatches::Initialize()
{
  // Initialize base class
  SurfaceFilter::Initialize();

  // Initialize patch labels
  vtkSmartPointer<vtkDataArray> labels = GetLabelsArray();
  if (labels == nullptr || labels->GetNumberOfComponents() != 1) {
    labels = NewCellArray("PatchLabel", 1, VTK_UNSIGNED_SHORT);
    _Output->GetCellData()->AddArray(labels);
  }
  labels->FillComponent(0, 0.);

  // Reset counters
  _NumberOfPatches = 0;
  _PatchSize.clear();
}

// -----------------------------------------------------------------------------
void SurfacePatches::Execute()
{
  vtkDataArray * const labels = GetLabelsArray();

  Queue<vtkIdType>  active;
  vtkIdType         seedId, curId, nbrId, npts, *pts, ncells;
  vtkNew<vtkIdList> cellIds;

  for (seedId = 0; seedId < _Output->GetNumberOfCells(); ++seedId) {
    if (labels->GetComponent(seedId, 0) == 0.) {
      ncells = 0;
      ++_NumberOfPatches;
      active.push(seedId);
      while (!active.empty()) {
        curId = active.front();
        active.pop();
        if (labels->GetComponent(curId, 0) == 0.) {
          ++ncells;
          labels->SetComponent(curId, 0, _NumberOfPatches);
          _Output->GetCellPoints(curId, npts, pts);
          for (vtkIdType i = 0; i < npts; ++i) {
            _Output->GetCellEdgeNeighbors(curId, pts[i], pts[(i+1)%npts], cellIds.GetPointer());
            if (cellIds->GetNumberOfIds() == 1) {
              nbrId = cellIds->GetId(0);
              if (labels->GetComponent(nbrId, 0) == 0.) {
                active.push(nbrId);
              }
            }
          }
        }
      }
      _PatchSize.push_back(ncells);
    }
  }
}

// -----------------------------------------------------------------------------
void SurfacePatches::Finalize()
{
  // Relabel components based on patch size order
  if (_Ordering != NoOrdering) {
    Array<int> order;
    if (_Ordering == SmallestFirst) {
      order = IncreasingOrder(_PatchSize);
    } else {
      order = DecreasingOrder(_PatchSize);
    }
    Array<int> new_labels(_NumberOfPatches);
    for (int i = 0; i < _NumberOfPatches; ++i) {
      new_labels[order[i]] = i + 1;
    }
    const Array<int> new_size = _PatchSize; // make copy
    int current, new_label;
    vtkDataArray * const labels = GetLabelsArray();
    for (vtkIdType cellId = 0; cellId < _Output->GetNumberOfCells(); ++cellId) {
      current = static_cast<int>(labels->GetComponent(cellId, 0));
      if (current != 0.) {
        new_label = new_labels[current-1];
        _PatchSize[new_label-1] = new_size[current-1];
        labels->SetComponent(cellId, 0, static_cast<double>(new_label));
      }
    }
  }

  // Finalize base class
  SurfaceFilter::Finalize();
}


} // namespace mirtk

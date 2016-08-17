/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
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

#include "mirtk/ConnectedComponents.h"

#include "mirtk/Assert.h"
#include "mirtk/Queue.h"
#include "mirtk/OrderedSet.h"
#include "mirtk/Algorithm.h"


namespace mirtk {


// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Marks all voxels as unvisited and collects set of unique segmentation labels
template <class TLabel>
void MarkUnvisited(const GenericImage<TLabel> &segmentation,
                   GenericImage<TLabel>       &components)
{
  for (int idx = 0; idx < segmentation.NumberOfVoxels(); ++idx) {
    if (segmentation(idx)) {
      components(idx) = voxel_limits<TLabel>::max_value();
    } else {
      components(idx) = TLabel(0);
    }
  }
}

// -----------------------------------------------------------------------------
/// Finds index of next unlabeled component voxel or -1 when done
template <class TLabel>
int NextSeed(const GenericImage<TLabel> &components)
{
  for (int idx = 0; idx < components.NumberOfVoxels(); ++idx) {
    if (components(idx) == voxel_limits<TLabel>::max_value()) return idx;
  }
  return -1;
}

// -----------------------------------------------------------------------------
/// Performs region growing for a specific input label and seed
template <class TLabel>
int LabelComponent(const NeighborhoodOffsets  &offsets,
                   const GenericImage<TLabel> &segmentation,
                   GenericImage<TLabel>       &components,
                   int                         seed,
                   TLabel                      component_label)
{
  const TLabel unvisited     = voxel_limits<TLabel>::max_value();
  const TLabel current_label = segmentation(seed);
  const TLabel * start_label = segmentation.Data();

  int idx, neighbor_idx, num = 0;
  Queue<int> active;
  active.push(seed);

  while (!active.empty()) {
    idx = active.front();
    active.pop();
    if (components(idx) == unvisited) {
      components(idx) = component_label, ++num;
      if (segmentation.IsBoundary(idx)) {
        const TLabel *neighbor_label, *label = segmentation.Data(idx);
        const TLabel *neighbor_comp,  *comp  = components  .Data(idx);
        for (int n = 0; n < offsets.Size(); ++n) {
          neighbor_label = label + offsets(n);
          neighbor_idx   = static_cast<int>(neighbor_label - start_label);
          if (0 <= neighbor_idx && neighbor_idx < segmentation.NumberOfVoxels()) {
            neighbor_comp = comp + offsets(n);
            if (*neighbor_label == current_label && *neighbor_comp == unvisited) {
              active.push(neighbor_idx);
            }
          }
        }
      } else {
        const TLabel *neighbor_label, *label = segmentation.Data(idx);
        const TLabel *neighbor_comp,  *comp  = components  .Data(idx);
        for (int n = 0; n < offsets.Size(); ++n) {
          neighbor_label = label + offsets(n);
          neighbor_comp  = comp  + offsets(n);
          if (*neighbor_label == current_label && *neighbor_comp == unvisited) {
            active.push(static_cast<int>(neighbor_label - start_label));
          }
        }
      }
    }
  }

  return num;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
ConnectedComponents<VoxelType>::ConnectedComponents(ConnectedComponentsOrdering ordering, ConnectivityType conn)
:
  _Ordering(ordering),
  _Connectivity(conn),
  _NumberOfComponents(0)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
ConnectedComponents<VoxelType>::~ConnectedComponents()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void ConnectedComponents<VoxelType>::Initialize()
{
  ImageToImage<VoxelType>::Initialize();

  _NumberOfComponents = 0;
  _ComponentSize.clear();

  MarkUnvisited(*this->Input(), *this->Output());
  _Offsets.Initialize(this->Input(), _Connectivity);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void ConnectedComponents<VoxelType>::Run()
{
  this->Initialize();

  int num, seed;
  VoxelType c(0);
  while ((seed = NextSeed(*this->Output())) != -1) {
    if (++c == voxel_limits<VoxelType>::max_value()) {
      cerr << "ConnectedComponents::Run: No. of components exceeded maximum label value!" << endl;
      exit(1);
    }
    num = LabelComponent(_Offsets, *this->Input(), *this->Output(), seed, c);
    mirtkAssert(num > 0, "found seed must be valid and no. of unvisited voxels decrease");
    _ComponentSize.push_back(num);
  }
  _NumberOfComponents = static_cast<int>(c);

  this->Finalize();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void ConnectedComponents<VoxelType>::Finalize()
{
  if (_Ordering != CC_NoOrdering) {
    Array<int> order;
    if (_Ordering == CC_SmallestFirst) {
      order = IncreasingOrder(_ComponentSize);
    } else {
      order = DecreasingOrder(_ComponentSize);
    }
    Array<VoxelType> new_label(_NumberOfComponents);
    for (int i = 0; i < _NumberOfComponents; ++i) {
      new_label[order[i]] = voxel_cast<VoxelType>(i + 1);
    }
    GenericImage<VoxelType> &output = *this->Output();
    const Array<int> new_size = _ComponentSize; // make copy
    VoxelType current;
    const VoxelType zero(0);
    for (int idx = 0; idx < output.NumberOfVoxels(); ++idx) {
      current = output(idx);
      if (current > zero) {
        output(idx)                   = new_label[current-1];
        _ComponentSize[output(idx)-1] = new_size [current-1];
      }
    }
  }

  ImageToImage<VoxelType>::Finalize();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void ConnectedComponents<VoxelType>::DeleteComponent(VoxelType c)
{
  ++c; // 1-indexed component labels
  const VoxelType zero(0);
  GenericImage<VoxelType> &output = *this->Output();
  for (int idx = 0; idx < output.NumberOfVoxels(); ++idx) {
    if (output(idx) == c) {
      output(idx) = zero;
    }
  }
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class ConnectedComponents<BytePixel>;
template class ConnectedComponents<GreyPixel>;


} // namespace mirtk

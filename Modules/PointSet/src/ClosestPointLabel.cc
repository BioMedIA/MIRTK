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

#include "mirtk/ClosestPointLabel.h"

#include "mirtk/Math.h"
#include "mirtk/Array.h"
#include "mirtk/Queue.h"
#include "mirtk/UnorderedMap.h"
#include "mirtk/UnorderedSet.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/VtkMath.h"

#include "vtkSmartPointer.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkIdTypeArray.h"
#include "vtkFloatArray.h"
#include "vtkIdList.h"
#include "vtkOctreePointLocator.h"


namespace mirtk {


// =============================================================================
// Auxiliaries
// =============================================================================

namespace ClosestPointLabelUtils {


// -----------------------------------------------------------------------------
/// Determine number of distinct labels
int NumberOfLabels(vtkDataArray *labels, Array<int> &labelSet, UnorderedMap<int, int> &labelPos)
{
  UnorderedSet<int> uniqueLabelSet;
  for (vtkIdType i = 0; i < labels->GetNumberOfTuples(); ++i) {
    uniqueLabelSet.insert(labels->GetComponent(i, 0));
  }
  int idx = 0;
  labelSet.reserve(uniqueLabelSet.size());
  labelPos.clear();
  for (UnorderedSet<int>::const_iterator i = uniqueLabelSet.begin(); i != uniqueLabelSet.end(); ++i, ++idx) {
    labelSet.push_back(*i);
    labelPos[*i] = idx;
  }
  return static_cast<int>(labelSet.size());
}

// -----------------------------------------------------------------------------
struct ComputeClosestLabeledPoints
{
  vtkPointSet      *_DataSet;
  EdgeTable        *_EdgeTable;
  const Array<int> *_LabelSet;
  vtkDataArray     *_Labels;
  vtkDataArray     *_MinDistance;
  vtkDataArray     *_ClosestPoint;

  bool IsBoundaryPoint(int ptId, int label) const
  {
    const int *adjPtId, *adjPtEnd;
    for (_EdgeTable->GetAdjacentPoints(ptId, adjPtId, adjPtEnd); adjPtId != adjPtEnd; ++adjPtId) {
      if (static_cast<int>(_Labels->GetComponent(*adjPtId, 0)) != label) return true;
    }
    return false;
  }

  void Initialize(int c, Queue<int> &active) const
  {
    const int label = _LabelSet->at(c);
    while (!active.empty()) active.pop();
    for (vtkIdType ptId = 0; ptId < _DataSet->GetNumberOfPoints(); ++ptId) {
      if (static_cast<int>(_Labels->GetComponent(ptId, 0)) == label) {
        _ClosestPoint->SetComponent(ptId, c, ptId);
        _MinDistance ->SetComponent(ptId, c, .0);
        if (IsBoundaryPoint(ptId, label)) active.push(static_cast<int>(ptId));
      } else {
        _ClosestPoint->SetComponent(ptId, c, -1.0);
        _MinDistance ->SetComponent(ptId, c, numeric_limits<double>::infinity());
      }
    }
  }

  void operator ()(const blocked_range<int> &re) const
  {
    double p1[3], p2[3], curDistance, minDistance, newDistance;
    int ptId, closestPtId;
    Queue<int> active;
    UnorderedMap<int, int>::const_iterator it;
    const int *adjPtId, *adjPtEnd;

    for (int c = re.begin(); c != re.end(); ++c) {
      Initialize(c, active);
      while (!active.empty()) {
        ptId = active.front(), active.pop();
        _DataSet->GetPoint(ptId, p1);
        minDistance = _MinDistance ->GetComponent(ptId, c);
        closestPtId = static_cast<int>(_ClosestPoint->GetComponent(ptId, c));
        for (_EdgeTable->GetAdjacentPoints(ptId, adjPtId, adjPtEnd); adjPtId != adjPtEnd; ++adjPtId) {
          curDistance = _MinDistance->GetComponent(*adjPtId, c);
          if (curDistance <= minDistance) continue;
          _DataSet->GetPoint(*adjPtId, p2);
          newDistance = minDistance + vtkMath::Distance2BetweenPoints(p1, p2);
          if (curDistance <= newDistance) continue;
          _ClosestPoint->SetComponent(*adjPtId, c, closestPtId);
          _MinDistance ->SetComponent(*adjPtId, c, newDistance);
          active.push(*adjPtId);
        }
      }
    }
  }
};

// -----------------------------------------------------------------------------
struct FindClosestLabeledPoints
{
  vtkPoints                    *_TargetPoints;
  vtkDataArray                 *_TargetLabels;
  vtkAbstractPointLocator      *_Source;
  const UnorderedMap<int, int> *_LabelIndex;
  vtkDataArray                 *_ClosestPoint;
  Array<int>                   *_SourceIndex;
  Array<double>                *_Distance;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    int label;
    double p1[3], p2[3];
    vtkIdType closestPtId;
    UnorderedMap<int, int>::const_iterator it;
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      label = static_cast<int>(_TargetLabels->GetComponent(ptId, 0));
      it = _LabelIndex->find(label);
      if (it == _LabelIndex->end()) {
        (*_SourceIndex)[ptId] = -1;
      } else {
        _TargetPoints->GetPoint(ptId, p1);
        closestPtId = _Source->FindClosestPoint(p1);
        closestPtId = static_cast<vtkIdType>(_ClosestPoint->GetComponent(closestPtId, it->second));
        (*_SourceIndex)[ptId] = static_cast<int>(closestPtId);
        _Source->GetDataSet()->GetPoint(closestPtId, p2);
        (*_Distance)[ptId] = vtkMath::Distance2BetweenPoints(p1, p2);
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Label points based on simple majority voting of neighboring cell labels
struct ConvertCellLabelsToPointLabels
{
  vtkPointSet  *_DataSet;
  vtkDataArray *_CellLabels;
  vtkDataArray *_PointLabels;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    UnorderedMap<int, int> count;
    UnorderedMap<int, int>::iterator it;
    int label, maxcount;

    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
    for (vtkIdType ptId = 0; ptId < _DataSet->GetNumberOfPoints(); ++ptId) {
      _DataSet->GetPointCells(ptId, cellIds);
      count.clear();
      for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
        label = static_cast<int>(_CellLabels->GetComponent(cellIds->GetId(i), 0));
        it = count.find(label);
        if (it != count.end()) ++it->second;
        else count[label] = 1;
      }
      maxcount = 0;
      for (it = count.begin(); it != count.end(); ++it) {
        if (it->second > maxcount) {
          _PointLabels->SetComponent(ptId, 0, it->first);
          maxcount = it->second;
        }
      }
    }
  }
};

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkDataArray> CellLabelsToPointLabels(vtkPointSet *dataset, vtkDataArray *cellLabels)
{
  vtkSmartPointer<vtkDataArray> pointLabels;
  pointLabels.TakeReference(cellLabels->NewInstance());
  pointLabels->SetName("Labels");
  pointLabels->SetNumberOfComponents(1);
  pointLabels->SetNumberOfTuples(dataset->GetNumberOfPoints());
  ConvertCellLabelsToPointLabels convert;
  convert._DataSet     = dataset;
  convert._CellLabels  = cellLabels;
  convert._PointLabels = pointLabels;
  parallel_for(blocked_range<vtkIdType>(0, dataset->GetNumberOfPoints()), convert);
  return pointLabels;
}

// -----------------------------------------------------------------------------
void UpdateClosestLabeledPoints(const RegisteredPointSet        *target,
                                const RegisteredPointSet        *source,
                                vtkSmartPointer<vtkIdTypeArray> &closest_points,
                                UnorderedMap<int, int>          &label_pos)
{
  if (closest_points) closest_points->Initialize();
  else closest_points = vtkSmartPointer<vtkIdTypeArray>::New();
  // Get source point labels
  vtkSmartPointer<vtkDataArray> targetLabels, sourceLabels;
  targetLabels = GetArrayByCaseInsensitiveName(target->PointSet()->GetPointData(), "labels");
  if (!targetLabels) {
    vtkDataArray *cellLabels = GetArrayByCaseInsensitiveName(target->PointSet()->GetCellData(), "labels");
    if (cellLabels) targetLabels = CellLabelsToPointLabels(target->PointSet(), cellLabels);
  }
  sourceLabels = GetArrayByCaseInsensitiveName(source->PointSet()->GetPointData(), "labels");
  if (!sourceLabels) {
    vtkDataArray *cellLabels = GetArrayByCaseInsensitiveName(source->PointSet()->GetCellData(), "labels");
    if (cellLabels) sourceLabels = CellLabelsToPointLabels(source->PointSet(), cellLabels);
  }
  if (targetLabels == NULL || sourceLabels == NULL) {
    cerr << "ClosestPointLabel::Update: Missing data array named \"Labels\" (case insensitive)" << endl;
    exit(1);
  }
  Array<int> labelSet;
  const int nlabels = NumberOfLabels(sourceLabels, labelSet, label_pos);
  if (nlabels == 0) return;
  // Build edge table from source mesh
  MIRTK_START_TIMING();
  EdgeTable sourceEdges(source->PointSet());
  MIRTK_DEBUG_TIMING(5, "building the edge table");
  // Determine closest points to each label
  MIRTK_RESET_TIMING();
  vtkSmartPointer<vtkFloatArray> minDistance = vtkSmartPointer<vtkFloatArray>::New();
  minDistance->SetName("MinDistance");
  minDistance->SetNumberOfComponents(nlabels);
  minDistance->SetNumberOfTuples(source->NumberOfPoints());
  closest_points->SetName("ClosestPointId");
  closest_points->SetNumberOfComponents(nlabels);
  closest_points->SetNumberOfTuples(source->NumberOfPoints());
  ComputeClosestLabeledPoints init;
  init._DataSet      = source->PointSet();
  init._Labels       = sourceLabels;
  init._EdgeTable    = &sourceEdges;
  init._LabelSet     = &labelSet;
  init._MinDistance  = minDistance;
  init._ClosestPoint = closest_points;
  parallel_for(blocked_range<int>(0, nlabels), init);
  MIRTK_DEBUG_TIMING(5, "distance transform for each label");
}

// -----------------------------------------------------------------------------
/// Find corresponding source points with label identical to the one of the target point
void UpdateCorrespondences(const RegisteredPointSet     *target,
                           const RegisteredPointSet     *source,
                           vtkIdTypeArray               *closest_points,
                           const UnorderedMap<int, int> &label_pos,
                           Array<int>                   &cor,
                           Array<double>                &distance)
{
  MIRTK_START_TIMING();
  cor     .resize(target->NumberOfPoints(), -1);
  distance.resize(target->NumberOfPoints(), .0);
  vtkSmartPointer<vtkOctreePointLocator> octree = vtkSmartPointer<vtkOctreePointLocator>::New();
  octree->SetDataSet(source->PointSet());
  octree->BuildLocator();
  FindClosestLabeledPoints find;
  find._TargetPoints = target->Points();
  find._TargetLabels = GetArrayByCaseInsensitiveName(target->PointSet()->GetPointData(), "labels");
  find._ClosestPoint = closest_points;
  find._LabelIndex   = &label_pos;
  find._Source       = octree;
  find._SourceIndex  = &cor;
  find._Distance     = &distance;
  parallel_for(blocked_range<vtkIdType>(0, target->NumberOfPoints()), find);
  MIRTK_DEBUG_TIMING(5, "update of point correspondences");
}


} // namespace ClosestPointLabelUtils
using namespace ClosestPointLabelUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
ClosestPointLabel::ClosestPointLabel()
{
}

// -----------------------------------------------------------------------------
ClosestPointLabel::ClosestPointLabel(const ClosestPointLabel &other)
:
  ClosestPoint(other),
  _TargetComponent(other._TargetComponent),
  _SourceComponent(other._SourceComponent)
{
  if (other._TargetIds) {
    _TargetIds = vtkSmartPointer<vtkIdTypeArray>::New();
    _TargetIds->DeepCopy(other._TargetIds);
  }
  if (other._SourceIds) {
    _SourceIds = vtkSmartPointer<vtkIdTypeArray>::New();
    _SourceIds->DeepCopy(other._SourceIds);
  }
}

// -----------------------------------------------------------------------------
PointCorrespondence *ClosestPointLabel::NewInstance() const
{
  return new ClosestPointLabel(*this);
}

// -----------------------------------------------------------------------------
ClosestPointLabel::~ClosestPointLabel()
{
}

// -----------------------------------------------------------------------------
ClosestPointLabel::TypeId ClosestPointLabel::Type() const
{
  return TypeId::ClosestPointLabel;
}

// =============================================================================
// Correspondences
// =============================================================================

// -----------------------------------------------------------------------------
void ClosestPointLabel::Init()
{
  if (_FromTargetToSource) {
    UpdateClosestLabeledPoints(_Target, _Source, _SourceIds, _SourceComponent);
  }
  if (_FromSourceToTarget) {
    UpdateClosestLabeledPoints(_Source, _Target, _TargetIds, _TargetComponent);
  }
}

// -----------------------------------------------------------------------------
void ClosestPointLabel::Initialize()
{
  // Initialize base class
  ClosestPoint::Initialize();

  // Initialize this class
  ClosestPointLabel::Init();
}

// -----------------------------------------------------------------------------
void ClosestPointLabel::Reinitialize()
{
  // Reinitialize base class
  ClosestPoint::Reinitialize();

  // Reinitialize this class
  ClosestPointLabel::Init();
}

// -----------------------------------------------------------------------------
void ClosestPointLabel::Update()
{
  // Update point features
  PointCorrespondence::Update();

  // Find closest points with identical label
  if (_FromTargetToSource) {
    UpdateCorrespondences(_Target, _Source, _SourceIds, _SourceComponent, _SourceIndex, _SourceDistance);
  }
  if (_FromSourceToTarget) {
    UpdateCorrespondences(_Source, _Target, _TargetIds, _TargetComponent, _TargetIndex, _TargetDistance);
  }
}


} // namespace mirtk

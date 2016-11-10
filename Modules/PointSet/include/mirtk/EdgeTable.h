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

#ifndef MIRTK_EdgeTable_H
#define MIRTK_EdgeTable_H

#include "mirtk/SparseMatrix.h"

#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"

#include "vtkSmartPointer.h"
#include "vtkDataSet.h"


namespace mirtk {


/**
 * Edge table / adjacency matrix
 *
 * This class represents the adjacency matrix of point set nodes. It provides
 * efficient access to the set of nodes adjacent to a given point. The non-zero
 * entries stored in the sparse matrix are the one-based edge IDs such that
 * sparse matrix entries are non-zero. To efficiently iterate all edges or a
 * subset of these, use the thread-safe EdgeIterator.
 */
class EdgeTable : public GenericSparseMatrix<int>
{
  mirtkObjectMacro(EdgeTable);

  /// Pointer to data set which this edge table was computed from
  mirtkReadOnlyAttributeMacro(vtkSmartPointer<vtkDataSet>, Mesh);

  /// Number of undirected edges
  mirtkReadOnlyAttributeMacro(int, NumberOfEdges);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const EdgeTable &);

public:

  /// Construct edge table for given dataset
  EdgeTable(vtkDataSet * = nullptr);

  /// Copy constructor
  EdgeTable(const EdgeTable &);

  /// Assignment operator
  EdgeTable &operator =(const EdgeTable &);

  /// Destructor
  virtual ~EdgeTable();

  /// Initialize edge table from given dataset
  void Initialize(vtkDataSet *);

  /// Clear edge table
  virtual void Clear();

  /// Number of nodes
  int NumberOfPoints() const;

  /// Determine whether two nodes are connected by an edge
  template <class IdType1, class IdType2>
  bool IsEdge(IdType1, IdType2) const;

  /// Get ID of undirected edge connecting two nodes
  ///
  /// \return Zero-based edge ID or -1 if edge does not exist.
  template <class IdType1, class IdType2>
  int EdgeId(IdType1, IdType2) const;

  /// Get IDs of edge nodes (ptId1 < ptId2)
  template <class EdgeIdType, class IdType1, class IdType2>
  bool GetEdge(EdgeIdType, IdType1 &ptId1, IdType2 &ptId2) const;

  /// Get number of adjacent points
  template <class IdType>
  int NumberOfAdjacentPoints(IdType) const;

  /// Get maximum number of adjacent points
  int MaxNumberOfAdjacentPoints() const;

  /// Access list of adjacent nodes (thread-safe)
  void GetAdjacentPoints(int, int &, const int *&) const;

  /// Get start and end pointer into list of adjacent nodes (thread-safe)
  void GetAdjacentPoints(int, const int *&, const int *&) const;
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int EdgeTable::NumberOfPoints() const
{
  return Rows();
}

// -----------------------------------------------------------------------------
template <class IdType1, class IdType2>
inline bool EdgeTable::IsEdge(IdType1 ptId1, IdType2 ptId2) const
{
  return static_cast<bool>(Get(static_cast<int>(ptId1), static_cast<int>(ptId2)));
}

// -----------------------------------------------------------------------------
template <class IdType1, class IdType2>
inline int EdgeTable::EdgeId(IdType1 ptId1, IdType2 ptId2) const
{
  return Get(static_cast<int>(ptId1), static_cast<int>(ptId2)) - 1;
}

// -----------------------------------------------------------------------------
template <class EdgeIdType, class IdType1, class IdType2>
inline bool EdgeTable::GetEdge(EdgeIdType edgeId, IdType1 &ptId1, IdType2 &ptId2) const
{
  const int eid = static_cast<int>(edgeId) + 1;
  if (eid <= 0) return false;

  int i, j;
  const int *row = _Row;
  const int *col = _Col;
  if (_Layout == CRS) swap(row, col);

  const IdType2 npoints = static_cast<IdType2>(NumberOfPoints());
  for (ptId2 = 0; ptId2 < npoints; ++ptId2) {
    for (i = col[ptId2], j = col[ptId2 + 1]; i < j; ++i) {
      ptId1 = static_cast<IdType1>(row[i]);
      if (_Data[i] == eid) return true;
      if (ptId1 > ptId2) break;
    }
  }

  ptId1 = static_cast<IdType1>(-1);
  ptId2 = static_cast<IdType2>(-1);
  return false;
}

// -----------------------------------------------------------------------------
template <class IdType>
inline int EdgeTable::NumberOfAdjacentPoints(IdType ptId) const
{
  if (_Layout == CCS) {
    return _Col[ptId+1] - _Col[ptId];
  } else {
    return _Row[ptId+1] - _Row[ptId];
  }
}

// -----------------------------------------------------------------------------
inline void EdgeTable::GetAdjacentPoints(int ptId, int &numAdjPts, const int *&adjPtIds) const
{
  if (_Layout == CCS) {
    numAdjPts = _Col[ptId+1] - _Col[ptId];
    adjPtIds  = _Row + _Col[ptId];
  } else {
    numAdjPts = _Row[ptId+1] - _Row[ptId];
    adjPtIds  = _Col + _Row[ptId];
  }
}

// -----------------------------------------------------------------------------
inline void EdgeTable::GetAdjacentPoints(int ptId, const int *&begin, const int *&end) const
{
  if (_Layout == CCS) {
    begin = _Row + _Col[ptId];
    end   = _Row + _Col[ptId + 1];
  } else {
    begin = _Col + _Row[ptId];
    end   = _Col + _Row[ptId + 1];
  }
}

////////////////////////////////////////////////////////////////////////////////
// EdgeIterator
////////////////////////////////////////////////////////////////////////////////

/**
 * Thread-safe helper class for iteration of edges
 */
class EdgeIterator
{
  const EdgeTable &_Table;
  int              _EdgeId;
  int              _EndId;
  const int       *_PointId1;
  const int       *_ListEnd;
  int              _PointId2;

public:

  /// Constructor
  EdgeIterator(const EdgeTable &table)
  :
    _Table(table), _EdgeId(-1), _EndId(-1),
    _PointId1(NULL), _ListEnd(NULL), _PointId2(-1)
  {}

  /// Constructor
  EdgeIterator(const EdgeTable &table, int begin, int end = -1)
  :
    EdgeIterator(table)
  {
    InitTraversal(begin, end);
  }

  /// Initialize traversal of edges
  ///
  /// @param[in] begin ID of first edge to iterate.
  /// @param[in] end   ID one behind last edge to iterate. If negative,
  ///                  all edges from @p begin until the end are iterated.
  ///
  /// @bug Not sure if this functions works correctly for begin > 0, cf.,
  ///      AverageEdgeLength implementation.
  inline void InitTraversal(int begin = 0, int end = -1)
  {
    _EdgeId = begin;
    _EndId  = ((end < 0 || end > _Table.NumberOfEdges()) ? _Table.NumberOfEdges() : end);
    if (_EndId > _EdgeId) {
      for (_PointId2 = 0; _PointId2 < _Table.NumberOfPoints(); ++_PointId2) {
        _Table.GetAdjacentPoints(_PointId2, _PointId1, _ListEnd);
        while (_PointId1 != _ListEnd) {
          if (*_PointId1 > _PointId2) {
            _PointId1 = _ListEnd;
            break;
          }
          if (begin == 0) break;
          ++_PointId1, --begin;
        }
        if (begin == 0 && _PointId1 != _ListEnd) break;
      }
    } else {
      _PointId1 = _ListEnd = NULL;
      _PointId2 = -1;
    }
  }

  /// Initialize traversal of edges in range
  template <class IdType>
  inline void InitTraversal(const blocked_range<IdType> &re)
  {
    InitTraversal(static_cast<int>(re.begin()), static_cast<int>(re.end()));
  }

  /// Get next edge
  ///
  /// @param[out] ptId1 ID of first edge point.
  /// @param[out] ptId2 ID of second edge point.
  ///
  /// @return ID of undirected edge (ptId1, ptId2) or -1 if end reached
  template <class IdType>
  int GetNextEdge(IdType &ptId1, IdType &ptId2)
  {
    if (_EdgeId >= _EndId) return -1;
    ptId1 = static_cast<IdType>(*_PointId1);
    ptId2 = static_cast<IdType>(_PointId2);
    int edgeId = _EdgeId;
    if (++_EdgeId < _EndId) {
      if (++_PointId1 == _ListEnd || (*_PointId1) > _PointId2) {
        do {
          _Table.GetAdjacentPoints(++_PointId2, _PointId1, _ListEnd);
        } while (_PointId1 == _ListEnd || (*_PointId1) > _PointId2);
      }
    }
    return edgeId;
  }
};


} // namespace mirtk

#endif // MIRTK_EdgeTable_H

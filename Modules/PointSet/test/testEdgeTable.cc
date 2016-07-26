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

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/EdgeTable.h"
#include "mirtk/PointSetIO.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"

#include "gtest/gtest.h"

using namespace mirtk;


// =============================================================================
// Global variables / test arguments
// =============================================================================

const char *input_name = nullptr;

// =============================================================================
// Tests
// =============================================================================

// -----------------------------------------------------------------------------
TEST(EdgeTable, Initialize)
{
  // Read input point set
  ASSERT_TRUE(input_name != nullptr) << "input point set given";
  vtkSmartPointer<vtkPolyData> input = ReadPolyData(input_name, false);
  const int npoints = static_cast<int>(input->GetNumberOfPoints());
  ASSERT_GT(npoints, 0) << "input point set contains points";

  // Construct edge table
  EdgeTable edgeTable(input);
  ASSERT_EQ(EdgeTable::CCS, edgeTable.Layout()) << "sparse matrix layout is CCS";
  ASSERT_EQ(npoints, edgeTable.Rows()) << "one row for each point";
  ASSERT_EQ(npoints, edgeTable.Cols()) << "one column for each point";
  ASSERT_TRUE(edgeTable.IsSymmetric()) << "edge table is a symmetric matrix";

  // All upper triangular entries are unique and in ascending order
  // and verify that CCS sparse matrix entries are sorted with unique row entries
  int edgeId = 0;
  OrderedSet<int> rows;
  EdgeTable::Entries entries;
  EdgeTable::Entries::const_iterator it;
  OrderedSet<int>::const_iterator rowsIt;
  for (int c = 0, r; c < npoints; ++c) {
    edgeTable.GetCol(c, entries);
    rows.clear();
    for (it = entries.begin(); it != entries.end(); ++it) {
      rows.insert(it->first);
    }
    ASSERT_EQ(rows.size(), entries.size()) << "row entries are unique";
    for (it = entries.begin(), rowsIt = rows.begin(); it != entries.end(); ++it, ++rowsIt) {
      ASSERT_EQ(*rowsIt, it->first) << "row entries are sorted";
    }
    for (it = entries.begin(); it != entries.end(); ++it) {
      r = it->first;
      if (r > c) break;
      ASSERT_EQ(++edgeId, it->second);
    }
  }
  ASSERT_EQ(edgeId, edgeTable.NumberOfEdges()) << "all edges have a unique ID";
}

// -----------------------------------------------------------------------------
TEST(EdgeTable, EdgeIterator)
{
  vtkIdType pid1, pid2;

  // Read input point set
  ASSERT_TRUE(input_name != nullptr) << "input point set given";
  vtkSmartPointer<vtkPolyData> input = ReadPolyData(input_name, false);
  ASSERT_TRUE(input->GetNumberOfPoints() > 0) << "input point set contains points";

  // Construct edge table and iterator
  EdgeTable edgeTable(input);
  EdgeIterator it(edgeTable);

  // Check that each edge is returned by iterator exactly one time
  // and that the order of edge IDs is the ascening sequence of positive integers
  it.InitTraversal();
  int edgeId, edgeCount = 0;
  while ((edgeId = it.GetNextEdge(pid1, pid2)) != -1) {
    ASSERT_EQ(edgeCount++, edgeId) << "edges are iterated by ascending edge ID";
    ASSERT_EQ(edgeId, edgeTable.EdgeId(pid1, pid2))
        << "EdgeTable::EdgeId consistent with ID returned by EdgeTable::GetNextEdge";
    ASSERT_LE(pid1, pid2) << "order of returned point IDs is such that ptId1 <= ptId2";
  }
  ASSERT_EQ(edgeCount, edgeTable.NumberOfEdges()) << "all edges are being iterated";
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  input_name = argv[1];
  return RUN_ALL_TESTS();
}

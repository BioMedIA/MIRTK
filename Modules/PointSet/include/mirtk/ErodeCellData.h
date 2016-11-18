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

#ifndef MIRTK_ErodeCellData_H
#define MIRTK_ErodeCellData_H

#include "mirtk/CellDataFilter.h"


namespace mirtk {


/**
 * Component-wise morphological erosion of mesh cell data
 */
class ErodeCellData : public CellDataFilter
{
  mirtkObjectMacro(ErodeCellData);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Number of iterations
  mirtkPublicAttributeMacro(int, Iterations);

  /// Copy attributes of this filter from another instance
  void CopyAttributes(const ErodeCellData &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  ErodeCellData();

  /// Copy constructor
  ErodeCellData(const ErodeCellData &);

  /// Assignment operator
  ErodeCellData &operator =(const ErodeCellData &);

  /// Destructor
  virtual ~ErodeCellData();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Execute filter
  virtual void Execute();

};


} // namespace mirtk

#endif // MIRTK_ErodeCellData_H

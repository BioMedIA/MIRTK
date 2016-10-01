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

#ifndef MIRTK_DataSelection_H
#define MIRTK_DataSelection_H

#include "mirtk/Object.h"

#include "mirtk/UnorderedSet.h"
#include "mirtk/Algorithm.h"


namespace mirtk { namespace data {

/// Set of selected data point IDs
typedef UnorderedSet<int> Selection;

// -----------------------------------------------------------------------------
/// Base class of data selectors
class Selector : public Object
{
  mirtkAbstractMacro(Selector);

public:

  /// Destructor
  virtual ~Selector() {};

  /// Run selection query and return IDs of selected data points
  ///
  /// \param[in] values Values of the n data points.
  ///
  /// \returns Set of IDs of the selected data points.
  virtual Selection Evaluate(const Array<double> &values) const = 0;
};

// -----------------------------------------------------------------------------
/// Combine one or more data selection criteria using a logical operator
class LogicalOp : public Selector
{
  mirtkAbstractMacro(LogicalOp);

  typedef SharedPtr<const Selector> SelectorPointer;
  typedef List<SelectorPointer>     SelectorList;

protected:

  /// Individual data selectors
  mirtkAttributeMacro(SelectorList, Criteria);

public:

  /// Number of data selection criteria
  int NumberOfCriteria() const
  {
    return static_cast<int>(_Criteria.size());
  }

  /// Get n-th data selection criterium
  SelectorPointer Criterium(int i) const
  {
    auto it = _Criteria.begin();
    for (int pos = 0; pos < i; ++pos) ++it;
    return *it;
  }

  /// Add data selection criterium
  void Push(const SelectorPointer &criterium)
  {
    _Criteria.push_back(criterium);
  }
};

// -----------------------------------------------------------------------------
/// Combine one or more data selection criteria using a logical AND
class LogicalAnd : public LogicalOp
{
  mirtkObjectMacro(LogicalAnd);

public:

  /// Run selection query and return IDs of selected data points
  virtual Selection Evaluate(const Array<double> &values) const
  {
    if (_Criteria.empty()) return Selection();
    auto criterium = _Criteria.begin();
    auto selection = (*criterium)->Evaluate(values);
    while (++criterium != _Criteria.end()) {
      selection = Intersection(selection, (*criterium)->Evaluate(values));
    }
    return selection;
  }
};

// -----------------------------------------------------------------------------
/// Combine one or more data selection criteria using a logical OR
class LogicalOr : public LogicalOp
{
  mirtkObjectMacro(LogicalOr);

public:

  /// Run selection query and return IDs of selected data points
  virtual Selection Evaluate(const Array<double> &values) const
  {
    Selection selection;
    for (auto criterium : _Criteria) {
      auto ids = criterium->Evaluate(values);
      for (auto id : ids) {
        selection.insert(id);
      }
    }
    return selection;
  }
};

// -----------------------------------------------------------------------------
/// Base class of data selection criteria
class SelectionCriterium : public Selector
{
  mirtkAbstractMacro(SelectionCriterium);

public:

  /// Run selection query and return IDs of selected data points
  virtual Selection Evaluate(const Array<double> &values) const
  {
    Selection selection;
    selection.reserve(values.size());
    const int n = static_cast<int>(values.size());
    for (int id = 0; id < n; ++id) {
      if (this->Select(values[id])) {
        selection.insert(id);
      }
    }
    return selection;
  }

  /// Evaluate criterium for given data value
  virtual bool Select(double) const = 0;
};

// =============================================================================
// Data selection criteria
// =============================================================================

namespace select {


// -----------------------------------------------------------------------------
/// Select data points with a value equal the specified value
class Equal : public SelectionCriterium
{
  mirtkObjectMacro(Equal);

  /// Upper data value threshold
  mirtkPublicAttributeMacro(double, Value);

public:

  /// Constructor
  Equal(double value) : _Value(value) {}

  /// Evaluate criterium for given data value
  virtual bool Select(double value) const
  {
    return fequal(value, _Value);
  }
};

// -----------------------------------------------------------------------------
/// Select data points a value different from the specified value
class NotEqual : public SelectionCriterium
{
  mirtkObjectMacro(NotEqual);

  /// Upper data value threshold
  mirtkPublicAttributeMacro(double, Value);

public:

  /// Constructor
  NotEqual(double value) : _Value(value) {}

  /// Evaluate criterium for given data value
  virtual bool Select(double value) const
  {
    return !fequal(value, _Value);
  }
};

// -----------------------------------------------------------------------------
/// Select data points with a value less than the specified threshold
class LessThan : public SelectionCriterium
{
  mirtkObjectMacro(LessThan);

  /// Upper data value threshold
  mirtkPublicAttributeMacro(double, Threshold);

public:

  /// Constructor
  LessThan(double threshold) : _Threshold(threshold) {}

  /// Evaluate criterium for given data value
  virtual bool Select(double value) const
  {
    return value < _Threshold;
  }
};

// -----------------------------------------------------------------------------
/// Select data points with a value less than or equal the specified threshold
class LessOrEqual : public SelectionCriterium
{
  mirtkObjectMacro(LessOrEqual);

  /// Upper data value threshold
  mirtkPublicAttributeMacro(double, Threshold);

public:

  /// Constructor
  LessOrEqual(double threshold) : _Threshold(threshold) {}

  /// Evaluate criterium for given data value
  virtual bool Select(double value) const
  {
    return value <= _Threshold;
  }
};

// -----------------------------------------------------------------------------
/// Select data points with a value greater than the specified threshold
class GreaterThan : public SelectionCriterium
{
  mirtkObjectMacro(GreaterThan);

  /// Lower data value threshold
  mirtkPublicAttributeMacro(double, Threshold);

public:

  /// Constructor
  GreaterThan(double threshold) : _Threshold(threshold) {}

  /// Evaluate criterium for given data value
  virtual bool Select(double value) const
  {
    return value > _Threshold;
  }
};

// -----------------------------------------------------------------------------
/// Select data points with a value greater than or equal the specified threshold
class GreaterOrEqual : public SelectionCriterium
{
  mirtkObjectMacro(GreaterOrEqual);

  /// Lower data value threshold
  mirtkPublicAttributeMacro(double, Threshold);

public:

  /// Constructor
  GreaterOrEqual(double threshold) : _Threshold(threshold) {}

  /// Evaluate criterium for given data value
  virtual bool Select(double value) const
  {
    return value >= _Threshold;
  }
};


} } } // namespace mirtk::data::select

#endif // MIRTK_DataSelection_H

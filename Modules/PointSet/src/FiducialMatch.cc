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

#include "mirtk/FiducialMatch.h"

#include "mirtk/Point.h"
#include "mirtk/Array.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
FiducialMatch::FiducialMatch()
{
}

// -----------------------------------------------------------------------------
FiducialMatch
::FiducialMatch(const FiducialMatch &other)
:
  PointCorrespondence(other)
{
}

// -----------------------------------------------------------------------------
PointCorrespondence *FiducialMatch::NewInstance() const
{
  return new FiducialMatch(*this);
}

// -----------------------------------------------------------------------------
FiducialMatch::~FiducialMatch()
{
}

// -----------------------------------------------------------------------------
FiducialMatch::TypeId FiducialMatch::Type() const
{
  return TypeId::FiducialMatch;
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool FiducialMatch::Set(const char *name, const char *value)
{
  if (strcmp(name, "Correspondence map") == 0) {
    _CorrespondenceMap = value;
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
ParameterList FiducialMatch::Parameter() const
{
  ParameterList params;
  Insert(params, "Correspondence map", _CorrespondenceMap);
  return params;
}

// =============================================================================
// Correspondences
// =============================================================================

// -----------------------------------------------------------------------------
void FiducialMatch::ValidateCorrespondenceMap(const RegisteredPointSet *target,
                                              const RegisteredPointSet *source,
                                              const Array<int>         &map,
                                              const char               *map_name) const
{
  const int m = static_cast<int>(target->NumberOfPoints());
  const int n = static_cast<int>(source->NumberOfPoints());
  if (map.size() < static_cast<size_t>(m)) {
    cerr << NameOfType() << "::Initialize: Correspondence map has no entries for target indices t >= " << map.size() << endl;
    if (map_name && map_name[0] != '\0') cerr << "  Correspondence map input file: " << map_name << endl;
    exit(1);
  }
  if (map.size() > static_cast<size_t>(m)) {
    cerr << NameOfType() << "::Initialize: Correspondence map has additional entries for unused indices t >= " << m << endl;
    if (map_name && map_name[0] != '\0') cerr << "  Correspondence map input file: " << map_name << endl;
  }
  for (int t = 0; t < m; ++t) {
    if (map[t] == -1) {
      cerr << NameOfType() << "::Initialize: Missing correspondence map entry for t=" << t << endl;
      if (map_name && map_name[0] != '\0') cerr << "  Correspondence map input file: " << map_name << endl;
      exit(1);
    } else if (map[t] < 0 || map[t] >= n) {
      cerr << NameOfType() << "::Initialize: Invalid correspondence map entry: t=" << t << ", s=" << map[t] << endl;
      if (map_name && map_name[0] != '\0') cerr << "  Correspondence map input file: " << map_name << endl;
      exit(1);
    }
  }
}

// -----------------------------------------------------------------------------
void FiducialMatch::InvertCorrespondenceMap(const RegisteredPointSet *target,
                                            const RegisteredPointSet *source,
                                            const Array<int>         &map,
                                            Array<int>               &inv) const
{
  const int m = static_cast<int>(target->NumberOfPoints());
  const int n = static_cast<int>(source->NumberOfPoints());
  inv.clear();
  inv.resize(n, -1);
  for (int t = 0; t < m; ++t) {
    inv[map[t]] = t;
  }
}

// -----------------------------------------------------------------------------
void FiducialMatch::Initialize()
{
  // Initialize base class
  PointCorrespondence::Initialize();

  const int m = static_cast<int>(_Target->NumberOfPoints());
  const int n = static_cast<int>(_Source->NumberOfPoints());

  // Either read corresponding indices from input file
  if (!_CorrespondenceMap.empty()) {
    ifstream ifs(_CorrespondenceMap);
    if (!ifs.is_open()) {
      cerr << NameOfType()
                << "::Initialize: Failed to open correspondence input file: "
                << _CorrespondenceMap << endl;
      exit(1);
    }
    string line;
    int          t, s;
    _SourceIndex.resize(m, -1);
    while (getline(ifs, line)) {
      istringstream is(line);
      if (!(is >> t >> s)) {
        cerr << NameOfType()
                  << "::Initialize: Failed to read correspondence map from file: "
                  << _CorrespondenceMap << endl;
        exit(1);
      }
      if (t < 0 || t >= m) {
        cerr << NameOfType()
                  << "::Initialize: Invalid target index in correspondence map: "
                  << _CorrespondenceMap << endl;
        exit(1);
      }
      _SourceIndex[t] = s;
    }
    ifs.close();
  // Or ensure that both data sets contain same number of fiducial points
  } else {
    if (m != n) {
      cerr << NameOfType() << "::Initialize: Data sets must have same number of fiducial points" << endl;
      exit(1);
    }
    _SourceIndex.clear();
    _TargetIndex.clear();
  }

  if (!_SourceIndex.empty()) {
    // Check correspondence map
    ValidateCorrespondenceMap(_Target, _Source, _SourceIndex, _CorrespondenceMap.c_str());
    // Compute inverse correspondence map
    InvertCorrespondenceMap(_Target, _Source, _SourceIndex, _TargetIndex);
    // Check inverse correspondence map
    ValidateCorrespondenceMap(_Source, _Target, _TargetIndex);
  }
}

// -----------------------------------------------------------------------------
bool FiducialMatch::GetInputTargetPoint(int i, Point &p) const
{
  if (_SourceSample) i = (*_SourceSample)[i];
  if (!_TargetIndex.empty()) i = _TargetIndex[i];
  _Target->GetInputPoint(i, p);
  return true;
}

// -----------------------------------------------------------------------------
bool FiducialMatch::GetInputSourcePoint(int i, Point &p) const
{
  if (_TargetSample) i = (*_TargetSample)[i];
  if (!_SourceIndex.empty()) i = _SourceIndex[i];
  _Source->GetInputPoint(i, p);
  return true;
}

// -----------------------------------------------------------------------------
bool FiducialMatch::GetTargetPoint(int i, Point &p) const
{
  if (_SourceSample) i = (*_SourceSample)[i];
  if (!_TargetIndex.empty()) i = _TargetIndex[i];
  _Target->GetPoint(i, p);
  return true;
}

// -----------------------------------------------------------------------------
bool FiducialMatch::GetSourcePoint(int i, Point &p) const
{
  if (_TargetSample) i = (*_TargetSample)[i];
  if (!_SourceIndex.empty()) i = _SourceIndex[i];
  _Source->GetPoint(i, p);
  return true;
}

// -----------------------------------------------------------------------------
int FiducialMatch::GetTargetIndex(int i) const
{
  if (_SourceSample) i = (*_SourceSample)[i];
  if (!_TargetIndex.empty()) i = _TargetIndex[i];
  return i;
}

// -----------------------------------------------------------------------------
int FiducialMatch::GetSourceIndex(int i) const
{
  if (_TargetSample) i = (*_TargetSample)[i];
  if (!_SourceIndex.empty()) i = _SourceIndex[i];
  return i;
}


} // namespace mirtk

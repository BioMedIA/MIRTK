/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2015 Imperial College London
 * Copyright 2015 Andreas Schuh
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

#include "mirtk/Configurable.h"

#include "mirtk/Array.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
Configurable::Configurable(const char *name)
:
  _Name(name)
{
}

// -----------------------------------------------------------------------------
void Configurable::CopyAttributes(const Configurable &other)
{
  _Name            = other._Name;
  _ParameterPrefix = other._ParameterPrefix;
}

// -----------------------------------------------------------------------------
Configurable::Configurable(const Configurable &other)
:
  Observable(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
Configurable &Configurable::operator =(const Configurable &other)
{
  if (this != &other) {
    Observable::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
Configurable::~Configurable()
{
}

// -----------------------------------------------------------------------------
string Configurable::DefaultName() const
{
  string prefix(this->NameOfClass());
  if (prefix.compare(0, 4, "irtk") == 0) prefix.erase(0, 4);
  for (string::size_type i = 1; i < prefix.length(); ++i) {
    if ('A' <= prefix[i] && prefix[i] <= 'Z') {
      prefix[i] = 'a' + (prefix[i] - 'A');
      prefix.insert(i, " ");
      ++i;
    }
  }
  return prefix;
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool Configurable::HasName() const
{
  return !_Name.empty();
}

// -----------------------------------------------------------------------------
bool Configurable::HasPrefix() const
{
  return !_Name.empty() || !_ParameterPrefix.empty();
}

// -----------------------------------------------------------------------------
string Configurable::DefaultPrefix() const
{
  if (!_Name.empty()) return _Name + " ";
  return _ParameterPrefix.empty() ? string() : _ParameterPrefix.front();
}

// -----------------------------------------------------------------------------
string Configurable::ParameterNameWithPrefix(const string &str) const
{
  string name = DefaultPrefix();
  if (name.empty()) return str;
  name += ::tolower(str[0]);
  name += str.substr(1);
  return name;
}

// -----------------------------------------------------------------------------
string Configurable::ParameterNameWithPrefix(const char *str) const
{
  string name = DefaultPrefix();
  if (name.empty()) return str;
  name += ::tolower(str[0]);
  name += (str + 1);
  return name;
}

// -----------------------------------------------------------------------------
string Configurable::ParameterNameWithoutPrefix(const char *str) const
{
  string        param_name;
  Array<string> prefixes = _ParameterPrefix;
  if (!_Name.empty()) prefixes.push_back(_Name + " ");

  const size_t len = strlen(str);
  Array<string>::const_iterator prefix;
  for (prefix = prefixes.begin(); prefix != prefixes.end(); ++prefix) {
    if (len > prefix->length() &&
        strncmp(str, prefix->c_str(), prefix->length()) == 0 &&
        *(str + prefix->length()) != '\0') {
      param_name = str + prefix->length();
      if (islower(param_name[0])) param_name[0] = toupper(param_name[0]);
      break;
    }
  }

  return param_name;
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool Configurable::SetWithPrefix(const char *, const char *)
{
  return false;
}

// -----------------------------------------------------------------------------
bool Configurable::SetWithoutPrefix(const char *, const char *)
{
  return false;
}

// -----------------------------------------------------------------------------
bool Configurable::Set(const char *param, const char *value)
{
  if (this->SetWithPrefix(param, value)) return true;
  const string name = ParameterNameWithoutPrefix(param);
  return !name.empty() && this->SetWithoutPrefix(name.c_str(), value);
}


} // namespace mirtk

/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
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

#include "mirtk/String.h"

#include "mirtk/Array.h"
#include "mirtk/Stream.h"
#include "mirtk/Algorithm.h"


namespace mirtk {


// ========================================================================
// Read from file
// ========================================================================

#define MAX_LINE 255

// ------------------------------------------------------------------------
int ReadInt(ifstream &in)
{
  char c;
  char *string, *s;
  int data;

  string = new char[MAX_LINE];
  s = string;
  do {
    in.get(string, MAX_LINE, '\n');
    in.get(c);
  } while ((strlen(string) == 0) || (string[0] == '#'));
  if ((string = strchr(string, '=')) == NULL) {
    cerr << "ReadInt: No valid line format" << endl;
    exit(1);
  }
  do {
    string++;
  } while ((*string == ' ') || (*string == '\t'));

  data = atoi(string);
  delete[] s;
  return (data);
}

// ------------------------------------------------------------------------
float ReadFloat(ifstream &in)
{
  char c;
  char *string, *s;
  float data;

  string = new char[MAX_LINE];
  s = string;
  do {
    in.get(string, MAX_LINE, '\n');
    in.get(c);
  } while ((strlen(string) == 0) || (string[0] == '#'));
  if ((string = strchr(string, '=')) == NULL) {
    cerr << "ReadFloat: No valid line format" << endl;
    exit(1);
  }
  do {
    string++;
  } while ((*string == ' ') || (*string == '\t'));

  data = static_cast<float>(atof(string));
  delete[] s;
  return (data);
}

// ------------------------------------------------------------------------
char *ReadString(ifstream &in)
{
  char c;
  char *string;

  string = new char[MAX_LINE];
  do {
    in.get(string, MAX_LINE, '\n');
    in.get(c);
  } while ((strlen(string) == 0) || (string[0] == '#'));
  if ((string = strchr(string, '=')) == NULL) {
    cerr << "ReadString: No valid line format" << endl;
    exit(1);
  }
  do {
    string++;
  } while ((*string == ' ') || (*string == '\t'));

  return (string);
}

// ------------------------------------------------------------------------
string ToLower(const string &s)
{
  string o(s);
  transform(s.begin(), s.end(), o.begin(), tolower);
  return o;
}

// ------------------------------------------------------------------------
string ToUpper(const string &s)
{
  string o(s);
  transform(s.begin(), s.end(), o.begin(), toupper);
  return o;
}

// ------------------------------------------------------------------------
string Trim(const string &str, const string &what)
{
  const auto pos = str.find_first_not_of(what);
  if (pos == string::npos) return "";
  const auto end = str.find_last_not_of(what);
  return str.substr(pos, end - pos + 1u);
}

// ------------------------------------------------------------------------
string TrimAll(const string &str, const string &what)
{
  string res;
  size_t pos = 0, end;
  while ((pos = str.find_first_not_of(what, pos)) != string::npos) {
    end  = str.find_first_of(what, pos);
    res += str.substr(pos, end - pos);
    pos  = end;
  }
  return res;
}

// ------------------------------------------------------------------------
Array<string> Split(string s, const char *d, int n, bool discard_empty, bool handle_quotes)
{
  size_t i, j, pos, start, end;
  Array<string> parts;
  const size_t delimiter_length = strlen(d);
  string whitespace(" \t\f\v\n\r");
  if (delimiter_length == 0) {
    parts.push_back(s);
    return parts;
  }
  pos = whitespace.find(d[0]);
  if (pos != string::npos) {
    whitespace.erase(pos, 1);
  }
  if (n < 0 && !handle_quotes) {
    n = abs(n);
    while (parts.size() < static_cast<size_t>(n)) {
      start = s.rfind(d);
      if (start == string::npos) {
        parts.push_back(s);
        break;
      }
      parts.push_back(s.substr(start + delimiter_length));
      if (discard_empty && parts.back().empty()) {
        parts.pop_back();
      }
      s.erase(start);
    }
    reverse(parts.begin(), parts.end());
  } else {
    start = 0;
    size_t m = static_cast<size_t>(n > 0 ? n : 0);
    while (start < s.length() && (m == 0 || parts.size() < m)) {
      i = j = string::npos;
      if (handle_quotes) {
        i = s.find_first_not_of(whitespace, start);
        if (i != string::npos && (s[i] == '\"' || s[i] == '\'')) {
          j = i;
          do {
            j = s.find(s[i], ++j);
            if (j == string::npos) break;
            pos = s.find_first_not_of(whitespace, j + 1);
            end = s.find(d, j + 1);
          } while (pos < end);
          ++i;
        }
      }
      if (i == string::npos || j == string::npos) {
        i = start;
        j = end = s.find(d, start);
      }
      if (j != string::npos) j -= i;
      parts.push_back(s.substr(i, j));
      if (discard_empty && parts.back().empty()) {
        parts.pop_back();
      }
      if (end == string::npos) break;
      start = end + delimiter_length;
    }
    if (n < 0) {
      parts.erase(parts.begin(), parts.end() + n);
    }
  }
  return parts;
}

// ------------------------------------------------------------------------
Array<string> Split(string s, char d, int n, bool discard_empty, bool handle_quotes)
{
  const char delim[2] = {d, '\0'};
  return Split(s, delim, n, discard_empty, handle_quotes);
}

// ------------------------------------------------------------------------
string CamelCaseToPrettyParameterName(const string &s)
{
  string param;
  param.reserve(s.length() + 10);
  for (auto c = s.begin(); c != s.end(); ++c) {
    if (param.empty()) {
      param += toupper(*c);
    } else {
      if (isupper(*c)) param += ' ';
      param += tolower(*c);
    }
  }
  param.shrink_to_fit();
  return param;
}

// ------------------------------------------------------------------------
string StandardUnits(const string &str)
{
  if (str.empty()) return str;
  string units = Trim(ToLower(str));
  if (units.size() > 2u && units.front() == '[' && units.back() == ']') {
    units = units.substr(1, units.size() - 2u);
    if (units.empty()) return units;
  }
  if (units == "voxel" || units == "voxels") {
    units = "vox";
  } else if (units == "millimeters" || units == "millimeter") {
    units = "mm";
  } else if (units == "absolute") {
    units = "abs";
  } else if (units == "relative") {
    units = "rel";
  } else if (units == "percentage" || units == "percentile" || units == "pct") {
    units = "%";
  } else if (units == "gauss" || units == "gaussian" || units == "sd" || units == "stddev" || units == "stdev") {
    units = "sigma";
  }
  return units;
}

// ------------------------------------------------------------------------
string ParameterUnits(const string &str, string *name, const char *dflt)
{
  const string trimmed = Trim(str);
  string units;
  Array<string> parts = Split(trimmed, " ", -1);
  if (parts.size() == 1u) {
    string last = parts.back();
    const auto pos = last.rfind('[');
    if (pos != string::npos) last = last.substr(pos);
    if (last.size() > 2u && last.front() == '[' && last.back() == ']') {
      units = last.substr(1u, last.size() - 2u);
    }
  }
  if (name) {
    if (units.empty()) {
      *name = trimmed;
    } else {
      *name = Trim(trimmed.substr(0u, trimmed.size() - units.size() - 2u));
    }
  }
  if (units.empty()) return dflt;
  return StandardUnits(units);
}

// ------------------------------------------------------------------------
string ValueUnits(const string &str, string *value, const char *dflt)
{
  const string trimmed = Trim(str);
  string units;
  Array<string> parts = Split(trimmed, " ", 0, true);
  if (!parts.empty()) {
    const string last = parts.back();
    const auto pos = last.find_first_of("0123456789");
    if (pos != string::npos) {
      const auto pos = last.find_last_of("0123456789");
      parts.back() = last.substr(0u, pos + 1u);
      parts.push_back(last.substr(pos + 1u));
    }
    if (parts.size() > 1u) {
      // When any of the parameter values is not numeric, the units specification
      // must be enclosed in brackets to distinguish it from the parameter value
      double v;
      for (size_t i = 0u; i < parts.size() - 1u; ++i) {
        if (!FromString(parts[i], v)) {
          return ParameterUnits(str, value);
        }
      }
      if (!FromString(parts.back(), v)) {
        units = parts.back();
      }
    }
  }
  if (value) {
    *value = Trim(trimmed.substr(0u, trimmed.size() - units.size()));
  }
  if (units.empty()) return dflt;
  return StandardUnits(units);
}


} // namespace mirtk

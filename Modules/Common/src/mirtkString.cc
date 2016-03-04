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

#include <mirtkString.h>

#include <mirtkArray.h>
#include <mirtkStream.h>
#include <mirtkAlgorithm.h>


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

  data = atof(string);
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
Array<string> Split(string s, const char *d, int n, bool discard_empty)
{
  const size_t delimiter_length = strlen(d);
  Array<string> parts;
  if (n <= 0) {
    n = abs(n);
    size_t pos;
    while (n == 0 || parts.size() < static_cast<size_t>(n)) {
      pos = s.rfind(d);
      if (pos == string::npos) {
        parts.push_back(s);
        break;
      }
      parts.push_back(s.substr(pos + delimiter_length));
      s.erase(pos);
    }
    reverse(parts.begin(), parts.end());
  } else {
    size_t start = 0;
    size_t end   = string::npos;
    while (start < s.length() && (n == 0 || parts.size() < static_cast<size_t>(n))) {
      end = s.find(d, start);
      parts.push_back(s.substr(start, end));
      if (end == string::npos) break;
      start = end + delimiter_length;
    }
  }
  if (discard_empty) {
    for (Array<string>::iterator part = parts.begin(); part != parts.end(); ++part) {
      if (part->empty()) parts.erase(part--);
    }
  }
  return parts;
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
  } else if (units == "percentage") {
    units = "%";
  }
  return units;
}

// ------------------------------------------------------------------------
string ParameterUnits(const string &str, string *name)
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
    *name = Trim(trimmed.substr(0u, trimmed.size() - units.size() - 2u));
  }
  return StandardUnits(units);
}

// ------------------------------------------------------------------------
string ValueUnits(const string &str, string *value)
{
  const string trimmed = Trim(str);
  string units;
  Array<string> parts = Split(trimmed, " ", 0, true);
  if (!parts.empty()) {
    const string last = parts.back();
    const auto pos = last.find_first_of("0123456789");
    if (pos != string::npos) {
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
  return StandardUnits(units);
}


} // namespace mirtk

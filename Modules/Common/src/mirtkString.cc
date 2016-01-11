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
Array<string> Split(string s, const char *d, int n)
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
  return parts;
}


} // namespace mirtk

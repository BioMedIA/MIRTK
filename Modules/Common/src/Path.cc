/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2017 Andreas Schuh
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

#include "mirtk/Path.h"
#include "mirtk/Config.h" // WINDOWS
#include "mirtk/Array.h"
#include "mirtk/Algorithm.h"


namespace mirtk {


// -----------------------------------------------------------------------------
#ifdef WINDOWS
const char PATHSEP = '\\';
#else
const char PATHSEP = '/';
#endif

// -----------------------------------------------------------------------------
string Extension(const char *path, ExtensionMode mode)
{
  // Get base name without directory path
  string s = BaseName(path);
  if (s.empty() || mode == EXT_None) return "";
  // Split file path, note that first part is the file name itself and that
  // an empty initial part is caused by the leading '.' of hidden files on Unix
  Array<string> parts = Split(path, ".");
#ifdef WINDOWS
  const bool hidden = false;
#else
  const bool hidden = (s[0] == '.');
#endif
  parts.erase(parts.begin(), parts.begin() + (hidden ? 2 : 1));
  // If no parts left, return empty string
  if (parts.empty()) return "";
  // Concatenate requested file extension parts only
  s.clear();
  switch (mode) {
    case EXT_Default:
    case EXT_LastWithGz:
      s  = '.';
      s += parts.back();
      s = ToLower(s);
      if (s == ".gz" && parts.size() > 1) {
        s = '.' + parts[parts.size()-2] + s;
      }
      break;
    case EXT_LastWithoutGz:
      s  = '.';
      s += parts.back();
      s = ToLower(s);
      if (s == ".gz" && parts.size() > 1) {
        s = '.' + parts[parts.size()-2];
      }
      break;
    case EXT_Last:
      s  = '.';
      s += ToLower(parts.back());
      break;
    case EXT_All:
      for (size_t i = 0; i < parts.size(); ++i) {
        s += '.';
        s += ToLower(parts[i]);
      }
      break;
    case EXT_None:
      // Dealt with before
      return "";
  }
  return s;
}

// -----------------------------------------------------------------------------
string Extension(const string &path, ExtensionMode mode)
{
  return Extension(path.c_str(), mode);
}

// -----------------------------------------------------------------------------
string Directory(const char *path)
{
  string dir(path);
  size_t n = dir.find_last_of(PATHSEP);
  if (n == string::npos) dir.clear();
  else                   dir.resize(n);
  return dir;
}

// -----------------------------------------------------------------------------
string Directory(const string &path)
{
  return Directory(path.c_str());
}

// -----------------------------------------------------------------------------
string BaseName(const char *path)
{
  const char *name = path;
  while (name[0] != '\0') ++name;
  while (name != path && name[0] != PATHSEP) --name;
  if (name[0] == PATHSEP) ++name;
  return name;
}

// -----------------------------------------------------------------------------
string BaseName(const string &path)
{
  return BaseName(path.c_str());
}

// -----------------------------------------------------------------------------
string FileName(const char *path, ExtensionMode mode)
{
  string name = BaseName(path);
  string ext  = Extension(path, mode);
  return name.substr(0, name.length() - ext.length());
}

// -----------------------------------------------------------------------------
string FileName(const string &path, ExtensionMode mode)
{
  return FileName(path.c_str(), mode);
}

// -----------------------------------------------------------------------------
string FilePrefix(const char *path, ExtensionMode mode)
{
  return FilePrefix(string(path), mode);
}

// -----------------------------------------------------------------------------
string FilePrefix(const string &path, ExtensionMode mode)
{
  string ext = Extension(path, mode);
  return path.substr(0, path.length() - ext.length());
}


} // namespace mirtk

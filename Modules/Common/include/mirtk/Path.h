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

#ifndef MIRTK_Path_H
#define MIRTK_Path_H

#include "mirtk/CommonExport.h"

#include "mirtk/String.h"


namespace mirtk {

// =============================================================================
// Global constants
// =============================================================================

/// Path separating character
MIRTK_Common_EXPORT extern const char PATHSEP;

// =============================================================================
// Split file path
// =============================================================================

/// Enumeration of file path extension retrival modes
enum ExtensionMode
{
  EXT_Default,       ///< Default extension mode
  EXT_None,          ///< No part of the file is considered to be an extension
  EXT_Last,          ///< Last file extension only
  EXT_LastWithGz,    ///< Last file extension possibly plus ".gz"
  EXT_LastWithoutGz, ///< Last file extension with possibly ".gz" removed first
  EXT_All            ///< All file extensions, i.e., everything after first dot (besides leading dot of hidden files on Unix)
};

/// Get file name extension in lower case incl. leading dot ('.')
string Extension(const char *, ExtensionMode = EXT_Default);

/// Get file name extension in lower case incl. leading dot ('.')
string Extension(const string &, ExtensionMode = EXT_Default);

/// Get directory part of file path
string Directory(const char *);

/// Get directory part of file path
string Directory(const string &);

/// Get file name of file path incl. file extension
string BaseName(const char *);

/// Get file name of file path incl. file extension
string BaseName(const string &);

/// Get file name of file path excl. file extension
string FileName(const char *, ExtensionMode = EXT_Default);

/// Get file name of file path excl. file extension
string FileName(const string &, ExtensionMode = EXT_Default);

/// Get file path excl. file extension
string FilePrefix(const char *, ExtensionMode = EXT_Default);

/// Get file path excl. file extension
string FilePrefix(const string &, ExtensionMode = EXT_Default);


} // namespace mirtk

#endif // MIRTK_Path_H

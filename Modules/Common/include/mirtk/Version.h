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

#ifndef MIRTK_Version_H
#define MIRTK_Version_H

#include "mirtk/CommonExport.h"
#include "mirtk/VersionInfo.h"

#include "mirtk/Object.h"


namespace mirtk {


// -----------------------------------------------------------------------------
/// Software version object
class Version : public Object
{
  mirtkObjectMacro(Version);

  mirtkPublicAttributeMacro(unsigned int, Major); ///< Major version number
  mirtkPublicAttributeMacro(unsigned int, Minor); ///< Minor version number
  mirtkPublicAttributeMacro(unsigned int, Patch); ///< Patch number
 
public:

  /// Constructor
  Version(unsigned int major = 0u,
          unsigned int minor = 0u,
          unsigned int patch = 0u);

  /// Constructor
  Version(int major, int minor = 0, int patch = 0);

  /// Constructor
  Version(const char *);
 
  /// Copy constructor
  Version(const Version &);
 
  /// Assignment operator
  Version &operator =(const Version &);
 
  /// Whether this version is valid
  operator bool() const;

  /// Equality operator
  bool operator ==(const Version &) const;
 
  /// Inequality operator
  bool operator !=(const Version &) const;

  /// Less operator
  bool operator <(const Version &) const;

  /// Greater operator
  bool operator >(const Version &) const;

  /// Less or equal operator
  bool operator <=(const Version &) const;

  /// Greater or equal operator
  bool operator >=(const Version &) const;

  /// Get version information as string
  string ToString() const;

};

////////////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////////////

/// Current software version
const Version current_version(MIRTK_VERSION_MAJOR,
                              MIRTK_VERSION_MINOR,
                              MIRTK_VERSION_PATCH);

/// Version to emulate
MIRTK_Common_EXPORT extern Version version;

/// Print software revision number (or version if not available) only
void PrintRevision(ostream &);

/// Print build time stamp as version string
void PrintVersion(ostream &, const char* = NULL);

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
/// Write version to output stream
inline ostream &operator <<(ostream &os, const Version &version)
{
  os << version.Major() << "." << version.Minor();
  if (version.Patch() != 0u) os << "." << version.Patch();
  return os;
}

// -----------------------------------------------------------------------------
/// Read version from input stream
inline istream &operator >>(istream &is, Version &version)
{
  string token;
  char c;
  while (is.get(c) && (isdigit(c) || c == '.')) token.push_back(c);
  if (is.bad()) return is;
  if (is.eof()) is.clear(ios::eofbit); // reset failbit
  istringstream ss(token);
  unsigned int number[3] = {0u, 0u, 0u};
  int          i         = -1;
  while (getline(ss, token, '.')) {
    if (++i == 3 || !FromString(token.c_str(), number[i])) {
      i = -1;
      break;
    }
  }
  if (i == -1) is.setstate(ios::failbit);
  else {
    version.Major(number[0]);
    version.Minor(number[1]);
    version.Patch(number[2]);
  }
  return is;
}

// -----------------------------------------------------------------------------
inline string Version::ToString() const
{
  ostringstream ss;
  ss << (*this);
  return ss.str();
}

// -----------------------------------------------------------------------------
inline Version::Version(unsigned int major,
                        unsigned int minor,
                        unsigned int patch)
:
  _Major(major), _Minor(minor), _Patch(patch)
{
}

// -----------------------------------------------------------------------------
inline Version::Version(int major, int minor, int patch)
:
  _Major(major >= 0 ? static_cast<unsigned int>(major) : 0u),
  _Minor(minor >= 0 ? static_cast<unsigned int>(minor) : 0u),
  _Patch(patch >= 0 ? static_cast<unsigned int>(patch) : 0u)
{
}

// -----------------------------------------------------------------------------
inline Version::Version(const char *str)
:
  _Major(0u), _Minor(0u), _Patch(0u)
{
  FromString(str, *this);
}

// -----------------------------------------------------------------------------
inline Version::Version(const Version &other)
:
  _Major(other._Major), _Minor(other._Minor), _Patch(other._Patch)
{
}

// -----------------------------------------------------------------------------
inline Version &Version::operator =(const Version &other)
{
  if (this != &other) {
    _Major = other._Major;
    _Minor = other._Minor;
    _Patch = other._Patch;
  }
  return *this;
}

// -----------------------------------------------------------------------------
inline Version::operator bool() const
{
  return (_Major != 0u || _Minor != 0u || _Patch != 0u);
}

// -----------------------------------------------------------------------------
inline bool Version::operator ==(const Version &rhs) const
{
  return _Major == rhs._Major && _Minor == rhs._Minor && _Patch == rhs._Patch;
}

// -----------------------------------------------------------------------------
inline bool Version::operator !=(const Version &rhs) const
{
  return !(*this == rhs);
}

// -----------------------------------------------------------------------------
// Note: Invalid/development branch version 0.0.0 is greater than any other
inline bool Version::operator <(const Version &rhs) const
{
  return bool(*this) && (!rhs || (_Major < rhs._Major || (_Major == rhs._Major && (_Minor < rhs._Minor || (_Minor == rhs._Minor && _Patch < rhs._Patch)))));
}

// -----------------------------------------------------------------------------
inline bool Version::operator <=(const Version &rhs) const
{
  return (*this == rhs) || (*this < rhs);
}

// -----------------------------------------------------------------------------
// Note: Invalid/development branch version 0.0.0 is greater than any other
inline bool Version::operator >(const Version &rhs) const
{
  return !(*this <= rhs);
}

// -----------------------------------------------------------------------------
inline bool Version::operator >=(const Version &rhs) const
{
  return (*this == rhs) || (*this > rhs);
}


} // namespace mirtk

#endif // MIRTK_Version_H

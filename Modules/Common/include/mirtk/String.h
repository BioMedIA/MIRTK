/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#ifndef MIRTK_String_H
#define MIRTK_String_H

#include "mirtk/CommonExport.h"

#include <cctype>
#include <cstring>
#include <string>
#include <sstream>
#include <iomanip>
#include <limits>

#include "mirtk/Array.h"


namespace mirtk {


// Basic C string functions
using ::isalnum;
using ::isdigit;
using ::islower;
using ::isupper;
using ::tolower;
using ::toupper;
using ::strstr;
using ::strcmp;
using ::strncmp;

// Import C++ library type into mirtk namespace
using std::string;


// =============================================================================
// String comparison
// =============================================================================

/// Case insensitive string comparison
inline bool iequal(char const *a, char const *b)
{
  while (*a && *b) {
    if (tolower(*a) != tolower(*b)) return false;
    ++a, ++b;
  }
  return !(*a || *b);
}

/// Case insensitive string comparison
inline bool iequal(const string &a, char const *b)
{
  return iequal(a.c_str(), b);
}

/// Case insensitive string comparison
inline bool iequal(const char *a, const string &b)
{
  return iequal(a, b.c_str());
}

/// Case insensitive string comparison
inline bool iequal(const string &a, const string &b)
{
  return iequal(a.c_str(), b.c_str());
}

// =============================================================================
// Custom string/stream functions
// =============================================================================

/// General routine to read float from a file stream
MIRTK_Common_DEPRECATED int ReadInt(std::ifstream &);

/// General routine to read float from a file stream
MIRTK_Common_DEPRECATED float ReadFloat(std::ifstream &);

/// General routine to read list of char (string) from a file stream
MIRTK_Common_DEPRECATED char *ReadString(std::ifstream &);

/// Convert string to numeric value
template <typename T>
bool FromString(const char *str, T &value)
{
  if (str == NULL || str[0] == '\0') return false;
  std::istringstream is(str);
  return !(is >> value).fail() && is.eof();
}

/// Convert string to numeric value
template <typename T>
bool FromString(const string &s, T &value)
{
  return FromString(s.c_str(), value);
}

/// Convert string to boolean value
template <>
inline bool FromString(const char *str, bool &value)
{
  if (strcmp(str, "yes") == 0 || strcmp(str, "Yes") == 0 || strcmp(str, "YES") == 0) {
    value = true;
    return true;
  } else if (strcmp(str, "no") == 0 || strcmp(str, "No") == 0 || strcmp(str, "NO") == 0) {
    value = false;
    return true;
  } else if (strcmp(str, "true") == 0 || strcmp(str, "True") == 0 || strcmp(str, "TRUE") == 0) {
    value = true;
    return true;
  } else if (strcmp(str, "false") == 0 || strcmp(str, "False") == 0 || strcmp(str, "FALSE") == 0) {
    value = false;
    return true;
  } else if (strcmp(str, "on") == 0 || strcmp(str, "On") == 0 || strcmp(str, "ON") == 0) {
    value = true;
    return true;
  } else if (strcmp(str, "off") == 0 || strcmp(str, "Off") == 0 || strcmp(str, "OFF") == 0) {
    value = false;
    return true;
  } else {
    std::istringstream is(str);
    return !(is >> value).fail() && is.eof();
  }
}

/// Convert string to float value
template <>
inline bool FromString(const char *str, float &value)
{
  if (strcmp(str, "nan") == 0 || strcmp(str, "NaN") == 0) {
    value = std::numeric_limits<float>::quiet_NaN();
    return true;
  } else if (strcmp(str, "-inf") == 0 || strcmp(str, "-Inf") == 0) {
    value = -std::numeric_limits<float>::infinity();
    return true;
  } else if (strcmp(str, "+inf") == 0 || strcmp(str, "inf") == 0 || strcmp(str, "+Inf") == 0 || strcmp(str, "Inf") == 0) {
    value = +std::numeric_limits<float>::infinity();
    return true;
  } else {
    std::istringstream is(str);
    return !(is >> value).fail() && is.eof();
  }
}

/// Convert string to double value
template <>
inline bool FromString(const char *str, double &value)
{
  if (strcmp(str, "nan") == 0 || strcmp(str, "NaN") == 0) {
    value = std::numeric_limits<double>::quiet_NaN();
    return true;
  } else if (strcmp(str, "-inf") == 0 || strcmp(str, "-Inf") == 0) {
    value = -std::numeric_limits<double>::infinity();
    return true;
  } else if (strcmp(str, "+inf") == 0 || strcmp(str, "inf") == 0 || strcmp(str, "+Inf") == 0 || strcmp(str, "Inf") == 0) {
    value = +std::numeric_limits<double>::infinity();
    return true;
  } else {
    std::istringstream is(str);
    return !(is >> value).fail() && is.eof();
  }
}

/// Check if given string represents a value of the specified template type
template <typename T> bool IsValueOfType(const char *str)
{
  T value;
  return FromString(str, value);
}

/// Check if given string is a (floating point) number
inline bool IsNumber(const char *str)
{
  return IsValueOfType<double>(str);
}

/// Check if given string is an integer value
inline bool IsInteger(const char *str)
{
  return IsValueOfType<int>(str);
}

/// Convert numeric value to string
template <typename T>
string ToString(const T &value, int w = 0, char c = ' ', bool left = false)
{
  std::ostringstream os;
  os.fill(c);
  if (left) os << std::left;
  else      os << std::right;
  os << std::setw(w) << value;
  return os.str();
}

/// Convert boolean value to string
template <>
inline string ToString(const bool &value, int w, char c, bool left)
{
  std::ostringstream os;
  os.fill(c);
  if (left) os << std::left;
  else      os << std::right;
  os << std::setw(w) << (value ? "Yes" : "No");
  return os.str();
}

/// Write "<name> = <value>" configuration entry to output stream
inline void PrintParameter(std::ostream &os, const char *name, const char *value)
{
  const std::streamsize w = os.width(40);
  os << std::left << name << std::setw(0) << " = " << value << "\n";
  os.width(w);
}

/// Write "<name> = <value>" configuration entry to output stream
inline void PrintParameter(std::ostream &os, const char *name, const string &value)
{
  PrintParameter(os, name, value.c_str());
}

/// Write "<name> = <value>" configuration entry to output stream
inline void PrintParameter(std::ostream &os, const string &name, const string &value)
{
  PrintParameter(os, name.c_str(), value.c_str());
}

/// Write "<name> = <value>" configuration entry to output stream
template <class TValue>
inline void PrintParameter(std::ostream &os, const char *name, const TValue &value)
{
  PrintParameter(os, name, ToString(value));
}

/// Write "<name> = <value>" configuration entry to output stream
template <class TValue>
inline void PrintParameter(std::ostream &os, const string &name, const TValue &value)
{
  PrintParameter(os, name, ToString(value));
}

/// Print single argument to output stream
template <typename T>
std::ostream &Print(std::ostream &os, T value)
{
  os << ToString(value);
  return os;
}

/// Print any number of arguments to output stream, i.e.,
/// concatenating the string representations of the arguments
template <typename T, typename... Args>
std::ostream &Print(std::ostream &os, T value, Args... args)
{
  os << ToString(value);
  return Print(os, args...);
}

/// Convert string to lowercase letters
string ToLower(const string &);

/// Convert string to uppercase letters
string ToUpper(const string &);

/// Trim leading and trailing (whitespace) characters from string
///
/// \param[in] str  Input string.
/// \param[in] what Set of characters to remove from start and end of string.
///
/// \returns String with leading and trailing (whitespace) characters removed.
string Trim(const string &str, const string &what = " \t\r\n");

/// Remove (whitespace) characters from string
///
/// \param[in] str  Input string.
/// \param[in] what Set of characters to remove from string.
///
/// \returns String with leading and trailing (whitespace) characters removed.
string TrimAll(const string &str, const string &what = " \t\r\n");

/// Split string into parts separated by specified delimiting sequence of characters
///
/// @param s String to be split.
/// @param d Delimiting sequence of characters.
/// @param n Maximum number of parts. If zero, all parts are returned,
///          if negative, the last n parts are returned, and if positive,
///          the first n parts are returned.
/// @param e Discard empty strings.
/// @param q Do not split quoted parts. Double quotes within quotes
///          have to be escaped with a preceding backslash character.
///
/// @returns Parts of the string.
Array<string> Split(string s, const char *d, int n = 0, bool e = false, bool q = false);

/// Split string into parts separated by specified delimiting character
///
/// @param s String to be split.
/// @param d Delimiting character.
/// @param n Maximum number of parts. If zero, all parts are returned,
///          if negative, the last n parts are returned, and if positive,
///          the first n parts are returned.
/// @param e Discard empty strings.
/// @param q Do not split quoted parts. Double quotes within quotes
///          have to be escaped with a preceding backslash character.
///
/// @returns Parts of the string.
Array<string> Split(string s, char d, int n = 0, bool e = false, bool q = false);

/// Convert (upper) camel case string to space separated string
///
/// \param[in] s Camel case string.
///
/// \return String starting with uppercase letter followed by lowercase letters
///         only and a space character before each uppercase letter in the
///         camel case string.
string CamelCaseToPrettyParameterName(const string &s);

/// Convert units specification to standard lowercase string
///
/// For example, this function returns "vox" for any units specification allowed
/// for voxel units including "voxel" and "VOxELS". Same for other units for
/// which alternative units specifications are allowed. Removes enclosing
/// brackets from units string if present, e.g., returns "mm" when the input
/// string is "[mm]".
///
/// \param[in] str Units specification string.
///
/// \returns Standard units specification in all lowercase.
string StandardUnits(const string &str);

/// Splits a parameter name string such as "Resolution [mm]" into prefix and units
///
/// \param[in]  str  Parameter name string with optional units specification.
/// \param[out] name Name of parameter without units specification.
/// \param[in]  dflt Default units if none specified.
///
/// \returns Standard units string or \p dflt string if units specification missing.
string ParameterUnits(const string &str, string *name = nullptr, const char *dflt = "");

/// Splits a parameter value string such as "1 mm", "1 [mm]", "1 2 3 mm", "foo [rel]"
///
/// \note When the parameter value is not numeric, the units specification must
///       be enclosed in square brackets ([]) and separated by at least one white
///       space character from the parameter value(s).
///
/// \param[in]  str   Numeric parameter value string with optional units specification.
/// \param[out] value Value(s) of parameter without units specification.
/// \param[in]  dflt  Default units if none specified.
///
/// \returns Standard units string or \p dflt string if units specification missing.
string ValueUnits(const string &str, string *value = nullptr, const char *dflt = "");


} // namespace mirtk

#endif // MIRTK_String_H

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

#include "mirtk/Version.h"
#include "mirtk/VersionInfo.h"


namespace mirtk {

// Default: "Emulate" latest version
MIRTK_Common_EXPORT Version version = current_version;

// -----------------------------------------------------------------------------
void PrintRevision(ostream &out)
{
  if (MIRTK_REVISION) out << MIRTK_REVISION;
  else                out << current_version << endl;
}

// -----------------------------------------------------------------------------
void PrintVersion(ostream &out, const char *name)
{
  if (name) out << name << " ";
  if (current_version) {
    out << MIRTK_VERSION_STRING;
    if (version < current_version) out << ", emulating version " << version;
  }
  out << " (";
  if (MIRTK_REVISION && MIRTK_REVISION[0] != '\0') {
    out << "rev " << MIRTK_REVISION << ", ";
  }
  out << "built on " << (__DATE__);
  out << ")";
  out << endl;
}


} // namespace mirtk

# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2016 Imperial College London
# Copyright 2013-2016 Andreas Schuh
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ============================================================================

# Build and install Python modules and scripts when Python is available
# without a WITH_Python user option. Just look for it and use it if available.
set(USE_Python TRUE)

# Initialize project
mirtk_project_begin()

# Process subdirectories
foreach (SUBDIR IN LISTS PROJECT_SUBDIRS)
  basis_add_subdirectory(${SUBDIR})
endforeach ()

# Dynamic linker search path for external dependencies
# added to this global property by mirtk_add_library calls
get_property(LDPATH GLOBAL PROPERTY MIRTK_LIBRARY_PATH)
if (LDPATH)
  list(REMOVE_DUPLICATES LDPATH)
endif ()
if (UNIX)
  list(REMOVE_ITEM LDPATH "/lib")
  list(REMOVE_ITEM LDPATH "/lib64")
  list(REMOVE_ITEM LDPATH "/usr/lib")
  list(REMOVE_ITEM LDPATH "/usr/lib64")
  list(REMOVE_ITEM LDPATH "/usr/local/lib")
  list(REMOVE_ITEM LDPATH "/opt/local/lib")
endif ()

# Set script configuration of Python modules
basis_set_target_properties(${BASIS_PYTHON_LIBRARY_TARGET} PROPERTIES
  SCRIPT_DEFINITIONS
    "basis_list_to_delimited_string(LDPATH \"', '\" NOAUTOQUOTE \"'${LDPATH}'\")
     basis_set_script_path(LIBEXEC_DIR \"${BINARY_LIBEXEC_DIR}/${MIRTK_TOOLS_SUBDIR}\" \"${MIRTK_TOOLS_DIR}\")"
)

# Finalize project
mirtk_project_end()

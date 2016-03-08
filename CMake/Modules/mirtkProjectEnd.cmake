# ==============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2016 Imperial College London
# Copyright 2016 Andreas Schuh
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
# ==============================================================================

if (COMMAND mirtk_project_end)
  return()
endif ()

# ------------------------------------------------------------------------------
## Finalize MIRTK project (see basis_project_end)
macro(mirtk_project_end)
  # Finalize BASIS project
  basis_project_end()
  # Hide unused BASIS installation directory settings
  if (NOT PROJECT_IS_MODULE)
    foreach (lang IN ITEMS BASH JYTHON PYTHON MATLAB PERL PYTHON)
      set_property(CACHE INSTALL_${lang}_SITE_DIR PROPERTY TYPE INTERNAL)
    endforeach ()
  endif ()
endmacro ()

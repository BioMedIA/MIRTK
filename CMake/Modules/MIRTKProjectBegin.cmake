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

if (COMMAND mirtk_project_begin)
  return()
endif ()

# ------------------------------------------------------------------------------
## Initialize MIRTK project (see basis_project_begin)
macro(mirtk_project_begin)
  # Override default BASIS settings and hide unused cache entries
  if (NOT PROJECT_IS_MODULE)
    set(BASIS_NAMESPACE_DELIMITER_CMAKE ::)
    set(BUILD_MODULES_BY_DEFAULT    OFF       CACHE INTERNAL "" FORCE)
    set(BASIS_BUILD_ONLY            OFF       CACHE INTERNAL "" FORCE)
    set(BASIS_ALL_DOC               OFF       CACHE INTERNAL "" FORCE)
    set(BASIS_COMPILE_MATLAB        ON        CACHE INTERNAL "" FORCE)
    set(BASIS_COMPILE_SCRIPTS       OFF       CACHE INTERNAL "" FORCE)
    set(BASIS_CREATE_EXPORTS_FILE   ON        CACHE INTERNAL "" FORCE)
    set(BASIS_INSTALL_RPATH         ON        CACHE INTERNAL "" FORCE)
    set(BASIS_INSTALL_SCHEME        "default" CACHE INTERNAL "" FORCE)
    set(BASIS_INSTALL_SITE_PACKAGES OFF       CACHE INTERNAL "" FORCE)
    set(BASIS_REGISTER              ON        CACHE INTERNAL "" FORCE)
    set(BASIS_REVISION_INFO         ON        CACHE INTERNAL "" FORCE)
    set(BASIS_SUPERBUILD_MODULES    OFF       CACHE INTERNAL "" FORCE)
    foreach (lang IN ITEMS BASH CXX PYTHON PERL)
      set(BUILD_BASIS_UTILITIES_FOR_${lang} OFF CACHE INTERNAL "" FORCE)
    endforeach ()
    unset(lang)
  endif ()
  # Start BASIS project
  basis_project_begin()
endmacro ()

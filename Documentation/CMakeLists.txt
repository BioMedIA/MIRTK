# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2015 Imperial College London
# Copyright 2013-2015 Andreas Schuh
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

##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of MIRTK Documentation.
#
# This build configuration builds and/or installs the documentation of the
# software package. The documentation can be divided into user documentation
# (e.g., software manual) and developer documentation (e.g., developer's guide).
# For developers, both those using a library provided by this package and the
# package maintainers, the documentation of the API is of particular interest.
# Developers who are responsible for maintaining this software may be provided
# with even more detailed documentation of the implementation in the so-called
# developer's guide.
#
# See the basis_add_doc() command for details on the documentation build tools
# supported by BASIS for the generation of documentation from plain text files
# and in-source code comments.
##############################################################################

# LaTeX required for Doxygen formulas
#basis_find_package(LATEX COMPONENTS PDFLATEX REQUIRED)

# ============================================================================
# API reference (Doxygen-style in-source code comments)
# ============================================================================

if (WIN32)
  set(WINDOWS_DEFINE WINDOWS)
else ()
  set(WINDOWS_DEFINE)
endif ()

basis_find_package(Doxygen REQUIRED)

basis_add_doc (
  apidoc                                                        # documentation build target name
  GENERATOR             Doxygen                                 # select Doxygen as the generator
  DOXYFILE              apidoc/doxygen_config.in                # template Doxygen configuration file
  HTML_DESTINATION      "${INSTALL_DOC_DIR}/html/apidoc/latest" # output location of html documentation
  HTML_HEADER	          "apidoc/doxygen_header.html.in"         # custom top of page with CMake substitutions to fit with Sphinx
  HTML_FOOTER           "apidoc/doxygen_footer.html.in"         # custom bottom of page with CMake substitutions to fit with Sphinx
  HTML_EXTRA_STYLESHEET "apidoc/doxygen_extra.css.in"           # custom formatting to fit in with the Sphinx docs
  OUTPUT                html xml                                # request output in HTML for websites and XML for reprocessing,
                                                                # such as with the breathe Sphinx extension in particular
  # compile definitions
  PREDEFINED
    ${WINDOWS_DEFINE}
  # exclude CMake BASIS Module groups and CMake files
  EXCLUDE_BASIS_MODULES
  EXCLUDE_BASIS_UTILITIES
  # exclude third-party code
  EXCLUDE_PATTERNS
    *.cc
    mirtkParallel.h
    mirtkWindows.h
    mirtkArpack.h
    mirtkUmfpack.h
    mirtkEigen.h
    mirtkVtk.h
    mirtkMath.h
    mirtkCuda.h
    mirtkCudaRuntime.h
    mirtkCutilMath.h
    mirtkForEach*VoxelFunction.h
    dfsurface.h
    nifti1.h nifti2*.* znzlib.*
    gifti_io.* gifti_xml.*
    */Eigen/*
    */Boost/*
    */LBFGS/*
)
 
# ============================================================================
# Write commands.rst and commands-overview.rst
# ============================================================================

if (BUILD_DOCUMENTATION_SOURCES)

  # Get list of MIRTK command executable targets added by mirtk_add_executable
  get_property(COMMANDS GLOBAL PROPERTY MIRTK_COMMANDS)

  # Remove commands which do not (yet) support the -h, -help options
  if (COMMANDS)
    list(REMOVE_ITEM COMMANDS
      convert-mris
      padding
      display
      view
    )
  endif ()

  if (COMMANDS)

    # Check if required auxiliary tools are available
    basis_get_target_uid(mirtk_target_uid    mirtk)
    basis_get_target_uid(help_rst_target_uid help-rst)
    if (TARGET ${mirtk_target_uid} AND TARGET ${help_rst_target_uid})
      basis_get_target_location(mirtk_path ${mirtk_target_uid} ABSOLUTE)

      # Phony target to generate all help pages
      basis_get_target_uid(commands_help_target_uid commands-help)
      add_custom_target(${commands_help_target_uid})

      # Generate help page and brief description section for each command
      foreach (target_name IN LISTS COMMANDS)
        basis_get_target_uid(target_uid ${target_name})
        add_custom_target(${target_uid}-description
          COMMAND "${mirtk_path}" help-rst ${target_name} -generated -description -noheaders
                  -output "${TOPLEVEL_PROJECT_DOC_DIR}/commands/_descriptions/${target_name}.rst"
                  -output-brief-description "${TOPLEVEL_PROJECT_DOC_DIR}/commands/_summaries/${target_name}.rst"
          DEPENDS ${mirtk_target_uid} ${help_rst_target_uid} ${target_uid}
          COMMENT "Extracting command description from mirtk ${target_name} -help"
        )
        add_custom_target(${target_uid}-help
          COMMAND "${mirtk_path}" help-rst ${target_name} -generated -orphan
                  "-include-description" "_descriptions/${target_name}.rst"
                  -output "${TOPLEVEL_PROJECT_DOC_DIR}/commands/${target_name}.rst"
          DEPENDS ${mirtk_target_uid} ${help_rst_target_uid} ${target_uid} ${target_uid}-description
          COMMENT "Generating documentation page for command ${target_name}"
        )
        add_dependencies(${commands_help_target_uid} ${target_uid}-help)
      endforeach ()

      # Write overview pages: commands.rst, commands-overview.rst
      if (COMMANDS)
        list(SORT COMMANDS)
        message(STATUS "Writing commands/_overview.rst...")
        file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/commands")
        execute_process(
          COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tools/write-commands-overview.py" ${COMMANDS}
          OUTPUT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/commands/_overview.rst"
        )
        message(STATUS "Writing commands/_overview.rst... done")
        message(STATUS "Writing commands/_summary.rst...")
        execute_process(
          COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tools/write-commands-summary.py" ${COMMANDS}
          OUTPUT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/commands/_summary.rst"
        )
        message(STATUS "Writing commands/_summary.rst... done")
      endif ()

    else ()
      message(WARNING "Target of mirtk or help-rst tool missing! Skipping auto-generation of command documentation.")
    endif ()
  endif ()
endif ()

# ============================================================================
# Online documentation
# ============================================================================

basis_find_package(Sphinx COMPONENTS build REQUIRED)

basis_add_doc (
  pages                                                 # documenation build target name
  GENERATOR        Sphinx                               # select Sphinx as the generator
  BUILDER          html                                 # selected Sphinx output formats
  HTML_THEME       bizstyle                             # select the Sphinx layout theme
  HTML_SIDEBARS    searchbox globaltoc                  # sidebar options to use
  HTML_DESTINATION "${INSTALL_DOC_DIR}/html"            # output directory for completed documentation
  SIDEBARWIDTH     300                                  # sidebar width in pixels
  MASTER_DOC       "sidebar"                            # .rst file to start with when generating HTML
  OUTPUT_NAME      "${PROJECT_NAME}_Software_Manual"    # general output file name, e.g., LaTeX/PDF files
  DOXYLINK         apidoc                               # Doxygen generator build target for integrated API reference
  DOXYLINK_URL     "${PROJECT_WEBSITE}/apidoc/latest"   # location of Doxygen output files
  EXTENSIONS       mathjax                              # additional Sphinx extensions
  NO_HTML_MODINDEX NO_HTML_INDEX                        # disable currently unused index page
  EXCLUDE_PATTERN  **/_descriptions
  EXCLUDE_PATTERN  **/_summaries
  EXCLUDE_PATTERN  commands/_overview.rst
  EXCLUDE_PATTERN  commands/_usage.rst
  EXCLUDE_PATTERN  commands/_options.rst
  EXCLUDE_PATTERN  commands/_summary.rst
  EXCLUDE_PATTERN  modules/_overview.rst
  EXCLUDE_PATTERN  modules/_refs.rst
  EXCLUDE_PATTERN  contents.rst
)

# (re-)generate commands/<command_name>.rst files from command -help output
if (BUILD_DOCUMENTATION_SOURCES AND TARGET commands-help)
  add_dependencies(pages_html commands-help)
endif ()

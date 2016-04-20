# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2016 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

##############################################################################
# @file  CMakeLists.txt
# @brief Common build configuration file to add uninstaller.
#
# CMake does not support the creation of an uninstaller. Therefore, we add
# code to the cmake_install.cmake file which writes the content of the
# install_manifest.txt file to a file in the installation tree. As this
# code has to be executed after all files were executed and CMake includes
# the cmake_install.cmake files of subdirectories at last, this CMakeLists.txt
# is placed into a shared directory of the BASIS installation which is added
# to every BASIS project using add_subdirctory() at the very end of the root
# CMakeLists.txt file.
##############################################################################

# install CMake-based uninstall script
if (NOT PROJECT_IS_SUBMODULE)
  install (
    FILES       "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake"
    DESTINATION "${INSTALL_CONFIG_DIR}"
    RENAME      "${PROJECT_PACKAGE_CONFIG_PREFIX}Uninstall.cmake"
  )
endif ()

# write executable uninstaller
if (NOT PROJECT_IS_SUBMODULE)
  set (UNINSTALLER_NAME "uninstall-${PROJECT_PACKAGE_NAME_L}")
  if (WIN32)
    set (UNINSTALLER_NAME "${UNINSTALLER_NAME}.cmd")
  endif ()
  
  install (
    CODE
      "
      set (UNINSTALLER \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${INSTALL_RUNTIME_DIR}/${UNINSTALLER_NAME}\")
      message (\"-- Installing: \${UNINSTALLER}\")
      file (WRITE \"\${UNINSTALLER}\" \"\")
      if (UNIX)
        file (APPEND \"\${UNINSTALLER}\" \"#! /bin/sh\nexec \")
      endif ()
      file (APPEND \"\${UNINSTALLER}\" \"\\\"${CMAKE_COMMAND}\\\" -P \\\"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${INSTALL_CONFIG_DIR}/${PROJECT_PACKAGE_CONFIG_PREFIX}Uninstall.cmake\\\"\")
      if (UNIX)
        execute_process (COMMAND /bin/chmod +x \"\${UNINSTALLER}\")
      endif ()
      list (APPEND CMAKE_INSTALL_MANIFEST_FILES \"\${UNINSTALLER}\")
      "
  )
endif ()

# write install manifest at the very end of the installation
install (
  CODE
    "
    set (INSTALL_MANIFEST_FILE \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${INSTALL_CONFIG_DIR}/${PROJECT_PACKAGE_CONFIG_PREFIX}InstallManifest.txt\")
    list (APPEND CMAKE_INSTALL_MANIFEST_FILES \"\${INSTALL_MANIFEST_FILE}\")
    message (\"-- Installing: \${INSTALL_MANIFEST_FILE}\")
    file (WRITE \"\${INSTALL_MANIFEST_FILE}\" \"\")
    foreach (F \${CMAKE_INSTALL_MANIFEST_FILES})
      file (APPEND \"\${INSTALL_MANIFEST_FILE}\" \"\${F}\n\")
    endforeach ()
    "
)

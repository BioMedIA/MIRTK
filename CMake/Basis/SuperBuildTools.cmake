# ============================================================================
# Copyright (c) 2014 Carnegie Mellon University
# All rights reserved.
#
# See COPYING file for license information or visit
# https://cmake-basis.github.io/download.html#license
# ============================================================================

if (__BASIS_SUPER_BUILD_INCLUDED)
  return ()
else ()
  set (__BASIS_SUPER_BUILD_INCLUDED TRUE)
endif ()

include(ExternalProject)
         
##
# @brief When enabled CMake will always reconfigure super build modules. Slows performance but won't ignore changes in external projects.
#
# @note The global variable BASIS_SUPER_BUILD_ARGS is passed to the CMAKE_ARGS 
#       parameter of ExternalProject_Add in case custom variables need to be supplied.
#
option(BASIS_ALWAYS_RECONFIGURE_SUPER_BUILD "Enable to always reconfigure super build modules. Slows performance but won't ignore changes." OFF)
mark_as_advanced(BASIS_ALWAYS_RECONFIGURE_SUPER_BUILD)

##
# @brief super build for BASIS modules
#
function(basis_super_build PACKAGE_NAME)
  set(options )
  set(singleValueArgs DIR CMAKE_MODULE_PATH BINARY_DIR CMAKE_INSTALL_PREFIX)
  set(multiValueArgs  DEPENDS)
  
  cmake_parse_arguments(${PACKAGE_NAME} ${options} ${singleValueArgs} ${multiValueArgs} ${ARGN})

  # TODO: consider combining this variable with MODULE_${PACKAGE_NAME} variable
  #option (USE_SYSTEM_${PACKAGE_NAME} "Skip build of ${PACKAGE_NAME} if already installed." OFF)
  
  if(NOT ${PACKAGE_NAME}_CMAKE_MODULE_PATH)
    set(${PACKAGE_NAME}_CMAKE_MODULE_PATH "-DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")
  endif()
  
    # TODO: make sure default install prefix does not typically trample the installation
  if(NOT ${PACKAGE_NAME}_CMAKE_INSTALL_PREFIX)
    set(${PACKAGE_NAME}_CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
  endif()
  
  # set directory where binaries will build if it was not already set by the arguments
  if(NOT ${PACKAGE_NAME}_BINARY_DIR)
    if(MODULE_${PACKAGE_NAME}_BINARY_DIR)
      set(${PACKAGE_NAME}_BINARY_DIR ${MODULE_${PACKAGE_NAME}_BINARY_DIR})
    elseif(NOT MODULE_${PACKAGE_NAME}_BINARY_DIR)
      set(MODULE_${PACKAGE_NAME}_BINARY_DIR ${PROJECT_BINARY_DIR})
    endif()
  endif()
  
  if(NOT ${PACKAGE_NAME}_DIR AND MODULE_${MODULE}_SOURCE_DIR)
    set(${PACKAGE_NAME}_DIR "${MODULE_${MODULE}_SOURCE_DIR}")
  endif()

  # TODO: Fix DEPENDS parameter. May need to separate basis module and regular dependencies so they can specified separately for the super build.
  # TODO: Consider using the EP_BASE path with SET(ep_base ${${PACKAGE_NAME}_BINARY_DIR}) instead. (ep stands for external project) 
  # TODO: Figure out why a few intermediate files are still being put in the ${CMAKE_BINARY_DIR}/${PACKAGE_NAME}-prefix/ directory
  # TODO: Check for additional useful -D parameters.

  # passing semicolons has odd side effects because they may be automatically
  # dereferenced, so substitute another character, in this case pipe |
  string(REPLACE ";" "|" CMAKE_PREFIX_PATH_PIPE "${CMAKE_PREFIX_PATH}")
  
  # only specifiy dependencies that are actual targets
  # otherwise there would be an error
  set(SUPER_BUILD_TARGET_DEPENDENCIES)
  foreach(DEPENDENCY IN ${DEPENDS})
    if(TARGET DEPENDENCY)
      list(APPEND SUPER_BUILD_TARGET_DEPENDENCIES ${DEPENDENCY})
    endif()
  endforeach()
    
  string(REPLACE ";" " " SUPER_BUILD_TARGET_DEPENDENCIES "${SUPER_BUILD_TARGET_DEPENDENCIES}")
  
  if(BASIS_DEBUG)
      message(STATUS 
    "basis_super_build() Module:
      ExternalProject_Add(${PACKAGE_NAME}
                          DEPENDS ${SUPER_BUILD_TARGET_DEPENDENCIES}
                          SOURCE_DIR ${${PACKAGE_NAME}_DIR}
                          CMAKE_ARGS 
                            -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> 
                            -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} 
                            -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS} 
                            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} 
                            ${${PACKAGE_NAME}_CMAKE_MODULE_PATH}
                            ${BASIS_SUPER_BUILD_ARGS}
                          CMAKE_CHACHE_ARGS
                            -DCMAKE_PREFIX_PATH:STRING=${CMAKE_PREFIX_PATH_PIPE}
                          CMAKE_GENERATOR
                            ${CMAKE_GENERATOR}
                          CMAKE_TOOLSET
                            ${CMAKE_TOOLSET}
                          BINARY_DIR
                            ${${PACKAGE_NAME}_BINARY_DIR}
                          INSTALL_DIR
                            ${${PACKAGE_NAME}_CMAKE_INSTALL_PREFIX}
                          )
    "  )
  endif()
  
    #if(USE_SYSTEM_${PACKAGE_NAME})
    #  find_package(${PACKAGE_NAME})
    #elseif()
    
      ExternalProject_Add(${PACKAGE_NAME}
                          DEPENDS ${SUPER_BUILD_TARGET_DEPENDENCIES}
                          SOURCE_DIR ${${PACKAGE_NAME}_DIR}
                          LIST_SEPARATOR "|"
                          CMAKE_ARGS 
                            -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> 
                            -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} 
                            -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS} 
                            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} 
                            ${${PACKAGE_NAME}_CMAKE_MODULE_PATH}
                            ${BASIS_SUPER_BUILD_ARGS}
                          CMAKE_CHACHE_ARGS
                            -DCMAKE_PREFIX_PATH:STRING=${CMAKE_PREFIX_PATH_PIPE}
                          CMAKE_GENERATOR
                            ${CMAKE_GENERATOR}
                          CMAKE_TOOLSET
                            ${CMAKE_TOOLSET}
                          BINARY_DIR
                            ${${PACKAGE_NAME}_BINARY_DIR}
                          INSTALL_DIR
                            ${${PACKAGE_NAME}_CMAKE_INSTALL_PREFIX}
                          )
                        
  
      if(BASIS_ALWAYS_RECONFIGURE_SUPER_BUILD)
        ExternalProject_Add_Step(${PACKAGE_NAME} reconfigure
          COMMAND ${CMAKE_COMMAND} -E echo "Force configure of ${PACKAGE_NAME}"
          DEPENDEES update
          DEPENDERS configure
          ALWAYS 1)
      endif()
        
  #endif()
endfunction()

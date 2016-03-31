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
# @brief Build configuration of (core) MIRTK Applications.
##############################################################################

# ----------------------------------------------------------------------------
# Include optional dependencies (cf. Config/Depends.cmake)
if (MODULE_Image)
  if (VTK_FOUND)
    basis_use_package(VTK)
  endif ()
endif ()
if (MODULE_PointSet)
  if (VTK_FOUND)
    basis_use_package(VTK)
  endif ()
  if (FLANN_FOUND)
    basis_use_package(FLANN)
  else ()
    unset(FLANN_CPP_LIBRARY)
  endif ()
  if (MATLAB_FOUND)
    basis_use_package(MATLAB)
  else ()
    set(MATLAB_mwmclmcrrt_LIBRARY)
  endif ()
endif ()

# Generic file info tool
mirtk_add_executable(info
  DEPENDS  LibCommon LibNumerics LibImage LibIO
  OPTIONAL LibTransformation LibPointSet ${VTK_LIBRARIES}
)

# ----------------------------------------------------------------------------
# Image tools
if (MODULE_Image)

  macro(add_image_command cmd)
    mirtk_add_executable(${cmd} DEPENDS LibCommon LibNumerics LibImage LibIO ${ARGN})
  endmacro()

  add_image_command(average-images OPTIONAL LibTransformation)
  add_image_command(calculate OPTIONAL ${VTK_LIBRARIES})
  add_image_command(calculate-distance-map)
  add_image_command(calculate-exponential-map)
  add_image_command(calculate-lie-bracket)
  add_image_command(calculate-logarithmic-map)
  add_image_command(close-image)
  add_image_command(convert-image)
  add_image_command(cut-brain)
  add_image_command(detect-edges)
  add_image_command(dilate-image)
  add_image_command(downsample-image)
  add_image_command(edit-image OPTIONAL LibTransformation)
  add_image_command(erode-image)
  add_image_command(evaluate-overlap)
  add_image_command(extract-connected-components)
  add_image_command(extract-image-region)
  add_image_command(flip-image)
  add_image_command(concatenate-images)
  add_image_command(open-image)
  add_image_command(reflect-image)
  add_image_command(resample-image)
  add_image_command(smooth-image)
  add_image_command(transform-image LibTransformation)

endif ()

# ----------------------------------------------------------------------------
# Registration tools
if (MODULE_Registration)

  macro(add_registration_command cmd)
    mirtk_add_executable(${cmd} DEPENDS LibCommon LibNumerics LibImage LibTransformation LibRegistration ${ARGN})
  endmacro()

  add_registration_command(register LibIO OPTIONAL LibLBFGS LibPointSet LibDeformable)

endif ()

# ----------------------------------------------------------------------------
# Transformation tools
if (MODULE_Transformation)

  macro(add_transformation_command cmd)
    mirtk_add_executable(${cmd} DEPENDS LibCommon LibNumerics LibImage LibTransformation ${ARGN})
  endmacro()

  add_transformation_command(edit-dof)
  add_transformation_command(invert-dof)

  if (MODULE_Image)

    macro(add_image_transformation_command cmd)
      add_transformation_command(${cmd} LibIO ${ARGN})
    endmacro()

    add_image_transformation_command(average-dofs)
    add_image_transformation_command(bisect-dof)
    add_image_transformation_command(concatenate-dofs)
    add_image_transformation_command(convert-dof)
    add_image_transformation_command(convert-dof2csv)
    add_image_transformation_command(convert-dof2velo)
    add_image_transformation_command(evaluate-dof)
    add_image_transformation_command(evaluate-jacobian)
    add_image_transformation_command(init-dof OPTIONAL LibPointSet)

  endif ()
endif ()

# ----------------------------------------------------------------------------
# Point set tools
if (MODULE_PointSet)

  macro(add_pointset_command cmd)
    mirtk_add_executable(${cmd} DEPENDS LibCommon LibNumerics LibImage LibPointSet ${VTK_LIBRARIES} ${ARGN})
  endmacro()

  add_pointset_command(calculate-surface-attributes)
  add_pointset_command(calculate-surface-spectrum)
  add_pointset_command(convert-mris) # requires Python VTK modules and FreeSurfer installation
  add_pointset_command(copy-pointset-attributes)
  add_pointset_command(decimate-surface)
  add_pointset_command(delete-pointset-attributes)
  add_pointset_command(evaluate-distance)
  add_pointset_command(evaluate-distortion)
  add_pointset_command(evaluate-surface-overlap ${FLANN_CPP_LIBRARY})
  add_pointset_command(extract-connected-points)
  add_pointset_command(offset-surface)
  add_pointset_command(project-onto-surface)
  add_pointset_command(remesh)
  add_pointset_command(smooth-surface)

  if (MODULE_Image)

    macro(add_image_pointset_command cmd)
      add_pointset_command(${cmd} LibIO ${ARGN})
    endmacro()

    add_image_pointset_command(convert-pointset)
    add_image_pointset_command(extract-pointset-surface)
    add_image_pointset_command(extract-surface)

    if (MODULE_Transformation)
      add_image_pointset_command(transform-points LibTransformation)
    endif ()

  endif ()

  if (MODULE_Transformation)
    add_pointset_command(match-points LibTransformation)
  endif ()

  if (MATLAB_FOUND)
    add_pointset_command(convert-pointset-to-mat ${MATLAB_mwmclmcrrt_LIBRARY})
    add_pointset_command(copy-pointset-attributes-from-mat ${MATLAB_mwmclmcrrt_LIBRARY})
  endif ()

endif ()

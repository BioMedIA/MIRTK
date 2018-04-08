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

##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of basic MIRTK Commands.
##############################################################################

# ----------------------------------------------------------------------------
# Help re-formatting tools
if (Python_FOUND)
  mirtk_add_executable(help-rst)
endif ()

# ----------------------------------------------------------------------------
# Generic file info tool
if (MIRTK_Image_FOUND)
  mirtk_add_executable(info)
  mirtk_target_dependencies(info LibCommon LibNumerics LibImage LibIO)
  if (MIRTK_Transformation_FOUND)
    mirtk_target_dependencies(info LibTransformation)
  endif ()
  if (MIRTK_PointSet_FOUND)
    mirtk_target_dependencies(info LibPointSet vtkCommonCore vtkCommonDataModel vtkFiltersCore vtkFiltersGeneral)
  endif ()
endif ()

# ----------------------------------------------------------------------------
# Image tools
if (MIRTK_Image_FOUND)

  macro(add_image_command cmd)
    mirtk_add_executable(${cmd})
    mirtk_target_dependencies(${cmd} LibCommon LibNumerics LibImage LibIO ${ARGN})
  endmacro()

  add_image_command(aggregate-images)
  add_image_command(average-images)
  add_image_command(calculate-element-wise)
  add_image_command(calculate-distance-map)
  add_image_command(calculate-exponential-map)
  add_image_command(calculate-lie-bracket)
  add_image_command(calculate-logarithmic-map)
  add_image_command(close-image)
  add_image_command(combine-images)
  add_image_command(convert-image)
  add_image_command(cut-brain)
  add_image_command(detect-edges)
  add_image_command(dilate-image)
  add_image_command(downsample-image)
  add_image_command(edit-image)
  add_image_command(erode-image)
  add_image_command(evaluate-overlap)
  add_image_command(extract-connected-components)
  add_image_command(extract-image-region)
  add_image_command(extract-image-slice)
  add_image_command(extract-image-volume)
  add_image_command(flip-image)
  add_image_command(match-histogram)
  add_image_command(open-image)
  add_image_command(reflect-image)
  add_image_command(resample-image)
  add_image_command(smooth-image)
  add_image_command(subdivide-brain-image)

  if (MIRTK_Transformation_FOUND)
    mirtk_target_dependencies(average-images LibTransformation)
    mirtk_target_dependencies(edit-image     LibTransformation)
  endif ()

  if (MIRTK_Image_WITH_VTK)
    mirtk_target_dependencies(calculate-element-wise vtkCommonCore)
  endif ()

  if (MIRTK_Transformation_FOUND AND MIRTK_Registration_FOUND)
    add_image_command(evaluate-similarity LibTransformation LibRegistration)
  endif ()

  if (Python_FOUND)
    mirtk_add_executable(average-measure)
    mirtk_add_executable(average-overlap)
  endif ()

endif ()

# ----------------------------------------------------------------------------
# Registration tools
if (MIRTK_Registration_FOUND)

  macro(add_registration_command cmd)
    mirtk_add_executable(${cmd})
    mirtk_target_dependencies(${cmd} LibCommon LibNumerics LibImage LibIO LibTransformation LibRegistration ${ARGN})
  endmacro()

  add_registration_command(register)

  if (MIRTK_Registration_WITH_PointSet)
    mirtk_target_dependencies(register LibPointSet)
  endif ()
  if (MIRTK_Registration_WITH_Deformable)
    mirtk_target_dependencies(register LibDeformable)
  endif ()
  if (TARGET TBB::tbb)
    mirtk_target_dependencies(register TBB::tbb)
  endif ()

  if (Python_FOUND)
    mirtk_add_executable(construct-atlas)
    mirtk_add_executable(evaluate-atlas)
  endif ()

endif ()

# ----------------------------------------------------------------------------
# Transformation tools
if (MIRTK_Transformation_FOUND)

  macro(add_transformation_command cmd)
    mirtk_add_executable(${cmd})
    mirtk_target_dependencies(${cmd} LibCommon LibNumerics LibImage LibTransformation ${ARGN})
  endmacro()

  add_transformation_command(edit-dof)
  add_transformation_command(invert-dof)

  if (MIRTK_Image_FOUND)

    macro(add_image_transformation_command cmd)
      add_transformation_command(${cmd} LibIO ${ARGN})
    endmacro()

    add_image_transformation_command(average-dofs)
    add_image_transformation_command(bisect-dof)
    add_image_transformation_command(compose-dofs)
    add_image_transformation_command(convert-dof)
    add_image_transformation_command(evaluate-dof)
    add_image_transformation_command(evaluate-jacobian)
    add_image_transformation_command(init-dof)
    add_image_transformation_command(transform-image)

    if (MIRTK_IO_WITH_VTK AND TARGET vtkCommonCore)
      mirtk_target_dependencies(init-dof vtkCommonCore)
    endif ()
    if (TARGET TBB::tbb)
      mirtk_target_dependencies(evaluate-jacobian TBB::tbb)
    endif ()

  endif ()
endif ()

# ----------------------------------------------------------------------------
# Point set tools
if (MIRTK_PointSet_FOUND)

  macro(add_pointset_command cmd)
    mirtk_add_executable(${cmd})
    mirtk_target_dependencies(${cmd} LibCommon LibNumerics LibImage LibPointSet LibIO ${VTK_LIBRARIES} ${ARGN})
  endmacro()

  add_pointset_command(calculate-surface-attributes)
  add_pointset_command(calculate-surface-spectrum)
  add_pointset_command(close-scalars)
  add_pointset_command(copy-pointset-attributes)
  add_pointset_command(decimate-surface)
  add_pointset_command(delete-pointset-attributes)
  add_pointset_command(dilate-scalars)
  add_pointset_command(erode-scalars)
  add_pointset_command(evaluate-cardiac-motion)
  add_pointset_command(evaluate-distance)
  add_pointset_command(evaluate-distortion)
  add_pointset_command(evaluate-surface-mesh)
  add_pointset_command(evaluate-surface-overlap)
  add_pointset_command(extract-connected-points)
  add_pointset_command(offset-surface)
  add_pointset_command(open-scalars)
  add_pointset_command(project-onto-surface)
  add_pointset_command(remesh-surface)
  add_pointset_command(smooth-surface)
  add_pointset_command(blend-surface)

  if (FLANN_FOUND)
    mirtk_target_dependencies(evaluate-surface-overlap ${FLANN_CPP_LIBRARY})
  endif ()

  if (Python_FOUND)
    # requires Python VTK modules and FreeSurfer installation
    mirtk_add_executable(convert-mris)
  endif ()

  if (MIRTK_Image_FOUND)

    add_pointset_command(convert-pointset)
    add_pointset_command(extract-pointset-cells)
    add_pointset_command(extract-pointset-surface)
    add_pointset_command(extract-surface)
    add_pointset_command(merge-surfaces)

    if (MIRTK_Transformation_FOUND)
      add_pointset_command(transform-points LibTransformation)
    endif ()

  endif ()

  if (MIRTK_Transformation_FOUND)
    add_pointset_command(match-points LibTransformation)
  endif ()

  if (MATLAB_FOUND)
    add_pointset_command(convert-pointset-to-mat ${MATLAB_mwmclmcrrt_LIBRARY})
    add_pointset_command(copy-pointset-attributes-from-mat ${MATLAB_mwmclmcrrt_LIBRARY})
  endif ()

endif ()

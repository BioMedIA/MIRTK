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

macro (add_image_test class_name)
  mirtk_add_test(${class_name} DEPENDS LibImage)
endmacro ()


# Parallel voxel functions
add_image_test(ConvolutionFunction) # TODO: Requires arguments
add_image_test(UnaryVoxelFunction)

# Image interpolation/extrapolation
add_image_test(InterpolateExtrapolateImageFunction)

# Core image filters
add_image_test(Downsampling) # TODO: Requires arguments

# Exponential/Logartihmic map of vector field
#add_image_test(DisplacementToVelocityField)

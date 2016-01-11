# ==============================================================================
# Medical Image Registration ToolKit (MIRTK) LBFGS Library
#
# Copyright 2013-2015 Imperial College London
# Copyright 2013-2015 Andreas Schuh
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ==============================================================================

################################################################################
# @file  BasisProject.cmake
# @brief Sets basic information about the MIRTK module and calls basis_project().
#
# This file defines basic information about a project by calling 
# the basis_project() function. This basic information, also known as metadata, 
# is used by CMake BASIS to setup the project. The dependencies to other modules
# have to be specified here such that the top-level IRTK project can analyze the
# inter-module dependencies, as well as dependencies on third-party libraries.
#
# @sa http://opensource.andreasschuh.com/cmake-basis/standard/modules.html
#
# @ingroup BasisSettings
################################################################################

# Note: The #<*> dependency patterns are required by the basisproject tool and
#       should be kept on a separate line as last commented argument of the
#       corresponding options of the basis_project() command. The TEMPLATE
#       option and set argument are also required by this tool and should not
#       be changed manually. The argument is updated by basisproject --update.

basis_project (

  # ----------------------------------------------------------------------------
  # meta-data
  NAME        "LBFGS"
  PACKAGE     "MIRTK"
  AUTHORS     "Andreas Schuh, Naoaki Okazaki"
  DESCRIPTION "LBFGS library of the Medical Image Registration ToolKit."
  COPYRIGHT   "2013-2015 Imperial College London, Andreas Schuh"
  LICENSE     "MIT License"
  CONTACT     "Andreas Schuh <andreas.schuh.84@gmail.com>"
  TEMPLATE    "mirtk-module/1.0"

  # ----------------------------------------------------------------------------
  # dependencies
  DEPENDS
    MIRTK{Common,Numerics}
    #<dependency>
  OPTIONAL_DEPENDS
    #<optional-dependency>
  TEST_DEPENDS
    #<test-dependency>
  OPTIONAL_TEST_DEPENDS
    #<optional-test-dependency>

)

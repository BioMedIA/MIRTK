# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2014 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# http://opensource.andreasschuh.com/cmake-basis/download.html#license
# ============================================================================

##############################################################################
# @file  FindSparseBayes.cmake
# @brief Find SparseBayes package from Vector Anomaly Limited.
#
# @sa http://www.vectoranomaly.com/downloads/downloads.htm
##############################################################################

include (FindPackageHandleStandardArgs)

find_path (
  SparseBayes_DIR SparseBayes.m
  DOC "The directory containing SparseBayes.m file of the SparseBayes package."
)

set (SparseBayes_INCLUDE_DIR "${SparseBayes_DIR}")

find_package_handle_standard_args (
  SparseBayes
  REQUIRED_ARGS
    SparseBayes_INCLUDE_DIR
)

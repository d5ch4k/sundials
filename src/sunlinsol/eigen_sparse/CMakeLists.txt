# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU, Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the eigen_sparse SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_EIGENSPARSE\n\")")

# Add the library
sundials_add_library(sundials_sunlinsol_eigen_sparse
  SOURCES
        sunlinsol_eigen_sparse.cpp
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunlinsol/sunlinsol_eigen_sparse.h
  INCLUDE_SUBDIR
    sunlinsol
  INCLUDE_DIRECTORIES
    PUBLIC
    ${EIGEN_INCLUDE_DIR}
  OBJECT_LIBRARIES
    sundials_generic_obj
  LINK_LIBRARIES
    PUBLIC
    sundials_sunmatrixsparse
  OUTPUT_NAME
        sunlinsol_eigen_sparse
  VERSION
    ${sunlinsollib_VERSION}
  SOVERSION
    ${sunlinsollib_VERSION}
)

message(STATUS "Added SUNLINSOL_EIGENSPARSE module")

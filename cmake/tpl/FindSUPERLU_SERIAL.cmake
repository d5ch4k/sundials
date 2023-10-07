# ---------------------------------------------------------------
# Programmer(s): Eddy Banks and David J. Gardner @ LLNL
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
# SuperLU_serial find module that creates an imported target for
# SuperLU_SERIAL. The target is SUNDIALS::SUPERLU_SERIAL.
#
# The variable SUPERLU_SERIAL_LIBRARY_DIR can be used to control
# where the module looks for the library.
#
# The variable SUPERLU_SERIAL_INCLUDE_DIR can be used to set the
# include path for the library.
#
# Additional libraries can be passed in SUPERLU_SERIAL_LIBRARIES.
#
# This module also defines variables, but it is best to use
# the defined target to ensure includes and compile/link
# options are correctly passed to consumers.
#
#   SUPERLU_SERIAL_FOUND       - system has SuperLU_serial library
#   SUPERLU_SERIAL_LIBRARY     - the SuperLU_serial library
# ---------------------------------------------------------------

# Set SuperLU_serial library name
set(SUPERLU_SERIAL_LIBRARY_NAME superlu)

if(MSVC)
  set(CMAKE_FIND_LIBRARY_PREFIXES lib ${CMAKE_FIND_LIBRARY_PREFIXES})
endif()

# Check if SUPERLU_SERIAL_LIBRARIES contains the superlu
# library as well as TPLs. If so, extract it into the
# SUPERLU_SERIAL_LIBRARY variable.
if(SUPERLU_SERIAL_LIBRARIES MATCHES "${SUPERLU_SERIAL_LIBRARY_NAME}")
  foreach(lib ${SUPERLU_SERIAL_LIBRARIES})
    if(lib MATCHES "${SUPERLU_SERIAL_LIBRARY_NAME}")
      set(SUPERLU_SERIAL_LIBRARY ${lib})
    endif()
  endforeach()
endif()

# find library
if(NOT SUPERLU_SERIAL_LIBRARY)
  # search user provided directory path
  find_library(SUPERLU_SERIAL_LIBRARY ${SUPERLU_SERIAL_LIBRARY_NAME}
          PATHS ${SUPERLU_SERIAL_LIBRARY_DIR} NO_DEFAULT_PATH)
  # if user didn't provide a path, search anywhere
  if(NOT (SUPERLU_SERIAL_LIBRARY_DIR OR SUPERLU_SERIAL_LIBRARY))
    find_library(SUPERLU_SERIAL_LIBRARY ${SUPERLU_SERIAL_LIBRARY_NAME})
  endif()
  mark_as_advanced(SUPERLU_SERIAL_LIBRARY)
endif()

# set the libraries, stripping out 'NOTFOUND' from previous attempts
string(REPLACE "SUPERLU_SERIAL_LIBRARY-NOTFOUND" "" SUPERLU_SERIAL_LIBRARIES "${SUPERLU_SERIAL_LIBRARIES}")
set(SUPERLU_SERIAL_LIBRARIES "${SUPERLU_SERIAL_LIBRARY};${SUPERLU_SERIAL_LIBRARIES}" CACHE STRING "" FORCE)

# set the library dir option if it wasn't preset
if(SUPERLU_SERIAL_LIBRARY AND (NOT SUPERLU_SERIAL_LIBRARY_DIR))
  get_filename_component(SUPERLU_SERIAL_LIBRARY_DIR ${SUPERLU_SERIAL_LIBRARY} DIRECTORY)
  set(SUPERLU_SERIAL_LIBRARY_DIR ${SUPERLU_SERIAL_LIBRARY_DIR} CACHE PATH "" FORCE)
endif()

# set the include dir option if it wasn't preset
if(SUPERLU_SERIAL_LIBRARY AND (NOT SUPERLU_SERIAL_INCLUDE_DIR))
  get_filename_component(SUPERLU_SERIAL_INCLUDE_DIR ${SUPERLU_SERIAL_LIBRARY_DIR} DIRECTORY)
  set(SUPERLU_SERIAL_INCLUDE_DIR "${SUPERLU_SERIAL_INCLUDE_DIR}/include" CACHE PATH "" FORCE)
endif()

# set a more informative error message in case the library was not found
set(SUPERLU_SERIAL_NOT_FOUND_MESSAGE "\
************************************************************************\n\
ERROR: Could not find SuperLU_serial. Please check the variables:\n\
       SUPERLU_SERIAL_INCLUDE_DIR and SUPERLU_SERIAL_LIBRARY_DIR\n\
************************************************************************")

# set package variables including SUPERLU_SERIAL_FOUND
find_package_handle_standard_args(SUPERLU_SERIAL
        REQUIRED_VARS
        SUPERLU_SERIAL_LIBRARY
        SUPERLU_SERIAL_LIBRARIES
        SUPERLU_SERIAL_INCLUDE_DIR
        FAIL_MESSAGE
        "${SUPERLU_SERIAL_NOT_FOUND_MESSAGE}"
)

# Create target for SuperLU_serial
if(SUPERLU_SERIAL_FOUND)

  if(NOT TARGET SUNDIALS::SUPERLU_SERIAL)
    add_library(SUNDIALS::SUPERLU_SERIAL UNKNOWN IMPORTED)
  endif()

  set_target_properties(SUNDIALS::SUPERLU_SERIAL PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${SUPERLU_SERIAL_INCLUDE_DIR}"
          INTERFACE_LINK_LIBRARIES "${SUPERLU_SERIAL_LIBRARIES}"
          IMPORTED_LOCATION "${SUPERLU_SERIAL_LIBRARY}")

endif()

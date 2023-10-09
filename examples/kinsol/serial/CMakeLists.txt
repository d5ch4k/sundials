# ---------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
#                David J. Gardner @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2023, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for KINSOL serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(KINSOL_examples
  "kinAnalytic_fp\;\;"
  "kinAnalytic_fp\;--damping_fp 0.5\;develop"
  "kinAnalytic_fp\;--m_aa 2\;"
  "kinAnalytic_fp\;--m_aa 2 --delay_aa 2\;"
  "kinAnalytic_fp\;--m_aa 2 --damping_aa 0.5\;develop"
  "kinAnalytic_fp\;--m_aa 2 --orth_aa 1\;"
  "kinAnalytic_fp\;--m_aa 2 --orth_aa 2\;"
  "kinAnalytic_fp\;--m_aa 2 --orth_aa 3\;"
  "kinFerTron_dns\;\;develop"
  "kinFoodWeb_kry\;\;exclude-single"
  "kinKrylovDemo_ls\;\;exclude-single"
  "kinLaplace_bnd\;\;exclude-single"
  "kinLaplace_picard_bnd\;\;exclude-single"
  "kinLaplace_picard_kry\;\;exclude-single"
  "kinRoberts_fp\;\;develop"
  "kinRoboKin_dns\;\;exclude-single"
  )

# Examples using LAPACK linear solvers
set(KINSOL_examples_BL
  )

# Examples using KLU linear solver
set(KINSOL_examples_KLU
  "kinFerTron_klu\;develop"
  )

# Examples using SuperLU_MT linear solver
set(KINSOL_examples_SUPERLUMT
  "kinRoboKin_slu\;develop"
  )

# Examples using SuperLU serial linear solver
set(KINSOL_examples_SUPERLU_SERIAL
        "kinRoboKin_slu_serial\;c\;develop"
        "kinRoboKin_slu_serial_cpp\;cpp\;develop"
        "kinFerTron_slu_serial\;c\;develop"
)

# Examples using Eigen Sparse serial linear solver
set(KINSOL_examples_EIGENSPARSE
        "kinRoboKin_eigen_sparse\;cpp\;develop"
        "kinRoboKin_more_eigen_sparse\;cpp\;develop"
)

# Auxiliary files to install
set(KINSOL_extras
  kinRoboKin_dns_stats.csv
  )

# Specify libraries to link against
set(KINSOL_LIB sundials_kinsol)
set(NVECS_LIB sundials_nvecserial)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${KINSOL_LIB} ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})


# Add the build and install targets for each example
foreach(example_tuple ${KINSOL_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})

    # example source files
    add_executable(${example} ${example}.c)

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})

  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    set(test_name ${example})
  else()
    string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

  # find all .out files for this example
  file(GLOB example_out ${example}*.out)

  # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.c ${example}.out
      DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/serial)
  endif()

endforeach(example_tuple ${KINSOL_examples})


# Add the build and install targets for each LAPACK example (if needed)
if(BUILD_SUNLINSOL_LAPACKBAND AND BUILD_SUNLINSOL_LAPACKDENSE)

  # Sundials LAPACK linear solver modules
  set(SUNLINSOLLAPACK_LIBS
    sundials_sunlinsollapackband
    sundials_sunlinsollapackdense)

  # LAPACK libraries
  list(APPEND SUNLINSOLLAPACK_LIBS ${LAPACK_LIBRARIES})

  foreach(example_tuple ${KINSOL_examples_BL})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_type)

    # example source files
    add_executable(${example} ${example}.c)

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # add example to regression tests
    sundials_add_test(${example} ${example}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLLAPACK_LIBS})

    # install example source and out files
    if(EXAMPLES_INSTALL)
      install(FILES ${example}.c ${example}.out
        DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/serial)
    endif()

  endforeach(example_tuple ${KINSOL_examples_BL})

endif()


# Add the build and install targets for each KLU example (if needed)
if(BUILD_SUNLINSOL_KLU)

  # Sundials KLU linear solver module
  set(SUNLINSOLKLU_LIBS sundials_sunlinsolklu)

  # KLU libraries
  list(APPEND SUNLINSOLKLU_LIBS)

  foreach(example_tuple ${KINSOL_examples_KLU})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_type)

    # example source files
    add_executable(${example} ${example}.c)

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # add example to regression tests
    sundials_add_test(${example} ${example}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLKLU_LIBS})

    # install example source and out files
    if(EXAMPLES_INSTALL)
      install(FILES ${example}.c ${example}.out
        DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/serial)
    endif()

  endforeach(example_tuple ${KINSOL_examples_KLU})

endif()


# Add the build and install targets for each SuperLU_MT example (if needed)
if(BUILD_SUNLINSOL_SUPERLUMT)

  # Sundials SuperLU_MT linear solver module
  set(SUNLINSOLSLUMT_LIBS sundials_sunlinsolsuperlumt)

  # SuperLU_MT libraries
  list(APPEND SUNLINSOLSLUMT_LIBS ${SUPERLUMT_LIBRARIES})

  foreach(example_tuple ${KINSOL_examples_SUPERLUMT})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_type)

    # example source files
    add_executable(${example} ${example}.c)

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # add example to regression tests
    sundials_add_test(${example} ${example}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLSLUMT_LIBS})

    # install example source and out files
    if(EXAMPLES_INSTALL)
      install(FILES ${example}.c ${example}.out
        DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/serial)
    endif()

  endforeach(example_tuple ${KINSOL_examples_SUPERLUMT})

endif()


# Add the build and install targets for each SuperLU serial example (if needed)
if(BUILD_SUNLINSOL_SUPERLU_SERIAL)

  # Sundials SuperLU serial linear solver module
  set(SUNLINSOLSLU_SERIAL_LIBS sundials_sunlinsolsuperlu_serial)

  # SuperLU serial libraries
  list(APPEND SUNLINSOLSLU_SERIAL_LIBS ${SUPERLU_SERIAL_LIBRARIES})

  foreach(example_tuple ${KINSOL_examples_SUPERLU_SERIAL})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 c_or_cpp)
    list(GET example_tuple 2 example_type)

    # example source files
    add_executable(${example} ${example}.${c_or_cpp})
#    if(${example} MATCHES "_cpp$")
#        add_executable(${example} ${example}.cpp)
#    else()
#        add_executable(${example} ${example}.c)
#    endif ()

    if(${c_or_cpp} STREQUAL "cpp")
        set_target_properties(${example} PROPERTIES LINKER_LANGUAGE CXX)
        target_include_directories(${example} PUBLIC /opt/local/include/eigen3)
    endif ()
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # add example to regression tests
    sundials_add_test(${example} ${example}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLSLU_SERIAL_LIBS} sundials_sunmatrix_eigensparse)
    target_link_directories(${example} PUBLIC ${SUPERLU_SERIAL_LIBRARY_DIR})

    # install example source and out files
    if(EXAMPLES_INSTALL)
        install(FILES ${example}.${c_or_cpp} ${example}.out
                DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/serial)
    endif()

  endforeach(example_tuple ${KINSOL_examples_SUPERLU_SERIAL})

endif()

# Add the build and install targets for each SuperLU serial example (if needed)
if(BUILD_SUNLINSOL_EIGENSPARSE)

    # Sundials Eigen Sparse serial linear solver module
    set(SUNLINSOL_EIGENSPARSE_LIBS sundials_sunlinsol_eigen_sparse)

    # Eigen Sparse serial libraries
    list(APPEND SUNLINSOL_EIGENSPARSE_LIBS ${EIGENSPARSE_LIBRARIES})

    foreach(example_tuple ${KINSOL_examples_EIGENSPARSE})

        # parse the example tuple
        list(GET example_tuple 0 example)
        list(GET example_tuple 1 c_or_cpp)
        list(GET example_tuple 2 example_type)

        # example source files
        add_executable(${example} ${example}.${c_or_cpp})
        #    if(${example} MATCHES "_cpp$")
        #        add_executable(${example} ${example}.cpp)
        #    else()
        #        add_executable(${example} ${example}.c)
        #    endif ()

        if(${c_or_cpp} STREQUAL "cpp")
            set_target_properties(${example} PROPERTIES LINKER_LANGUAGE CXX)
        endif ()
        set_target_properties(${example} PROPERTIES FOLDER "Examples")

        # add example to regression tests
        sundials_add_test(${example} ${example}
                ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                ANSWER_FILE ${example}.out
                EXAMPLE_TYPE ${example_type})

        # libraries to link against
        target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOL_EIGENSPARSE_LIBS})
        target_link_directories(${example} PUBLIC ${EIGENSPARSE_LIBRARY_DIR})

        # install example source and out files
        if(EXAMPLES_INSTALL)
            install(FILES ${example}.${c_or_cpp} ${example}.out
                    DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/serial)
        endif()

    endforeach(example_tuple ${KINSOL_examples_EIGENSPARSE})

endif()


# create Makfile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/serial)

  # Install the extra files
  foreach(extrafile ${KINSOL_extras})
    install(FILES ${extrafile} DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/serial)
  endforeach()

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "KINSOL")
  set(SOLVER_LIB "sundials_kinsol")

  examples2string(KINSOL_examples EXAMPLES)

  if(BUILD_SUNLINSOL_LAPACKBAND AND BUILD_SUNLINSOL_LAPACKDENSE)
    examples2string(KINSOL_examples_BL EXAMPLES_BL)
  else()
    set(EXAMPLES_BL "")
  endif()

  if(BUILD_SUNLINSOL_KLU)
    examples2string(KINSOL_examples_KLU EXAMPLES_KLU)
  else()
    set(EXAMPLES_KLU "")
  endif()

  if(BUILD_SUNLINSOL_SUPERLUMT)
    examples2string(KINSOL_examples_SUPERLUMT EXAMPLES_SLUMT)
    if (SUNDIALS_SUPERLUMT_THREAD_TYPE STREQUAL "PTHREAD")
      set(THREAD_LIBRARY_SLUMT ${CMAKE_THREAD_LIBS_INIT})
    else()
      set(THREAD_LIBRARY_SLUMT "")
    endif()
  else()
    set(EXAMPLES_SLUMT "")
    set(THREAD_LIBRARY_SLUMT "")
  endif()
  if(BUILD_SUNLINSOL_SUPERLU_SERIAL)
    examples2string(KINSOL_examples_SUPERLU_SERIAL EXAMPLES_SLU_SERIAL)
  else()
    set(EXAMPLES_SLU_SERIAL "")
  endif()
  if(BUILD_SUNLINSOL_EIGENSPARSE)
      examples2string(KINSOL_examples_EIGENSPARSE EXAMPLES_EIGENSPARSE)
  else()
      set(EXAMPLES_EIGENSPARSE "")
  endif()

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
    ${PROJECT_BINARY_DIR}/examples/kinsol/serial/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/kinsol/serial/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/serial
    )

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/kinsol/serial/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/kinsol/serial/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/kinsol/serial
      RENAME Makefile
      )
  endif()

  # add test_install target
  sundials_add_test_install(kinsol serial EXECUTABLE kinAnalytic_fp)

endif()

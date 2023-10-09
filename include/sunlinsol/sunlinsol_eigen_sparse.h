/*
 * -----------------------------------------------------------------
 * Programmer(s): Martin Posch
 * Based on codes sundials_superlu_impl.h and <solver>_superlu.h
 *     written by Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the SUNLINSOL_EIGENSPARSE implementation of the
 * SUNLINSOL module, SUNLINSOL_EIGENSPARSE.
 *
 * Note:
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file sundials_linearsolver.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_EIGEN_SPARSE_H
#define _SUNLINSOL_EIGEN_SPARSE_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

using SpMat = Eigen::SparseMatrix<realtype, Eigen::ColMajor, sunindextype>;
using Vec = Eigen::Matrix<realtype , Eigen::Dynamic, 1>;

/* --------------------------------------------
 * Eigen Sparse Implementation of SUNLinearSolver
 * -------------------------------------------- */

struct _SUNLinearSolverContent_EigenSparse {
  int          last_flag;
  int          first_factorize;
  //SpMat *A;
  Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<sunindextype>> * LU;
  //sunindextype N;
  //int          ordering;
};

typedef struct _SUNLinearSolverContent_EigenSparse *SUNLinearSolverContent_EigenSparse;


/* -------------------------------------------
 * Exported Functions for SUNLINSOL_EIGENSPARSE
 * ------------------------------------------- */

SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_EigenSparse(N_Vector y, SUNMatrix A, SUNContext sunctx);
SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_EigenSparse(SUNLinearSolver S);
SUNDIALS_EXPORT SUNLinearSolver_ID SUNLinSolGetID_EigenSparse(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_EigenSparse(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetup_EigenSparse(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_EigenSparse(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_EigenSparse(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSpace_EigenSparse(SUNLinearSolver S, long int *lenrwLS, long int *leniwLS);
SUNDIALS_EXPORT int SUNLinSolFree_EigenSparse(SUNLinearSolver S);


#ifdef __cplusplus
}
#endif

#endif

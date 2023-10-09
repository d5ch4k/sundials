/* -----------------------------------------------------------------
 * Programmer(s): Martin Posch
 * -----------------------------------------------------------------
 * Based on codes <solver>_superlu.c, written by
 * Carol S. Woodward @ LLNL
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
 * This is the implementation file for the SUNLINSOL_EIGEN_SPARSE serial implementation
 * of the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_eigen_sparse.h>
#include <sundials/sundials_math.h>

#define ZERO      RCONST(0.0)
#define ONE       RCONST(1.0)
#define TWO       RCONST(2.0)

/*
 * -----------------------------------------------------------------
 * SUNLINSOL_EIGEN_SPARSE solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define SLU_CONTENT(S)    ( (SUNLinearSolverContent_EigenSparse)(S->content) )
#define LASTFLAG(S)         (  SLU_CONTENT(S)->last_flag )
#define FIRSTFACTORIZE(S)   (  SLU_CONTENT(S)->first_factorize )
//#define SM_A(S)             (  SLU_CONTENT(S)->A )
#define SM_LU(S)             (  SLU_CONTENT(S)->LU )
//#define SIZE(S)             (  SLU_CONTENT(S)->N )
//#define ORDERING(S)         (  SLU_CONTENT(S)->ordering )
//#define OPTIONS(S)          (&(SLU_CONTENT(S)->options) )

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new SUNLINSOL_EIGEN_SPARSE linear solver
 */

SUNLinearSolver SUNLinSol_EigenSparse(N_Vector y, SUNMatrix A, SUNContext sunctx)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_EigenSparse content;
  sunindextype MatrixRows;

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE) return(NULL);

  if (SUNSparseMatrix_Rows(A) != SUNSparseMatrix_Columns(A)) return(NULL);

  if (N_VGetVectorID(y) != SUNDIALS_NVEC_SERIAL) return(NULL);

  MatrixRows = SUNSparseMatrix_Rows(A);
  if (MatrixRows != N_VGetLength(y)) return(NULL);

  /* Create an empty linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty(sunctx);
  if (S == NULL) return(NULL);

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_EigenSparse;
  S->ops->getid      = SUNLinSolGetID_EigenSparse;
  S->ops->initialize = SUNLinSolInitialize_EigenSparse;
  S->ops->setup      = SUNLinSolSetup_EigenSparse;
  S->ops->solve      = SUNLinSolSolve_EigenSparse;
  S->ops->lastflag   = SUNLinSolLastFlag_EigenSparse;
  S->ops->space      = SUNLinSolSpace_EigenSparse;
  S->ops->free       = SUNLinSolFree_EigenSparse;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_EigenSparse) malloc(sizeof *content);
  if (content == NULL) { SUNLinSolFree(S); return(NULL); }

  /* Attach content */
  S->content = content;

  /* Fill content */
  //content->N                 = MatrixRows;
  content->last_flag         = 0;
  //content->ordering          = 0;
  //content->A                 = NULL;
  content->LU                = NULL;

  /* Allocate content */
  //content->A = new SpMat(MatrixRows, MatrixRows);
  //if (content->A == NULL) { SUNLinSolFree(S); return(NULL); }
  return(S);
}


/* ----------------------------------------------------------------------------
 * Function to set the ordering type for a EigenSparse linear solver
 */

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */


SUNLinearSolver_Type SUNLinSolGetType_EigenSparse(SUNLinearSolver S) { return (SUNLINEARSOLVER_DIRECT); }

SUNLinearSolver_ID SUNLinSolGetID_EigenSparse(SUNLinearSolver S) { return (SUNLINEARSOLVER_EIGEN_SPARSE); }

int SUNLinSolInitialize_EigenSparse(SUNLinearSolver S)
{
  /* force a first factorization */
  FIRSTFACTORIZE(S) = 1;

  /* Initialize statistics variables */
  // StatInit(SLUSTAT(S));

  if (SM_LU(S)) delete (SM_LU(S));
  (SM_LU(S)) = new Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<sunindextype>>();

  LASTFLAG(S) = SUNLS_SUCCESS;
  return (LASTFLAG(S));
}

int SUNLinSolSetup_EigenSparse(SUNLinearSolver S, SUNMatrix A)
{
  // int n = SIZE(S);
  int nrows             = SUNSparseMatrix_Rows(A);
  int ncols             = SUNSparseMatrix_Columns(A);
  int nnz               = SUNSparseMatrix_NNZ(A);
  sunindextype* colptrs = SUNSparseMatrix_IndexPointers(A);
  sunindextype* rowvals = SUNSparseMatrix_IndexValues(A);
  realtype* data        = SUNSparseMatrix_Data(A);

  Eigen::Map<Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>> SpA(nrows, ncols, nnz, colptrs, rowvals, data);
  /* On first decomposition, set up reusable pieces */
  if (FIRSTFACTORIZE(S))
  {
    SM_LU(S)->compute(SpA); // analyze and factorize
    FIRSTFACTORIZE(S) = 0;
  }
  else { SM_LU(S)->factorize(SpA); }
  if (SM_LU(S)->info() == Eigen::Success) { LASTFLAG(S) = SUNLS_SUCCESS; }
  else { LASTFLAG(S) = SUNLS_LUFACT_FAIL; }

  return (LASTFLAG(S));
}

int SUNLinSolSolve_EigenSparse(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, realtype tol)
{
  /* copy b into x */
  N_VScale(ONE, b, x);

  /* access x data array */
  realtype* xdata = N_VGetArrayPointer(x);
  realtype* bdata = N_VGetArrayPointer(b);
  sunindextype nx = N_VGetLength(x);
  sunindextype nb = N_VGetLength(b);
  if (xdata == NULL)
  {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return (LASTFLAG(S));
  }
  if (bdata == NULL)
  {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return (LASTFLAG(S));
  }

  int nrows             = SUNSparseMatrix_Rows(A);
  int ncols             = SUNSparseMatrix_Columns(A);
  int nnz               = SUNSparseMatrix_NNZ(A);
  sunindextype* colptrs = SUNSparseMatrix_IndexPointers(A);
  sunindextype* rowvals = SUNSparseMatrix_IndexValues(A);
  realtype* data        = SUNSparseMatrix_Data(A);

  Eigen::Map<Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>> SpA(nrows, ncols, nnz, colptrs, rowvals, data);

  // factorize already done in Setup
  /*
  SM_LU(S)->factorize(SpA);
  if (SM_LU(S)->info() == Eigen::Success) {
    LASTFLAG(S) = SUNLS_SUCCESS;
  }
  else {
    LASTFLAG(S) = SUNLS_LUFACT_FAIL;
  }
  */

  Eigen::Map<Vec> bb(bdata, nb);
  Eigen::Map<Vec> xx(xdata, nx);

  xx = SM_LU(S)->solve(bb);

  if (SM_LU(S)->info() == Eigen::Success) { LASTFLAG(S) = SUNLS_SUCCESS; }
  else { LASTFLAG(S) = SUNLS_LUFACT_FAIL; }

  return (LASTFLAG(S));
}

sunindextype SUNLinSolLastFlag_EigenSparse(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  if (S == NULL) return (-1);
  return (LASTFLAG(S));
}

int SUNLinSolSpace_EigenSparse(SUNLinearSolver S, long int* lenrwLS, long int* leniwLS)
{
  /* TODO */
  *leniwLS = 5;
  *lenrwLS = 1;
  return (SUNLS_SUCCESS);
}

int SUNLinSolFree_EigenSparse(SUNLinearSolver S)
{
  /* return with success if already freed */
  if (S == NULL) return (SUNLS_SUCCESS);

  /* delete items from the contents structure (if it exists) */
  if (S->content)
  {
    if (SM_LU(S))
    {
      free(SM_LU(S));
      SM_LU(S) = NULL;
    }

    free(S->content);
    S->content = NULL;
  }

  /* delete generic structures */
  if (S->ops)
  {
    free(S->ops);
    S->ops = NULL;
  }
  free(S);
  S = NULL;
  return (SUNLS_SUCCESS);
}


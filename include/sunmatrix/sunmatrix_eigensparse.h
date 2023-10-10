/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner @ LLNL
 * Based on code sundials_eigensparse.h by: Carol Woodward and
 *     Slaven Peles @ LLNL, and Daniel R. Reynolds @ SMU
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
 * This is the header file for the sparse implementation of the
 * SUNMATRIX module, SUNMATRIX_EIGENSPARSE.
 *
 * Notes:
 *   - The definition of the generic SUNMatrix structure can be found
 *     in the header file sundials_matrix.h.
 *   - The definition of the type 'realtype' can be found in the
 *     header file sundials_types.h, and it may be changed (at the
 *     configuration stage) according to the user's needs.
 *     The sundials_types.h file also contains the definition
 *     for the type 'booleantype' and 'indextype'.
 * -----------------------------------------------------------------
 */

#ifndef _SUNMATRIX_EIGENSPARSE_H
#define _SUNMATRIX_EIGENSPARSE_H

#include <stdio.h>
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_band.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ------------------------
 * Matrix Type Definitions
 * ------------------------ */

#define CSC_MAT 0
#define CSR_MAT 1


/* ------------------------------------------
 * Sparse Implementation of SUNMATRIX_EIGENSPARSE
 * ------------------------------------------ */

struct _SUNMatrixContent_EigenSparse {
  Eigen::SparseMatrix<realtype, Eigen::ColMajor, sunindextype>* A_CSC = nullptr;
  Eigen::SparseMatrix<realtype, Eigen::RowMajor, sunindextype>* A_CSR = nullptr;
//  sunindextype M;
//  sunindextype N;
  sunindextype NNZ;
//  sunindextype NP;
//  realtype *data;
//  int sparsetype;
//  sunindextype *indexvals;
//  sunindextype *indexptrs;
//  /* CSC indices */
//  sunindextype **rowvals;
//  sunindextype **colptrs;
//  /* CSR indices */
//  sunindextype **colvals;
//  sunindextype **rowptrs;
};

typedef struct _SUNMatrixContent_EigenSparse *SUNMatrixContent_EigenSparse;


/* ---------------------------------------
 * Macros for access to SUNMATRIX_EigenSparse
 * --------------------------------------- */

#define SM_CONTENT_S(A)     ( (SUNMatrixContent_EigenSparse)(A->content) )

#define SM_ROWS_S(A)        ( SM_CONTENT_S(A)->M )

#define SM_COLUMNS_S(A)     ( SM_CONTENT_S(A)->N )

#define SM_NNZ_S(A)         ( SM_CONTENT_S(A)->NNZ )

#define SM_NP_S(A)          ( SM_CONTENT_S(A)->NP )

#define SM_SPARSETYPE_S(A)  ( SM_CONTENT_S(A)->sparsetype )

#define SM_DATA_S(A)        ( SM_CONTENT_S(A)->data )

#define SM_INDEXVALS_S(A)   ( SM_CONTENT_S(A)->indexvals )

#define SM_INDEXPTRS_S(A)   ( SM_CONTENT_S(A)->indexptrs )

/* ----------------------------------------
 * Exported Functions for SUNMATRIX_EigenSparse
 * ---------------------------------------- */

SUNDIALS_EXPORT SUNMatrix SUNEigenSparseMatrix(sunindextype M, sunindextype N, sunindextype NNZ, int sparsetype, SUNContext sunctx);

SUNDIALS_EXPORT SUNMatrix SUNEigenSparseFromDenseMatrix(SUNMatrix A, realtype droptol, int sparsetype);

//SUNDIALS_EXPORT SUNMatrix SUNEigenSparseFromBandMatrix(SUNMatrix A, realtype droptol, int sparsetype);

//SUNDIALS_EXPORT int SUNEigenSparseMatrix_ToCSR(const SUNMatrix A, SUNMatrix* Bout);
//SUNDIALS_EXPORT int SUNEigenSparseMatrix_ToCSC(const SUNMatrix A, SUNMatrix* Bout);

//SUNDIALS_EXPORT int SUNEigenSparseMatrix_Realloc(SUNMatrix A);

//SUNDIALS_EXPORT int SUNEigenSparseMatrix_Reallocate(SUNMatrix A, sunindextype NNZ);

//SUNDIALS_EXPORT void SUNEigenSparseMatrix_Print(SUNMatrix A, FILE* outfile);

//SUNDIALS_EXPORT sunindextype SUNEigenSparseMatrix_Rows(SUNMatrix A);
//SUNDIALS_EXPORT sunindextype SUNEigenSparseMatrix_Columns(SUNMatrix A);
//SUNDIALS_EXPORT sunindextype SUNEigenSparseMatrix_NNZ(SUNMatrix A);
//SUNDIALS_EXPORT sunindextype SUNEigenSparseMatrix_NP(SUNMatrix A);
//SUNDIALS_EXPORT int SUNEigenSparseMatrix_SparseType(SUNMatrix A);
//SUNDIALS_EXPORT realtype* SUNEigenSparseMatrix_Data(SUNMatrix A);
//SUNDIALS_EXPORT sunindextype* SUNEigenSparseMatrix_IndexValues(SUNMatrix A);
//SUNDIALS_EXPORT sunindextype* SUNEigenSparseMatrix_IndexPointers(SUNMatrix A);

SUNDIALS_EXPORT SUNMatrix_ID SUNMatGetID_EigenSparse(SUNMatrix A);
//SUNDIALS_EXPORT SUNMatrix SUNMatClone_EigenSparse(SUNMatrix A);
SUNDIALS_EXPORT void SUNMatDestroy_EigenSparse(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatZero_EigenSparse(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatCopy_EigenSparse(SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatEigenCSCCopy_EigenSparse(Eigen::SparseMatrix<realtype, Eigen::ColMajor, sunindextype> const & A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatEigenCSRCopy_EigenSparse(Eigen::SparseMatrix<realtype, Eigen::RowMajor, sunindextype>const & A, SUNMatrix B);
//SUNDIALS_EXPORT int SUNMatScaleAdd_EigenSparse(realtype c, SUNMatrix A, SUNMatrix B);
//SUNDIALS_EXPORT int SUNMatScaleAddI_EigenSparse(realtype c, SUNMatrix A);
//SUNDIALS_EXPORT int SUNMatMatvec_EigenSparse(SUNMatrix A, N_Vector x, N_Vector y);
//SUNDIALS_EXPORT int SUNMatSpace_EigenSparse(SUNMatrix A, long int *lenrw, long int *leniw);


#ifdef __cplusplus
}
#endif

#endif

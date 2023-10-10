/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner @ LLNL
 * Based on code sundials_sparse.c by: Carol Woodward and
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
 * This is the implementation file for the sparse implementation of
 * the SUNMATRIX package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sunmatrix/sunmatrix_eigensparse.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* Private function prototypes */
static booleantype SMCompatible_EigenSparse(SUNMatrix A, SUNMatrix B);
static booleantype SMCompatible2_EigenSparse(SUNMatrix A, N_Vector x, N_Vector y);
static int Matvec_EigenSparseCSC(SUNMatrix A, N_Vector x, N_Vector y);
static int Matvec_EigenSparseCSR(SUNMatrix A, N_Vector x, N_Vector y);
static int format_convert(const SUNMatrix A, SUNMatrix B);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/*
 * ==================================================================
 * Private function prototypes (functions working on SlsMat)
 * ==================================================================
 */

/* ----------------------------------------------------------------------------
 * Function to create a new sparse matrix
 */

SUNMatrix SUNEigenSparseMatrix(sunindextype M, sunindextype N,
                          sunindextype NNZ, int sparsetype,
                          SUNContext sunctx)
{

  /* return with NULL matrix on illegal input */
  if ( (M <= 0) || (N <= 0) || (NNZ < 0) ) return(NULL);
  if ( (sparsetype != CSC_MAT) && (sparsetype != CSR_MAT) ) return(NULL);

  /* Create an empty matrix object */
  SUNMatrix A;
  A = NULL;
  A = SUNMatNewEmpty(sunctx);
  if (A == NULL) return(NULL);

  /* Attach operations */
  A->ops->getid     = SUNMatGetID_EigenSparse;
  A->ops->clone     = 0;//SUNMatClone_EigenSparse;
  A->ops->destroy   = SUNMatDestroy_EigenSparse;
  A->ops->zero      = SUNMatZero_EigenSparse;
  A->ops->copy      = 0;//SUNMatCopy_EigenSparse;
  A->ops->scaleadd  = 0;//SUNMatScaleAdd_EigenSparse;
  A->ops->scaleaddi = 0;//SUNMatScaleAddI_EigenSparse;
  A->ops->matvec    = 0;//SUNMatMatvec_EigenSparse;
  A->ops->space     = 0;//SUNMatSpace_EigenSparse;

  /* Create content */
  SUNMatrixContent_EigenSparse content;
  content = NULL;
  content = (SUNMatrixContent_EigenSparse) malloc(sizeof *content);
  if (content == NULL) { SUNMatDestroy(A); return(NULL); }

  /* Attach content */
  A->content = content;

  /* Fill content */
  content->A_CSC = nullptr;
  content->A_CSR = nullptr;

  if (sparsetype == CSC_MAT)
  {
    content->A_CSC = new Eigen::SparseMatrix<realtype, Eigen::ColMajor, sunindextype>(M, N);
  }
  else
  {
    content->A_CSR = new Eigen::SparseMatrix<realtype, Eigen::RowMajor, sunindextype>(M, N);
  }

  content->NNZ = NNZ;

  return(A);
}




/* ----------------------------------------------------------------------------
 * Function to create a new sparse matrix from an existing dense matrix
 * by copying all nonzero values into the sparse matrix structure.  Returns NULL
 * if the request for matrix storage cannot be satisfied.
 */

SUNMatrix SUNEigenSparseFromDenseMatrix(SUNMatrix Ad, realtype droptol,
                                   int sparsetype)
{
  sunindextype i, j, nnz;
  sunindextype M, N;
  SUNMatrix As;

  /* check for legal sparsetype, droptol and input matrix type */
  if ( (sparsetype != CSR_MAT) && (sparsetype != CSC_MAT) )
    return NULL;
  if ( droptol < ZERO )
    return NULL;
  if (SUNMatGetID(Ad) != SUNMATRIX_DENSE)
    return NULL;

  /* set size of new matrix */
  M = SM_ROWS_D(Ad);
  N = SM_COLUMNS_D(Ad);

  Eigen::Map<Eigen::Matrix<realtype,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> aa(SM_DATA_D(Ad),M,N);

  /* determine total number of nonzeros */
  nnz = 0;
//  for (j=0; j<N; j++)
//    for (i=0; i<M; i++)
//      nnz += (SUNRabs(SM_ELEMENT_D(Ad,i,j)) > droptol);

  /* allocate sparse matrix */
  As = NULL;
  As = SUNEigenSparseMatrix(M, N, nnz, sparsetype, Ad->sunctx);
  if (As == NULL)  return NULL;

  /* Create content */
  SUNMatrixContent_EigenSparse content;
  content = NULL;
  content = (SUNMatrixContent_EigenSparse) malloc(sizeof *content);
  if (content == NULL) { SUNMatDestroy(As); return(NULL); }

  /* Attach content */
  As->content = content;

  /* Fill content */
  content->A_CSC = nullptr;
  content->A_CSR = nullptr;

  if (sparsetype == CSC_MAT)
  {
    content->A_CSC = new Eigen::SparseMatrix<realtype, Eigen::ColMajor, sunindextype>(M, N);
    *(content->A_CSC) = aa.sparseView(droptol);
    content->NNZ = content->A_CSC->nonZeros();
  }
  else
  {
    content->A_CSR = new Eigen::SparseMatrix<realtype, Eigen::RowMajor, sunindextype>(M, N);
    *(content->A_CSR) = aa.sparseView(droptol);
    content->NNZ = content->A_CSR->nonZeros();
  }

  return(As);
}


/* ----------------------------------------------------------------------------
 * Function to create a new sparse matrix from an existing band matrix
 * by copying all nonzero values into the sparse matrix structure.  Returns NULL
 * if the request for matrix storage cannot be satisfied.
 */

//SUNMatrix SUNEigenSparseFromBandMatrix(SUNMatrix Ad, realtype droptol, int sparsetype)
//{
//  sunindextype i, j, nnz;
//  sunindextype M, N;
//  SUNMatrix As;
//
//  /* check for legal sparsetype, droptol and input matrix type */
//  if ( (sparsetype != CSR_MAT) && (sparsetype != CSC_MAT) )
//    return NULL;
//  if ( droptol < ZERO )
//    return NULL;
//  if (SUNMatGetID(Ad) != SUNMATRIX_BAND)
//    return NULL;
//
//  /* set size of new matrix */
//  M = SM_ROWS_B(Ad);
//  N = SM_COLUMNS_B(Ad);
//
//  /* determine total number of nonzeros */
//  nnz = 0;
//  for (j=0; j<N; j++)
//    for (i=SUNMAX(0,j-SM_UBAND_B(Ad)); i<=SUNMIN(M-1,j+SM_LBAND_B(Ad)); i++)
//      nnz += (SUNRabs(SM_ELEMENT_B(Ad,i,j)) > droptol);
//
//  /* allocate sparse matrix */
//  As = SUNEigenSparseMatrix(M, N, nnz, sparsetype, Ad->sunctx);
//  if (As == NULL)  return NULL;
//
//  /* copy nonzeros from Ad into As, based on CSR/CSC type */
//  nnz = 0;
//  if (sparsetype == CSC_MAT) {
//    for (j=0; j<N; j++) {
//      (SM_INDEXPTRS_S(As))[j] = nnz;
//      for (i=SUNMAX(0,j-SM_UBAND_B(Ad)); i<=SUNMIN(M-1,j+SM_LBAND_B(Ad)); i++) {
//        if ( SUNRabs(SM_ELEMENT_B(Ad,i,j)) > droptol ) {
//          (SM_INDEXVALS_S(As))[nnz] = i;
//          (SM_DATA_S(As))[nnz++] = SM_ELEMENT_B(Ad,i,j);
//        }
//      }
//    }
//    (SM_INDEXPTRS_S(As))[N] = nnz;
//  } else {       /* CSR_MAT */
//    for (i=0; i<M; i++) {
//      (SM_INDEXPTRS_S(As))[i] = nnz;
//      for (j=SUNMAX(0,i-SM_LBAND_B(Ad)); j<=SUNMIN(N-1,i+SM_UBAND_B(Ad)); j++) {
//        if ( SUNRabs(SM_ELEMENT_B(Ad,i,j)) > droptol ) {
//          (SM_INDEXVALS_S(As))[nnz] = j;
//          (SM_DATA_S(As))[nnz++] = SM_ELEMENT_B(Ad,i,j);
//        }
//      }
//    }
//    (SM_INDEXPTRS_S(As))[M] = nnz;
//  }
//
//  return(As);
//}


/* ----------------------------------------------------------------------------
 * Functions to access the contents of the sparse matrix structure
 */

sunindextype SUNEigenSparseMatrix_Rows(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_EIGENSPARSE)
  {
  auto pc = (SUNMatrixContent_EigenSparse)A->content;
  if (!pc->A_CSC && !pc->A_CSR)
      return SUNMAT_ILL_INPUT;
  if (pc->A_CSC)
      return pc->A_CSC->rows();
  else
      return pc->A_CSR->rows();
  }
  else
    return SUNMAT_ILL_INPUT;
}

sunindextype SUNEigenSparseMatrix_Columns(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_EIGENSPARSE)
  {
    auto pc = (SUNMatrixContent_EigenSparse)A->content;
    if (!pc->A_CSC && !pc->A_CSR)
      return SUNMAT_ILL_INPUT;
    if (pc->A_CSC)
      return pc->A_CSC->cols();
    else
      return pc->A_CSR->cols();
  }
  else
    return SUNMAT_ILL_INPUT;
}

sunindextype SUNEigenSparseMatrix_NNZ(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_EIGENSPARSE)
  {
    auto pc = (SUNMatrixContent_EigenSparse)A->content;
    if (!pc->A_CSC && !pc->A_CSR)
      return SUNMAT_ILL_INPUT;
    if (pc->A_CSC)
      return pc->A_CSC->nonZeros();
    else
      return pc->A_CSR->nonZeros();
  }
  else
    return SUNMAT_ILL_INPUT;
}

sunindextype SUNEigenSparseMatrix_NP(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_EIGENSPARSE)
  {
    auto pc = (SUNMatrixContent_EigenSparse)A->content;
    if (!pc->A_CSC && !pc->A_CSR)
      return SUNMAT_ILL_INPUT;
    if (pc->A_CSC)
      return pc->A_CSC->outerSize()+1;
    else
      return pc->A_CSR->outerSize()+1;
  }
  else
    return SUNMAT_ILL_INPUT;
}

int SUNEigenSparseMatrix_EigenSparseType(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_EIGENSPARSE)
  {
    auto pc = (SUNMatrixContent_EigenSparse)A->content;
    if (!pc->A_CSC && !pc->A_CSR)
      return SUNMAT_ILL_INPUT;
    if (pc->A_CSC)
      return CSC_MAT;
    else
      return CSR_MAT;
  }
  else
    return SUNMAT_ILL_INPUT;
}

realtype* SUNEigenSparseMatrix_Data(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_EIGENSPARSE)
  {
    auto pc = (SUNMatrixContent_EigenSparse)A->content;
    if (!pc->A_CSC && !pc->A_CSR)
      return nullptr;
    if (pc->A_CSC)
      return pc->A_CSC->valuePtr();
    else
      return pc->A_CSR->valuePtr();
  }
  else
    return nullptr;
}

sunindextype* SUNEigenSparseMatrix_IndexValues(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_EIGENSPARSE)
  {
    auto pc = (SUNMatrixContent_EigenSparse)A->content;
    if (!pc->A_CSC && !pc->A_CSR)
      return nullptr;
    if (pc->A_CSC)
      return pc->A_CSC->innerIndexPtr();
    else
      return pc->A_CSR->innerIndexPtr();
  }
  else
    return nullptr;
}

sunindextype* SUNEigenSparseMatrix_IndexPointers(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_EIGENSPARSE)
  {
    auto pc = (SUNMatrixContent_EigenSparse)A->content;
    if (!pc->A_CSC && !pc->A_CSR)
      return nullptr;
    if (pc->A_CSC)
      return pc->A_CSC->outerIndexPtr();
    else
      return pc->A_CSR->outerIndexPtr();
  }
  else
    return nullptr;
}


/*
 * -----------------------------------------------------------------
 * implementation of matrix operations
 * -----------------------------------------------------------------
 */

SUNMatrix_ID SUNMatGetID_EigenSparse(SUNMatrix A)
{
  return SUNMATRIX_EIGENSPARSE;
}

void SUNMatDestroy_EigenSparse(SUNMatrix A)
{
  if (A == NULL) return;

  /* free content */
  if (A->content != NULL) {

    /* free data array */
    auto pc = (SUNMatrixContent_EigenSparse)A->content;
    if (pc->A_CSC) {
      delete pc->A_CSC;
      pc->A_CSC = nullptr;
    }
    if (pc->A_CSR)
    {
      delete pc->A_CSR;
      pc->A_CSR = nullptr;
    }
    /* free content struct */
    free(A->content);
    A->content = NULL;
  }


  /* free ops and matrix */
  if (A->ops) { free(A->ops); A->ops = NULL; }
  free(A); A = NULL;

  return;
}

int SUNMatZero_EigenSparse(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_EIGENSPARSE)
  {
    auto pc = (SUNMatrixContent_EigenSparse)A->content;
    if (!pc->A_CSC && !pc->A_CSR)
      return SUNMAT_ILL_INPUT;
    if (pc->A_CSC)
       pc->A_CSC->setZero();
    else
       pc->A_CSR->setZero();
  }
  else
    return SUNMAT_ILL_INPUT;

  return SUNMAT_SUCCESS;
}

int SUNMatCopy_EigenSparse(SUNMatrix A, SUNMatrix B)
{
  /* Verify that A and B are compatible */
  if (!SMCompatible_EigenSparse(A, B))
    return SUNMAT_ILL_INPUT;

  /* Perform operation */
  auto ca = (SUNMatrixContent_EigenSparse)A->content;
  auto cb = (SUNMatrixContent_EigenSparse)B->content;

  if (ca->A_CSC){
    *cb->A_CSC = *ca->A_CSC;
  }

  if (ca->A_CSR){
    *cb->A_CSR = *ca->A_CSR;
  }

  return SUNMAT_SUCCESS;
}

int SUNMatEigenCSCCopy_EigenSparse(Eigen::SparseMatrix<realtype, Eigen::ColMajor, sunindextype> const & A, SUNMatrix B)
{
  auto cb = (SUNMatrixContent_EigenSparse)B->content;

  if (ca->A_CSC){
    *cb->A_CSC = A;
  }
  else
    return SUNMAT_ILL_INPUT;
  return SUNMAT_SUCCESS;
}

int SUNMatEigenCSRCopy_EigenSparse(Eigen::SparseMatrix<realtype, Eigen::RowMajor, sunindextype> const & A, SUNMatrix B)
{
  auto cb = (SUNMatrixContent_EigenSparse)B->content;

  if (ca->A_CSR){
    *cb->A_CSR = A;
  }
  else
    return SUNMAT_ILL_INPUT;
  return SUNMAT_SUCCESS;
}

//int SUNMatScaleAddI_EigenSparse(realtype c, SUNMatrix A)
//{
//  sunindextype j, i, p, nz, newvals, M, N, cend, nw;
//  booleantype newmat, found;
//  sunindextype *w, *Ap, *Ai, *Cp, *Ci;
//  realtype *x, *Ax, *Cx;
//  SUNMatrix C;
//
//  /* store shortcuts to matrix dimensions (M is inner dimension, N is outer) */
//  if (SM_SPARSETYPE_S(A) == CSC_MAT) {
//    M = SM_ROWS_S(A);
//    N = SM_COLUMNS_S(A);
//  }
//  else {
//    M = SM_COLUMNS_S(A);
//    N = SM_ROWS_S(A);
//  }
//
//  /* access data arrays from A (return if failure) */
//  Ap = Ai = NULL;
//  Ax = NULL;
//  if (SM_INDEXPTRS_S(A))  Ap = SM_INDEXPTRS_S(A);
//  else  return (SUNMAT_MEM_FAIL);
//  if (SM_INDEXVALS_S(A))  Ai = SM_INDEXVALS_S(A);
//  else  return (SUNMAT_MEM_FAIL);
//  if (SM_DATA_S(A))       Ax = SM_DATA_S(A);
//  else  return (SUNMAT_MEM_FAIL);
//
//
//  /* determine if A: contains values on the diagonal (so I can just be added in);
//     if not, then increment counter for extra storage that should be required. */
//  newvals = 0;
//  for (j=0; j < SUNMIN(M,N); j++) {
//    /* scan column (row if CSR) of A, searching for diagonal value */
//    found = SUNFALSE;
//    for (i=Ap[j]; i<Ap[j+1]; i++) {
//      if (Ai[i] == j) {
//        found = SUNTRUE;
//        break;
//      }
//    }
//    /* if no diagonal found, increment necessary storage counter */
//    if (!found)  newvals += 1;
//  }
//
//  /* If extra nonzeros required, check whether matrix has sufficient storage space
//     for new nonzero entries  (so I can be inserted into existing storage) */
//  newmat = SUNFALSE;   /* no reallocation needed */
//  if (newvals > (SM_NNZ_S(A) - Ap[N]))
//    newmat = SUNTRUE;
//
//
//  /* perform operation based on existing/necessary structure */
//
//  /*   case 1: A already contains a diagonal */
//  if (newvals == 0) {
//
//    /* iterate through columns, adding 1.0 to diagonal */
//    for (j=0; j < SUNMIN(M,N); j++)
//      for (i=Ap[j]; i<Ap[j+1]; i++)
//        if (Ai[i] == j) {
//          Ax[i] = ONE + c*Ax[i];
//        } else {
//          Ax[i] = c*Ax[i];
//        }
//
//
//  /*   case 2: A has sufficient storage, but does not already contain a diagonal */
//  } else if (!newmat) {
//
//    /* create work arrays for nonzero row (column) indices and values in a single column (row) */
//    w = (sunindextype *) malloc(M * sizeof(sunindextype));
//    x = (realtype *) malloc(M * sizeof(realtype));
//
//    /* determine storage location where last column (row) should end */
//    nz = Ap[N] + newvals;
//
//    /* store pointer past last column (row) from original A,
//       and store updated value in revised A */
//    cend = Ap[N];
//    Ap[N] = nz;
//
//    /* iterate through columns (rows) backwards */
//    for (j=N-1; j>=0; j--) {
//
//      /* reset diagonal entry, in case it's not in A */
//      x[j] = ZERO;
//
//      /* iterate down column (row) of A, collecting nonzeros */
//      for (p=Ap[j], i=0; p<cend; p++, i++) {
//        w[i] = Ai[p];        /* collect row (column) index */
//        x[Ai[p]] = c*Ax[p];    /* collect/scale value */
//      }
//
//      /* NNZ in this column (row) */
//      nw = cend - Ap[j];
//
//      /* add identity to this column (row) */
//      if (j < M) {
//        x[j] += ONE;   /* update value */
//      }
//
//      /* fill entries of A with this column's (row's) data */
//      /* fill entries past diagonal */
//      for (i=nw-1; i>=0 && w[i]>j; i--) {
//        Ai[--nz] = w[i];
//        Ax[nz] = x[w[i]];
//      }
//      /* fill diagonal if applicable */
//      if (i < 0 /* empty or insert at front */ || w[i] != j /* insert behind front */) {
//        Ai[--nz] = j;
//        Ax[nz] = x[j];
//      }
//      /* fill entries before diagonal */
//      for (; i>=0; i--) {
//        Ai[--nz] = w[i];
//        Ax[nz] = x[w[i]];
//      }
//
//      /* store ptr past this col (row) from orig A, update value for new A */
//      cend = Ap[j];
//      Ap[j] = nz;
//
//    }
//
//    /* clean up */
//    free(w);
//    free(x);
//
//
//  /*   case 3: A must be reallocated with sufficient storage */
//  } else {
//
//    /* create work array for nonzero values in a single column (row) */
//    x = (realtype *) malloc(M * sizeof(realtype));
//
//    /* create new matrix for sum */
//    C = SUNEigenSparseMatrix(SM_ROWS_S(A), SM_COLUMNS_S(A),
//                        Ap[N] + newvals,
//                        SM_SPARSETYPE_S(A), A->sunctx);
//
//    /* access data from CSR structures (return if failure) */
//    Cp = Ci = NULL;
//    Cx = NULL;
//    if (SM_INDEXPTRS_S(C))  Cp = SM_INDEXPTRS_S(C);
//    else  return (SUNMAT_MEM_FAIL);
//    if (SM_INDEXVALS_S(C))  Ci = SM_INDEXVALS_S(C);
//    else  return (SUNMAT_MEM_FAIL);
//    if (SM_DATA_S(C))       Cx = SM_DATA_S(C);
//    else  return (SUNMAT_MEM_FAIL);
//
//    /* initialize total nonzero count */
//    nz = 0;
//
//    /* iterate through columns (rows for CSR) */
//    for (j=0; j<N; j++) {
//
//      /* set current column (row) pointer to current # nonzeros */
//      Cp[j] = nz;
//
//      /* reset diagonal entry, in case it's not in A */
//      x[j] = ZERO;
//
//      /* iterate down column (along row) of A, collecting nonzeros */
//      for (p=Ap[j]; p<Ap[j+1]; p++) {
//        x[Ai[p]] = c*Ax[p];    /* collect/scale value */
//      }
//
//      /* add identity to this column (row) */
//      if (j < M) {
//        x[j] += ONE;   /* update value */
//      }
//
//      /* fill entries of C with this column's (row's) data */
//      /* fill entries before diagonal */
//      for (p=Ap[j]; p<Ap[j+1] && Ai[p]<j; p++) {
//        Ci[nz] = Ai[p];
//        Cx[nz++] = x[Ai[p]];
//      }
//      /* fill diagonal if applicable */
//      if (p >= Ap[j+1] /* empty or insert at end */ ||  Ai[p] != j /* insert before end */) {
//        Ci[nz] = j;
//        Cx[nz++] = x[j];
//      }
//      /* fill entries past diagonal */
//      for (; p<Ap[j+1]; p++) {
//        Ci[nz] = Ai[p];
//        Cx[nz++] = x[Ai[p]];
//      }
//    }
//
//    /* indicate end of data */
//    Cp[N] = nz;
//
//    /* update A's structure with C's values; nullify C's pointers */
//    SM_NNZ_S(A) = SM_NNZ_S(C);
//
//    if (SM_DATA_S(A))
//      free(SM_DATA_S(A));
//    SM_DATA_S(A) = SM_DATA_S(C);
//    SM_DATA_S(C) = NULL;
//
//    if (SM_INDEXVALS_S(A))
//      free(SM_INDEXVALS_S(A));
//    SM_INDEXVALS_S(A) = SM_INDEXVALS_S(C);
//    SM_INDEXVALS_S(C) = NULL;
//
//    if (SM_INDEXPTRS_S(A))
//      free(SM_INDEXPTRS_S(A));
//    SM_INDEXPTRS_S(A) = SM_INDEXPTRS_S(C);
//    SM_INDEXPTRS_S(C) = NULL;
//
//    /* clean up */
//    SUNMatDestroy_EigenSparse(C);
//    free(x);
//
//  }
//  return SUNMAT_SUCCESS;
//
//}

//int SUNMatScaleAdd_EigenSparse(realtype c, SUNMatrix A, SUNMatrix B)
//{
//  sunindextype j, i, p, nz, newvals, M, N, cend;
//  booleantype newmat;
//  sunindextype *w, *Ap, *Ai, *Bp, *Bi, *Cp, *Ci;
//  realtype *x, *Ax, *Bx, *Cx;
//  SUNMatrix C;
//
//  /* Verify that A and B are compatible */
//  if (!SMCompatible_EigenSparse(A, B))
//    return SUNMAT_ILL_INPUT;
//
//  /* store shortcuts to matrix dimensions (M is inner dimension, N is outer) */
//  if (SM_SPARSETYPE_S(A) == CSC_MAT) {
//    M = SM_ROWS_S(A);
//    N = SM_COLUMNS_S(A);
//  }
//  else {
//    M = SM_COLUMNS_S(A);
//    N = SM_ROWS_S(A);
//  }
//
//  /* access data arrays from A and B (return if failure) */
//  Ap = Ai = Bp = Bi = NULL;
//  Ax = Bx = NULL;
//  if (SM_INDEXPTRS_S(A))  Ap = SM_INDEXPTRS_S(A);
//  else  return(SUNMAT_MEM_FAIL);
//  if (SM_INDEXVALS_S(A))  Ai = SM_INDEXVALS_S(A);
//  else  return(SUNMAT_MEM_FAIL);
//  if (SM_DATA_S(A))       Ax = SM_DATA_S(A);
//  else  return(SUNMAT_MEM_FAIL);
//  if (SM_INDEXPTRS_S(B))  Bp = SM_INDEXPTRS_S(B);
//  else  return(SUNMAT_MEM_FAIL);
//  if (SM_INDEXVALS_S(B))  Bi = SM_INDEXVALS_S(B);
//  else  return(SUNMAT_MEM_FAIL);
//  if (SM_DATA_S(B))       Bx = SM_DATA_S(B);
//  else  return(SUNMAT_MEM_FAIL);
//
//  /* create work arrays for row indices and nonzero column values */
//  w = (sunindextype *) malloc(M * sizeof(sunindextype));
//  x = (realtype *) malloc(M * sizeof(realtype));
//
//  /* determine if A already contains the sparsity pattern of B */
//  newvals = 0;
//  for (j=0; j<N; j++) {
//
//    /* clear work array */
//    for (i=0; i<M; i++)  w[i] = 0;
//
//    /* scan column of A, incrementing w by one */
//    for (i=Ap[j]; i<Ap[j+1]; i++)
//      w[Ai[i]] += 1;
//
//    /* scan column of B, decrementing w by one */
//    for (i=Bp[j]; i<Bp[j+1]; i++)
//      w[Bi[i]] -= 1;
//
//    /* if any entry of w is negative, A doesn't contain B's sparsity,
//       so increment necessary storage counter */
//    for (i=0; i<M; i++)
//      if (w[i] < 0)  newvals += 1;
//  }
//
//  /* If extra nonzeros required, check whether A has sufficient storage space
//     for new nonzero entries (so B can be inserted into existing storage) */
//  newmat = SUNFALSE;   /* no reallocation needed */
//  if (newvals > (SM_NNZ_S(A) - Ap[N]))
//    newmat = SUNTRUE;
//
//  /* perform operation based on existing/necessary structure */
//
//  /*   case 1: A already contains sparsity pattern of B */
//  if (newvals == 0) {
//
//    /* iterate through columns, adding matrices */
//    for (j=0; j<N; j++) {
//
//      /* clear work array */
//      for (i=0; i<M; i++)
//        x[i] = ZERO;
//
//      /* scan column of B, updating work array */
//      for (i = Bp[j]; i < Bp[j+1]; i++)
//        x[Bi[i]] = Bx[i];
//
//      /* scan column of A, updating array entries appropriately */
//      for (i = Ap[j]; i < Ap[j+1]; i++)
//        Ax[i] = c*Ax[i] + x[Ai[i]];
//
//    }
//
//
//  /*   case 2: A has sufficient storage, but does not already contain B's sparsity */
//  } else if (!newmat) {
//
//    /* determine storage location where last column (row) should end */
//    nz = Ap[N] + newvals;
//
//    /* store pointer past last column (row) from original A,
//       and store updated value in revised A */
//    cend = Ap[N];
//    Ap[N] = nz;
//
//    /* iterate through columns (rows) backwards */
//    for (j=N-1; j>=0; j--) {
//
//      /* clear out temporary arrays for this column (row) */
//      for (i=0; i<M; i++) {
//        w[i] = 0;
//        x[i] = RCONST(0.0);
//      }
//
//      /* iterate down column (row) of A, collecting nonzeros */
//      for (p=Ap[j]; p<cend; p++) {
//        w[Ai[p]] += 1;         /* indicate that row (column) is filled */
//        x[Ai[p]] = c*Ax[p];    /* collect/scale value */
//      }
//
//      /* iterate down column of B, collecting nonzeros */
//      for (p=Bp[j]; p<Bp[j+1]; p++) {
//        w[Bi[p]] += 1;       /* indicate that row is filled */
//        x[Bi[p]] += Bx[p];   /* collect value */
//      }
//
//      /* fill entries of A with this column's (row's) data */
//      for (i=M-1; i>=0; i--) {
//        if ( w[i] > 0 ) {
//          Ai[--nz] = i;
//          Ax[nz] = x[i];
//        }
//      }
//
//      /* store ptr past this col (row) from orig A, update value for new A */
//      cend = Ap[j];
//      Ap[j] = nz;
//
//    }
//
//
//  /*   case 3: A must be reallocated with sufficient storage */
//  } else {
//
//
//    /* create new matrix for sum */
//    C = SUNEigenSparseMatrix(SM_ROWS_S(A), SM_COLUMNS_S(A),
//                        Ap[N] + newvals, SM_SPARSETYPE_S(A), A->sunctx);
//
//    /* access data from CSR structures (return if failure) */
//    Cp = Ci = NULL;
//    Cx = NULL;
//    if (SM_INDEXPTRS_S(C))  Cp = SM_INDEXPTRS_S(C);
//    else  return(SUNMAT_MEM_FAIL);
//    if (SM_INDEXVALS_S(C))  Ci = SM_INDEXVALS_S(C);
//    else  return(SUNMAT_MEM_FAIL);
//    if (SM_DATA_S(C))       Cx = SM_DATA_S(C);
//    else  return(SUNMAT_MEM_FAIL);
//
//    /* initialize total nonzero count */
//    nz = 0;
//
//    /* iterate through columns (rows) */
//    for (j=0; j<N; j++) {
//
//      /* set current column (row) pointer to current # nonzeros */
//      Cp[j] = nz;
//
//      /* clear out temporary arrays for this column (row) */
//      for (i=0; i<M; i++) {
//        w[i] = 0;
//        x[i] = RCONST(0.0);
//      }
//
//      /* iterate down column of A, collecting nonzeros */
//      for (p=Ap[j]; p<Ap[j+1]; p++) {
//        w[Ai[p]] += 1;         /* indicate that row is filled */
//        x[Ai[p]] = c*Ax[p];    /* collect/scale value */
//      }
//
//      /* iterate down column of B, collecting nonzeros */
//      for (p=Bp[j]; p<Bp[j+1]; p++) {
//        w[Bi[p]] += 1;       /* indicate that row is filled */
//        x[Bi[p]] += Bx[p];   /* collect value */
//      }
//
//      /* fill entries of C with this column's data */
//      for (i=0; i<M; i++) {
//        if ( w[i] > 0 ) {
//          Ci[nz] = i;
//          Cx[nz++] = x[i];
//        }
//      }
//    }
//
//    /* indicate end of data */
//    Cp[N] = nz;
//
//    /* update A's structure with C's values; nullify C's pointers */
//    SM_NNZ_S(A) = SM_NNZ_S(C);
//
//    free(SM_DATA_S(A));
//    SM_DATA_S(A) = SM_DATA_S(C);
//    SM_DATA_S(C) = NULL;
//
//    free(SM_INDEXVALS_S(A));
//    SM_INDEXVALS_S(A) = SM_INDEXVALS_S(C);
//    SM_INDEXVALS_S(C) = NULL;
//
//    free(SM_INDEXPTRS_S(A));
//    SM_INDEXPTRS_S(A) = SM_INDEXPTRS_S(C);
//    SM_INDEXPTRS_S(C) = NULL;
//
//    /* clean up */
//    SUNMatDestroy_EigenSparse(C);
//
//  }
//
//  /* clean up */
//  free(w);
//  free(x);
//
//  /* return success */
//  return(0);
//
//}

//int SUNMatMatvec_EigenSparse(SUNMatrix A, N_Vector x, N_Vector y)
//{
//  /* Verify that A, x and y are compatible */
//  if (!SMCompatible2_EigenSparse(A, x, y))
//    return SUNMAT_ILL_INPUT;
//
//  /* Perform operation */
//  if(SM_SPARSETYPE_S(A) == CSC_MAT)
//    return Matvec_EigenSparseCSC(A, x, y);
//  else
//    return Matvec_EigenSparseCSR(A, x, y);
//}

int SUNMatSpace_EigenSparse(SUNMatrix A, long int *lenrw, long int *leniw)
{
  *lenrw = SM_NNZ_S(A);
  *leniw = 10 /*+ SM_NP_S(A) + SM_NNZ_S(A)*/; // TODO
  return SUNMAT_SUCCESS;
}


/*
 * =================================================================
 * private functions
 * =================================================================
 */

/* -----------------------------------------------------------------
 * Function to check compatibility of two sparse SUNMatrix objects
 */

static booleantype SMCompatible_EigenSparse(SUNMatrix A, SUNMatrix B)
{
  /* both matrices must be sparse */
  if ( (SUNMatGetID(A) != SUNMATRIX_EIGENSPARSE) ||
       (SUNMatGetID(B) != SUNMATRIX_EIGENSPARSE) )
    return SUNFALSE;

  /* both matrices must have the same shape and sparsity type */
  if (SUNEigenSparseMatrix_Rows(A) != SUNEigenSparseMatrix_Rows(B))
    return SUNFALSE;
  if (SUNEigenSparseMatrix_Columns(A) != SUNEigenSparseMatrix_Columns(B))
    return SUNFALSE;
  if (SUNEigenSparseMatrix_EigenSparseType(A) != SUNEigenSparseMatrix_EigenSparseType(B))
    return SUNFALSE;

  return SUNTRUE;
}


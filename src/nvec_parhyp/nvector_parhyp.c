/*
 * -----------------------------------------------------------------
 * $Revision: 4396 $
 * $Date: 2015-02-26 16:59:39 -0800 (Thu, 26 Feb 2015) $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                and Aaron Collier @ LLNL
 *                Jean M. Sexton
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for a parhyp MPI implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include <nvector/nvector_parhyp.h>
#include <sundials/sundials_math.h>
/* hypre header files */
#include <seq_mv.h>
#include <_hypre_parcsr_mv.h>
#include <HYPRE_parcsr_mv.h>
#include <_hypre_parcsr_ls.h>

/*
#include "_hypre_utilities.h"
#include <HYPRE_config.h>
*/
#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/* Error Message */

#define BAD_N1 "N_VNew_ParHyp -- Sum of local vector lengths differs from "
#define BAD_N2 "input global length. \n\n"
#define BAD_N   BAD_N1 BAD_N2

/* Private function prototypes */

/* Reduction operations add/max/min over the processor group */
static realtype VAllReduce_ParHyp(realtype d, int op, MPI_Comm comm);
/* z=x */
static void VCopy_ParHyp(N_Vector x, N_Vector z);
/* z=x+y */
static void VSum_ParHyp(N_Vector x, N_Vector y, N_Vector z);
/* z=x-y */
static void VDiff_ParHyp(N_Vector x, N_Vector y, N_Vector z);
/* z=-x */
static void VNeg_ParHyp(N_Vector x, N_Vector z);
/* z=c(x+y) */
static void VScaleSum_ParHyp(realtype c, N_Vector x, N_Vector y, N_Vector z);
/* z=c(x-y) */
static void VScaleDiff_ParHyp(realtype c, N_Vector x, N_Vector y, N_Vector z); 
/* z=ax+y */
static void VLin1_ParHyp(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* z=ax-y */
static void VLin2_ParHyp(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* y <- ax+y */
static void Vaxpy_ParHyp(realtype a, N_Vector x, N_Vector y);
/* x <- ax */
static void VScaleBy_ParHyp(realtype a, N_Vector x);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------
 * Function to create a new parhyp vector with empty data array
 */

N_Vector N_VNewEmpty_ParHyp(MPI_Comm comm, 
                              long int local_length,
                              long int global_length)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_ParHyp content;
  long int n, Nsum;
  /* Added variables for hypre_parhyp intialization */
  int nprocs,myid;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myid);
  long int local_lengths[nprocs];
  HYPRE_Int local_jlower[nprocs+1];
  HYPRE_Int partitioning[2];
  int i;

  /* Compute global length as sum of local lengths */
  n = local_length;
  MPI_Allreduce(&n, &Nsum, 1, PVEC_INTEGER_MPI_TYPE, MPI_SUM, comm);
  if (Nsum != global_length) {
    printf(BAD_N);
    return(NULL);
  }
  
  /* Compute local lengths for proper partitioning of the hypre_parhyp vector*/
  local_lengths[myid]=local_length;
  MPI_Allgather(&local_length,1,PVEC_INTEGER_MPI_TYPE, &local_lengths,1,PVEC_INTEGER_MPI_TYPE,comm);
  local_jlower[0]=0;
  for(i=0;i<nprocs; i++)
  {
    local_jlower[i+1]=local_jlower[i]+local_lengths[i];
  }
  partitioning[0]=local_jlower[myid];
  partitioning[1]=local_jlower[myid+1]-1;

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvclone           = N_VClone_ParHyp;
  ops->nvcloneempty      = N_VCloneEmpty_ParHyp;
  ops->nvdestroy         = N_VDestroy_ParHyp;
  ops->nvspace           = N_VSpace_ParHyp;
  ops->nvgetarraypointer = N_VGetArrayPointer_ParHyp;
  ops->nvsetarraypointer = N_VSetArrayPointer_ParHyp;
  ops->nvlinearsum       = N_VLinearSum_ParHyp;
  ops->nvconst           = N_VConst_ParHyp;
  ops->nvprod            = N_VProd_ParHyp;
  ops->nvdiv             = N_VDiv_ParHyp;
  ops->nvscale           = N_VScale_ParHyp;
  ops->nvabs             = N_VAbs_ParHyp;
  ops->nvinv             = N_VInv_ParHyp;
  ops->nvaddconst        = N_VAddConst_ParHyp;
  ops->nvdotprod         = N_VDotProd_ParHyp;
  ops->nvmaxnorm         = N_VMaxNorm_ParHyp;
  ops->nvwrmsnormmask    = N_VWrmsNormMask_ParHyp;
  ops->nvwrmsnorm        = N_VWrmsNorm_ParHyp;
  ops->nvmin             = N_VMin_ParHyp;
  ops->nvwl2norm         = N_VWL2Norm_ParHyp;
  ops->nvl1norm          = N_VL1Norm_ParHyp;
  ops->nvcompare         = N_VCompare_ParHyp;
  ops->nvinvtest         = N_VInvTest_ParHyp;
  ops->nvconstrmask      = N_VConstrMask_ParHyp;
  ops->nvminquotient     = N_VMinQuotient_ParHyp;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_ParHyp) malloc(sizeof(struct _N_VectorContent_ParHyp));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  /* Attach lengths and communicator */
  content->local_length  = local_length;
  content->global_length = global_length;
  content->comm          = comm;
  content->own_data      = FALSE;
  content->own_parvector = TRUE;
  content->data          = NULL;
  
  /* Create ParVector with partitioning info for whole system*/
  content->x             = hypre_ParVectorCreate(comm,global_length, local_jlower);
  /* Set up hypre_ParVector, reserve ParVectorInitialize for when there is data to allocate */
  hypre_VectorData(hypre_ParVectorLocalVector(content->x))=NULL;
  hypre_ParVectorSetPartitioningOwner(content->x,0);

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}

/* ---------------------------------------------------------------- 
 * Function to create a new parhyp vector
 */

N_Vector N_VNew_ParHyp(MPI_Comm comm, 
                         long int local_length,
                         long int global_length)
{
  N_Vector v;
  realtype *data;

  v = NULL;
  v = N_VNewEmpty_ParHyp(comm, local_length, global_length);
  if (v == NULL) return(NULL);

  /* Create data */
  if(local_length > 0) {

    /* Allocate memory */
    data = NULL;
    data = (realtype *) malloc(local_length * sizeof(realtype));
    if(data == NULL) { N_VDestroy_ParHyp(v); return(NULL); }

    /* Attach data */
    NV_OWN_DATA_PH(v) = TRUE;
    NV_DATA_PH(v)     = data;

  }

  return(v);
}

/* ---------------------------------------------------------------- 
 * Function to create a parhyp N_Vector with user data component 
 */

N_Vector N_VMake_ParHyp(MPI_Comm comm, 
                          long int local_length,
                          long int global_length,
                          realtype *v_data)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_ParHyp(comm, local_length, global_length);
  if (v == NULL) return(NULL);

  if (local_length > 0) {
    /* Attach data */
    NV_OWN_DATA_PH(v) = FALSE;
    NV_DATA_PH(v)     = v_data;
  }

  return(v);
}

/* ---------------------------------------------------------------- 
 * Function to create an array of new parhyp vectors. 
 */

N_Vector *N_VCloneVectorArray_ParHyp(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VClone_ParHyp(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_ParHyp(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ---------------------------------------------------------------- 
 * Function to create an array of new parhyp vectors with empty
 * (NULL) data array.
 */

N_Vector *N_VCloneVectorArrayEmpty_ParHyp(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VCloneEmpty_ParHyp(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_ParHyp(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------
 * Function to free an array created with N_VCloneVectorArray_ParHyp
 */

void N_VDestroyVectorArray_ParHyp(N_Vector *vs, int count)
{
  int j;

  for (j = 0; j < count; j++) N_VDestroy_ParHyp(vs[j]);

  free(vs); vs = NULL;

  return;
}

/* ---------------------------------------------------------------- 
 * Function to print a parhyp vector 
 */

void N_VPrint_ParHyp(N_Vector x)
{
  long int i, N;
  realtype *xd;

  xd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);

  for (i = 0; i < N; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%Lg\n", xd[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%g\n", xd[i]);
#else
    printf("%g\n", xd[i]);
#endif
  }
  printf("\n");

  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_ParHyp(N_Vector w)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_ParHyp content;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  /* Added variables for hypre_parhyp intialization */
  int nprocs,myid;
  MPI_Comm_size(NV_COMM_PH(w), &nprocs);
  MPI_Comm_rank(NV_COMM_PH(w), &myid);
  
  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }
  
  ops->nvclone           = w->ops->nvclone;
  ops->nvcloneempty      = w->ops->nvcloneempty;
  ops->nvdestroy         = w->ops->nvdestroy;
  ops->nvspace           = w->ops->nvspace;
  ops->nvgetarraypointer = w->ops->nvgetarraypointer;
  ops->nvsetarraypointer = w->ops->nvsetarraypointer;
  ops->nvlinearsum       = w->ops->nvlinearsum;
  ops->nvconst           = w->ops->nvconst;  
  ops->nvprod            = w->ops->nvprod;   
  ops->nvdiv             = w->ops->nvdiv;
  ops->nvscale           = w->ops->nvscale; 
  ops->nvabs             = w->ops->nvabs;
  ops->nvinv             = w->ops->nvinv;
  ops->nvaddconst        = w->ops->nvaddconst;
  ops->nvdotprod         = w->ops->nvdotprod;
  ops->nvmaxnorm         = w->ops->nvmaxnorm;
  ops->nvwrmsnormmask    = w->ops->nvwrmsnormmask;
  ops->nvwrmsnorm        = w->ops->nvwrmsnorm;
  ops->nvmin             = w->ops->nvmin;
  ops->nvwl2norm         = w->ops->nvwl2norm;
  ops->nvl1norm          = w->ops->nvl1norm;
  ops->nvcompare         = w->ops->nvcompare;    
  ops->nvinvtest         = w->ops->nvinvtest;
  ops->nvconstrmask      = w->ops->nvconstrmask;
  ops->nvminquotient     = w->ops->nvminquotient;

  /* Create content */  
  content = NULL;
  content = (N_VectorContent_ParHyp) malloc(sizeof(struct _N_VectorContent_ParHyp));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  /* Attach lengths and communicator */
  content->local_length  = NV_LOCLENGTH_PH(w);
  content->global_length = NV_GLOBLENGTH_PH(w);
  content->comm          = NV_COMM_PH(w);
  content->own_data      = FALSE;
  content->own_parvector = TRUE;
  content->data          = NULL;
  /* Using CloneShallow to get the proper partitioning and data sizes:
   * "returns a complete copy of a hypre_ParVector x - a shallow copy, re-using
   * the partitioning and data arrays of x" */
  content->x             = hypre_ParVectorCloneShallow(NV_HYPRE_PARVEC_PH(w));
  hypre_ParVectorSetPartitioningOwner(content->x,0);
  /* Changes the data array pointer in the local vector to NULL*/
  hypre_VectorData(hypre_ParVectorLocalVector(content->x))=NULL;
  
  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}

N_Vector N_VClone_ParHyp(N_Vector w)
{
  N_Vector v;
  realtype *data;
  long int local_length;

  v = NULL;
  v = N_VCloneEmpty_ParHyp(w);
  if (v == NULL) return(NULL);

  local_length  = NV_LOCLENGTH_PH(w);

  /* Create data */
  if(local_length > 0) {

    /* Allocate memory */
    data = NULL;
    data = (realtype *) malloc(local_length * sizeof(realtype));
    if(data == NULL) { N_VDestroy_ParHyp(v); return(NULL); }

    /* Attach data */
    NV_OWN_DATA_PH(v) = TRUE;
    hypre_ParVectorLocalVector(NV_HYPRE_PARVEC_PH(v))=hypre_SeqVectorCloneDeep(hypre_ParVectorLocalVector(NV_HYPRE_PARVEC_PH(w)));
  }

  return(v);
}

void N_VDestroy_ParHyp(N_Vector v)
{
  if ((NV_OWN_PARVEC_PH(v) == TRUE)) {
    hypre_ParVectorDestroy(NV_HYPRE_PARVEC_PH(v));
  }
  free(v->content); v->content = NULL;
  free(v->ops); v->ops = NULL;
  free(v); v = NULL;

  return;
}


void N_VSpace_ParHyp(N_Vector v, long int *lrw, long int *liw)
{
  MPI_Comm comm;
  int npes;

  comm = NV_COMM_PH(v);
  MPI_Comm_size(comm, &npes);
  
  *lrw = NV_GLOBLENGTH_PH(v);
  *liw = 2*npes;

  return;
}

realtype *N_VGetArrayPointer_ParHyp(N_Vector v)
{
  return((realtype *) NV_DATA_PH(v));
}

void N_VSetArrayPointer_ParHyp(realtype *v_data, N_Vector v)
{
  if (NV_LOCLENGTH_PH(v) > 0) NV_DATA_PH(v) = v_data;

  return;
}

void N_VLinearSum_ParHyp(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype c, *xd, *yd, *zd;
  N_Vector v1, v2;
  booleantype test;

  xd = yd = zd = NULL;

  if ((b == ONE) && (z == y)) {    /* BLAS usage: axpy y <- ax+y */
    HYPRE_Complex   alpha=a;
    HYPRE_ParVectorAxpy( alpha, (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(x),
                     (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(y));
    return;
  }

  if ((a == ONE) && (z == x)) {    /* BLAS usage: axpy x <- by+x */
    HYPRE_Complex   beta=b;
    HYPRE_ParVectorAxpy( beta, (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(y),
                     (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(x));
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE)) {
    VSum_ParHyp(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_ParHyp(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_ParHyp(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_ParHyp(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b) {
    VScaleSum_ParHyp(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b) {
    VScaleDiff_ParHyp(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */
  
  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])+(b*yd[i]);

  return;
}

void N_VConst_ParHyp(realtype c, N_Vector z)
{
  long int i, N;
  realtype *zd;

HYPRE_Complex    value=c;
HYPRE_ParVectorSetConstantValues( (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(z), value);
  return;
}

/* ----------------------------------------------------------------------------
 * Compute componentwise product z[i] = x[i]*y[i]
 */

void N_VProd_ParHyp(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

// \   num_threads(NV_NUM_THREADS_OMP(x))
//#pragma omp parallel for default(none) private(i) shared(N,xd,yd,zd) schedule(static) 
  for (i = 0; i < N; i++)
    zd[i] = xd[i]*yd[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute componentwise division z[i] = x[i]/y[i]
 */

void N_VDiv_ParHyp(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

// \   num_threads(NV_NUM_THREADS_OMP(x))
//#pragma omp parallel for default(none) private(i) shared(N,xd,yd,zd) schedule(static) 
  for (i = 0; i < N; i++)
    zd[i] = xd[i]/yd[i];

  return;
}

void N_VScale_ParHyp(realtype c, N_Vector x, N_Vector z)
{
  HYPRE_Complex   value=c;
  HYPRE_ParVectorCopy((HYPRE_ParVector) NV_HYPRE_PARVEC_PH(x),  (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(z));
  HYPRE_ParVectorScale(value, (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(z));
  return;
}

void N_VAbs_ParHyp(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = SUNRabs(xd[i]);

  return;
}

void N_VInv_ParHyp(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  zd = NV_DATA_PH(z);

//#pragma omp parallel for default(none) private(i) shared(N,xd,zd) schedule(static)
  for (i = 0; i < N; i++)
    zd[i] = ONE/xd[i];

  return;
}

void N_VAddConst_ParHyp(N_Vector x, realtype b, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  zd = NV_DATA_PH(z);
  
  for (i = 0; i < N; i++) zd[i] = xd[i]+b;

  return;
}

realtype N_VDotProd_ParHyp(N_Vector x, N_Vector y)
{

  HYPRE_Real gsum;
  HYPRE_ParVectorInnerProd( (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(x),
                          (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(y), &gsum);

  return(gsum);
}

realtype N_VMaxNorm_ParHyp(N_Vector x)
{
  long int i, N;
  realtype max, *xd, gmax;
  MPI_Comm comm;

  xd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  comm = NV_COMM_PH(x);

  max = ZERO;

  for (i = 0; i < N; i++) {
    if (SUNRabs(xd[i]) > max) max = SUNRabs(xd[i]);
  }
   
  gmax = VAllReduce_ParHyp(max, 2, comm);

  return(gmax);
}

realtype N_VWrmsNorm_ParHyp(N_Vector x, N_Vector w)
{
  long int i, N, N_global;
  realtype sum, prodi, *xd, *wd, gsum;
  MPI_Comm comm;

  sum = ZERO;
  xd = wd = NULL;

  N        = NV_LOCLENGTH_PH(x);
  N_global = NV_GLOBLENGTH_PH(x);
  xd       = NV_DATA_PH(x);
  wd       = NV_DATA_PH(w);
  comm = NV_COMM_PH(x);

  for (i = 0; i < N; i++) {
    prodi = xd[i]*wd[i];
    sum += SUNSQR(prodi);
  }

  gsum = VAllReduce_ParHyp(sum, 1, comm);

  return(SUNRsqrt(gsum/N_global));
}

realtype N_VWrmsNormMask_ParHyp(N_Vector x, N_Vector w, N_Vector id)
{
  long int i, N, N_global;
  realtype sum, prodi, *xd, *wd, *idd, gsum;
  MPI_Comm comm;

  sum = ZERO;
  xd = wd = idd = NULL;

  N        = NV_LOCLENGTH_PH(x);
  N_global = NV_GLOBLENGTH_PH(x);
  xd       = NV_DATA_PH(x);
  wd       = NV_DATA_PH(w);
  idd      = NV_DATA_PH(id);
  comm = NV_COMM_PH(x);

  for (i = 0; i < N; i++) {
    if (idd[i] > ZERO) {
      prodi = xd[i]*wd[i];
      sum += SUNSQR(prodi);
    }
  }

  gsum = VAllReduce_ParHyp(sum, 1, comm);

  return(SUNRsqrt(gsum/N_global));
}

realtype N_VMin_ParHyp(N_Vector x)
{
  long int i, N;
  realtype min, *xd, gmin;
  MPI_Comm comm;

  xd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  comm = NV_COMM_PH(x);

  min = BIG_REAL;

  if (N > 0) {

    xd = NV_DATA_PH(x);

    min = xd[0];

    for (i = 1; i < N; i++) {
      if (xd[i] < min) min = xd[i];
    }

  }

  gmin = VAllReduce_ParHyp(min, 3, comm);

  return(gmin);
}

realtype N_VWL2Norm_ParHyp(N_Vector x, N_Vector w)
{
  long int i, N;
  realtype sum, prodi, *xd, *wd, gsum;
  MPI_Comm comm;

  sum = ZERO;
  xd = wd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  wd = NV_DATA_PH(w);
  comm = NV_COMM_PH(x);

  for (i = 0; i < N; i++) {
    prodi = xd[i]*wd[i];
    sum += SUNSQR(prodi);
  }

  gsum = VAllReduce_ParHyp(sum, 1, comm);

  return(SUNRsqrt(gsum));
}

realtype N_VL1Norm_ParHyp(N_Vector x)
{
  long int i, N;
  realtype sum, gsum, *xd;
  MPI_Comm comm;

  sum = ZERO;
  xd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  comm = NV_COMM_PH(x);

  for (i = 0; i<N; i++) 
    sum += SUNRabs(xd[i]);

  gsum = VAllReduce_ParHyp(sum, 1, comm);

  return(gsum);
}

void N_VCompare_ParHyp(realtype c, N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++) {
    zd[i] = (SUNRabs(xd[i]) >= c) ? ONE : ZERO;
  }

  return;
}

booleantype N_VInvTest_ParHyp(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd, val, gval;
  MPI_Comm comm;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  zd = NV_DATA_PH(z);
  comm = NV_COMM_PH(x);

  val = ONE;
  for (i = 0; i < N; i++) {
    if (xd[i] == ZERO) 
      val = ZERO;
    else
      zd[i] = ONE/xd[i];
  }

  gval = VAllReduce_ParHyp(val, 3, comm);

  if (gval == ZERO)
    return(FALSE);
  else
    return(TRUE);
}

booleantype N_VConstrMask_ParHyp(N_Vector c, N_Vector x, N_Vector m)
{
  long int i, N;
  realtype temp;
  realtype *cd, *xd, *md;
  MPI_Comm comm;

  cd = xd = md = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  cd = NV_DATA_PH(c);
  md = NV_DATA_PH(m);
  comm = NV_COMM_PH(x);

  temp = ONE;

  for (i = 0; i < N; i++) {
    md[i] = ZERO;
    if (cd[i] == ZERO) continue;
    if (cd[i] > ONEPT5 || cd[i] < -ONEPT5) {
      if (xd[i]*cd[i] <= ZERO) { temp = ZERO; md[i] = ONE; }
      continue;
    }
    if (cd[i] > HALF || cd[i] < -HALF) {
      if (xd[i]*cd[i] < ZERO ) { temp = ZERO; md[i] = ONE; }
    }
  }

  temp = VAllReduce_ParHyp(temp, 3, comm);

  if (temp == ONE) return(TRUE);
  else return(FALSE);
}

realtype N_VMinQuotient_ParHyp(N_Vector num, N_Vector denom)
{
  booleantype notEvenOnce;
  long int i, N;
  realtype *nd, *dd, min;
  MPI_Comm comm;

  nd = dd = NULL;

  N  = NV_LOCLENGTH_PH(num);
  nd = NV_DATA_PH(num);
  dd = NV_DATA_PH(denom);
  comm = NV_COMM_PH(num);

  notEvenOnce = TRUE;
  min = BIG_REAL;

  for (i = 0; i < N; i++) {
    if (dd[i] == ZERO) continue;
    else {
      if (!notEvenOnce) min = SUNMIN(min, nd[i]/dd[i]);
      else {
        min = nd[i]/dd[i];
        notEvenOnce = FALSE;
      }
    }
  }

  return(VAllReduce_ParHyp(min, 3, comm));
}

/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

static realtype VAllReduce_ParHyp(realtype d, int op, MPI_Comm comm)
{
  /* 
   * This function does a global reduction.  The operation is
   *   sum if op = 1,
   *   max if op = 2,
   *   min if op = 3.
   * The operation is over all processors in the communicator 
   */

  realtype out;

  switch (op) {
   case 1: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_SUM, comm);
           break;

   case 2: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MAX, comm);
           break;

   case 3: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MIN, comm);
           break;

   default: break;
  }

  return(out);
}

static void VCopy_ParHyp(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]; 

  return;
}

static void VSum_ParHyp(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]+yd[i];

  return;
}

static void VDiff_ParHyp(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]-yd[i];

  return;
}

static void VNeg_ParHyp(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = -xd[i];

  return;
}

static void VScaleSum_ParHyp(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = c*(xd[i]+yd[i]);

  return;
}

static void VScaleDiff_ParHyp(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = c*(xd[i]-yd[i]);

  return;
}

static void VLin1_ParHyp(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])+yd[i];

  return;
}

static void VLin2_ParHyp(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])-yd[i];

  return;
}

static void Vaxpy_ParHyp(realtype a, N_Vector x, N_Vector y)
{
  long int i, N;
  realtype *xd, *yd;

  xd = yd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);

  if (a == ONE) {
    for (i = 0; i < N; i++)
      yd[i] += xd[i];
    return;
  }
  
  if (a == -ONE) {
    for (i = 0; i < N; i++)
      yd[i] -= xd[i];
    return;
  }    
  
  for (i = 0; i < N; i++)
    yd[i] += a*xd[i];

  return;
}

static void VScaleBy_ParHyp(realtype a, N_Vector x)
{
  long int i, N;
  realtype *xd;

  xd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);

  for (i = 0; i < N; i++)
    xd[i] *= a;

  return;
}
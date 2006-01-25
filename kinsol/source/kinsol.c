/*
 * -----------------------------------------------------------------
 * $Revision: 1.47 $
 * $Date: 2006-01-25 23:08:10 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the main KINSol solver.
 * It is independent of the KINSol linear solver in use.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include "kinsol_impl.h"
#include "sundials_math.h"

/*
 * -----------------------------------------------------------------
 * private constants
 * -----------------------------------------------------------------
 */

#define HALF      RCONST(0.5)
#define ZERO      RCONST(0.0)
#define ONE       RCONST(1.0)
#define ONEPT5    RCONST(1.5)
#define TWO       RCONST(2.0)
#define THREE     RCONST(3.0)
#define FIVE      RCONST(5.0)
#define POINT1    RCONST(0.1)
#define POINT01   RCONST(0.01)
#define POINT99   RCONST(0.99)
#define THOUSAND  RCONST(1000.0)
#define ONETHIRD  RCONST(0.3333333333333333)
#define TWOTHIRDS RCONST(0.6666666666666667)
#define POINT9    RCONST(0.9)
#define POINT0001 RCONST(0.0001)

/* KINStop return value requesting more iterations */

#define RETRY_ITERATION     -998
#define CONTINUE_ITERATIONS -999

/*
 * -----------------------------------------------------------------
 * keys for KINPrintInfo
 * -----------------------------------------------------------------
 */

#define PRNT_RETVAL     1
#define PRNT_NNI        2
#define PRNT_TOL        3
#define PRNT_FMAX       4
#define PRNT_PNORM      5
#define PRNT_PNORM1     6
#define PRNT_FNORM      7
#define PRNT_LAM        8
#define PRNT_ALPHA      9
#define PRNT_BETA      10
#define PRNT_ALPHABETA 11
#define PRNT_ADJ       12

/*
 * -----------------------------------------------------------------
 * private helper function prototypes
 * -----------------------------------------------------------------
 */

static booleantype KINCheckNvector(N_Vector tmpl);
static booleantype KINAllocVectors(KINMem kin_mem, N_Vector tmpl);
static int KINSolInit(KINMem kin_mem);
static int KINConstraint(KINMem kin_mem );
static void KINForcingTerm(KINMem kin_mem, realtype fnormp);
static void KINFreeVectors(KINMem kin_mem);

static int  KINFullNewton(KINMem kin_mem, realtype *fnormp, 
                          realtype *f1normp, booleantype *maxStepTaken);
static int  KINLineSearch(KINMem kin_mem, realtype *fnormp, 
                          realtype *f1normp, booleantype *maxSteptaken);

static int  KINLinSolDrv(KINMem kinmem);
static realtype KINScFNorm(KINMem kin_mem, N_Vector v, N_Vector scale);
static realtype KINScSNorm(KINMem kin_mem, N_Vector v, N_Vector u);
static int KINStop(KINMem kin_mem, booleantype maxStepTaken, int globalstratret);
static void KINPrintInfo(KINMem kin_mem, char *funcname, int key,...);

/* loop macro */
#define loop for(;;)

/*
 * -----------------------------------------------------------------
 * user-callable functions
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : KINCreate
 * -----------------------------------------------------------------
 * KINCreate creates an internal memory block for a problem to 
 * be solved by KINSOL. If successful, KINCreate returns a pointer
 * to the problem memory. This pointer should be passed to
 * KINMalloc. If an initialization error occurs, KINCreate prints
 * an error message to standard error and returns NULL. 
 * -----------------------------------------------------------------
 */

void *KINCreate(void)
{
  KINMem kin_mem;
  realtype uround;

  kin_mem = NULL;
  kin_mem = (KINMem) malloc(sizeof(struct KINMemRec));
  if (kin_mem == NULL) {
    fprintf(stderr, MSG_KINMEM_FAIL);
    return(NULL);
  }

  /* set uround (unit roundoff) */

  kin_mem->kin_uround = uround = UNIT_ROUNDOFF;
  
  /* set default values for solver optional inputs */

  kin_mem->kin_func           = NULL;
  kin_mem->kin_f_data         = NULL;
  kin_mem->kin_constraints    = NULL;
  kin_mem->kin_uscale         = NULL;
  kin_mem->kin_fscale         = NULL;
  kin_mem->kin_constraintsSet = FALSE;
  kin_mem->kin_errfp          = stderr;
  kin_mem->kin_infofp         = stdout;
  kin_mem->kin_printfl        = PRINTFL_DEFAULT;
  kin_mem->kin_mxiter         = MXITER_DEFAULT;
  kin_mem->kin_noInitSetup    = FALSE;
  kin_mem->kin_msbset         = MSBSET_DEFAULT;
  kin_mem->kin_noResMon       = FALSE;
  kin_mem->kin_msbset_sub     = MSBSET_SUB_DEFAULT;
  kin_mem->kin_mxnbcf         = MXNBCF_DEFAULT;
  kin_mem->kin_sthrsh         = TWO;
  kin_mem->kin_noMinEps       = FALSE;
  kin_mem->kin_mxnewtstep     = ZERO;
  kin_mem->kin_sqrt_relfunc   = RSqrt(uround);
  kin_mem->kin_scsteptol      = RPowerR(uround,TWOTHIRDS);
  kin_mem->kin_fnormtol       = RPowerR(uround,ONETHIRD);
  kin_mem->kin_etaflag        = KIN_ETACHOICE1;
  kin_mem->kin_eta            = POINT1;  /* default for KIN_ETACONSTANT */
  kin_mem->kin_eta_alpha      = TWO;     /* default for KIN_ETACHOICE2  */
  kin_mem->kin_eta_gamma      = POINT9;  /* default for KIN_ETACHOICE2  */
  kin_mem->kin_MallocDone     = FALSE;
  kin_mem->kin_setupNonNull   = FALSE;
  kin_mem->kin_omega_min      = OMEGA_MIN;
  kin_mem->kin_omega_max      = OMEGA_MAX;

  /* initialize lrw and liw */

  kin_mem->kin_lrw = 17;
  kin_mem->kin_liw = 22;

  return((void *) kin_mem);
}

#define errfp (kin_mem->kin_errfp)
#define liw   (kin_mem->kin_liw)
#define lrw   (kin_mem->kin_lrw)

/*
 * -----------------------------------------------------------------
 * Function : KINMalloc
 * -----------------------------------------------------------------
 * KINMalloc allocates memory for a problem or execution of KINSol. 
 * If memory is successfully allocated, KIN_SUCCESS is returned.
 * Otherwise, an error message is printed and an error flag
 * returned.
 * -----------------------------------------------------------------
 */

int KINMalloc(void *kinmem, KINSysFn func, N_Vector tmpl)
{
  long int liw1, lrw1;
  KINMem kin_mem;
  booleantype allocOK, nvectorOK;
  
  /* check kinmem */

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINM_NO_MEM);
    return(KIN_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (func == NULL) {
    fprintf(errfp, MSG_FUNC_NULL);
    return(KIN_ILL_INPUT);
  }

  /* check if all required vector operations are implemented */

  nvectorOK = KINCheckNvector(tmpl);
  if (!nvectorOK) {
    if (errfp != NULL) fprintf(errfp, MSG_BAD_NVECTOR);
    return(KIN_ILL_INPUT);
  }

  /* set space requirements for one N_Vector */

  if (tmpl->ops->nvspace != NULL) {
    N_VSpace(tmpl, &lrw1, &liw1);
    kin_mem->kin_lrw1 = lrw1;
    kin_mem->kin_liw1 = liw1;
  }
  else {
    kin_mem->kin_lrw1 = 0;
    kin_mem->kin_liw1 = 0;
  }

  /* allocate necessary vectors */

  allocOK = KINAllocVectors(kin_mem, tmpl);
  if (!allocOK) {
    fprintf(errfp, MSG_MEM_FAIL);
    free(kin_mem); kin_mem = NULL;
    return(KIN_MEM_FAIL);
  }

  /* copy the input parameter into KINSol state */

  kin_mem->kin_func = func;

  /* set the linear solver addresses to NULL */

  kin_mem->kin_linit  = NULL;
  kin_mem->kin_lsetup = NULL;
  kin_mem->kin_lsolve = NULL;
  kin_mem->kin_lfree  = NULL;
  kin_mem->kin_lmem   = NULL;
  
  /* problem memory has been successfully allocated */

  kin_mem->kin_MallocDone = TRUE;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * readability constants
 * -----------------------------------------------------------------
 */

#define func             (kin_mem->kin_func)
#define f_data           (kin_mem->kin_f_data)
#define infofp           (kin_mem->kin_infofp)
#define printfl          (kin_mem->kin_printfl)
#define mxiter           (kin_mem->kin_mxiter)
#define noInitSetup      (kin_mem->kin_noInitSetup)
#define noResMon         (kin_mem->kin_noResMon)
#define retry_nni        (kin_mem->kin_retry_nni)
#define update_fnorm_sub (kin_mem->kin_update_fnorm_sub)
#define msbset           (kin_mem->kin_msbset)
#define msbset_sub       (kin_mem->kin_msbset_sub)
#define etaflag          (kin_mem->kin_etaflag)
#define eta              (kin_mem->kin_eta)
#define ealpha           (kin_mem->kin_eta_alpha)
#define egamma           (kin_mem->kin_eta_gamma)
#define noMinEps         (kin_mem->kin_noMinEps)
#define mxnewtstep       (kin_mem->kin_mxnewtstep)
#define mxnbcf           (kin_mem->kin_mxnbcf)
#define relfunc          (kin_mem->kin_sqrt_relfunc)
#define fnormtol         (kin_mem->kin_fnormtol)
#define scsteptol        (kin_mem->kin_scsteptol)
#define constraints      (kin_mem->kin_constraints)

#define uround          (kin_mem->kin_uround)
#define nni             (kin_mem->kin_nni)
#define nfe             (kin_mem->kin_nfe)
#define nbcf            (kin_mem->kin_nbcf)  
#define nbktrk          (kin_mem->kin_nbktrk)
#define ncscmx          (kin_mem->kin_ncscmx)
#define stepl           (kin_mem->kin_stepl)
#define stepmul         (kin_mem->kin_stepmul)
#define sthrsh          (kin_mem->kin_sthrsh)
#define linit           (kin_mem->kin_linit)
#define lsetup          (kin_mem->kin_lsetup)
#define lsolve          (kin_mem->kin_lsolve) 
#define lfree           (kin_mem->kin_lfree)
#define constraintsSet  (kin_mem->kin_constraintsSet) 
#define jacCurrent      (kin_mem->kin_jacCurrent)          
#define nnilset         (kin_mem->kin_nnilset)
#define nnilset_sub     (kin_mem->kin_nnilset_sub)
#define lmem            (kin_mem->kin_lmem)        
#define inexact_ls      (kin_mem->kin_inexact_ls)
#define setupNonNull    (kin_mem->kin_setupNonNull)
#define fval            (kin_mem->kin_fval)      
#define fnorm           (kin_mem->kin_fnorm)
#define f1norm          (kin_mem->kin_f1norm)
#define etaflag         (kin_mem->kin_etaflag)
#define callForcingTerm (kin_mem->kin_callForcingTerm)
#define uu              (kin_mem->kin_uu)
#define uscale          (kin_mem->kin_uscale)
#define fscale          (kin_mem->kin_fscale)
#define globalstrategy  (kin_mem->kin_globalstrategy)     
#define sJpnorm         (kin_mem->kin_sJpnorm)
#define sfdotJp         (kin_mem->kin_sfdotJp)
#define unew            (kin_mem->kin_unew)
#define pp              (kin_mem->kin_pp)
#define vtemp1          (kin_mem->kin_vtemp1)
#define vtemp2          (kin_mem->kin_vtemp2)
#define eps             (kin_mem->kin_eps)
#define res_norm        (kin_mem->kin_res_norm)
#define fnorm_sub       (kin_mem->kin_fnorm_sub)
#define liw1            (kin_mem->kin_liw1)
#define lrw1            (kin_mem->kin_lrw1)

/*
 * -----------------------------------------------------------------
 * Function : KINSol
 * -----------------------------------------------------------------
 * KINSol (main KINSOL driver routine) manages the computational
 * process of computing an approximate solution of the nonlinear
 * system F(uu) = 0. The KINSol routine calls the following
 * subroutines:
 *
 *  KINSolInit  checks if initial guess satisfies user-supplied
 *              constraints and initializes linear solver
 *
 *  KINLinSolDrv  interfaces with linear solver to find a
 *                solution of the system J(uu)*x = b (calculate
 *                Newton step)
 *
 *  KINFullNewton/KINLineSearch  implement the global strategy
 *
 *  KINForcingTerm  computes the forcing term (eta)
 *
 *  KINStop  determines if an approximate solution has been found
 * -----------------------------------------------------------------
 */

int KINSol(void *kinmem, N_Vector u, int strategy,  
           N_Vector u_scale, N_Vector f_scale)
{
  realtype fnormp, f1normp, epsmin;
  KINMem kin_mem;
  int ret, globalstratret;
  booleantype maxStepTaken;

  globalstratret = 0;

  /* initialize epsmin to avoid compiler warning message */

  epsmin = ZERO;

  /* check for kinmem non-NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINSOL_NO_MEM);
    return(KIN_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if(kin_mem->kin_MallocDone == FALSE) {
    if (errfp != NULL) fprintf(errfp, MSG_KINSOL_NO_MALLOC);
    else fprintf(stderr, MSG_KINSOL_NO_MALLOC);
    return(KIN_NO_MALLOC);
  }

  /* load input arguments */

  globalstrategy = strategy;
  uu = u;
  uscale = u_scale;
  fscale = f_scale;

  /* initialize solver */

  ret = KINSolInit(kin_mem);
  if (ret != KIN_SUCCESS) return(ret);
  
  ncscmx = 0;

  /* Note: The following logic allows the choice of whether or not
     to force a call to the linear solver setup upon a given call to
     KINSol */

  if (noInitSetup) sthrsh = ONE;
  else             sthrsh = TWO;

  /* if eps is to be bounded from below, set the bound */

  if (inexact_ls && !noMinEps) epsmin = POINT01 * fnormtol;

  loop{

    retry_nni = FALSE;

    nni++;

    /* calculate the epsilon (stopping criteria for iterative linear solver)
       for this iteration based on eta from the routine KINForcingTerm */

    if (inexact_ls) {
      eps = (eta + uround) * fnorm;
      if(!noMinEps) eps = MAX(epsmin, eps);
    }

    repeat_nni:

    /* call KINLinSolDrv to calculate the (approximate) Newton step, pp */ 

    ret = KINLinSolDrv(kin_mem);
    if (ret != 0) break;

    /* call the appropriate routine to calculate an acceptable step pp */

    if (globalstrategy == KIN_NONE) {

      /* Full Newton Step*/
      globalstratret = KINFullNewton(kin_mem, &fnormp, &f1normp, &maxStepTaken);

    } else if (globalstrategy == KIN_LINESEARCH) {

      /* Line Search */
      globalstratret = KINLineSearch(kin_mem, &fnormp, &f1normp, &maxStepTaken);

      /* if too many beta condition failures, then stop iteration */
      if (nbcf > mxnbcf) {
        ret = KIN_LINESEARCH_BCFAIL;
        break;
      }

    }

    /* evaluate eta by calling the forcing term routine */

    if (callForcingTerm) KINForcingTerm(kin_mem, fnormp);

    fnorm = fnormp;

    /* call KINStop to check if tolerances where met by this iteration */

    ret = KINStop(kin_mem, maxStepTaken, globalstratret); 

    if (ret == RETRY_ITERATION) {
      retry_nni = TRUE;
      goto repeat_nni;
    }

    /* update uu after the iteration */

    N_VScale(ONE, unew, uu);

    f1norm = f1normp;

    /* print the current nni, fnorm, and nfe values if printfl > 0 */

    if (printfl>0) 
      KINPrintInfo(kin_mem, "KINSol", PRNT_NNI, &nni, &nfe, &fnorm);

    if (ret != CONTINUE_ITERATIONS) break; 

    fflush(errfp);
    
  }  /* end of loop; return */

  if (printfl > 0)
    KINPrintInfo(kin_mem, "KINSol", PRNT_RETVAL, &ret);
  
  fflush(infofp);

  return(ret);
}

/*
 * -----------------------------------------------------------------
 * Function : KINFree
 * -----------------------------------------------------------------
 * This routine frees the problem memory allocated by KINMalloc.
 * Such memory includes all the vectors allocated by
 * KINAllocVectors, and the memory lmem for the linear solver
 * (deallocated by a call to lfree).
 * -----------------------------------------------------------------
 */

void KINFree(void **kinmem)
{
  KINMem kin_mem;

  if (*kinmem == NULL) return;

  kin_mem = (KINMem) (*kinmem);
  KINFreeVectors(kin_mem);

  /* call lfree if non-NULL */

  if (lfree != NULL) lfree(kin_mem);

  free(*kinmem);
  *kinmem = NULL;
}

/*
 * -----------------------------------------------------------------
 * private helper functions
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : KINCheckNvector
 * -----------------------------------------------------------------
 * This routine checks if all required vector operations are
 * implemented (excluding those required by KINConstraint). If all
 * necessary operations are present, then KINCheckNvector returns
 * TRUE. Otherwise, FALSE is returned.
 * -----------------------------------------------------------------
 */

static booleantype KINCheckNvector(N_Vector tmpl)
{
  if ((tmpl->ops->nvclone     == NULL) ||
      (tmpl->ops->nvdestroy   == NULL) ||
      (tmpl->ops->nvlinearsum == NULL) ||
      (tmpl->ops->nvprod      == NULL) ||
      (tmpl->ops->nvdiv       == NULL) ||
      (tmpl->ops->nvscale     == NULL) ||
      (tmpl->ops->nvabs       == NULL) ||
      (tmpl->ops->nvinv       == NULL) ||
      (tmpl->ops->nvmaxnorm   == NULL) ||
      (tmpl->ops->nvmin       == NULL) ||
      (tmpl->ops->nvwl2norm   == NULL)) return(FALSE);
  else return(TRUE);
}

/*
 * -----------------------------------------------------------------
 * Function : KINAllocVectors
 * -----------------------------------------------------------------
 * This routine allocates the KINSol vectors. If all memory
 * allocations are successful, KINAllocVectors returns TRUE.
 * Otherwise all allocated memory is freed and KINAllocVectors
 * returns FALSE.
 * -----------------------------------------------------------------
 */

static booleantype KINAllocVectors(KINMem kin_mem, N_Vector tmpl)
{
  /* allocate unew, fval, pp, vtemp1 and vtemp2 */
  
  unew = NULL;
  unew = N_VClone(tmpl);
  if (unew == NULL) return(FALSE);

  fval = NULL;
  fval = N_VClone(tmpl);
  if (fval == NULL) {
    N_VDestroy(unew);
    return(FALSE);
  }

  pp = NULL;
  pp = N_VClone(tmpl);
  if (pp == NULL) {
    N_VDestroy(unew);
    N_VDestroy(fval);
    return(FALSE);
  }

  vtemp1 = NULL;
  vtemp1 = N_VClone(tmpl);
  if (vtemp1 == NULL) {
    N_VDestroy(unew);
    N_VDestroy(fval);
    N_VDestroy(pp);
    return(FALSE);
  }

  vtemp2 = NULL;
  vtemp2 = N_VClone(tmpl);
  if (vtemp2 == NULL) {
    N_VDestroy(unew);
    N_VDestroy(fval);
    N_VDestroy(pp);
    N_VDestroy(vtemp1);
    return(FALSE);
  }

  /* rpdate solver workspace lengths */

  liw += 5*liw1;
  lrw += 5*lrw1;

  return(TRUE);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSolInit
 * -----------------------------------------------------------------
 * KINSolInit initializes the problem for the specific input
 * received in this call to KINSol (which calls KINSolInit). All
 * problem specification inputs are checked for errors. If any error
 * occurs during initialization, it is reported to the file whose
 * file pointer is errfp.
 *
 * The possible return values for KINSolInit are:
 *   KIN_SUCCESS : indicates a normal initialization
 *
 *   KINS_ILL_INPUT : indicates that an input error has been found
 *
 *   KIN_INITIAL_GUESS_OK : indicates that the guess uu
 *                          satisfied the system func(uu) = 0
 *                          within the tolerances specified
 * -----------------------------------------------------------------
 */

static int KINSolInit(KINMem kin_mem)
{
  int ret;
  realtype fmax;
  
  /* check for illegal input parameters */

  if (uu == NULL) {
    fprintf(errfp, MSG_UU_NULL);   
    return(KIN_ILL_INPUT);
  }

  if ((globalstrategy != KIN_NONE) && (globalstrategy != KIN_LINESEARCH)) {
    fprintf(errfp, MSG_BAD_GLSTRAT);
    return(KIN_ILL_INPUT);
  }

  if (uscale == NULL)  {
    fprintf(errfp, MSG_BAD_USCALE);
    return(KIN_ILL_INPUT);
  }

  if (N_VMin(uscale) <= ZERO){
    fprintf(errfp, MSG_USCALE_NONPOSITIVE); 
    return(KIN_ILL_INPUT);
  }

  if (fscale == NULL)  {
    fprintf(errfp, MSG_BAD_FSCALE);
    return(KIN_ILL_INPUT);
  }

  if (N_VMin(fscale) <= ZERO){
    fprintf(errfp, MSG_FSCALE_NONPOSITIVE);
    return(KIN_ILL_INPUT);
  }

  /* set the constraints flag */

  if (constraints == NULL) 
    constraintsSet = FALSE;
  else {
    constraintsSet = TRUE;
    if ((constraints->ops->nvconstrmask  == NULL) ||
	(constraints->ops->nvminquotient == NULL)) {
      if (errfp != NULL) fprintf(errfp, MSG_BAD_NVECTOR);
      return(KIN_ILL_INPUT);
    }
  }

  /* check the initial guess uu against the constraints */

  if (constraintsSet) {
    if (!N_VConstrMask(constraints, uu, vtemp1)) {
      fprintf(errfp, MSG_INITIAL_CNSTRNT);
      return(KIN_ILL_INPUT);
    }
  }
  
  /* all error checking is complete at this point */

  if (printfl > 0)
    KINPrintInfo(kin_mem, "KINSolInit", PRNT_TOL, &scsteptol, &fnormtol);

  /* calculate the default value for mxnewtstep (maximum Newton step) */

  if (mxnewtstep == ZERO)
    mxnewtstep = THOUSAND * N_VWL2Norm(uu, uscale);
  if (mxnewtstep < ONE) mxnewtstep = ONE;


  /* additional set-up for inexact linear solvers */

  if (inexact_ls) {

    /* set up the coefficients for the eta calculation */

    callForcingTerm = (etaflag != KIN_ETACONSTANT);

    /* this value is always used for choice #1 */

    if (etaflag == KIN_ETACHOICE1) ealpha = (ONE + RSqrt(FIVE)) * HALF;

    /* initial value for eta set to 0.5 for other than the KIN_ETACONSTANT option */

    if (etaflag != KIN_ETACONSTANT) eta = HALF;

    /* disable residual monitoring if using an inexact linear solver */

    noResMon = TRUE;

  } else {

    callForcingTerm = FALSE;

  }

  /* initialize counters */

  nfe = nnilset = nnilset_sub = nni = nbcf = nbktrk = 0;

  /* see if the system func(uu) = 0 is satisfied by the initial guess uu */

  func(uu, fval, f_data); nfe++;
  fmax = KINScFNorm(kin_mem, fval, fscale);

  if (printfl > 1)
    KINPrintInfo(kin_mem, "KINSolInit", PRNT_FMAX, &fmax);

  if (fmax <= (POINT01 * fnormtol)) return(KIN_INITIAL_GUESS_OK);

  /* initialize the linear solver if linit != NULL */

  if (linit != NULL) {
    ret = linit(kin_mem);
    if (ret != 0) {
      fprintf(errfp, MSG_LINIT_FAIL);
      return(KIN_LINIT_FAIL);
    }
  }

  /* initialize the L2 (Euclidean) norms of f for the linear iteration steps */

  fnorm = N_VWL2Norm(fval, fscale);
  f1norm = HALF * fnorm * fnorm;

  fnorm_sub = fnorm;

  if (printfl > 0)
    KINPrintInfo(kin_mem, "KINSolInit", PRNT_NNI, &nni, &nfe, &fnorm);

  /* problem has now been successfully initialized */

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINConstraint
 * -----------------------------------------------------------------
 * This routine checks if the proposed solution vector uu + pp
 * violates any constraints. If a constraint is violated, then the
 * scalar stepmul is determined such that uu + stepmul * pp does
 * not violate any constraints.
 *
 * Note: This routine is called by the functions
 *       KINLineSearch and KINFullNewton.
 * -----------------------------------------------------------------
 */

static int KINConstraint(KINMem kin_mem)
{
  N_VLinearSum(ONE, uu, ONE, pp, vtemp1);

  /* if vtemp1[i] violates constraint[i] then vtemp2[i] = 1
     else vtemp2[i] = 0 (vtemp2 is the mask vector) */

  if(N_VConstrMask(constraints, vtemp1, vtemp2)) return(0);

  /* vtemp1[i] = ABS(pp[i]) */

  N_VAbs(pp, vtemp1);

  /* consider vtemp1[i] only if vtemp2[i] = 1 (constraint violated) */

  N_VProd(vtemp2, vtemp1, vtemp1);

  N_VAbs(uu, vtemp2);
  stepmul = POINT9 * N_VMinQuotient(vtemp2, vtemp1);

  return(1);
}

/*
 * -----------------------------------------------------------------
 * Function : KINForcingTerm
 * -----------------------------------------------------------------
 * This routine computes eta, the scaling factor in the linear
 * convergence stopping tolerance eps when choice #1 or choice #2
 * forcing terms are used. Eta is computed here for all but the
 * first iterative step, which is set to the default in routine
 * KINSolInit.
 *
 * This routine was written by Homer Walker of Utah State
 * University with subsequent modifications by Allan Taylor @ LLNL.
 *
 * It is based on the concepts of the paper 'Choosing the forcing
 * terms in an inexact Newton method', SIAM J Sci Comput, 17
 * (1996), pp 16 - 32, or Utah State University Research Report
 * 6/94/75 of the same title.
 * -----------------------------------------------------------------
 */

static void KINForcingTerm(KINMem kin_mem, realtype fnormp)
{
  realtype eta_max, eta_min, eta_safe, linmodel_norm;

  eta_max  = POINT9;
  eta_min  = POINT0001;
  eta_safe = HALF;

  /* choice #1 forcing term */

  if (etaflag == KIN_ETACHOICE1) {

    /* compute the norm of f + Jp , scaled L2 norm */

    linmodel_norm = RSqrt((fnorm * fnorm) + (TWO * sfdotJp) + (sJpnorm * sJpnorm));

    /* form the safeguarded for choice #1 */ 

    eta_safe = RPowerR(eta, ealpha); 
    eta = ABS(fnormp - linmodel_norm) / fnorm; 
  }

  /* choice #2 forcing term */

  if (etaflag == KIN_ETACHOICE2) {
    eta_safe = egamma * RPowerR(eta, ealpha); 
    eta = egamma * RPowerR((fnormp / fnorm), ealpha); 
  }

  /* apply safeguards */
 
  if(eta_safe < POINT1) eta_safe = ZERO;
  eta = MAX(eta, eta_safe); 
  eta = MAX(eta, eta_min); 
  eta = MIN(eta, eta_max); 

  return; 
}

/*
 * -----------------------------------------------------------------
 * Function : KINFreeVectors
 * -----------------------------------------------------------------
 * This routine frees the KINSol vectors allocated by
 * KINAllocVectors.
 * -----------------------------------------------------------------
 */

static void KINFreeVectors(KINMem kin_mem)
{
  if (unew != NULL)   N_VDestroy(unew);
  if (fval != NULL)   N_VDestroy(fval);
  if (pp != NULL)     N_VDestroy(pp);
  if (vtemp1 != NULL) N_VDestroy(vtemp1);
  if (vtemp2 != NULL) N_VDestroy(vtemp2);

  lrw -= 5*lrw1;
  liw -= 5*liw1;

  if (kin_mem->kin_constraintsSet) {
    if (constraints != NULL) N_VDestroy(constraints);
    lrw -= lrw1;
    liw -= liw1;
  }

  return;
}

/*
 * -----------------------------------------------------------------
 * Function : KINFullNewton
 * -----------------------------------------------------------------
 * This routine is the main driver for the Full Newton
 * algorithm. Its purpose is to compute unew = uu + pp in the
 * direction pp from uu, taking the full Newton step. The
 * step may be constrained if the constraint conditions are
 * violated, or if the norm of pp is greater than mxnewtstep. 
 * -----------------------------------------------------------------
 */

static int KINFullNewton(KINMem kin_mem, realtype *fnormp, realtype *f1normp,
                         booleantype *maxStepTaken)
{
  realtype pnorm, ratio;
  int ret;

  *maxStepTaken = FALSE;
  pnorm = N_VWL2Norm(pp, uscale);
  ratio = ONE;
  if (pnorm > mxnewtstep) {
    ratio = mxnewtstep / pnorm;
    N_VScale(ratio, pp, pp);
    pnorm = mxnewtstep;
  }

  if (printfl > 0)
    KINPrintInfo(kin_mem, "KINFullNewton", PRNT_PNORM, &pnorm);

  /* if constraints are active, then constrain the step accordingly */

  stepl = pnorm;
  stepmul = ONE;
  if (constraintsSet) {
    ret = KINConstraint(kin_mem);  /* Note: This routine changes the step pp */
    if (ret == 1) {
      ratio *= stepmul;
      N_VScale(stepmul, pp, pp);
      pnorm *= stepmul;
      stepl = pnorm;
      if (printfl > 0)
        KINPrintInfo(kin_mem, "KINFullNewton", PRNT_PNORM, &pnorm);
      if (pnorm <= scsteptol) return(1);
    }
  }
 
  /* scale these two expressions by ratio for later use in KINForcingTerm */

  sfdotJp *= ratio;
  sJpnorm *= ratio;
 
  /* compute the iterate unew = uu + pp */

  N_VLinearSum(ONE, uu, ONE, pp, unew);

  /* evaluate func(unew) and its norm, and return */
 
  func(unew, fval, f_data); nfe++;
  *fnormp = N_VWL2Norm(fval,fscale);
  *f1normp = HALF * (*fnormp) * (*fnormp);

  if (printfl > 1) 
    KINPrintInfo(kin_mem, "KINFullNewton", PRNT_FNORM, fnormp);

  if (pnorm > (POINT99 * mxnewtstep)) *maxStepTaken = TRUE; 

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : KINLineSearch
 * -----------------------------------------------------------------
 * The routine KINLineSearch implements the LineSearch algorithm.
 * Its purpose is to find unew = uu + rl * pp in the direction pp
 * from uu so that:
 *                                    t
 *  func(unew) <= func(uu) + alpha * g  (unew - uu) (alpha = 1.e-4)
 *
 *    and
 *                                   t
 *  func(unew) >= func(uu) + beta * g  (unew - uu) (beta = 0.9)
 *
 * where 0 < rlmin <= rl <= rlmax.
 *
 * Note:
 *             mxnewtstep
 *  rlmax = ----------------   if uu+pp does not violateany constraints
 *          ||uscale*pp||_L2
 *
 *  rlmax = 1   otherwise
 *
 *    and
 *
 *                 scsteptol
 *  rlmin = --------------------------
 *          ||           pp         ||
 *          || -------------------- ||_L-infinity
 *          || (1/uscale + ABS(uu)) ||
 *
 * -----------------------------------------------------------------
 */

static int KINLineSearch(KINMem kin_mem, realtype *fnormp, realtype *f1normp,
                         booleantype *maxStepTaken)
{
  realtype pnorm, ratio, slpi, rlmin, rlength, rl, rlmax, rldiff;
  realtype rltmp, rlprev, pt1trl, f1nprv, rllo, rlinc, alpha, beta;
  realtype alpha_cond, beta_cond, rl_a, tmp1, rl_b, tmp2, disc;
  long int nfesave, rladjust;
  int ret, ivio;
  booleantype useCubicBktr;

  rladjust = 0;
  *maxStepTaken = FALSE;
  ratio = ONE;
  alpha = POINT0001;
  beta = POINT9;
  useCubicBktr = FALSE;

  /* initialize rlprev and f1nprv to avoid compiler warning message */

  rlprev = f1nprv = ZERO;

  pnorm = N_VWL2Norm(pp, uscale);
  if(pnorm > mxnewtstep ) {
    ratio = mxnewtstep / pnorm;
    N_VScale(ratio, pp, pp);
    pnorm = mxnewtstep;
  }

  ivio = 0;

  /* check if constraints are active and if so constrain the step by the constraints */

  stepl = pnorm;
  stepmul = ONE;

  if(constraintsSet){
    ret = KINConstraint(kin_mem);
    if(ret == 1){
      ivio = 1;
      ratio *= stepmul;
      N_VScale(stepmul, pp, pp);
      pnorm *= stepmul;
      stepl = pnorm;
      if (printfl > 0)
        KINPrintInfo(kin_mem, "KINLineSearch", PRNT_PNORM1, &pnorm);
    }
  }

  slpi = sfdotJp * ratio;
  rlength = KINScSNorm(kin_mem, pp, uu);
  rlmin = scsteptol / rlength;
  rl = ONE;

  if (printfl > 2)
    KINPrintInfo(kin_mem, "KINLineSearch", PRNT_LAM, &rlmin, &f1norm, &pnorm);

  /* now begin the iteration to find an rl value which satisfies both the
     alpha and beta conditions */

  /* if rl < rlmin, then terminate and return 1 */

  nfesave = nfe;

  loop {

    N_VLinearSum(ONE, uu, rl, pp, unew);
    func(unew, fval, f_data); nfe++;
    *fnormp = N_VWL2Norm(fval, fscale);
    *f1normp = HALF * (*fnormp) * (*fnormp) ;
    alpha_cond = f1norm + (alpha * slpi * rl);

    if (printfl > 2 && rladjust > 0)
      KINPrintInfo(kin_mem, "KINLinesearch", PRNT_ALPHA, fnormp, f1normp, &alpha_cond, &rl);

    if ((*f1normp) <= alpha_cond) break;

    /* alpha condition not satisfied so perform backtracking to compute a new rl value */

    if (rl < rlmin) {

      /* no satisfactory unew can be found sufficiently distinct from uu
         so copy uu into unew and return */

      /* step remains unchanged */

      N_VScale(ONE, uu, unew);
      func(unew, fval, f_data); nfe++;
      *fnormp = N_VWL2Norm(fval, fscale);
      *f1normp = HALF * (*fnormp) * (*fnormp);
      nbktrk = nfe - nfesave - 1;
      return(1);
    }

    /* use cubic fit for all remaining backtracks */

    if (useCubicBktr) {
      tmp1 = (*f1normp) - f1norm - (rl * slpi);
      tmp2 = f1nprv - f1norm - (rlprev * slpi);
      rl_a = ((ONE / (rl * rl)) * tmp1) - ((ONE / (rlprev * rlprev)) * tmp2);
      rl_b = ((-rlprev / (rl * rl)) * tmp1) + ((rl / (rlprev * rlprev)) * tmp2);
      tmp1 = ONE / (rl - rlprev);
      rl_a *= tmp1;
      rl_b *= tmp1;
      disc = (rl_b * rl_b) - (THREE * rl_a * slpi);

      /* cubic is actually just a quadratic (rl_a ~ 0) */

      if (ABS(rl_a) < uround) rltmp = -slpi / (TWO * rl_b);

      /* real cubic */

      else rltmp = (-rl_b + RSqrt(disc)) / (THREE * rl_a);

      if (rltmp > (HALF * rl)) rltmp = HALF * rl;
    }

    /* use quadratic fit only for initial backtrack */

    else if (!useCubicBktr) {
      rltmp = -slpi / (TWO * ((*f1normp) - f1norm - slpi));
      useCubicBktr = TRUE;
    }

    rlprev = rl;
    f1nprv = (*f1normp);
    pt1trl = POINT1 * rl;
    rl = MAX(pt1trl, rltmp);
    rladjust++;
  }

  /* alpha condition is satisfied so now check the beta condition */

  beta_cond = f1norm + (beta * slpi * rl);
  if ((*f1normp) < beta_cond) {

    if ((rl == ONE) && (pnorm < mxnewtstep)) {
      rlmax = mxnewtstep / pnorm;
      if (ivio == 1) rlmax = ONE;
      do {
        rlprev = rl;
        f1nprv = *f1normp;
        rl = MIN((TWO * rl), rlmax);
        rladjust++;
        N_VLinearSum(ONE, uu, rl, pp, unew);
        func(unew, fval, f_data); nfe++;
        *fnormp = N_VWL2Norm(fval, fscale);
        *f1normp = HALF * (*fnormp) * (*fnormp);
        alpha_cond = f1norm + (alpha * slpi * rl);
        beta_cond = f1norm + (beta * slpi * rl);
        if (printfl > 2)
          KINPrintInfo(kin_mem, "KINLineSearch", PRNT_BETA, f1normp, &beta_cond, &rl);
      } while (((*f1normp) <= alpha_cond) && 
	       ((*f1normp) < beta_cond) && (rl < rlmax));
    }
    if ((rl < ONE) || ((rl > ONE) && (*f1normp > alpha_cond))) {
      rllo = MIN(rl, rlprev);
      rldiff = ABS(rlprev - rl);

      do {
        rlinc = HALF * rldiff;
        rl = rllo + rlinc;
        rladjust++;
        N_VLinearSum(ONE, uu, rl, pp, unew);
        func(unew, fval, f_data); nfe++;
        *fnormp = N_VWL2Norm(fval, fscale);
        *f1normp = HALF * (*fnormp) * (*fnormp);
        alpha_cond = f1norm + (alpha * slpi * rl);
        beta_cond = f1norm + (beta * slpi * rl);
        if ((printfl > 2) && (rladjust > 0))
          KINPrintInfo(kin_mem, "KINLineSearch", PRNT_ALPHABETA, 
                       f1normp, &alpha_cond, &beta_cond, &rl);
        if ((*f1normp) > alpha_cond) rldiff = rlinc;
        else if (*f1normp < beta_cond) {
          rllo = rl;
          rldiff = rldiff - rlinc;
        }

      } while ((*f1normp > alpha_cond) ||
	       ((*f1normp < beta_cond) && (rldiff >= rlmin)));

      if ((*f1normp) < beta_cond) {

	/* beta condition could not be satisfied so set unew to last u value
	   that satisfied the alpha condition and continue */

	/* increment counter on number of steps beta condition not satisfied */

        N_VLinearSum(ONE, uu, rllo, pp, unew);
        func(unew, fval, f_data); nfe++;
        *fnormp = N_VWL2Norm(fval, fscale);
        *f1normp = HALF * (*fnormp) * (*fnormp);   
        nbcf++;
      }
    }  /* end of if (rl < ONE) block */
  }  /* end of (f1normp < beta_cond) loop */

  nbktrk= nfe - nfesave;
  if ((printfl > 1) && (rladjust > 0))
    KINPrintInfo(kin_mem, "KINLineSearch", PRNT_ADJ, &rladjust);

  /* scale these two expressions by rl * ratio for later use in KINForcingTerm */

  sfdotJp = sfdotJp * rl * ratio;
  sJpnorm = sJpnorm * rl * ratio;

  if ((rl * pnorm) > (POINT99 * mxnewtstep)) *maxStepTaken = TRUE;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : KINLinSolvDrv
 * -----------------------------------------------------------------
 * This routine handles the process of solving for the approximate
 * solution of the Newton equations in the Newton iteration.
 * Subsequent routines handle the nonlinear aspects of its
 * application. 
 * -----------------------------------------------------------------
 */

static int KINLinSolDrv(KINMem kin_mem)
{
  N_Vector x, b;
  int ret;

  if ((nni - nnilset) >= msbset) {
    sthrsh = TWO;
    update_fnorm_sub = TRUE;
  }

  loop{

    jacCurrent = FALSE;

    if ((sthrsh > ONEPT5) && setupNonNull) {
      ret = lsetup(kin_mem);
      jacCurrent = TRUE;
      nnilset = nni;
      nnilset_sub = nni;
      if (ret != 0) return(KIN_LSETUP_FAIL);
    }

    /* rename vectors for readability */

    b = unew;
    x = pp;

    /* load b with the current value of -fval */

    N_VScale(-ONE, fval, b);

    /* call the generic 'lsolve' routine to solve the system Jx = b */

    ret = lsolve(kin_mem, x, b, &res_norm);

    if (ret == 0) return(0);
    else if (ret < 0) return(KIN_LSOLVE_FAIL);
    else if ((!setupNonNull) || (jacCurrent)) return(KIN_LINSOLV_NO_RECOVERY);

    /* loop back only if the linear solver setup is in use and Jacobian information
       is not current */

    sthrsh = TWO;

  }
}

/*
 * -----------------------------------------------------------------
 * Function : KINScFNorm
 * -----------------------------------------------------------------
 * This routine computes the max norm for scaled vectors. The
 * scaling vector is scale, and the vector of which the norm is to
 * be determined is vv. The returned value, fnormval, is the
 * resulting scaled vector norm. This routine uses N_Vector
 * functions from the vector module.
 * -----------------------------------------------------------------
 */

static realtype KINScFNorm(KINMem kin_mem, N_Vector v, N_Vector scale)
{
  N_VProd(scale, v, vtemp1);
  return(N_VMaxNorm(vtemp1));
}

/*
 * -----------------------------------------------------------------
 * Function : KINScSNorm
 * -----------------------------------------------------------------
 * This routine computes the max norm of the scaled steplength, ss.
 * Here ucur is the current step and usc is the u scale factor.
 * -----------------------------------------------------------------
 */

static realtype KINScSNorm(KINMem kin_mem, N_Vector v, N_Vector u)
{
  realtype length;

  N_VInv(uscale, vtemp1);
  N_VAbs(u, vtemp2);
  N_VLinearSum(ONE, vtemp1, ONE, vtemp2, vtemp1);
  N_VDiv(v, vtemp1, vtemp1);

  length = N_VMaxNorm(vtemp1);

  return(length);
}

#define omega_min (kin_mem->kin_omega_min)
#define omega_max (kin_mem->kin_omega_max)

/*
 * -----------------------------------------------------------------
 * Function : KINStop
 * -----------------------------------------------------------------
 * This routine checks the current iterate unew to see if the
 * system func(unew) = 0 is satisfied by a variety of tests.
 * -----------------------------------------------------------------
 */

static int KINStop(KINMem kin_mem, booleantype maxStepTaken, int globalstratret)
{
  realtype fmax, rlength;
  realtype omega;
  N_Vector delta;

  if (globalstratret == 1) {

    if (setupNonNull && !jacCurrent) {

      /* if the globalstratret was caused (potentially) by the 
	 Jacobian info being out of date, then update it */

      sthrsh = TWO;

      return(CONTINUE_ITERATIONS);
    }

    else return((globalstrategy == KIN_NONE) ?
		KIN_STEP_LT_STPTOL : KIN_LINESEARCH_NONCONV);
  }

  /* check tolerance on scaled norm of func at the current iterate */

  fmax = KINScFNorm(kin_mem, fval, fscale);

  if (printfl > 1) 
    KINPrintInfo(kin_mem, "KINStop", PRNT_FMAX, &fmax);

  if (fmax <= fnormtol) return(KIN_SUCCESS);

  /* check if the scaled distance between the last two steps is too small */
  /* NOTE: pp used as work space to store this distance */

  delta = pp;
  N_VLinearSum(ONE, unew, -ONE, uu, delta);
  rlength = KINScSNorm(kin_mem, delta, unew);

  if (rlength <= scsteptol) {

    if (setupNonNull && !jacCurrent) {

      /* for rlength too small and the Jacobian info not current,
         try again with current Jacobian info */

      sthrsh = TWO;

      return(CONTINUE_ITERATIONS);
    }

    else

      /* otherwise return failure flag */

      return(KIN_STEP_LT_STPTOL);
  }

  /* if the maximum number of iterations is reached, return failure flag */

   if (nni >= mxiter) return(KIN_MAXITER_REACHED);

   /* check for consecutive number of steps taken of size mxnewtstep
      and if not maxStepTaken, then set ncscmx to 0 */
 
  if (maxStepTaken) ncscmx++;
  else ncscmx = 0;
 
  if (ncscmx == 5) return(KIN_MXNEWT_5X_EXCEEDED);

  /* load threshold for reevaluating the Jacobian info (lsetup) */

  if (inexact_ls) sthrsh = rlength;

  /* residual monitoring scheme for modified Newton method */

  else if (!inexact_ls && !noResMon) {
    if ((nni - nnilset_sub) >= msbset_sub) {
      nnilset_sub = nni;

      if (omega_min == -ONE) omega = omega_max;
      else omega = MIN(omega_min*exp(MAX(ZERO, (fnorm/fnormtol)-ONE)), omega_max);

      /* check if making satisfactory progress */

      if (fnorm > omega*fnorm_sub) {

	/* update jacobian and retry iteration */

	if (setupNonNull && !jacCurrent) {
	  sthrsh = TWO;
	  return(RETRY_ITERATION);
	}

	/* otherwise, we cannot do anything, so just return */

      }
      else {
	fnorm_sub = fnorm;
	sthrsh = ONE;
      }
    }

    /* must reset sthrsh */

    else {
      if (retry_nni || update_fnorm_sub) fnorm_sub = fnorm;
      if (update_fnorm_sub) update_fnorm_sub = FALSE;
      sthrsh = ONE;
    }
  }

  /* if made it to here, then the iteration process is not finished
     so return CONTINUE_ITERATIONS flag */

  return(CONTINUE_ITERATIONS);
}

/*
 * -----------------------------------------------------------------
 * KINPrintInfo
 * -----------------------------------------------------------------
 */

static void KINPrintInfo(KINMem kin_mem, char *funcname, int key,...)
{
  va_list ap;
  realtype rnum1, rnum2, rnum3, rnum4;
  long int inum1, inum2;
  int ret;

  fprintf(infofp, "---%s\n   ", funcname);

  /* initialize argument processing */

  va_start(ap, key); 

  switch(key) {

  case PRNT_RETVAL:
    ret = *(va_arg(ap, int *));
    fprintf(infofp,"return value: %d ", ret);
    switch(ret) {
    case KIN_SUCCESS:
      fprintf(infofp, "(KIN_SUCCESS)\n");
      break;
    case KIN_STEP_LT_STPTOL:
      fprintf(infofp, "(KIN_STEP_LT_STPTOL)\n");
      break;
    case KIN_LINESEARCH_NONCONV:
      fprintf(infofp, "(KIN_LINESEARCH_NONCONV)\n");
      break;
    case KIN_LINESEARCH_BCFAIL:
      fprintf(infofp, "(KIN_LINESEARCH_BCFAIL)\n");
      break;
    case KIN_MAXITER_REACHED:
      fprintf(infofp, "(KIN_MAXITER_REACHED)\n");
      break;
    case KIN_MXNEWT_5X_EXCEEDED:
      fprintf(infofp, "(KIN_MXNEWT_5X_EXCEEDED)\n");
      break;
    case KIN_LINSOLV_NO_RECOVERY:
      fprintf(infofp, "(KIN_LINSOLV_NO_RECOVERY)\n");
      break;
    case KIN_LSETUP_FAIL:
      fprintf(infofp, "(KIN_PRECONDSET_FAILURE)\n");
      break;
    case KIN_LSOLVE_FAIL:
      fprintf(infofp, "(KIN_PRECONDSOLVE_FAILURE)\n");
      break;
    }
    break;

  case PRNT_NNI:
    inum1 = *(va_arg(ap, long int *));
    inum2 = *(va_arg(ap, long int *));
    rnum1 = *(va_arg(ap, realtype *));
    fprintf(infofp, "nni = %4ld  ", inum1);
    fprintf(infofp, "nfe = %6ld  ", inum2);

#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "fnorm = %26.16Lg\n", rnum1);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(infofp, "fnorm = %26.16lg\n", rnum1);
#else
    fprintf(infofp, "fnorm = %26.16g\n", rnum1);
#endif
    break;

  case PRNT_TOL:
    rnum1 = *(va_arg(ap, realtype *));
    rnum2 = *(va_arg(ap, realtype *));

#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "scsteptol = %12.3Lg  fnormtol = %12.3Lg\n", rnum1, rnum2);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(infofp, "scsteptol = %12.3lg  fnormtol = %12.3lg\n", rnum1, rnum2);
#else
    fprintf(infofp, "scsteptol = %12.3g  fnormtol = %12.3g\n", rnum1, rnum2);
#endif
    break;
    
  case PRNT_FMAX:
    rnum1 = *(va_arg(ap, realtype *));

#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "scaled f norm (for stopping) = %12.3Lg\n", rnum1);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(infofp, "scaled f norm (for stopping) = %12.3lg\n", rnum1);
#else
    fprintf(infofp, "scaled f norm (for stopping) = %12.3g\n", rnum1);
#endif
    break;
    
  case PRNT_PNORM:
    rnum1 = *(va_arg(ap, realtype *));

#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "pnorm = %12.4Le\n", rnum1);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(infofp, "pnorm = %12.4le\n", rnum1);
#else
    fprintf(infofp, "pnorm = %12.4e\n", rnum1);
#endif
    break;

  case PRNT_PNORM1:
    rnum1 = *(va_arg(ap, realtype *));

#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "(ivio=1) pnorm = %12.4Le\n", rnum1);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(infofp, "(ivio=1) pnorm = %12.4le\n", rnum1);
#else
    fprintf(infofp, "(ivio=1) pnorm = %12.4e\n", rnum1);
#endif
    break;
     
  case PRNT_FNORM:
    rnum1 = *(va_arg(ap, realtype *));

#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "fnorm(L2) = %20.8Le\n", rnum1);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(infofp, "fnorm(L2) = %20.8le\n", rnum1);
#else
    fprintf(infofp, "fnorm(L2) = %20.8e\n", rnum1);
#endif
    break;

  case PRNT_LAM:
    rnum1 = *(va_arg(ap, realtype *));
    rnum2 = *(va_arg(ap, realtype *));
    rnum3 = *(va_arg(ap, realtype *));

#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "min_lam = %11.4Le  ", rnum1);
    fprintf(infofp, "f1norm = %11.4Le  ", rnum2);
    fprintf(infofp, "pnorm = %11.4Le\n", rnum3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(infofp, "min_lam = %11.4le  ", rnum1);
    fprintf(infofp, "f1norm = %11.4le  ", rnum2);
    fprintf(infofp, "pnorm = %11.4le\n", rnum3);
#else
    fprintf(infofp, "min_lam = %11.4e  ", rnum1);
    fprintf(infofp, "f1norm = %11.4e  ", rnum2);
    fprintf(infofp, "pnorm = %11.4e\n", rnum3);
#endif
    break;

  case PRNT_ALPHA:
    rnum1 = *(va_arg(ap, realtype *));
    rnum2 = *(va_arg(ap, realtype *));
    rnum3 = *(va_arg(ap, realtype *));
    rnum4 = *(va_arg(ap, realtype *));

#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "fnorm = %15.8Le  ", rnum1);
    fprintf(infofp, "f1norm = %15.8Le  ", rnum2);
    fprintf(infofp, "alpha_cond = %15.8Le  ", rnum3);
    fprintf(infofp, "pnorm = %15.8Le\n", rnum4);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(infofp, "fnorm = %15.8le  ", rnum1);
    fprintf(infofp, "f1norm = %15.8le  ", rnum2);
    fprintf(infofp, "alpha_cond = %15.8le  ", rnum3);
    fprintf(infofp, "pnorm = %15.8le\n", rnum4);
#else
    fprintf(infofp, "fnorm = %15.8e  ", rnum1);
    fprintf(infofp, "f1norm = %15.8e  ", rnum2);
    fprintf(infofp, "alpha_cond = %15.8e  ", rnum3);
    fprintf(infofp, "pnorm = %15.8e\n", rnum4);
#endif
    break;

  case PRNT_BETA:
    rnum1 = *(va_arg(ap, realtype *));
    rnum2 = *(va_arg(ap, realtype *));
    rnum3 = *(va_arg(ap, realtype *));

#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "f1norm = %15.8Le  ", rnum1);
    fprintf(infofp, "beta_cond = %15.8Le  ", rnum2);
    fprintf(infofp, "lam = %15.8Le\n", rnum3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(infofp, "f1norm = %15.8le  ", rnum1);
    fprintf(infofp, "beta_cond = %15.8le  ", rnum2);
    fprintf(infofp, "lam = %15.8le\n", rnum3);
#else
    fprintf(infofp, "f1norm = %15.8e  ", rnum1);
    fprintf(infofp, "beta_cond = %15.8e  ", rnum2);
    fprintf(infofp, "lam = %15.8e\n", rnum3);
#endif
    break;

  case PRNT_ALPHABETA:
    rnum1 = *(va_arg(ap, realtype *));
    rnum2 = *(va_arg(ap, realtype *));
    rnum3 = *(va_arg(ap, realtype *));
    rnum4 = *(va_arg(ap, realtype *));

#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "f1norm = %15.8Le  ", rnum1);
    fprintf(infofp, "alpha_cond = %15.8Le  ", rnum2);
    fprintf(infofp, "beta_cond = %15.8Le  ", rnum3);
    fprintf(infofp, "lam = %15.8Le\n", rnum4);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(infofp, "f1norm = %15.8le  ", rnum1);
    fprintf(infofp, "alpha_cond = %15.8le  ", rnum2);
    fprintf(infofp, "beta_cond = %15.8le  ", rnum3);
    fprintf(infofp, "lam = %15.8le\n", rnum4);
#else
    fprintf(infofp, "f1norm = %15.8e  ", rnum1);
    fprintf(infofp, "alpha_cond = %15.8e  ", rnum2);
    fprintf(infofp, "beta_cond = %15.8e  ", rnum3);
    fprintf(infofp, "lam = %15.8e\n", rnum4);
#endif
    break;

  case PRNT_ADJ:
    inum1 = *(va_arg(ap, long int *));
    fprintf(infofp, "no. of lambda adjustments = %ld\n", inum1);
    break;

  }

  /* finalize argument processing */

  va_end(ap);

  return;
}

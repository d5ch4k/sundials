/*******************************************************************
 * File          : kinbbdpre.h                                     *
 * Programmers   : Allan Grant Taylor, Alan C Hindmarsh, and       *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 07 February 2004                                *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/kinsol/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the header file for the KINBBDPRE module, for a         *
 * band-block-diagonal preconditioner, i.e. a block-diagonal       *
 * matrix with banded blocks, for use with KINSol, KINSpgmr, and   *
 * the parallel implementaion of the NVECTOR module.               *
 *                                                                 *
 * Summary:                                                        *
 *                                                                 *
 * These routines provide a preconditioner matrix for KINSol that  *
 * is block-diagonal with banded blocks.  The blocking corresponds *
 * to the distribution of the dependent variable vector u among    *
 * the processors.  Each preconditioner block is generated from    *
 * the Jacobian of the local part (on the current processor) of a  *
 * given function g(u) approximating f(u).  The blocks are         *
 * generated by a difference quotient scheme on each processor     *
 * independently, utilizing the assumed banded structure with      *
 * given half-bandwidths.                                          *
 *                                                                 *
 * The user's calling program should have the following form:      *
 *                                                                 *
 *   #include "nvector_parallel.h"                                 *
 *   #include "kinbbdpre.h"                                        *
 *   ...                                                           *
 *   void *p_data;                                                 *
 *   ...                                                           *
 *   nvSpec = NV_SpecInit_Parallel(...);                           *
 *   ...                                                           *
 *   kin_mem = KINCreate(...);                                     *
 *   KINMalloc(kin_mem,...);                                       *
 *   ...                                                           *
 *   p_data = KBBDPrecAlloc(Nlocal, mu, ml, ...);                  *
 *   ...                                                           *
 *   flag = KBBDSpgmr(kin_mem, maxl, p_data);                      *
 *   ...                                                           *
 *   KINSol(...);                                                  *
 *   ...                                                           *
 *   KBBDPrecFree(p_data);                                         *
 *   ...                                                           *
 *   KINFree(...);                                                 *
 *   ...                                                           *
 *   NV_SpecFree_Parallel(nvSpec);                                 *
 *                                                                 *
 *                                                                 *
 * The user-supplied routines required are:                        *
 *                                                                 *
 *  func   is the function f(u) defining the system to be solved:  *
 *                   f(u) = 0                                      *
 *                                                                 *
 *  glocal is the function defining the approximation g(u) to f(u) *
 *                                                                 *
 *  gcomm  is the function to do communication needed for glocal   *
 *                                                                 *
 *                                                                 *
 * Notes:                                                          *
 *                                                                 *
 * 1) This header file is included by the user for the definition  *
 *    of the KBBDData type and for needed function prototypes.     *
 *                                                                 *
 * 2) The KBBDPrecAlloc call includes half-bandwiths mu and ml to  *
 *    be used in the approximate Jacobian.  They need not be the   *
 *    true half-bandwidths of the Jacobian of the local block of g,*
 *    when smaller values may provide a greater efficiency.        *
 *    Also, mu and ml need not be the same on every processor.     *
 *                                                                 *
 * 3) The actual name of the user's f function is passed to        *
 *    KINMalloc, and the names of the user's glocal and gcomm      *
 *    functions are passed to KBBDPrecAlloc.                       *
 *                                                                 *
 * 4) The pointer to the user-defined data block f_data, which is  *
 *    passed to KINMalloc, is also passed to KBBDPrecAlloc, and    *
 *    is available to the user in glocal and gcomm.                *
 *                                                                 *
 * 5) Optional outputs specific to this module are available by    *
 *    way of macros listed below.  These include work space sizes  *
 *    and the cumulative number of glocal calls.  The costs        *
 *    associated with this module also include nsetups banded LU   *
 *    factorizations, nsetups gcomm calls, and nps banded          *
 *    backsolve calls, where nsetups and nps are KINSol optional   *
 *    outputs.                                                     *
 *******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif
#ifndef _kbbdpre_h
#define _kbbdpre_h

#include "kinsol.h"
#include "sundialstypes.h"
#include "nvector.h"
#include "band.h"


/******************************************************************
 * Type : KINCommFn                                               *
 *----------------------------------------------------------------*        
 * The user must supply a function of type KINCommFn which        *
 * performs all inter-processor communication necessary to        *
 * evaluate the approximate system function described above.      *
 *                                                                *
 * This function takes as input the local vector size Nlocal,     *
 * the solution vector u, and a pointer to the user-defined       *
 * data block f_data.                                             *
 * The f_data parameter is the same as that passed by the user to *
 * the KINMalloc routine.  The KINCommFn gcomm is expected to save*
 * communicated data in space defined with the structure *f_data. *
 * A KINCommFn gcomm does not have a return value.                *
 *                                                                *
 * Each call to the KINCommFn is preceded by a call to the system *
 * function func at the current iterate uu. Thus functions of the * 
 * type KINCommFn can omit any communications done by f if        *
 * relevant to the evaluation of the local function gloc.         *
 ******************************************************************/

typedef void (*KINCommFn)(long int Nlocal, N_Vector u, void *f_data);

/******************************************************************
 * Type : KINLocalFn                                              *
 *----------------------------------------------------------------*        
 * The user must supply a function g(u) which approximates the    *
 * function f for the system f(u) = 0. , and which  is computed   *
 * locally (without inter-processor communication).               *
 * (The case where g is mathematically identical to f is allowed.)*
 * The implementation of this function must have type KINLocalFn. *
 * It takes as input the local vector size Nlocal, the local      *
 * solution vector uu, the returned local g values vector, and a  *
 * pointer to the user-defined  data block f_data.  It is to      *
 * compute the local part of g(u) and store in the vector gval    *
 * Providing memory for uu and gval is handled within  the        *
 * preconditioner module.) It is expected that this routine will  *
 * save communicated data in work space defined by the user, and  *
 * made available to the preconditioner function for the problem. *
 * The f_data parameter is the same as that passed by the user    *
 * to the KINMalloc routine.                                      *
 * A KINLocalFn gloc does not have a return value.                *
 ******************************************************************/

typedef void (*KINLocalFn)(long int Nlocal, N_Vector uu,
                           N_Vector gval, void *f_data);
 
/*********************** Definition of KBBDData *****************/

typedef struct {

  /* passed by user to KBBDPrecAlloc, used by Precond/Psolve functions: */
  long int ml, mu;
  KINLocalFn gloc;
  KINCommFn gcomm;

  /* relative error for the Jacobian DQ routine */
  realtype rel_uu;

  /* allocated for use by KBBDPrecSetup */
  N_Vector vtemp3;

  /* set by KBBDPrecSetup and used by KBBDPrecSolve: */
  BandMat PP;
  long int *pivots;

  /* set by KBBDPrecAlloc and used by KBBDPrecSetup: */
  long int n_local;

  /* available for optional output: */
  long int rpwsize;
  long int ipwsize;
  long int nge;

  /* Pointer to KINSOL memory */
  KINMem kin_mem;

} *KBBDPrecData;


/******************************************************************
 * Function : KBBDPrecAlloc                                       *
 *----------------------------------------------------------------*
 * KBBDPrecAlloc allocates and initializes a KBBDData structure   *
 * to be passed to KINSpgmr (and then used by KBBDPrecSetup and   *
 * KBBDPrecSolve).                                                *
 *                                                                *
 * The parameters of KBBDPrecAlloc are as follows:                *
 *                                                                *
 * Nlocal  is the length of the local block of the vectors u etc. *
 *         on the current processor.                              *
 *                                                                *
 *                                                                *
 * mu, ml  are the upper and lower half-bandwidths to be used     *
 *         in the computation of the local Jacobian blocks.       *
 *                                                                *
 * dq_rel_uu is the relative error to be used in the difference   *
 *           quotient Jacobian calculation in the preconditioner. *
 *           The default is sqrt(unit roundoff), obtained by      *
 *           passing in 0.                                        *
 *                                                                *
 * gloc    is the name of the user-supplied function g(u) that    *
 *         approximates f and whose local Jacobian blocks are     *
 *         to form the preconditioner.                            *
 *                                                                *
 * gcomm   is the name of the user-defined function that performs *
 *         necessary inter-processor communication for the        *
 *         execution of gloc.                                     *
 *                                                                *
 * f_data  is a pointer to the optional user data block.          *
 *                                                                *
 * KBBDPrecAlloc returns the storage allocated,                   *
 * or NULL if the request for storage cannot be satisfied.        *
 ******************************************************************/

void *KBBDPrecAlloc(void *kinmem, long int Nlocal, 
                    long int mu, long int ml,
                    realtype dq_rel_uu, 
                    KINLocalFn gloc, KINCommFn gcomm);

/******************************************************************
 * Function : KBBDSpgmr                                           *
 *----------------------------------------------------------------*
 * KBBDSpgmr links the KINBBDPRE preconditioner to the KINSPGMR   *
 * linear solver. It performs the following actions:              *
 *  1) Calls the KINSPGMR specification routine and attaches the  *
 *     KINSPGMR linear solver to the KINSOL solver;               *
 *  2) Sets the preconditioner data structure for KINSPGMR        *
 *  3) Sets the preconditioner setup routine for KINSPGMR         *
 *  4) Sets the preconditioner solve routine for KINSPGMR         *
 *                                                                *
 * Its first 3 arguments are the same as for KINSpgmr (see        *
 * kinspgmr.h). The last argument is the pointer to the KBBDPRE   *
 * memory block returned by KBBDPrecAlloc.                        * 
 * Note that the user need not call KINSpgmr anymore.             *
 *                                                                *
 * Possible return values are:                                    *
 *   (from kinsol.h) SUCCESS                                      *
 *                   LIN_NO_MEM                                   *
 *                   LMEM_FAIL                                    *
 *                   LIN_NO_LMEM                                  *
 *                   LIN_ILL_INPUT                                *
 *   Additionaly, if KBBDPrecAlloc was not previously called,     *
 *   KBBDSpgmr returns BBDP_NO_PDATA (defined below).             *
 *                                                                *
 ******************************************************************/

int KBBDSpgmr(void *kinmem, int maxl, void *p_data);

/******************************************************************
 * Function : KBBDPrecFree                                        *
 *----------------------------------------------------------------*
 * KBBDPrecFree frees the memory block p_data allocated by the    *
 * call to KBBDPrecAlloc.                                         *
 ******************************************************************/

void KBBDPrecFree(void *p_data);

/******************************************************************
 * Function : KBBDPrecGet*                                        *
 *----------------------------------------------------------------*
 *                                                                *
 ******************************************************************/

int KBBDPrecGetIntWorkSpace(void *p_data, long int *leniwBBDP);
int KBBDPrecGetRealWorkSpace(void *p_data, long int *lenrwBBDP);
int KBBDPrecGetNumGfnEvals(void *p_data, long int *ngevalsBBDP);

/* Return values for KBBDPrecGet* functions */
/* OKAY = 0 */
enum { BBDP_NO_PDATA = -1 };

/****** Prototypes of functions KBBDPrecSetup and KBBDPrecSolve *********/

int KBBDPrecSetup(N_Vector uu, N_Vector uscale,
                  N_Vector fval, N_Vector fscale, 
                  void *p_data,
                  N_Vector vtemp1, N_Vector vtemp2);

int KBBDPrecSolve(N_Vector uu, N_Vector uscale,
                  N_Vector fval, N_Vector fscale, 
                  N_Vector vv, void *p_data,
                  N_Vector vtemp);

#endif
#ifdef __cplusplus
}
#endif

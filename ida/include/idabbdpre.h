/******************************************************************
 * File          : idabbdpre.h                                    *
 * Programmers   : Allan G. Taylor, Alan C Hindmarsh, and         *
 *                 Radu Serban @ LLNL                             *
 * Version of    : 11 July 2002                                   *
 *----------------------------------------------------------------*
 * This is the header file for the IDABBDPRE module, for a        *
 * band-block-diagonal preconditioner, i.e. a block-diagonal      *
 * matrix with banded blocks, for use with IDA and IDASpgmr.      *
 *                                                                *
 * Summary:                                                       *
 *                                                                *
 * These routines provide a preconditioner matrix for IDA that    *
 * is block-diagonal with banded blocks.  The blocking corresponds*
 * to the distribution of the dependent variable vector y among   *
 * the processors.  Each preconditioner block is generated from   *
 * the Jacobian of the local part (on the current processor) of a *
 * given function G(t,y,y') approximating F(t,y,y').  The blocks  *
 * are generated by a difference quotient scheme on each processor*
 * independently.  This scheme utilizes an assumed banded         *
 * structure with given half-bandwidths, mudq and mldq.           *
 * However, the banded Jacobian block kept by the scheme has      *
 * half-bandwiths mukeep and mlkeep, which may be smaller.        *
 *                                                                *
 * The user's calling program should have the following form:     *
 *                                                                *
 *   #include "idabbdpre.h"                                       *
 *   #include "nvector_parallel.h"                                *
 *   ...                                                          *
 *   IBBDData p_data;                                             *
 *   ...                                                          *
 *   machEnv = M_EnvInit_Parallel(...);                           *
 *   ...                                                          *
 *   ida_mem = IDAMalloc(...);                                    *
 *   ...                                                          *
 *   p_data = IBBDAlloc(Nlocal, mudq, mldq, mukeep,mlkeep,        *
 *              dq_rel_yy, glocal, gcomm, ida_mem, res_data);     *
 *   ...                                                          *
 *   ier = IDASpgmr(ida_mem, IBBDPrecon, IBBDPSol,                *
 *               gstype, maxl, maxrs, eplifac, dqincfac, p_data); *
 *   ...                                                          *
 *   ier = IDASolve(...);                                         *
 *   ...                                                          *
 *   IBBDFree(p_data);                                            *
 *   ...                                                          *
 *   IDAFree(...);                                                *
 *                                                                *
 *   M_EnvFree_Parallel(machEnv);                                 *
 *                                                                *
 * The user-supplied routines required are:                       *
 *                                                                *
 *   res  is the function F(t,y,y') defining the DAE system to    *
 *   be solved:  F(t,y,y') = 0.                                   *
 *                                                                *
 *   glocal  is the function defining a local approximation       *
 *   G(t,y,y') to F, for the purposes of the preconditioner.      *
 *                                                                *
 *   gcomm  is the function performing communication needed       *
 *                   for glocal.                                  *
 *                                                                *
 *                                                                *
 * Notes:                                                         *
 *                                                                *
 * 1) This header file is included by the user for the definition *
 *    of the IBBDData type and for needed function prototypes.    *
 *                                                                *
 * 2) The IBBDAlloc call includes half-bandwidths mudq and        *
 *    mldq to be used in the approximate Jacobian.  They need     *
 *    not be the true half-bandwidths of the Jacobian of the      *
 *    local block of G, when smaller values may provide a greater *
 *    efficiency. Similarly, mukeep and mlkeep, specifying the    *
 *    bandwidth kept for the approximate Jacobian, need not be    *
 *    the true half-bandwidths. Also, mukeep, mlkeep, mudq, and   *
 *    mldq need not be the same on every processor.               *
 *                                                                *
 * 3) The actual name of the user's res function is passed to     *
 *    IDAMalloc, and the names of the user's glocal and gcomm     *
 *    functions are passed to IBBDAlloc.                          *
 *                                                                *
 * 4) The pointer to the user-defined data block res_data, which  *
 *    is passed to IDAMalloc, is also passed to IBBDAlloc, and    *
 *    is available to the user in glocal and gcomm.               *
 *                                                                *
 * 5) The two functions IBBDPrecon and IBBDPSol are never called  *
 *    by the user explicitly; their names are simply passed to    *
 *    IDASpgmr as in the above.                                   *
 *                                                                *
 * 6) Optional outputs specific to this module are available by   *
 *    way of macros listed below.  These include work space sizes *
 *    and the cumulative number of glocal calls.  The costs       *
 *    associated with this module also include nsetups banded LU  *
 *    factorizations, nsetups gcomm calls, and nps banded         *
 *    backsolve calls, where nsetups and nps are IDA optional     *
 *    outputs.                                                    *
 ******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif
#ifndef _ibbdpre_h
#define _ibbdpre_h

#include "ida.h"
#include "sundialstypes.h"
#include "nvector.h"
#include "band.h"


/******************************************************************
 * Type : IDALocalFn                                              *
 *----------------------------------------------------------------*        
 * The user must supply a function G(t,y,y') which approximates   *
 * the function F for the system F(t,y,y') = 0, and which is      *
 * computed locally (without inter-processor communication).      *
 * (The case where G is mathematically identical to F is allowed.)*
 * The implementation of this function must have type IDALocalFn. *
 *                                                                *
 * This function takes as input the independent variable value tt,*
 * the current solution vector yy, the current solution           *
 * derivative vector yp, and a pointer to the user-defined data   *
 * block res_data.  It is to compute the local part of G(t,y,y')  *
 * and store it in the vector gval. (Providing memory for yy and  *
 * gval is handled within this preconditioner module.) It is      *
 * expected that this routine will save communicated data in work *
 * space defined by the user, and made available to the           *
 * preconditioner function for the problem. The res_data          *
 * parameter is the same as that passed by the user to the        *
 * IDAMalloc routine.                                             *
 *                                                                *
 * A IDALocalFn glocal is to return an int, defined in the same   *
 * way as for the residual function: 0 (success), +1 or -1 (fail).*
 ******************************************************************/

typedef int (*IDALocalFn)(realtype tt, N_Vector yy, 
                          N_Vector yp, N_Vector gval, void *res_data);
 

/******************************************************************
 * Type : IDACommFn                                               *
 *----------------------------------------------------------------*        
 * The user must supply a function of type IDACommFn which        *
 * performs all inter-processor communication necessary to        *
 * evaluate the approximate system function described above.      *
 *                                                                *
 * This function takes as input the solution vectors yy and yp,   *
 * and a pointer to the user-defined data block res_data. The     *
 * res_data parameter is the same as that passed by the user to   *
 * the IDAMalloc routine.                                         *
 *                                                                *
 * The IDACommFn gcomm is expected to save communicated data in   *
 * space defined with the structure *res_data.                    *
 *                                                                *
 * A IDACommFn gcomm returns an int value equal to 0 (success),   *
 * > 0 (recoverable error), or < 0 (unrecoverable error).         *
 *                                                                *
 * Each call to the IDACommFn is preceded by a call to the system *
 * function res with the same vectors yy and yp. Thus the         *
 * IDACommFn gcomm can omit any communications done by res if     *
 * relevant to the evaluation of the local function glocal.       *
 ******************************************************************/

typedef int (*IDACommFn)(N_Vector yy, N_Vector yp, void *res_data);


/*********************** Definition of IBBDData *****************/

typedef struct {

  /* passed by user to IBBDAlloc, used by Precond/Psolve functions: */
  void *res_data;
  integertype mudq, mldq, mukeep, mlkeep;
  realtype rel_yy;
  IDALocalFn glocal;
  IDACommFn gcomm;

 /* allocated for use by IBBDPrecon */
  N_Vector tempv4;

  /* set by IBBDPrecon and used by IBBDPSol: */
  BandMat PP;
  integertype *pivots;

  /* set by IBBDAlloc and used by IBBDPrecond */
  integertype n_local;

  /* available for optional output: */
  integertype rpwsize;
  integertype ipwsize;
  integertype nge;

} *IBBDData;


/*************** Macros for optional outputs **********************
 *                                                                *
 * IBBD_RPWSIZE(pdata) returns the size of the real work space,   *
 * in realtype words, used by this preconditioner module.         *
 * This size is local to the current processor.                   *
 *                                                                *
 * IBBD_IPWSIZE(pdata) returns the size of the integer work space,*
 * in integertype words, used by this preconditioner module.      *
 * This size is local to the current processor.                   *
 *                                                                *
 * IBBD_NGE(pdata) returns the number of G(t,y,y') evaluations,   *
 * i.e. the number of calls to the glocal function, so far.       *
 ******************************************************************/

#define IBBD_RPWSIZE(pdata) (pdata->rpwsize)
#define IBBD_IPWSIZE(pdata) (pdata->ipwsize)
#define IBBD_NGE(pdata) (pdata->nge)


/******************************************************************
 * Function : IBBDAlloc                                           *
 *----------------------------------------------------------------*
 * IBBDAlloc allocates and initializes an IBBDData structure      *
 * to be passed to IDASpgmr (and subsequently used by IBBDPrecon  *
 * and IBBDPSol.                                                  *
 *                                                                *
 * The parameters of IBBDAlloc are as follows:                    *
 *                                                                *
 * Nlocal  is the length of the local block of the vectors yy etc.*
 *         on the current processor.                              *
 *                                                                *
 * mudq, mldq  are the upper and lower half-bandwidths to be used *
 *         in the computation of the local Jacobian blocks.       *
 *                                                                *
 * mukeep, mlkeep are the upper and lower half-bandwidths to be   *
 *         used in saving the Jacobian elements in the local      *
 *         block of the preconditioner matrix PP.                 *
 *                                                                *
 * dq_rel_yy is an optional input.  It is the relative increment  *
 *         to be used in the difference quotient routine for      *
 *         Jacobian calculation in the preconditioner.  The       *
 *         default is sqrt(unit roundoff), and specified by       *
 *         passing dq_rel_yy = 0.                                 *
 *                                                                *
 * glocal    is the name of the user-supplied function G(t,y,y')  *
 *         that approximates F and whose local Jacobian blocks    *
 *         are to form the preconditioner.                        *
 *                                                                *
 * gcomm   is the name of the user-defined function that performs *
 *         necessary inter-processor communication for the        *
 *         execution of glocal.                                   *
 *                                                                *
 * idamem  is the pointer to the IDA memory returned by IDAMalloc.*
 *                                                                *
 * res_data is a pointer to the optional user data block, as      *
 *         passed to IDAMalloc.                                   *
 *                                                                *
 * IBBDAlloc returns the storage allocated (type IBBDData),       *
 * or NULL if the request for storage cannot be satisfied.        *
 ******************************************************************/

IBBDData IBBDAlloc(integertype Nlocal, integertype mudq, integertype mldq, 
                   integertype mukeep, integertype mlkeep, realtype dq_rel_yy, 
                   IDALocalFn glocal, IDACommFn gcomm, 
                   void *idamem, void *res_data);

/******************************************************************
 * Function : IDAReInitBBD                                        *
 *----------------------------------------------------------------*
 * IDAReInitBBD re-initializes the IDABBDPRE module when solving a*
 * sequence of problems of the same size with IDASPGMR/IDABBDPRE, *
 * provided there is no change in Nlocal, mukeep, or mlkeep.      *
 * After solving one problem, and after calling IDAReInit to      *
 * re-initialize IDA for a subsequent problem, call IDAReInitBBD. *
 * Then call IDAReInitSpgmr or IDASpgmr, if necessary, to         *
 * re-initialize the Spgmr linear solver, depending on changes    *
 * made in its input parameters, before calling IDASolve.         *
 *                                                                *
 * The first argument to IDAReInitBBD must be the pointer p_data  *
 * that was returned by IBBDAlloc.  All other arguments have      *
 * the same names and meanings as those of IBBDAlloc.             *
 *                                                                *
 * The return value of IDAReInitBBD is 0, indicating success.     *
 ******************************************************************/

int IDAReInitBBD(IBBDData p_data, integertype Nlocal, integertype mudq,
                 integertype mldq, integertype mukeep, integertype mlkeep,
                 realtype dq_rel_yy, IDALocalFn glocal, IDACommFn gcomm, 
                 void *idamem, void *res_data);


/******************************************************************
 * Function : IBBDFree                                            *
 *----------------------------------------------------------------*
 * IBBDFree frees the memory block p_data allocated by the call   *
 * to IBBDAlloc.                                                  *
 ******************************************************************/

void IBBDFree(IBBDData p_data);

/* Prototypes of IBDPrecon and IBBDPSol */

int IBBDPrecon(integertype Neq, realtype tt, N_Vector yy,
               N_Vector yp, N_Vector rr, realtype cj, ResFn res,
               void *res_data, void *P_data, N_Vector ewt, N_Vector constraints,
               realtype hh, realtype uround, long int *nrePtr,
               N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);
 
int IBBDPSol(integertype Neq, realtype tt, N_Vector yy, N_Vector yp, 
             N_Vector rr, realtype cj, ResFn res, void *res_data, 
             void *P_data, N_Vector ewt, realtype delta, N_Vector rvec, 
             N_Vector zvec, long int *nrePtr, N_Vector tempv);

#endif
#ifdef __cplusplus
}
#endif


#ifndef ITSOL_FGMRES_H__
#define ITSOL_FGMRES_H__

#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------
|                 *** Preconditioned FGMRES ***                  
+-----------------------------------------------------------------------
| on entry:
|---------- 
|
|(Amat)   = matrix struct. the matvec operation is Amat->matvec.
|(lu)     = preconditioner struct.. the preconditioner is lu->precon
|           if (lu == NULL) the no-preconditioning option is invoked.
| rhs     = real vector of length n containing the right hand side.
| sol     = real vector of length n containing an initial guess to the
|           solution on input.
|
| on return:
|---------- 
| fgmr      int =  0 --> successful return.
|           int =  1 --> convergence not achieved in itmax iterations.
| sol     = contains an approximate solution (upon successful return).
| nits    = has changed. It now contains the number of steps required
|           to converge -- 
| res     = relative residual.
+-----------------------------------------------------------------------
| internal work arrays:
|----------       
| vv      = work array of length [im+1][n] (used to store the Arnoldi
|           basis)
| hh      = work array of length [im][im+1] (Arnoldi matrix)
| z       = work array of length [im][n] to store preconditioned vectors
+-----------------------------------------------------------------------
| subroutines called :
|     matvec and
|     preconditionning operation 
+---------------------------------------------------------------------*/
int itsol_solver_fgmres(ITS_SMat *Amat, ITS_PC *lu, double *rhs, double *sol, ITS_PARS io,
        int *nits, double *res);

#ifdef __cplusplus
}
#endif
#endif

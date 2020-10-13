
#ifndef ITSOL_ILUT_H__
#define ITSOL_ILUT_H__

#include <math.h>
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------------
 * ILUT preconditioner
 * incomplete LU factorization with dual truncation mechanism
 * NOTE : no pivoting implemented as yet in GE for diagonal elements
 *----------------------------------------------------------------------------
 * Parameters
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * csmat    = matrix stored in SpaFmt format -- see heads.h for details
 * lu       = pointer to a ILUSpar struct -- see heads.h for details
 * lfil     = integer. The fill-in parameter. Each column of L and
 *            each column of U will have a maximum of lfil elements.
 *            WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
 *            EARLIER VERSIONS. 
 *            lfil must be .ge. 0.
 * tol      = real*8. Sets the threshold for dropping small terms in the
 *            factorization. See below for details on dropping strategy.
 * fp       = file pointer for error log (might be stdout)
 *
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr  = -1  --> Illegal value for lfil
 *            ierr  = -2  --> zero diagonal or zero col encountered
 * lu->n    = dimension of the matrix
 *   ->L    = L part -- stored in SpaFmt format
 *   ->D    = Diagonals
 *   ->U    = U part -- stored in SpaFmt format
 *----------------------------------------------------------------------------
 * Notes:
 * ======
 * All the diagonals of the input matrix must not be zero
 *----------------------------------------------------------------------------
 * Dual drop-off strategy works as follows. 
 *
 * 1) Theresholding in L and U as set by tol. Any element whose size
 *    is less than some tolerance (relative to the norm of current
 *    row in u) is dropped.
 *
 * 2) Keeping only the largest lfil elements in the i-th column of L
 *    and the largest lfil elements in the i-th column of U.
 *
 * Flexibility: one can use tol=0 to get a strategy based on keeping the
 * largest elements in each column of L and U. Taking tol .ne. 0 but lfil=n
 * will give the usual threshold strategy (however, fill-in is then
 * impredictible).
 *--------------------------------------------------------------------------*/
int itsol_pc_ilut(ITS_SparMat *csmat, ITS_ILUSpar *lu, int lfil, double tol, FILE *fp);

int itsol_pc_lutsolC(double *y, double *x, ITS_ILUSpar *lu);

#ifdef __cplusplus
}
#endif
#endif

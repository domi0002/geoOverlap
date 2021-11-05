
#ifndef ITSOL_VBILUT_H__
#define ITSOL_VBILUT_H__

#include "utils.h"
#include "mat-utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------------
 * Block ILUT (BILUT) preconditioner
 * Block incomplete LU factorization with dual truncation mechanism
 * NOTE : no pivoting implemented as yet in GE for diagonal blocks 
 *----------------------------------------------------------------------------
 * Parameters
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * vbmat    = block matrix stored in VBSpaFmt format -- see type-defs.h for
 *            details on format, the block sizes might be different
 * lu       = pointer to a VBILUSpar struct -- see type-defs.h for details
 *            on format
 * lfil     = integer. The fill-in parameter. Each row of L and
 *            each row of U will have a maximum of lfil elements.
 *            WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
 *            EARLIER VERSIONS. 
 *            lfil must be .ge. 0.
 * tol      = real*8. Sets the threshold for dropping small terms in the
 *            factorization. See below for details on dropping strategy.
 * w        = working array
 * fp       = file pointer for error log (might be stdout)
 *
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr  = -1  --> Illegal value for lfil
 *            ierr  = -2  --> singular block or zero row encountered
 * lu->n    = dimension of the block matrix
 *   ->bsz  = the row/col of the first element of each diagonal block
 *            the size of the i-th row block should be bsz[i+1] - bsz[i]
 *   ->L    = L part -- stored in VBSpaFmt format
 *   ->D    = Diagonals
 *   ->U    = U part -- stored in VBSpaFmt format
 *----------------------------------------------------------------------------
 * Notes:
 * ======
 * All the diagonal blocks of the input block matrix must not be singular
 *----------------------------------------------------------------------------
 * Dual drop-off strategy works as follows. 
 *
 * 1) Theresholding in L and U as set by tol. Any element whose size
 *    is less than some tolerance (relative to the norm of current
 *    row in u) is dropped.
 *
 * 2) Keeping only the largest lfil elements in the i-th row of L
 *    and the largest lfil elements in the i-th row of U.
 *
 * Flexibility: one can use tol=0 to get a strategy based on keeping the
 * largest elements in each row of L and U. Taking tol .ne. 0 but lfil=n
 * will give the usual threshold strategy (however, fill-in is then
 * impredictible).
 *--------------------------------------------------------------------------*/
int itsol_pc_vbilutC(ITS_VBSparMat *vbmat, ITS_VBILUSpar *lu, int lfil, double tol, ITS_BData *w, FILE *fp);

#ifdef __cplusplus
}
#endif
#endif

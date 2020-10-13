
#ifndef ITSOL_VBILUK_H__
#define ITSOL_VBILUK_H__

#include "mat-utils.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------------
 * Block ILUK preconditioner
 * Block incomplete LU factorization with level of fill dropping
 * This version uses svd to invert diagonal blocks
 *----------------------------------------------------------------------------
 * Parameters
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * lofM     = level of fill: all entries with level of fill > lofM are
 *            dropped. Setting lofM = 0 gives BILU(0).
 * vbmat    = block matrix stored in VBSpaFmt format -- see type-defs.h for
 *            details on format, the block sizes might be different
 * lu       = pointer to a VBILUSpar struct -- see type-defs.h for details
 *            on format
 * fp       = file pointer for error log (might be stderr)
 *
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr  = -1  --> error in lofC
 *            ierr  = -2  --> singular diagonal block
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
 *--------------------------------------------------------------------------*/
int itsol_pc_vbilukC(int lofM, ITS_VBSparMat *vbmat, ITS_VBILUSpar *lu, FILE *fp);

/*--------------------------------------------------------------------
 * symbolic ilu factorization to calculate structure of ilu matrix
 * for specified level of fill
 *--------------------------------------------------------------------
 * on entry:
 * =========
 * lofM     = level of fill, lofM >= 0
 * vbmat    = block matrix stored in VBSpaFmt format -- see type-defs.h for
 *            details on format, size of blocks might be different
 * lu       = pointer to a VBILUSpar struct -- see type-defs.h for details
 *            on format
 * fp       = file pointer for error log (might be stderr)
 *--------------------------------------------------------------------
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr != 0   --> error
 * lu->n    = dimension of the block matrix
 *   ->L    = L part -- stored in BSpaFmt format, patterns only in lofC
 *   ->U    = U part -- stored in BSpaFmt format, patterns only in lofC
 *------------------------------------------------------------------*/
int itsol_pc_vblofC(int lofM, ITS_VBSparMat *vbmat, ITS_VBILUSpar *lu, FILE *fp);

#ifdef __cplusplus
}
#endif
#endif

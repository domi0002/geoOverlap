
#ifndef ITSOL_ILUTC_H__
#define ITSOL_ILUTC_H__

#include "utils.h"
#include "mat-utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------------------------------------------------------
 * Column-based ILUT (ILUTC) preconditioner
 * incomplete LU factorization with dropping strategy specified by input
 * paramter drop.
 * NOTE : no pivoting implemented as yet in GE for diagonal elements
 *---------------------------------------------------------------------
 * Parameters
 *---------------------------------------------------------------------
 * on entry:
 * =========
 * mt       = matrix stored in LUSparMat format -- see heads.h for details
 * lu       = pointer to a ILUSpar struct -- see heads.h for details
 * lfil     = integer. The fill-in parameter. Each COLUMN of L and
 *            each ROW of U will have a maximum of lfil elements.
 *            lfil must be .ge. 0.
 * tol      = Sets the threshold for dropping small terms in the
 *            factorization. See below for details on dropping strategy.
 * 
 * drop     = specifies the dropping strategy: 
 *            0 : standard dropping strategy (default) ;
 *            1 : drops a term if corresponding perturbation relative
 *                to diagonal entry is small. For example for U: 
 *                || L(:,i) || * |u(i,j)/D(i,i)| < tol * | D(j,j) | 
 *                
 *            2 : condition number estimator based on: 
 *                rhs = (1, 1, ..., 1)^T ;
 *            3 : condition number estimator based on: 
 *                rhs = ((+/-)1, ..., (+/-)1)^T 
 *                + maximizing |v_k| ;
 *            4 : condition number estimator based on: 
 *                rhs = ((+/-)1, ..., (+/-)1)^T based on ;
 *                + maximizing the norm of (v_{k+1}, ..., v_n)^T);
 * 
 * fp       = file pointer for error log (might be stdout)
 *
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  =  0  --> successful return.
 *            ierr  = -1  --> Illegal value for lfil
 *            ierr  = -2  --> zero diagonal or zero col encountered
 *            if ierr = -2, try to use diagonal scaling
 * lu->n    = dimension of the matrix
 *   ->L    = L part -- stored in SparCol format
 *   ->D    = Diagonals
 *   ->U    = U part -- stored in SparRow format
 *
 * DETAILS on the data structure : 
 * ============================== 
 * 
 *     |  1   0  0 . . . 0 |     | U11 U12 U13 . . . U1n |
 *     | L21  1  0 . . . 0 |     |  0  U22 U23 . . . U2n |
 *     | L31 L32 1 . . . 0 |     |  0   0  U33 . . . U3n |
 * A = |  .   .  . .     . |  *  |  .   .   .  . . .  .  |
 *     |  .   .  .   .   . |     |  .   .   .    .    .  |
 *     |  .   .  .     . . |     |  .   .   .      .  .  |
 *     | Ln1 Ln2 . . . . 1 |     |  0   0   0  . . . Unn |
 * 
 * L - D - U preconditioner : 
 * 
 *     |  0                 |    | 0 U12 U13 . . . U1n |
 *     | L21  0             |    |    0  U23 . . . U2n |
 *     | L31 L32 0          |    |        0  . . . U3n |
 * L = |  .   .  .  .       | U= |           . . .  .  |
 *     |  .   .  .    .     |    |             .    .  |
 *     |  .   .  .      .   |    |               .  .  |
 *     | Ln1 Ln2 . . .  . 0 |    |                   0 | 
 * D = diag(1/U11, 1/U22, 1/U33, ..., 1/Unn)
 * 
 *   \ - - . . . - - >
 *   | \ - . . . - - >
 *   | | \ - . . . - > U            L: SparCol format
 *   . . | \ . . . . .              U: SparRow format
 *   . . . . \ . . . .
 *   . . . . . \ . . .
 *   | | . . . . \ . .
 *   | | | . . . . \ .
 *   V V V . . . . . \
 *     L               D
 *---------------------------------------------------------------------
 * Workspace
 * =========
 * wL(n)      nonzero values in current column
 * Lid(n)     row numbers of nonzeros in current column
 * Lfirst(n)  At column k
 *            Lfirst(1:k-1), Lfirst(j) points to next value in column j to use
 *            Lfirst(k+1:n), Lfirst(j) indicates if nonzero value exists in
 *                           row j of current column
 * Llist(n)   Llist(j) points to a linked list of columns that will update the
 *            j-th row in U part
 * 
 * wU(n)      nonzero values in current row
 * Uid(n)     column numbers of nonzeros in current row
 * Ufirst(n)  At row k
 *            Ufirst(1:k-1), Ufirst(j) points to next value in row j to use
 *            Ufirst(k+1:n), Ufirst(j) indicates if nonzero value exists in
 *                           column j of current row
 * Ulist(n)   Ulist(j) points to a linked list of rows that will update the
 *            j-th column in L part
 *----------------------------------------------------------------------*/
int itsol_pc_ilutc(ITS_ILUSpar *mt, ITS_ILUSpar *lu, int lfil, double tol, int drop, FILE *fp);

#ifdef __cplusplus
}
#endif
#endif

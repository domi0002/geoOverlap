#ifndef ITSOL_MatOps_H__
#define ITSOL_MatOps_H__

#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/* *-------------------- inversion by svd
   This calls lapack routines for inverting a dense matrix.
   dgetrf and dgetri

   ON ENTRY
 ** A = square matrix of size n x n -- dimensioned with
 ** leading dimension =  n

 ON RETURN A contains the inverse  of the input matrix.
 */
int itsol_invGauss(int nn, double *A);

/* *-------------------- inversion by svd
   This calls lapack routine dgesvd --
   ON ENTRY
 ** A = square matrix of size n x n -- dimensioned with
 ** leading dimension =  n
 ON RETURN A contains the truncated SVD inverse of input matrix.
 ** tolerance set for truncation is TOL and can be changed in
 ** above define statement
 **--------------------
 */
int itsol_invSVD(int nn, double *A);

void itsol_matvecC(ITS_SparMat *mat, double *x, double *y);
void itsol_matvecCSC(ITS_SMat *mat, double *x, double *y);
int itsol_CondestC(ITS_ILUSpar *lu, FILE * fp);
int itsol_diag_scal(ITS_VBSparMat *vbmat);
int itsol_diagvec(ITS_VBSparMat *vbmat, ITS_BData x, ITS_BData y);

/* y = Ax */
void itsol_matvec(ITS_SparMat *A, double *x, double *y); 

/* y = a * Ax + b * y*/
void itsol_amxpby(double a, ITS_SparMat *A, double *x, double b, double *y); 

/* z = a * Ax + b * y*/
void itsol_amxpbyz(double a, ITS_SparMat *A, double *x, double b, double *y, double *z); 

void itsol_matvecCSR(ITS_SMat *mat, double *x, double *y);
void itsol_matvecz(ITS_SparMat *mata, double *x, double *y, double *z);

void itsol_vbmatvec(ITS_VBSparMat *vbmat, double *x, double *y);
void itsol_luinv(int n, double *a, double *x, double *y); 
int itsol_vblusolC(double *y, double *x, ITS_VBILUSpar *lu); 
int itsol_lusolC(double *y, double *x, ITS_ILUSpar *lu); 
int itsol_rpermC(ITS_SparMat *mat, int *perm); 
int itsol_cpermC(ITS_SparMat *mat, int *perm) ; 
int itsol_dpermC(ITS_SparMat *mat, int *perm) ; 
int itsol_CSparTran(ITS_SparMat *amat, ITS_SparMat *bmat, ITS_CompressType *compress);
double itsol_vbnorm2(int sz, double *a);
void itsol_Lsol(ITS_SparMat *mata, double *b, double *x);
void itsol_Usol(ITS_SparMat *mata, double *b, double *x);
int itsol_ascend (ITS_Per4Mat *levmat, double *x, double *wk);
int itsol_descend(ITS_Per4Mat *levmat, double *x, double *wk);
int itsol_armsol2(double *x, ITS_ARMSpar *Prec);
int itsol_condestArms(ITS_ARMSpar *armspre, double *y, FILE *fp);
int itsol_VBcondestC(ITS_VBILUSpar *, FILE *fp); 
int itsol_CondestLUM(ITS_ILUSpar *lu, double *y, double *x, FILE *fp);
void itsol_matvecVBR(ITS_SMat *mat, double *x, double *y);
void itsol_matvecLDU(ITS_SMat *mat, double *x, double *y);
int itsol_preconILU(double *x, double *y, ITS_PC *mat);
int itsol_preconVBR(double *x, double *y, ITS_PC *mat);
int itsol_preconLDU(double *x, double *y, ITS_PC *mat);
int itsol_preconARMS(double *x, double *y, ITS_PC *mat);
ITS_Per4Mat *itsol_Lvsol2(double *x, int nlev, ITS_Per4Mat *levmat, ITS_ILUTSpar * ilusch) ;
int itsol_Uvsol2(double *x, int nlev, int n, ITS_Per4Mat *levmat, ITS_ILUTSpar * ilusch); 
void itsol_SchLsol(ITS_ILUTSpar * ilusch, double *y) ;
void itsol_SchUsol(ITS_ILUTSpar * ilusch, double *y) ;
int itsol_lumsolC(double *y, double *x, ITS_ILUSpar *lu);
int itsol_condestLU(ITS_ILUSpar *lu, FILE *fp);
int itsol_invGauss(int nn, double *A); 

/*----------------------------------------------------------------------------
 * Setup Blocks (rows and columns might be permuted to get better results)
 *----------------------------------------------------------------------------
 * Na Li, Aug 2001
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * csmat   = a matrix stored in SpaFmt format
 * eps     = parameter for deciding when to do a union of two rows
 *           into the same group.  Two rows u and v are merged into a 
 *           block  when cos(<u,v>) == (u,v)/(|u|*|v|), is > eps. 
 *           eps should be <= 1. 
 *----------------------------------------------------------------------------
 * on return:
 * ==========
 * csmat   = matrix stored in SpaFmt format after permutation
 * pnBlock = dimension of the block matrix
 * pnB     = dimension of each block
 *
 *----------------------------------------------------------------------------
 * Combination of hash method and angle method:
 *----------------------------------------------------------------------------
 * Designed for the matrices with symmetric patterns
 * (1) Hash method
 *     a. Calculate hash values
 *     b. qsort rows according to their hash values
 *     c. Get compressed graph as the following format:
 * (2) Angle method
 *     a. Calculate A^T
 *     b. for i-th row, calculate dot product (row_i, row_j) using A*A^T
 *        algorithm where j = i+1, ..., n-1 and group[j] == -1
 *        if cos(<row_i, row_j>) = (row_i,row_j)/|row_i||row_j| is > eps,
 *        we merge row_i and row_j by resetting
 *        group[j] = i and size[i] = size[i]+size[j]
 *--------------------------------------------------------------------------*/
int itsol_init_blocks(ITS_SparMat *csmat, int *pnBlock, int **pnB, int **pperm, double eps);

#ifdef __cplusplus
}
#endif

#endif 

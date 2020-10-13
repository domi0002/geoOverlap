#ifndef ITSOL_UTILS_H__
#define ITSOL_UTILS_H__

#include <stdarg.h>
#include "data-types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* sets.c */
int itsol_nnz_arms(ITS_ARMSpar *PreSt,  FILE *ft);
void itsol_errexit(char *f_str, ...);
void * itsol_malloc(int nbytes, char *msg); 
int itsol_setupCS(ITS_SparMat *amat, int len, int job); 
int itsol_cleanCS(ITS_SparMat *amat);
int itsol_cleanCOO(ITS_CooMat *amat);
int itsol_nnz_cs (ITS_SparMat *A) ;
int itsol_cscpy(ITS_SparMat *amat, ITS_SparMat *bmat);
int itsol_setupP4 (ITS_Per4Mat *amat, int Bn, int Cn,  ITS_SparMat *F,  ITS_SparMat *E);
int itsol_setupVBMat(ITS_VBSparMat *vbmat, int n, int *nB);
int itsol_setupILUT(ITS_ILUTSpar * amat, int len);
int itsol_cleanVBMat(ITS_VBSparMat *vbmat); 
int itsol_nnzVBMat(ITS_VBSparMat *vbmat) ;
int itsol_memVBMat(ITS_VBSparMat *vbmat); 
int itsol_setupVBILU(ITS_VBILUSpar *lu, int n, int *bsz);
int itsol_cleanVBILU(ITS_VBILUSpar *lu); 
int itsol_cleanILU(ITS_ILUSpar *lu);
int itsol_cleanILUT(ITS_ILUTSpar * amat, int indic);
int itsol_cleanP4(ITS_Per4Mat *amat);
int itsol_mallocVBRow(ITS_VBILUSpar *lu, int nrow); 
int itsol_mallocRow(ITS_ILUSpar *lu, int nrow);
void itsol_zrmC(int m, int n, ITS_BData data); 
void itsol_copyBData(int m, int n, ITS_BData dst, ITS_BData src, int isig); 
int itsol_CSRcs(int n, double *a, int *ja, int *ia, ITS_SparMat *mat, int rsa); 
int itsol_csrvbsrC(int job, int nBlk, int *nB, ITS_SparMat *csmat, ITS_VBSparMat *vbmat);  
int itsol_col2vbcol(int col, ITS_VBSparMat *vbmat);
int itsol_nnz_vbilu(ITS_VBILUSpar *lu); 
int itsol_nnz_lev4(ITS_Per4Mat *levmat, int *lev, FILE *ft);
int itsol_setupILU(ITS_ILUSpar *lu, int n);
int itsol_CS2lum(int n, ITS_SparMat *Amat, ITS_ILUSpar *mat, int typ);
int itsol_COOcs(int n, int nnz,  double *a, int *ja, int *ia, ITS_SparMat *bmat);
void itsol_coocsr_(int*, int*, double*, int*, int*, double*, int*, int*);

int itsol_csSplit4(ITS_SparMat *amat, int bsize, int csize, ITS_SparMat *B, ITS_SparMat *F, ITS_SparMat *E, ITS_SparMat *C);
void itsol_setup_arms (ITS_ARMSpar *Levmat);
int itsol_cleanARMS(ITS_ARMSpar *ArmsPre);
void itsol_coocsc(int n, int nnz, double *val, int *col, int *row, double **a, int **ja, int **ia, int job);
int itsol_nnz_ilu(ITS_ILUSpar *lu);
int itsol_CSClum(int n, double *a, int *ja, int *ia, ITS_ILUSpar *mat, int rsa);
int itsol_CSClumC(ITS_SparMat *amat, ITS_ILUSpar *mat, int rsa);

/* misc.c */
int itsol_SparTran(ITS_SparMat *amat, ITS_SparMat *bmat, int job, int flag); 
int itsol_coscalC(ITS_SparMat *mata, double *diag, int nrm);
void itsol_dscale(int n, double *dd, double *x, double * y);
void itsol_hilosort(ITS_SparMat *mat, int abval, int hilo);
void itsol_printmat(FILE *ft, ITS_SparMat *A, int i0, int i1);
void itsol_qqsort(int *ja, double *ma, int left, int right);
void itsol_qsort2C(int *ja, double *ma, int left, int right, int abval); 
void itsol_qsort3i(int *wa, int *cor1, int *cor2, int left, int right); 
void itsol_qsortC(int *ja, double *ma, int left, int right, int abval); 
void itsol_qsortR2I(double *wa, int *cor1, int *cor2, int left, int right); 
int itsol_qsplitC(double *a, int *ind, int n, int ncut);
int itsol_roscalC(ITS_SparMat *mata, double *diag, int nrm);
void itsol_swapj(int v[], int i, int j);
void itsol_swapm(double v[], int i, int j);
double itsol_get_time(void);
int itsol_dumpArmsMat(ITS_ARMSpar *PreSt, FILE *ft);
int itsol_outputLU(ITS_ILUSpar *lu, char *filename);
int itsol_checkperm(int *p, int n);
void itsol_qsortR1I(double *wa, int *cor1, int left, int right);

ITS_CooMat itsol_read_coo(char *Fname);
double itsol_norm(double *x, int n);
double itsol_dot(double *x, double *y, int n);
double itsol_norm2(double *x, int n);
void itsol_copy(double *d, double *s, int n);
void itsol_axpby(double a, double *x, double b, double *y, int n);

#ifdef __cplusplus
}
#endif
#endif 

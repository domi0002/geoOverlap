
#ifndef ITSOL_INDSETC_H__
#define ITSOL_INDSETC_H__

#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

int itsol_add2is(int *last, int nod, int *iord, int *riord);
int itsol_add2com(int *nback, int nod, int *iord, int *riord);

int itsol_PQperm(ITS_SparMat *mat, int bsize, int *Pord, int *Qord, int *nnod, double tol);

/*--------------------------------------------------------------------- 
| greedy algorithm for independent set ordering -- 
|----------------------------------------------------------------------
|     Input parameters:
|     -----------------
|     (mat)  =  matrix in SpaFmt format
|     
|     bsize  =  integer (input) the target size of each block.
|               each block is of size >= bsize. 
|
|     w      =  weight factors for the selection of the elements in the
|               independent set. If w(i) is small i will be left for the
|               vertex cover set. 
|
|     tol    =  a tolerance for excluding a row from independent set.
|
|     Output parameters:
|     ------------------ 
|     iord   = permutation array corresponding to the independent set 
|     ordering.  Row number i will become row number iord[i] in 
|     permuted matrix.
|     
|     nnod   = (output) number of elements in the independent set. 
|     
|----------------------------------------------------------------------- 
|     the algorithm searches nodes in lexicographic order and groups
|     the (BSIZE-1) nearest nodes of the current to form a block of
|     size BSIZE. The current algorithm does not use values of the matrix.
|---------------------------------------------------------------------*/ 
int itsol_indsetC(ITS_SparMat *mat, int bsize, int *iord, int *nnod, double tol);
int itsol_weightsC(ITS_SparMat *mat, double *w);

/*---------------------------------------------------------------------
| does a preselection of possible diagonal entries. will return a list
| of "count" bi-indices representing "good" entries to be selected as 
| diagonal elements -- the order is important (from best to
| to worst). The list is in the form (icor(ii), jcor(ii)) 
|
|      ON ENTRY: 
|       mat   = matrix in ITS_SparMat *format 
|       tol   = tolerance used for selecting best rows -|
|       job   = indicates whether or not to permute the max entry in 
|               each row to first position 
|        NOTE: CAN RECODE BY HAVING JCOR CARRY THE INDEX IN ROW[I] WHERE
|              MAX IS LOCATED.. 
|       
|      ON RETURN: 
|       icor  = list of row indices of entries selected 
|       jcor  = list of column indices of entries selected 
|       count = number of entries selected (size of B block) 
|--------------------------------------------------------------------*/
int itsol_preSel(ITS_SparMat *mat, int *icor, int *jcor, int job, double tol, int *count);

#ifdef __cplusplus
}
#endif
#endif

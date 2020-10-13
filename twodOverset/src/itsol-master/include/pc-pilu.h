
#ifndef ITSOL_PILUNEW_H__
#define ITSOL_PILUNEW_H__

#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------------------------------------------------------- 
| PARTIAL ILUT -
| Converted to C so that dynamic memory allocation may be implememted
| in order to have no dropping in block LU factors.
|----------------------------------------------------------------------
| Partial block ILU factorization with dual truncation. 
|                                                                      
| |  B   F  |        |    L      0  |   |  U   L^{-1} F |
| |         |   =    |              | * |               |
| |  E   C  |        | E U^{-1}  I  |   |  0       S    |                   
|                                                                      
| where B is a sub-matrix of dimension B->n.
| 
|----------------------------------------------------------------------
|
| on entry:
|========== 
| ( amat ) = Permuted matrix stored in a PerMat4 struct on entry -- 
|            Individual matrices stored in SpaFmt structs.
|            On entry matrices have C (0) indexing.
|            on return contains also L and U factors.
|            Individual matrices stored in SpaFmt structs.
|            On return matrices have C (0) indexing.
|
| lfil[0]  =  number nonzeros in L-part
| lfil[1]  =  number nonzeros in U-part
| lfil[2]  =  number nonzeros in L^{-1} F
| lfil[3]  =  not used
| lfil[4]  =  number nonzeros in Schur complement
|
| droptol[0] = threshold for dropping small terms in L during
|              factorization.
| droptol[1] = threshold for dropping small terms in U.
| droptol[2] = threshold for dropping small terms in L^{-1} F during
|              factorization.
| droptol[3] = threshold for dropping small terms in E U^{-1} during
|              factorization.
| droptol[4] = threshold for dropping small terms in Schur complement
|              after factorization is completed.
|
| On return:
|===========
|
| (schur)  = contains the Schur complement matrix (S in above diagram)
|            stored in SpaFmt struct with C (0) indexing.
|
|
|       integer value returned:
|
|             0   --> successful return.
|             1   --> Error.  Input matrix may be wrong.  (The 
|                         elimination process has generated a
|                         row in L or U whose length is > n.)
|             2   --> Memory allocation error.
|             5   --> Illegal value for lfil or last.
|             6   --> zero row in B block encountered.
|             7   --> zero row in [E C] encountered.
|             8   --> zero row in new Schur complement
|----------------------------------------------------------------------- 
| work arrays:
|=============
| jw, jwrev = integer work arrays of length B->n.
| w         = real work array of length B->n. 
| jw2, jwrev2 = integer work arrays of length C->n.
| w2          = real work array of length C->n. 
|----------------------------------------------------------------------- 
|     All processing is done using C indexing.
|--------------------------------------------------------------------*/
int itsol_pc_pilu(ITS_Per4Mat *amat, ITS_SparMat *B, ITS_SparMat *C, double *droptol, int *lfil, ITS_SparMat *schur);

#ifdef __cplusplus
}
#endif
#endif

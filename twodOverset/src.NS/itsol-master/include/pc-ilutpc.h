
#ifndef ITSOL_ILUTPC_H__
#define ITSOL_ILUTPC_H__

#include "utils.h"
#include "mat-utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------------------------------------------------------- 
| ILUTP -- ILUT with column pivoting -- adapted from ILUTP [Sparskit] 
| Converted to C so that dynamic memory allocation may be implememted.
| All indexing is in C format.
|----------------------------------------------------------------------
| ILUT factorization with dual truncation. 
|----------------------------------------------------------------------
|
| on entry:
|========== 
| (amat) =  Matrix stored in SpaFmt struct.
|
| lfil[5]  =  number nonzeros in L-part
| lfil[6]  =  number nonzeros in U-part     (lfil >= 0)
|
| droptol[5] = threshold for dropping small terms in L during
|              factorization.
| droptol[6] = threshold for dropping small terms in U.
|
| permtol  =  tolerance ratio used to  determine whether or not to permute
|             two columns.  At step i columns i and j are permuted when
|
|                     abs(a(i,j))*permtol > abs(a(i,i))
|
|           [0 --> never permute; good values 0.1 to 0.01]
|
| mband    =  permuting is done within a band extending to mband
|	      diagonals only. 
|             mband = 0 --> no pivoting. 
|	      mband = n --> pivot is searched in whole column 
|
|
| On return:
|===========
|
| (ilusch) =  Contains L and U factors in an LUfact struct.
|             Individual matrices stored in SpaFmt structs.
|             On return matrices have C (0) indexing.
|
| iperm    =  reverse permutation array.
|
|       integer value returned:
|
|             0   --> successful return.
|             1   --> Error.  Input matrix may be wrong.  (The 
|                         elimination process has generated a
|                         row in L or U whose length is > n.)
|             2   --> Memory allocation error.
|             5   --> Illegal value for lfil.
|             6   --> zero row encountered.
|----------------------------------------------------------------------- 
| work arrays:
|=============
| jw, jwrev = integer work arrays of length n
| w         = real work array of length n. 
|----------------------------------------------------------------------- 
|     All processing is done using C indexing.
|--------------------------------------------------------------------*/
int itsol_pc_ilutpC(ITS_SparMat *amat, double *droptol, int *lfil, double permtol,
		int mband, ITS_ILUTSpar *ilusch);

/*---------------------------------------------------------------------- 
| ILUT -
| Converted to C so that dynamic memory allocation may be implememted.
| All indexing is in C format.
|----------------------------------------------------------------------
| ILUT factorization with dual truncation. 
|
| March 1, 2000 - dropping in U: keep entries > tau * diagonal entry
|----------------------------------------------------------------------
|
| on entry:
|========== 
| (amat) =  Matrix stored in SpaFmt struct.
| (ilusch) =  Pointer to ILUTfac struct
|
| lfil[5]  =  number nonzeros in L-part 
| lfil[6]  =  number nonzeros in U-part     (lfil >= 0)
|
| droptol[5] = threshold for dropping small terms in L during 
|              factorization. 
| droptol[6] = threshold for dropping small terms in U.
|
| On return:
|===========
|
| (ilusch) =  Contains L and U factors in an LUfact struct.
|             Individual matrices stored in SpaFmt structs.
|             On return matrices have C (0) indexing.
|
|       integer value returned:
|
|             0   --> successful return.
|             1   --> Error.  Input matrix may be wrong.  (The 
|                         elimination process has generated a
|                         row in L or U whose length is > n.)
|             2   --> Memory allocation error.
|             5   --> Illegal value for lfil or last.
|             6   --> zero row encountered.
|----------------------------------------------------------------------- 
| work arrays:
|=============
| jw, jwrev = integer work arrays of length n
| w         = real work array of length n. 
|----------------------------------------------------------------------- 
|     All processing is done using C indexing.
|--------------------------------------------------------------------*/
int itsol_pc_ilutD(ITS_SparMat *amat, double *droptol, int *lfil, ITS_ILUTSpar *ilusch);

#ifdef __cplusplus
}
#endif
#endif

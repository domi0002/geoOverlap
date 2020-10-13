#ifndef VBLOCK_DATATYPE_DEFS_H__
#define VBLOCK_DATATYPE_DEFS_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h> 

#include "config.h"
#include "protos-deps.h"

#define ITS_MAX_BLOCK_SIZE   100
#define ITS_TOL_DD           0.7  /* diagonal dominance tolerance for arms */

/* FORTRAN style vblock format, compatible for many FORTRAN routines */
#define ITS_DATA(a,row,i,j)  (a[(j)*(row)+(i)])

/* the dimension of ith Block */
#define ITS_B_DIM(bs,i)      (bs[i+1]-bs[i])

#define ITS_MAX_LINE        256
#define ITS_MAX_HBNAME      64
#define ITS_MAX_MAT	        100
#define ITS_MaxNamLen       64

/*--------------------------------------------- 
  | C-style CSR format - used internally
  | for all matrices in CSR format 
  |---------------------------------------------*/
typedef struct ITS_SparMat_
{
    int n;
    int *nzcount;  /* length of each row */
    int **ja;      /* pointer-to-pointer to store column indices  */
    double **ma;   /* pointer-to-pointer to store nonzero entries */

} ITS_SparMat;

typedef struct ITS_CooMat_
{
    int n;
    int nnz;      /* length of each row */
    int *ia;      /* pointer-to-pointer to store column indices  */
    int *ja;      /* pointer-to-pointer to store column indices  */
    double *ma;   /* pointer-to-pointer to store nonzero entries */

} ITS_CooMat;

typedef double *ITS_BData;

typedef struct ITS_VBSparMat_
{
    int n;        /* the block row dimension of the matrix      */
    int *bsz;     /* the row/col of the first element of each   */

    /* diagonal block                             */
    int *nzcount;     /* length of each row                         */
    int **ja;         /* pointer-to-pointer to store column indices */
    ITS_BData **ba;   /* pointer-to-pointer to store nonzero blocks */
    ITS_BData *D;     /* to store inversion of diagonals            */

} ITS_VBSparMat;

typedef struct ITS_VBILUSpar_
{
    int n;
    int *bsz;     /* the row/col of the first element of each   */

    /* diagonal block                             */
    ITS_BData *D;     /* diagonal blocks                            */
    ITS_VBSparMat *L; /* L part blocks                              */
    ITS_VBSparMat *U; /* U part blocks                              */
    int *work;        /* working buffer                             */
    ITS_BData bf;     /* buffer of a temp block                     */
    int DiagOpt;  /* Option for diagonal inversion/solutiob     *
                   * opt =  1 -->> call luinv 
                   * opt == 2 -->> block inverted call dgemv    */

} ITS_VBILUSpar; 

typedef struct ITS_ILUSpar_
{
    int n;
    ITS_SparMat *L;   /* L part elements                            */
    double *D;        /* diagonal elements                          */
    ITS_SparMat *U;   /* U part elements                            */
    int *work;        /* working buffer */

} ITS_ILUSpar;

/*------------------------------------------------------------
  | struct for storing the block LU factorization 
  | contains all the block factors except the 
  | data related to the last block. 
  | n       = size of current block
  | symperm = whether or not permutations are symmetric.
  |           used only in cleanP4..
  | nB      = size of B-block
  | L, U    = ILU factors of B-block
  | F, E    = sparse matrices in (1,2) and (2,1) 
  |           parts of matrix. 
  | perm    = (symmetric) permutation used in factorization
  |           comes from the independent set ordering
  | rperm   = unsymmetric permutation (rows) not used in this
  |           version -- but left here for compatibility..
  | D1, D2  = diagonal matrices (left, right) used for scaling
  |           if scaling option is turned on. Note that the 
  |           method works by scaling the whole matrix first
  |           (at any level) before anything else is done. 
  | wk     = a work vector of length n needed for various tasks
  |            [reduces number of calls to malloc]           
  |----------------------------------------------------------*/ 
typedef struct ITS_Per4Mat_
{
    int n;                  
    int nB; 
    int symperm;

    /*   LU factors  */
    ITS_SparMat *L;
    ITS_SparMat *U;

    /* E, F blocks   */
    ITS_SparMat *E;
    ITS_SparMat *F;
    int *rperm;       /* row permutation         */ 
    int *perm;        /* col. permutation        */ 
    double *D1 ;      /* diagonal scaling row    */  
    double *D2 ;      /* diagonal scaling columns*/  
    double *wk;       /* work array              */

    /* pointer to next and previous struct         */
    struct ITS_Per4Mat_ *prev; 
    struct ITS_Per4Mat_ *next;

} ITS_Per4Mat; 

/*------------------------------------------------------------
  | struct for storing data related to the last schur complement 
  | we need to store the C matrix associated with the last block
  | and the ILUT factorization of the related Schur complement.
  | 
  | n       = size of C block = size of Schur complement
  | C       = C block of last level matrix. 
  | L, U    = ILU factors of last schur complement. 
  |
  | meth[4] = parameters for defining variants in factorization 
  |           - see function readin for details
  | rperm    = row permutation used for very nonsymmetric matrices 
  |            [such as bottleneck transversal] -- NOT IN THIS VERSION
  | perm2     = unsymmetric permutation (columns) - used primarily
  |           for the ILUTP version of ILUT/.. 
  | D1, D2  = diagonal matrices (left, right) used for scaling
  |           if scaling option is turned on. Note that the 
  |           method works by scaling the whole matrix first
  |           (at any level) before anything else is done. 
  | wk     = a work vector of length n needed for various tasks
  |            [reduces number of calls to malloc]           
  |-----------------------------------------------------------*/
typedef struct ITS_ILUTSpar_
{
    int n;                  

    /*-------------------- C matrix of last block */
    ITS_SparMat *C;

    /* LU factorization       */
    ITS_SparMat *L;
    ITS_SparMat *U;

    int *rperm;   /* row single-sinded permutation */
    int *perm;    /* column perm .                */
    int *perm2;   /* column permutation coming from pivoting in ILU */ 
    double *D1;
    double *D2;
    double *wk;

} ITS_ILUTSpar;

/* this is the arms preconditioner struct 
   | it consists of a linked list of p4mat structs
   | and the ILUt factorization (in the  form of an 
   | ITS_ILUTSpar struct  
   |---------------------------------------------- */
typedef struct ITS_ARMSpar_
{
    int n;                   /* dimension of matrix */
    int nlev;                /* number of levels    */
    ITS_ILUTSpar * ilus;     /* ILU for last level  */
    ITS_Per4Mat *levmat;     /* level structure     */

} ITS_ARMSpar;

typedef struct ITS_CompressType
{
    int grp;   /* -1: begin new group, >=0: belong to grp-th row */
    int count; /* block size, valid only if grp = -1 */

} ITS_CompressType;

/*-------------------- 3 types of matrices so far */
typedef struct ITS_SMat
{
    int n; 
    int Mtype;             /*--  type 1 = CSR, 2 = VBCSR, 3 = LDU    */
    ITS_SparMat *CS;       /* place holder for a CSR/CSC type matrix */
    ITS_ILUSpar *LDU;      /* struct for an LDU type matrix          */
    ITS_VBSparMat *VBCSR;  /* place holder for a block matrix        */
    void (*matvec)(struct ITS_SMat*, double *, double *);

} ITS_SMat;

/* types of pc */
typedef enum ITS_PC_TYPE_
{
    ITS_PC_NONE,
    ITS_PC_ARMS,
    ITS_PC_ILUK,
    ITS_PC_ILUT,
    ITS_PC_ILUC,
    ITS_PC_VBILUK,
    ITS_PC_VBILUT,

} ITS_PC_TYPE;

/* types of solver */
typedef enum ITS_SOLVER_TYPE_
{
    ITS_SOLVER_FGMRES,
    ITS_SOLVER_BICGSTAB,
    ITS_SOLVER_BICGSTABL,

} ITS_SOLVER_TYPE;

typedef struct ITS_PC
{
    ITS_PC_TYPE pc_type;

    ITS_ILUSpar *ILU;      /* struct for an ILU type preconditioner */
    ITS_ARMSpar *ARMS;     /* struct for a block preconditioner */

    ITS_VBILUSpar *VBILU;  /* struct for a block preconditioner */
    int *perm;

    int (*precon) (double *, double *, struct ITS_PC *); 
    FILE *log;

} ITS_PC;

typedef struct ITS_PARS_
{
    /* parameters from inputs -----------------------------------------*/
    int bgsl;                   /* parameter for BiCGSTAB(l)       */
    int restart;                /* Dim of Krylov subspace [fgmr]   */
    int maxits;                 /* maximum number of fgmres iters  */
    double tol;                 /* tolerance for stopping fgmres   */
    double eps;                 /* not available in Hash-based algorithm.  <= 1.  indicating
                                   how close are two rows or columns which can be grouped in
                                   the same block. */

    int ilut_p;                 /* ilut fill                    */
    double ilut_tol;            /* ilut drop tolerance          */
    int iluk_level;             /* level of fill for ILUK  */

    /* vbilu value always set to 1           */
    int perm_type;               /* indset perms (0) or PQ perms (1)*/
    int Bsize;                   /* block size - dual role */

    int diagscal;
    double tolind;
    int lfil_arr[7];
    double droptol[7], dropcoef[7];
    int ipar[18];

    FILE *fp;
    int verb;

} ITS_PARS;

typedef struct ITS_SOLVER_
{
    ITS_SOLVER_TYPE s_type;
    ITS_CooMat *A;

    /* internal mat */
    ITS_SMat smat;           /* Matrix structure for matvecs    */
    ITS_SparMat *csmat;

    ITS_PC_TYPE pc_type;
    ITS_PC pc;               /* general precond structure       */

    ITS_PARS pars;

    FILE *log;
    int nits;
    double res;
    int assembled;

} ITS_SOLVER;

#endif

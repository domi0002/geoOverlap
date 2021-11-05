/*-----------------------------------------------------------------------*
 * main test driver for VBILUK                                           *
 *-----------------------------------------------------------------------*
 * Na Li, Aug 26, 2001  -- Y.Saad 07/05                                  *
 *                                                                       *
 * Report bugs / send comments to: saad@cs.umn.edu, nli@cs.umn.edu       *
 *-----------------------------------------------------------------------*/
#include "itsol.h"

int main(void)
{
    int ierr = 0;
    /*-------------------------------------------------------------------
     * options
     *-----------------------------------------------------------------*/
    int lfil, nBlock, *nB = NULL, *perm = NULL;

    /*-------------------- main structs and wraper structs.   */
    ITS_SparMat *csmat = NULL;         /* matrix in csr formt             */
    ITS_VBSparMat *vbmat = NULL;
    ITS_VBILUSpar *lu = NULL;          /* vbilu preconditioner structure  */
    ITS_SMat *MAT;                     /* Matrix structure for matvecs    */
    ITS_PC *PRE;                       /* general precond structure       */
    double *sol = NULL, *x = NULL, *prhs = NULL, *rhs = NULL;

    int n, nnz;

    /*-------------------- IO */
    ITS_PARS io;
    int i;
    double terr, norm;
    ITS_CooMat A;
    int its;

    MAT = (ITS_SMat *) itsol_malloc(sizeof(ITS_SMat), "main:MAT");
    PRE = (ITS_PC *) itsol_malloc(sizeof(ITS_PC), "main:PRE");


    /*------------------ set parameters and other inputs  */
    itsol_solver_init_pars(&io);

    /* ------------------- Read in matrix and allocate memory */
    csmat = (ITS_SparMat *) itsol_malloc(sizeof(ITS_SparMat), "main");

    /*-------------------- case: COO formats */
    A = itsol_read_coo("pores3.coo");
    n = A.n;
    nnz = A.nnz;

    /*-------------------- conversion from COO to CSR format */
    if ((ierr = itsol_COOcs(n, nnz, A.ma, A.ja, A.ia, csmat)) != 0) {
        printf("mainARMS: COOcs error\n");
        return ierr;
    }

    /*------------------------------------------------------------*/
    x = (double *)itsol_malloc(n * sizeof(double), "main");
    rhs = (double *)itsol_malloc(n * sizeof(double), "main");
    sol = (double *)itsol_malloc(n * sizeof(double), "main");
    prhs = (double *)itsol_malloc(n * sizeof(double), "main");

    ierr = itsol_init_blocks(csmat, &nBlock, &nB, &perm, io.eps);

    /*--------------- permutes the rows and columns of the matrix */
    if (itsol_dpermC(csmat, perm) != 0) {
        printf("*** dpermC error ***\n");
        exit(9);
    }

    for( i = 0; i < n; i++ ) {
        rhs[i] = i;
    }

    /*-------------------- permute the right-hand-side  */
    for (i = 0; i < n; i++) prhs[perm[i]] = rhs[i];

    /*-------------------- convert to block matrix. */
    vbmat = (ITS_VBSparMat *) itsol_malloc(sizeof(ITS_VBSparMat), "main");
    ierr = itsol_csrvbsrC(1, nBlock, nB, csmat, vbmat);

    if (ierr != 0) {
        printf("*** in csrvbsr ierr != 0 ***\n");
        exit(10);
    }

    /*---------------------------*/
    lfil = io.iluk_level;

    lu = (ITS_VBILUSpar *) itsol_malloc(sizeof(ITS_VBILUSpar), "main");
    printf("begin vbiluk\n");

    /*-------------------- call VBILUK preconditioner set-up  */
    ierr = itsol_pc_vbilukC(lfil, vbmat, lu, stdout);

    /*-------------------- initial guess */
    for( i = 0; i < n; i++ ) {
        x[i] = 0.0;
    }

    /*-------------------- set up the structs before calling itsol_solver_fgmres */
    MAT->n = n;
    MAT->CS = csmat;
    MAT->matvec = itsol_matvecCSR;
    PRE->VBILU = lu;
    PRE->precon = itsol_preconVBR;

    /*-------------------- call itsol_solver_fgmres */
    itsol_solver_fgmres(MAT, PRE, prhs, x, io, &its, NULL);
    printf("solver converged in %d steps...\n\n", its);

    /*---------------------- calculate residual norm */
    itsol_matvec(csmat, x, sol);

    /* error */
    terr = 0.0;
    norm = 0.;
    for (i = 0; i < A.n; i++) {
        terr += (prhs[i] - sol[i]) * (prhs[i] - sol[i]);

        norm += prhs[i] * prhs[i];
    }

    printf("residual: %e, relative residual: %e\n\n", sqrt(terr), sqrt(terr / norm));

    itsol_cleanVBILU(lu);

    itsol_cleanCS(csmat);
    itsol_cleanCOO(&A);
    itsol_cleanVBMat(vbmat);
    free(nB);
    free(perm);
    free(sol);
    free(x);
    free(prhs);
    free(rhs);

    free(MAT);
    free(PRE);

    return 0;
}

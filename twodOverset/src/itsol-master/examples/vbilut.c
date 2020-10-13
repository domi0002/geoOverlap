/*-----------------------------------------------------------------*
 * main test driver for VBILUT                                     *
 *-----------------------------------------------------------------*
 * Na Li, Aug 26, 2001 -- YS 2005                                  *
 *                                                                 *
 * Report bugs / send comments to: saad@cs.umn.edu, nli@cs.umn.edu *
 *-----------------------------------------------------------------*/
#include "itsol.h"

int main(void)
{
    int ierr = 0;

    ITS_SparMat *csmat = NULL;        /* matrix in csr formt             */
    ITS_VBSparMat *vbmat = NULL;
    ITS_VBILUSpar *lu = NULL;         /* vbilu preconditioner structure  */
    ITS_SMat *MAT;                    /* Matrix structure for matvecs    */
    ITS_PC *PRE;                      /* general precond structure       */
    double *sol = NULL, *x = NULL, *prhs = NULL, *rhs = NULL;

    /*---------------------------------------------------------*/
    int n, nnz;
    ITS_BData *w = NULL;
    int lfil, max_blk_sz = ITS_MAX_BLOCK_SIZE * ITS_MAX_BLOCK_SIZE * sizeof(double);
    int nBlock, *nB = NULL, *perm = NULL;
    double tol;

    int i;
    double terr, norm;
    ITS_PARS io;
    ITS_CooMat A;
    int its;

    MAT = (ITS_SMat *) itsol_malloc(sizeof(ITS_SMat), "main:MAT");
    PRE = (ITS_PC *) itsol_malloc(sizeof(ITS_PC), "main:PRE");

    /*------------------ set parameters and other inputs  */
    itsol_solver_init_pars(&io);

    /* ------------------- Read in matrix and allocate memory-------- */
    csmat = (ITS_SparMat *) itsol_malloc(sizeof(ITS_SparMat), "main");

    /*-------------------- case: COO formats */
    A = itsol_read_coo("pores3.coo");
    n = A.n;
    nnz = A.nnz;

    /*-------------------- conversion from COO to CSR format */
    if ((ierr = itsol_COOcs(n, nnz, A.ma, A.ja, A.ia, csmat)) != 0) {
        fprintf(stderr, "mainARMS: COOcs error\n");
        return ierr;
    }

    /*------------------------------------------------------------*/
    sol = (double *)itsol_malloc(A.n * sizeof(double), "main");

    ierr = itsol_init_blocks(csmat, &nBlock, &nB, &perm, io.eps);

    if (ierr != 0) {
        printf("*** in init_blocks ierr != 0 ***\n");
        exit(8);
    }

    /*------------- permutes the rows and columns of the matrix */
    if (itsol_dpermC(csmat, perm) != 0) {
        printf("*** dpermC error ***\n");
        exit(9);
    }

    /*------------- permutes right hand side */
    prhs = (double *)itsol_malloc(n * sizeof(double), "main");
    rhs = (double *)itsol_malloc(n * sizeof(double), "main");
    x = (double *)itsol_malloc(n * sizeof(double), "main");
    for (i = 0; i < n; i++) rhs[i] = i;

    for (i = 0; i < n; i++) prhs[perm[i]] = rhs[i];

    /*-------------------- convert to block matrix. */
    vbmat = (ITS_VBSparMat *) itsol_malloc(sizeof(ITS_VBSparMat), "main");
    ierr = itsol_csrvbsrC(1, nBlock, nB, csmat, vbmat);

    if (ierr != 0) {
        printf("*** in csrvbsr ierr != 0 ***\n");
        exit(10);
    }

    /*---------------------------*/
    lfil = io.ilut_p;
    tol = io.ilut_tol;
    w = (ITS_BData *) itsol_malloc(vbmat->n * sizeof(ITS_BData), "main");

    for (i = 0; i < vbmat->n; i++) {
        w[i] = (double *)itsol_malloc(max_blk_sz, "main");
    }

    lu = (ITS_VBILUSpar *) itsol_malloc(sizeof(ITS_VBILUSpar), "main");

    /*-------------------- call VBILUT preconditioner set-up  */
    ierr = itsol_pc_vbilutC(vbmat, lu, lfil, tol, w, stdout);

    /*-------------------- initial guess */
    for( i = 0; i < A.n; i++ ) x[i] = 0.0;


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

    for (i = 0; i < vbmat->n; i++) free(w[i]);
    free(w);

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

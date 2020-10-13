/*-------------------------------------------------------------------*
 * main test driver for ILUT                                         *
 *-------------------------------------------------------------------*
 * Na Li, Oct 31, 2001                                               *
 *                                                                   *
 * Report bugs / send comments to: saad@cs.umn.edu, nli@cs.umn.edu   *
 *-------------------------------------------------------------------*/
#include "itsol.h"

int main(void)
{
    int ierr = 0;

    int lfil;
    double tol;

    /*-------------------- main structs and wraper structs.   */
    ITS_SparMat *csmat = NULL;         /* matrix in csr formt             */
    ITS_SMat *MAT;                /* Matrix structure for matvecs    */
    ITS_PC *PRE;                /* General precond structure       */
    ITS_ILUSpar *lu = NULL;           /* ilu preconditioner structure    */
    double *sol = NULL, *x = NULL, *rhs = NULL;

    int n, nnz;
    ITS_PARS io;
    int i;
    double terr, norm;
    ITS_CooMat A;
    int its;

    MAT = (ITS_SMat *) itsol_malloc(sizeof(ITS_SMat), "main:MAT");
    PRE = (ITS_PC *) itsol_malloc(sizeof(ITS_PC), "main:PRE");

    /*------------------ set parameters and other inputs  */
    itsol_solver_init_pars(&io);

    /*-------------------- Read matrix */
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

    /*---------------------------------------------------------*/
    x = (double *)itsol_malloc(n * sizeof(double), "main");
    rhs = (double *)itsol_malloc(n * sizeof(double), "main");
    sol = (double *)itsol_malloc(n * sizeof(double), "main");

    /*-------------------- set initial lfil and tol */
    lfil = io.ilut_p;
    tol = io.ilut_tol;

    lu = (ITS_ILUSpar *) itsol_malloc(sizeof(ITS_ILUSpar), "main");

    /*-------------------- call ILUT preconditioner set-up  */
    ierr = itsol_pc_ilut(csmat, lu, lfil, tol, stdout);

    /*-------------------- initial guess */
    for (i = 0; i < n; i++) {
        rhs[i] = i;
        x[i] = 0.0;
    }

    /*-------------------- set up the structs before calling itsol_solver_fgmres */
    MAT->n = n;
    MAT->CS = csmat;
    MAT->matvec = itsol_matvecCSR;
    PRE->ILU = lu;
    PRE->precon = itsol_preconILU;

    /*-------------------- call itsol_solver_fgmres */
    itsol_solver_fgmres(MAT, PRE, rhs, x, io, &its, NULL);

    printf("solver converged in %d steps...\n\n", its);

    /*-------------------- calculate residual norm */
    itsol_matvec(csmat, x, sol);

    /* error */
    terr = 0.0;
    norm = 0.;
    for (i = 0; i < A.n; i++) {
        terr += (rhs[i] - sol[i]) * (rhs[i] - sol[i]);

        norm += rhs[i] * rhs[i];
    }

    printf("residual: %e, relative residual: %e\n\n", sqrt(terr), sqrt(terr / norm));

    itsol_cleanILU(lu);

    /*-------------------- Test with next matrix   */
    itsol_cleanCS(csmat);
    itsol_cleanCOO(&A);

    free(sol);
    free(x);
    free(rhs);

    free(MAT);
    free(PRE);

    return 0;
}

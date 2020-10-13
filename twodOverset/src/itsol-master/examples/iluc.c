/*-------------------------------------------------------------------*
 * main test driver for Crout version of ILU
 *-------------------------------------------------------------------*
 * Na Li, Mar 25, 2002                                               *
 *                                                                   *
 * Report bugs / send comments to: saad@cs.umn.edu, nli@cs.umn.edu   *
 *-------------------------------------------------------------------*/
#include "itsol.h"

#define DRP_MTH 0               /* drop method see ilutc code */

int main(void)
{
    int ierr = 0;

    int dropmthd = DRP_MTH, nnz;
    int pattern_symm = 0;

    /*-------------------- main structs and wraper structs.     */
    ITS_SparMat *csmat = NULL, *Ac;     /* matrix in csr formt             */
    ITS_SMat *MAT = NULL;               /* Matrix structure for matvecs    */
    ITS_PC *PRE = NULL;                 /* General precond structure       */
    ITS_ILUSpar *lumat = NULL;          /* ilu preconditioner structure    */
    ITS_ILUSpar *lu = NULL;              /* a temporary lu matrix           */
    double *sol = NULL, *x = NULL, *rhs = NULL;

    /*-------------------- method for incrementing lfil is set here */
    int lfil;
    double tol;
    /*-------------------- harwell boeing temp. arrays */
    int rsa = 0;
    int n;

    ITS_PARS io;
    ITS_CooMat A;
    int i;
    double terr, norm;
    int its;

    MAT = (ITS_SMat *) itsol_malloc(sizeof(ITS_SMat), "main:MAT");
    PRE = (ITS_PC *) itsol_malloc(sizeof(ITS_PC), "main:PRE");

    /*------------------ set parameters and other inputs  */
    itsol_solver_init_pars(&io);

    /*-------------------- Read  matrix */
    lumat = (ITS_ILUSpar *) itsol_malloc(sizeof(ITS_ILUSpar), "main:lumat");
    csmat = (ITS_SparMat *) itsol_malloc(sizeof(ITS_SparMat), "main:csmat");

    /*-------------------- case: COO formats */
    A = itsol_read_coo("pores3.coo");
    n = A.n;
    nnz = A.nnz;

    if ((ierr = itsol_COOcs(n, nnz, A.ma, A.ia, A.ja, csmat)) != 0) {
        printf("mainILUC: COOcs error\n");
        return ierr;
    }

    /*-------------------- convert to lum format for iluc + symmetrize */
    if (rsa == 0 && pattern_symm) rsa = 2;

    if ((ierr = itsol_CSClumC(csmat, lumat, rsa)) != 0) {
        printf(" error: CSClum error\n");
        return (ierr);
    }

    /*---------------------------------------------------------*/
    x = (double *)itsol_malloc(n * sizeof(double), "main");
    rhs = (double *)itsol_malloc(n * sizeof(double), "main");
    sol = (double *)itsol_malloc(n * sizeof(double), "main");

    /*-------------------- set initial lfil and tol */
    lfil = io.ilut_p;
    tol = io.ilut_tol;

    lu = (ITS_ILUSpar *) itsol_malloc(sizeof(ITS_ILUSpar), "main");

    /*-------------------- call ILUC preconditioner set-up  */
    ierr = itsol_pc_ilutc(lumat, lu, lfil, tol, dropmthd, stdout);

    for (i = 0; i < n; i++) {
        rhs[i] = i;
        x[i] = 0.0;
    }

    /*-------------------- set up the structs before calling itsol_solver_fgmres */
    MAT->n = n;
    MAT->CS = csmat;    /* in column format */
    MAT->matvec = itsol_matvecCSC;    /* column matvec */
    PRE->ILU = lu;
    PRE->precon = itsol_preconLDU;

    /*-------------------- call itsol_solver_fgmres */
    itsol_solver_fgmres(MAT, PRE, rhs, x, io, &its, NULL);

    printf("solver converged in %d steps...\n\n", its);

    /*-------------------- calculate residual norm */
    {
        Ac = (ITS_SparMat *) itsol_malloc(sizeof(ITS_SparMat), "main");

        if ((ierr = itsol_COOcs(n, nnz, A.ma, A.ja, A.ia, Ac)) != 0) {
            fprintf(stderr, "mainARMS: COOcs error\n");
            return ierr;
        }

        itsol_matvec(Ac, x, sol);
        itsol_cleanCS(Ac);

        /* error */
        terr = 0.0;
        norm = 0.;
        for (i = 0; i < A.n; i++) {
            terr += (rhs[i] - sol[i]) * (rhs[i] - sol[i]);

            norm += rhs[i] * rhs[i];
        }

        printf("residual: %e, relative residual: %e\n\n", sqrt(terr), sqrt(terr / norm));
    }

    itsol_cleanILU(lu);

    /*-------------------- NEXT_MAT: */
    itsol_cleanCS(csmat);
    itsol_cleanCOO(&A);
    itsol_cleanILU(lumat);
    free(sol);
    free(x);
    free(rhs);

    free(MAT);
    free(PRE);
    return 0;
}

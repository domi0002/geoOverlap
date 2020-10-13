/*-----------------------------------------------------------------*
 * main test driver for the ARMS2 preconditioner for 
 * Matrices in the COO/Harwell Boeing format
 *-----------------------------------------------------------------*
 * Yousef Saad - Aug. 2005.                                        *
 *                                                                 *
 * Report bugs / send comments to: saad@cs.umn.edu                 *
 *-----------------------------------------------------------------*/
#include "itsol.h"

int main(void)
{
    int ierr = 0;

    int nnz = 0;

    /*-------------------- main structs and wraper structs.     */
    ITS_SparMat *csmat = NULL;         /* matrix in csr formt             */
    ITS_ARMSpar *ArmsSt = NULL;         /* arms preconditioner structure   */
    ITS_SMat *MAT = NULL;         /* Matrix structure for matvecs    */
    ITS_PC *PRE = NULL;         /* general precond structure       */
    double *sol = NULL, *x = NULL, *rhs = NULL;

    int n;

    ITS_PARS io;                /* structure for handling I/O 
                                   functions + a few other things */
    int i;
    double terr, norm;
    ITS_CooMat A;
    int its;

    MAT = (ITS_SMat *) itsol_malloc(sizeof(ITS_SMat), "main:MAT");
    PRE = (ITS_PC *) itsol_malloc(sizeof(ITS_PC), "main:PRE");

    /*------------------ set parameters and other inputs  */
    itsol_solver_init_pars(&io);

    csmat = (ITS_SparMat *) itsol_malloc(sizeof(ITS_SparMat), "main:csmat");

    /*-------------------- case: COO formats */
    A = itsol_read_coo("pores3.coo");
    n = A.n;
    nnz = A.nnz;

    /*-------------------- conversion from COO to CSR format */
    if ((ierr = itsol_COOcs(n, nnz, A.ma, A.ja, A.ia, csmat)) != 0) {
        printf("mainARMS: COOcs error\n");
        return ierr;
    }

    /*----------------------------------------------------------*/
    x = (double *)itsol_malloc(n * sizeof(double), "main:x");
    rhs = (double *)itsol_malloc(n * sizeof(double), "main");
    sol = (double *)itsol_malloc(n * sizeof(double), "main");

    ArmsSt = (ITS_ARMSpar *) itsol_malloc(sizeof(ITS_ARMSpar), "main:ArmsSt");
    itsol_setup_arms(ArmsSt);
    printf("begin arms\n");

    /*-------------------- call ARMS preconditioner set-up  */
    ierr = itsol_pc_arms2(csmat, io.ipar, io.droptol, io.lfil_arr, io.tolind, ArmsSt, stdout);

    /*-------------------- initial guess */
    for(i=0; i < n; i++) {
        rhs[i] = i;
        x[i] = 0.0;
    }

    /*-------------------- set up the structs before calling itsol_solver_fgmres */
    MAT->n = n;
    MAT->CS = csmat;
    MAT->matvec = itsol_matvecCSR;
    PRE->ARMS = ArmsSt;
    PRE->precon = itsol_preconARMS;

    /*-------------------- call itsol_solver_fgmres */
    itsol_solver_fgmres(MAT, PRE, rhs, x, io, &its, NULL);
    printf("solver converged in %d steps...\n", its);

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

    itsol_cleanARMS(ArmsSt);
    itsol_cleanCS(csmat);
    itsol_cleanCOO(&A);

    free(sol);
    free(x);
    free(rhs);
    free(MAT);
    free(PRE);

    return 0;
}

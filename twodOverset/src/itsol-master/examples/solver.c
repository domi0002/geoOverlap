/*-------------------------------------------------------------------*
 * Report bugs / send comments to: saad@cs.umn.edu                   *
 *-------------------------------------------------------------------*/

#include "itsol.h"

int main(void)
{

    double *sol = NULL, *x = NULL, *rhs = NULL;
    int n, nnz;
    int i, ierr;
    double terr, norm;
    ITS_CooMat A;
    int its;
    ITS_SparMat *csmat = NULL;
    ITS_SOLVER s;
    ITS_PARS io;

    /* case: COO formats */
    A = itsol_read_coo("pores3.coo");
    n = A.n;
    nnz = A.nnz;

    /* solution vectors */
    x = (double *)itsol_malloc(n * sizeof(double), "main");
    rhs = (double *)itsol_malloc(n * sizeof(double), "main");
    sol = (double *)itsol_malloc(n * sizeof(double), "main");

    for( i = 0; i < n; i++ ) x[i] = 0.0;
    for (i = 0; i < n; i++) rhs[i] = i;

    /* create */
    itsol_solver_initialize(&s, ITS_SOLVER_FGMRES, ITS_PC_ILUK, &A);

    /* tune parameters, optional */
    itsol_solver_init_pars(&io);
    io.tol = 1e-9;
    io.iluk_level = 3;

    /* override default parameters */
    itsol_solver_set_pars(&s, io);

    /* assemble, optional */
#if 0
    itsol_solver_assemble(&s);
#endif

    /* call itsol_solver_fgmres */
    itsol_solver_solve(&s, x, rhs);

    /* get results */
    its = s.nits;
    printf("solver converged in %d steps...\n\n", its);

    /* calculate residual norm, optional */
    {
        csmat = (ITS_SparMat *) itsol_malloc(sizeof(ITS_SparMat), "main");
        if ((ierr = itsol_COOcs(n, nnz, A.ma, A.ja, A.ia, csmat)) != 0) {
            fprintf(stderr, "mainARMS: COOcs error\n");
            return ierr;
        }

        itsol_matvec(csmat, x, sol);

        /* error */
        terr = 0.0;
        norm = 0.;
        for (i = 0; i < A.n; i++) {
            terr += (rhs[i] - sol[i]) * (rhs[i] - sol[i]);

            norm += rhs[i] * rhs[i];
        }

        printf("residual: %e, relative residual: %e\n\n", sqrt(terr), sqrt(terr / norm));
    }

    itsol_solver_finalize(&s);
    itsol_cleanCOO(&A);
    itsol_cleanCS(csmat);

    free(sol);
    free(x);
    free(rhs);

    return 0;
}

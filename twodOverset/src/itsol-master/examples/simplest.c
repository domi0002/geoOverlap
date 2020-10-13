#include "itsol.h"

int main(void)
{
    double *x = NULL, *rhs = NULL;
    int n, i;
    ITS_CooMat A;
    ITS_SOLVER s;

    /* case: COO formats */
    A = itsol_read_coo("pores3.coo");
    n = A.n;

    /* solution vectors */
    x = (double *)itsol_malloc(n * sizeof(double), "main");
    rhs = (double *)itsol_malloc(n * sizeof(double), "main");

    for( i = 0; i < n; i++ ) {
        x[i] = 0.0;
        rhs[i] = i;
    }

    /* init */
    itsol_solver_initialize(&s, ITS_SOLVER_BICGSTAB, ITS_PC_ILUK, &A);

    /* call solver */
    itsol_solver_solve(&s, x, rhs);

    /* get results */
    printf("solver converged in %d steps...\n\n", s.nits);

    /* cleanup */
    itsol_solver_finalize(&s);

    itsol_cleanCOO(&A);
    free(x);
    free(rhs);

    return 0;
}

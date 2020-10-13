
#include "solver-bicgstabl.h"

int itsol_solver_bicgstabl(ITS_SMat *Amat, ITS_PC *lu, double *rhs, double *x, ITS_PARS io,
        int *nits, double *res)
{
    int iter;

    double *rtld, *bp, *t, *tp, *xp, **r, **u;
    double *tau, *gamma, *gamma1, *gamma2;
    double *sigma;
    double alpha, beta, omega, rho0, rho1;
    double ires, nrm2;
    double nu;
    int z_dim;
    int l, i, j, n, k;
    int end_solve;
    int bgsl = io.bgsl;
    double tol = io.tol;
    int maxits = io.maxits;
    FILE * fp = io.fp;

    assert(x != NULL);
    assert(rhs != NULL);
    assert(Amat != NULL);

    l = bgsl;
    if (l <= 0) l = 4;

    z_dim = l + 1;
    n = Amat->n;
    rtld = itsol_malloc(n * sizeof(double), "bicgstabl");
    xp = itsol_malloc(n * sizeof(double), "bicgstabl");
    bp = itsol_malloc(n * sizeof(double), "bicgstabl");
    t = itsol_malloc(n * sizeof(double), "bicgstabl");
    tp = itsol_malloc(n * sizeof(double), "bicgstabl");

    r = itsol_malloc(sizeof(*r) * (l + 1), "bicgstabl");
    for (i = 0; i <= l; i++) r[i] = itsol_malloc(n * sizeof(double), "bicgstabl");

    u = itsol_malloc(sizeof(*u) * (l + 1), "bicgstabl");
    for (i = 0; i <= l; i++) u[i] = itsol_malloc(n * sizeof(double), "bicgstabl");

    tau = itsol_malloc(sizeof(*tau) * z_dim * (4 + l + 1), "bicgstabl");
    gamma = &tau[z_dim * z_dim];
    gamma1 = &gamma[z_dim];
    gamma2 = &gamma1[z_dim];
    sigma = &gamma2[z_dim];

    /* set terminate tol */
    Amat->matvec(Amat, x, tp);
    for (i = 0; i < n; i++) r[0][i] = rhs[i] - tp[i];

    itsol_copy(rtld, r[0], n);
    itsol_copy(bp, r[0], n);
    itsol_copy(xp, x, n);

    for (i = 0; i < n; i++) u[0][i] = 0.;

    nrm2 = ires = itsol_norm2(r[0], n);

    end_solve = 0;
    if (nrm2 == 0.) {
        iter = 0;
        end_solve = 1;
    }

    if (end_solve) goto end;

    /* init */
    tol = fabs(tol) * nrm2;
    alpha = 0.0;
    omega = 1.0;
    rho0 = 1.0;

    /* mail loop */
    for (iter = 0; iter < maxits; iter++) {
        /* BiCG PART, rho0 = -w*rho0 */
        rho0 = -omega * rho0;

        for (j = 0; j < l; j++) {
            iter++;

            /* rho1 = <rtld,r[j]> */
            rho1 = itsol_dot(rtld, r[j], n);

            /* test breakdown */
            if (fabs(rho1) == 0.0) {
                for (i = 0; i < n; i++) t[i] = 0.;

                /*  pc */
                if (lu == NULL) {
                    memcpy(t, x, n * sizeof(double));
                }
                else {
                    lu->precon(x, t, lu);
                }

                itsol_copy(x, t, n);

                for (i = 0; i < n; i++) x[i] += xp[i];

                end_solve = 1;
            }

            if (end_solve) goto end;

            /* beta = alpha * (rho1/rho0) */
            /* rho0 = rho1                */
            beta = alpha * (rho1 / rho0);
            rho0 = rho1;

            /* u[i] = r[i] - beta*u[i] (i=0,j) */
            for (i = 0; i <= j; i++) {
                for (k = 0; k < n; k++) u[i][k] = -beta * u[i][k] + r[i][k];
            }

            /* u[j+1] = A    * u[j]   */
            /* u[j+1] = M^-1 * u[j+1] */
            for (k = 0; k < n; k++) t[k] = 0.;

            if (lu == NULL) {
                memcpy(t, u[j], n * sizeof(double));
            }
            else {
                lu->precon(u[j], t, lu);
            }

            Amat->matvec(Amat, t, u[j + 1]);

            /* nu = <rtld, u[j+1]> */
            nu = itsol_dot(rtld, u[j + 1], n);

            /* test breakdown */
            if (fabs(nu) == 0.0) {
                for (k = 0; k < n; k++) t[k] = 0.;

                /*  pc */
                if (lu == NULL) {
                    memcpy(t, x, n * sizeof(double));
                }
                else {
                    lu->precon(x, t, lu);
                }

                itsol_copy(x, t, n);
                for (k = 0; k < n; k++) x[k] = x[k] + xp[k];
                end_solve = 1;
            }

            /* alpha = rho1 / nu */
            alpha = rho1 / nu;

            /* x = x + alpha*u[0] */
            for (k = 0; k < n; k++) x[k] = x[k] + alpha * u[0][k];

            /* r[i] = r[i] - alpha*u[i+1] (i=0,j) */
            for (i = 0; i <= j; i++) {
                for (k = 0; k < n; k++) r[i][k] = r[i][k] - alpha * u[i + 1][k];
            }

            nrm2 = itsol_norm2(r[0], n);
            if (io.verb > 0 && fp != NULL) fprintf(fp, "%8d   %10.2e\n", iter, nrm2 / ires);

            if (nrm2 <= tol) {
                for (k = 0; k < n; k++) t[k] = 0.;

                /*  pc */
                if (lu == NULL) {
                    memcpy(t, x, n * sizeof(double));
                }
                else {
                    lu->precon(x, t, lu);
                }

                itsol_copy(x, t, n);
                for (k = 0; k < n; k++) x[k] += xp[k];

                end_solve = 1;
            }

            if (end_solve) goto end;

            /* r[j+1] = A    * r[j]   */
            /* r[j+1] = M^-1 * r[j+1] */
            for (k = 0; k < n; k++) t[k] = 0.;

            /*  pc */
            if (lu == NULL) {
                memcpy(t, r[j], n * sizeof(double));
            }
            else {
                lu->precon(r[j], t, lu);
            }

            Amat->matvec(Amat, t, r[j + 1]);
        }

        /* MR PART */
        for (j = 1; j <= l; j++) {
            for (i = 1; i <= j - 1; i++) {
                nu = itsol_dot(r[j], r[i], n);
                nu = nu / sigma[i];
                tau[i * z_dim + j] = nu;

                for (k = 0; k < n; k++) r[j][k] = r[j][k] - nu * r[i][k];
            }

            sigma[j] = itsol_dot(r[j], r[j], n);
            nu = itsol_dot(r[0], r[j], n);
            gamma1[j] = nu / sigma[j];
        }

        gamma[l] = gamma1[l];
        omega = gamma[l];
        for (j = l - 1; j >= 1; j--) {
            nu = 0.0;
            for (i = j + 1; i <= l; i++) {
                nu += tau[j * z_dim + i] * gamma[i];
            }
            gamma[j] = gamma1[j] - nu;
        }
        for (j = 1; j <= l - 1; j++) {
            nu = 0.0;
            for (i = j + 1; i <= l - 1; i++) {
                nu += tau[j * z_dim + i] * gamma[i + 1];
            }
            gamma2[j] = gamma[j + 1] + nu;
        }

        /* UPDATE */
        itsol_axpby(gamma[1], r[0], 1, x, n);
        itsol_axpby(-gamma1[l], r[l], 1, r[0], n);
        itsol_axpby(-gamma[l], u[l], 1, u[0], n);

        for (j = 1; j <= l - 1; j++) {
            itsol_axpby(-gamma[j], u[j], 1, u[0], n);
            itsol_axpby(gamma2[j], r[j], 1, x, n);
            itsol_axpby(-gamma1[j], r[j], 1, r[0], n);
        }

        if (io.verb > 0 && fp != NULL) fprintf(fp, "%8d   %10.2e\n", iter, nrm2 / ires);

        if (nrm2 < tol) {
            for (k = 0; k < n; k++) t[k] = 0.;

            /*  pc */
            if (lu == NULL) {
                memcpy(t, x, n * sizeof(double));
            }
            else {
                lu->precon(x, t, lu);
            }

            itsol_copy(x, t, n);
            for (k = 0; k < n; k++) x[k] += xp[k];
            end_solve = 1;
        }

        if (end_solve) break;
    }

end:
    if (iter < maxits) iter += 1;

    if (nits != NULL) *nits = iter;
    if (res != NULL) *res = nrm2 / ires;

    free(rtld);
    free(xp);
    free(bp);
    free(tp);
    free(t);

    for (i = 0; i <= l; i++) free(r[i]);
    free(r);

    for (i = 0; i <= l; i++) free(u[i]);
    free(u);
    free(tau);

    return 0;
}

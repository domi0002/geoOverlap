
#include "itsol.h"

void itsol_solver_initialize(ITS_SOLVER *s, ITS_SOLVER_TYPE stype, ITS_PC_TYPE pctype, ITS_CooMat *A)
{
    assert(s != NULL);
    assert(A != NULL);

    /* init */
    bzero(s, sizeof(*s));

    s->s_type = stype;
    s->A = A;
    s->log = stdout;
    
    /* pc */
    s->pc_type = pctype;
    s->pc.log = s->log;
    itsol_pc_initialize(&s->pc, pctype);

    /* init parameters */
    itsol_solver_init_pars(&s->pars);
}

void itsol_solver_finalize(ITS_SOLVER *s)
{
    if (s == NULL) return;

    /* cleanup */
    if (s->csmat != NULL) itsol_cleanCS(s->csmat);
    s->csmat = NULL;

    itsol_pc_finalize(&s->pc);

    bzero(s, sizeof(*s));
}

int itsol_solver_assemble(ITS_SOLVER *s)
{
    ITS_PC_TYPE pctype;
    ITS_CooMat A;
    int ierr;
    FILE *log;

    assert(s != NULL);

    if (s->assembled) return 0;

    /* log */
    if (s->log == NULL) {
        log = stdout;
    }
    else {
        log = s->log;
    }

    /* assemble */
    pctype = s->pc_type;

    s->csmat = (ITS_SparMat *) itsol_malloc(sizeof(ITS_SparMat), "solver assemble");
    A = *s->A;

    if (pctype == ITS_PC_ILUC) {
        if ((ierr = itsol_COOcs(A.n, A.nnz, A.ma, A.ia, A.ja, s->csmat)) != 0) {
            fprintf(log, "solver assemble, COOcs error\n");
            return ierr;
        }

        /* smat */
        s->smat.n = A.n;
        s->smat.CS = s->csmat;               /* in column format */
        s->smat.matvec = itsol_matvecCSC;    /* column matvec */
    }
    else if(pctype == ITS_PC_ILUK || pctype == ITS_PC_ILUT || pctype == ITS_PC_VBILUK || pctype == ITS_PC_VBILUT
            || pctype == ITS_PC_ARMS) {
        if ((ierr = itsol_COOcs(A.n, A.nnz, A.ma, A.ja, A.ia, s->csmat)) != 0) {
            fprintf(log, "mainARMS: COOcs error\n");
            return ierr;
        }

        /* smat */
        s->smat.n = A.n;
        s->smat.CS = s->csmat;               /* in row format */
        s->smat.matvec = itsol_matvecCSR;    /* row matvec */
    }
    else {
        fprintf(log, "solver assemble, wrong preconditioner type\n");
        exit(-1);
    }

    /* pc assemble */
    itsol_pc_assemble(s);

    s->assembled = 1;
    return 0;
}

int itsol_solver_solve(ITS_SOLVER *s, double *x, double *rhs)
{
    ITS_PARS io;
    ITS_PC_TYPE pctype;
    ITS_SOLVER_TYPE stype;

    assert(s != NULL);
    assert(x != NULL);
    assert(rhs != NULL);

    /* assemble */
    itsol_solver_assemble(s);

    io = s->pars;
    pctype = s->pc_type;
    stype = s->s_type;

    if (stype == ITS_SOLVER_FGMRES) {
        if (pctype == ITS_PC_ILUC || pctype == ITS_PC_ILUK || pctype == ITS_PC_ILUT || pctype == ITS_PC_ARMS) {
            return itsol_solver_fgmres(&s->smat, &s->pc, rhs, x, io, &s->nits, &s->res);
        }
        else if (pctype == ITS_PC_VBILUK || pctype == ITS_PC_VBILUT) {
            if (s->pc.perm == NULL) {
                return itsol_solver_fgmres(&s->smat, &s->pc, rhs, x, io, &s->nits, &s->res);
            }
            else {
                double *px = NULL, *prhs = NULL;
                int i, rt;

                px = (double *)itsol_malloc(s->csmat->n * sizeof(double), "main");
                prhs = (double *)itsol_malloc(s->csmat->n * sizeof(double), "main");

                for (i = 0; i < s->csmat->n; i++) {
                    prhs[s->pc.perm[i]] = rhs[i];
                    px[s->pc.perm[i]] = x[i];
                }

                rt = itsol_solver_fgmres(&s->smat, &s->pc, prhs, px, io, &s->nits, &s->res);

                for (i = 0; i < s->csmat->n; i++) {
                    rhs[i] = prhs[s->pc.perm[i]];
                    x[i] = px[s->pc.perm[i]];
                }

                free(px);
                free(prhs);

                return rt;
            }
        }
        else if (pctype == ITS_PC_NONE) {
            return itsol_solver_fgmres(&s->smat, NULL, rhs, x, io, &s->nits, &s->res);
        }
        else {
            fprintf(s->pc.log, "wrong preconditioner type\n");
            exit(-1);
        }
    }
    else if (stype == ITS_SOLVER_BICGSTAB) {
        if (pctype == ITS_PC_ILUC || pctype == ITS_PC_ILUK || pctype == ITS_PC_ILUT || pctype == ITS_PC_ARMS) {
            return itsol_solver_bicgstab(&s->smat, &s->pc, rhs, x, io, &s->nits, &s->res);
        }
        else if (pctype == ITS_PC_VBILUK || pctype == ITS_PC_VBILUT) {
            if (s->pc.perm == NULL) {
                return itsol_solver_bicgstab(&s->smat, &s->pc, rhs, x, io, &s->nits, &s->res);
            }
            else {
                double *px = NULL, *prhs = NULL;
                int i, rt;

                px = (double *)itsol_malloc(s->csmat->n * sizeof(double), "main");
                prhs = (double *)itsol_malloc(s->csmat->n * sizeof(double), "main");

                for (i = 0; i < s->csmat->n; i++) {
                    prhs[s->pc.perm[i]] = rhs[i];
                    px[s->pc.perm[i]] = x[i];
                }

                rt = itsol_solver_bicgstab(&s->smat, &s->pc, prhs, px, io, &s->nits, &s->res);

                for (i = 0; i < s->csmat->n; i++) {
                    rhs[i] = prhs[s->pc.perm[i]];
                    x[i] = px[s->pc.perm[i]];
                }

                free(px);
                free(prhs);

                return rt;
            }
        }
        else if (pctype == ITS_PC_NONE) {
            return itsol_solver_bicgstab(&s->smat, NULL, rhs, x, io, &s->nits, &s->res);
        }
        else {
            fprintf(s->pc.log, "wrong preconditioner type\n");
            exit(-1);
        }
    }
    else if (stype == ITS_SOLVER_BICGSTABL) {
        if (pctype == ITS_PC_ILUC || pctype == ITS_PC_ILUK || pctype == ITS_PC_ILUT || pctype == ITS_PC_ARMS) {
            return itsol_solver_bicgstabl(&s->smat, &s->pc, rhs, x, io, &s->nits, &s->res);
        }
        else if (pctype == ITS_PC_VBILUK || pctype == ITS_PC_VBILUT) {
            if (s->pc.perm == NULL) {
                return itsol_solver_bicgstabl(&s->smat, &s->pc, rhs, x, io, &s->nits, &s->res);
            }
            else {
                double *px = NULL, *prhs = NULL;
                int i, rt;

                px = (double *)itsol_malloc(s->csmat->n * sizeof(double), "main");
                prhs = (double *)itsol_malloc(s->csmat->n * sizeof(double), "main");

                for (i = 0; i < s->csmat->n; i++) {
                    prhs[s->pc.perm[i]] = rhs[i];
                    px[s->pc.perm[i]] = x[i];
                }

                rt = itsol_solver_bicgstabl(&s->smat, &s->pc, prhs, px, io, &s->nits, &s->res);

                for (i = 0; i < s->csmat->n; i++) {
                    rhs[i] = prhs[s->pc.perm[i]];
                    x[i] = px[s->pc.perm[i]];
                }

                free(px);
                free(prhs);

                return rt;
            }
        }
        else if (pctype == ITS_PC_NONE) {
            return itsol_solver_bicgstabl(&s->smat, NULL, rhs, x, io, &s->nits, &s->res);
        }
        else {
            fprintf(s->pc.log, "wrong preconditioner type\n");
            exit(-1);
        }
    }
    else {
        fprintf(s->log, "wrong solver type\n");
        exit(-1);
    }


    return 0;
}

void itsol_pc_initialize(ITS_PC *pc, ITS_PC_TYPE pctype)
{
    assert(pc != NULL);

    pc->pc_type = pctype;

    if (pctype == ITS_PC_ILUC || pctype == ITS_PC_ILUK || pctype == ITS_PC_ILUT) {
        pc->ILU = (ITS_ILUSpar *) itsol_malloc(sizeof(ITS_ILUSpar), "pc init");
    }
    else if (pctype == ITS_PC_VBILUK || pctype == ITS_PC_VBILUT) {
        pc->VBILU = (ITS_VBILUSpar *) itsol_malloc(sizeof(ITS_VBILUSpar), "pc init");
    }
    else if (pctype == ITS_PC_ARMS) {
        pc->ARMS = (ITS_ARMSpar *) itsol_malloc(sizeof(ITS_ARMSpar), "pc init");
    }
    else {
        fprintf(pc->log, "wrong preconditioner type\n");
        exit(-1);
    }
}

void itsol_pc_finalize(ITS_PC *pc)
{
    ITS_PC_TYPE pctype;

    if (pc == NULL) return;

    pctype = pc->pc_type;
    if (pctype == ITS_PC_ILUC || pctype == ITS_PC_ILUK || pctype == ITS_PC_ILUT) {
        itsol_cleanILU(pc->ILU);
        pc->ILU = NULL;
    }
    else if (pctype == ITS_PC_VBILUK || pctype == ITS_PC_VBILUT) {
        itsol_cleanVBILU(pc->VBILU);
        pc->VBILU = NULL;

        if (pc->perm != NULL) free(pc->perm);
        pc->perm = NULL;
    }
    else if (pctype == ITS_PC_ARMS) {
        itsol_cleanARMS(pc->ARMS);
        pc->ARMS = NULL;
    }
    else {
        fprintf(pc->log, "wrong preconditioner type\n");
        exit(-1);
    }
}

int itsol_pc_assemble(ITS_SOLVER *s)
{
    ITS_PC_TYPE pctype;
    int ierr;
    ITS_PARS p;
    ITS_PC *pc;

    assert(s != NULL);
    pc = &s->pc;

    /* type */
    pctype = pc->pc_type;
    p = s->pars;

    if (pctype == ITS_PC_ILUC) {
        pc->precon = itsol_preconLDU;
    }
    else if (pctype == ITS_PC_ILUK) {
        ierr = itsol_pc_ilukC(p.iluk_level, s->csmat, pc->ILU, pc->log);

        if (ierr != 0) {
            fprintf(pc->log, "pc assemble, ILUK error\n");
            return ierr;
        }

        pc->precon = itsol_preconILU;
    }
    else if (pctype == ITS_PC_ILUT) {
        ierr = itsol_pc_ilut(s->csmat, pc->ILU, p.ilut_p, p.ilut_tol, pc->log);

        if (ierr != 0) {
            fprintf(pc->log, "pc assemble, ILUK error\n");
            return ierr;
        }

        pc->precon = itsol_preconILU;
    }
    else if (pctype == ITS_PC_VBILUK) {
        int nBlock, *nB = NULL, *perm = NULL;
        ITS_VBSparMat *vbmat = NULL;

        /* init */
        itsol_init_blocks(s->csmat, &nBlock, &nB, &perm, p.eps);

        /* save perm */
        pc->perm = perm;

        /* permutes the rows and columns of the matrix */
        if (itsol_dpermC(s->csmat, perm) != 0) {
            fprintf(pc->log, "*** dpermC error ***\n");
            exit(9);
        }

        /*-------------------- convert to block matrix. */
        vbmat = (ITS_VBSparMat *) itsol_malloc(sizeof(ITS_VBSparMat), "main");

        ierr = itsol_csrvbsrC(1, nBlock, nB, s->csmat, vbmat);
        if (ierr != 0) {
            fprintf(pc->log, "pc assemble in csrvbsr ierr != 0 ***\n");
            exit(10);
        }

        /* fac */
        ierr = itsol_pc_vbilukC(p.iluk_level, vbmat, pc->VBILU, pc->log);
        if (ierr != 0) {
            fprintf(pc->log, "pc assemble in vbilukC ierr != 0 ***\n");
            exit(10);
        }

        pc->precon = itsol_preconVBR;

        /* cleanup */
        itsol_cleanVBMat(vbmat);
        free(nB);
    }
    else if (pctype == ITS_PC_VBILUT) {
        int nBlock, *nB = NULL, *perm = NULL, i;
        ITS_VBSparMat *vbmat = NULL;
        int lfil, max_blk_sz = ITS_MAX_BLOCK_SIZE * ITS_MAX_BLOCK_SIZE * sizeof(double);
        double tol;
        ITS_BData *w = NULL;

        /* init */
        itsol_init_blocks(s->csmat, &nBlock, &nB, &perm, p.eps);

        /* save perm */
        pc->perm = perm;

        /* permutes the rows and columns of the matrix */
        if (itsol_dpermC(s->csmat, perm) != 0) {
            fprintf(pc->log, "*** dpermC error ***\n");
            exit(9);
        }

        /*-------------------- convert to block matrix. */
        vbmat = (ITS_VBSparMat *) itsol_malloc(sizeof(ITS_VBSparMat), "main");

        ierr = itsol_csrvbsrC(1, nBlock, nB, s->csmat, vbmat);
        if (ierr != 0) {
            fprintf(pc->log, "pc assemble in csrvbsr ierr != 0 ***\n");
            exit(10);
        }

        /* fac */
        lfil = p.ilut_p;
        tol = p.ilut_tol;
        w = (ITS_BData *) itsol_malloc(vbmat->n * sizeof(ITS_BData), "main");

        for (i = 0; i < vbmat->n; i++) {
            w[i] = (double *)itsol_malloc(max_blk_sz, "main");
        }

        ierr = itsol_pc_vbilutC(vbmat, pc->VBILU, lfil, tol, w, pc->log);
        if (ierr != 0) {
            fprintf(pc->log, "pc assemble in vbilutC ierr != 0 ***\n");
            exit(10);
        }

        pc->precon = itsol_preconVBR;

        /* cleanup */
        for (i = 0; i < vbmat->n; i++) free(w[i]);
        free(w);

        itsol_cleanVBMat(vbmat);
        free(nB);
    }
    else if (pctype == ITS_PC_ARMS) {
        /* setup */
        itsol_setup_arms(pc->ARMS);

        /* assemble */
        ierr = itsol_pc_arms2(s->csmat, p.ipar, p.droptol, p.lfil_arr, p.tolind, pc->ARMS, pc->log);

        if (ierr != 0) {
            fprintf(pc->log, "pc assemble, arms error\n");
            return ierr;
        }

        pc->precon = itsol_preconARMS;
    }
    else {
        fprintf(pc->log, "wrong preconditioner type\n");
        exit(-1);
    }

    return 0;
}


void itsol_solver_set_pars(ITS_SOLVER *s, ITS_PARS par)
{
    ITS_PARS *p;

    assert(s != NULL);

    memcpy(&s->pars, &par, sizeof(par));

    /* update arms pars */
    p = &s->pars;
    itsol_set_arms_pars(p, p->diagscal, p->ipar, p->dropcoef, p->lfil_arr);
}

void itsol_solver_init_pars(ITS_PARS *p)
{
    assert(p != NULL);

    p->fp = stdout;
    p->verb = 2;

    /* parameters from inputs -----------------------------------------*/
    p->bgsl = 4;
    p->restart = 30;               /* Dim of Krylov subspace [fgmr]   */
    p->maxits = 1000;              /* maximum number of fgmres iters  */
    p->tol = 1e-6;                 /* tolerance for stopping fgmres   */

    p->eps = 0.8;
    p->ilut_p = 50;                /* initial lfil                    */
    p->ilut_tol = 1e-3;            /* initial drop tolerance          */
    p->iluk_level = 1;             /* initial level of fill for ILUK  */

    /* value always set to 1           */
    p->perm_type = 0;              /* indset perms (0) or PQ perms (1)*/
    p->Bsize = 30;                 /* block size - dual role. see input file */

    /* arms */
    p->diagscal = 1;
    p->tolind = ITS_TOL_DD;

    /* init arms pars */
    itsol_set_arms_pars(p, p->diagscal, p->ipar, p->dropcoef, p->lfil_arr);
}


#include "pc-vbilut.h"

/*----------------------------------------------------------------------------
 * Block ILUT (BILUT) preconditioner
 * Block incomplete LU factorization with dual truncation mechanism
 * NOTE : no pivoting implemented as yet in GE for diagonal blocks 
 *----------------------------------------------------------------------------
 * Parameters
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * vbmat    = block matrix stored in VBSpaFmt format -- see type-defs.h for
 *            details on format, the block sizes might be different
 * lu       = pointer to a VBILUSpar struct -- see type-defs.h for details
 *            on format
 * lfil     = integer. The fill-in parameter. Each row of L and
 *            each row of U will have a maximum of lfil elements.
 *            WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
 *            EARLIER VERSIONS. 
 *            lfil must be .ge. 0.
 * tol      = real*8. Sets the threshold for dropping small terms in the
 *            factorization. See below for details on dropping strategy.
 * w        = working array
 * fp       = file pointer for error log ( might be stdout )
 *
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr  = -1  --> Illegal value for lfil
 *            ierr  = -2  --> singular block or zero row encountered
 * lu->n    = dimension of the block matrix
 *   ->bsz  = the row/col of the first element of each diagonal block
 *            the size of the i-th row block should be bsz[i+1] - bsz[i]
 *   ->L    = L part -- stored in VBSpaFmt format
 *   ->D    = Diagonals
 *   ->U    = U part -- stored in VBSpaFmt format
 *----------------------------------------------------------------------------
 * Notes:
 * ======
 * All the diagonal blocks of the input block matrix must not be singular
 *----------------------------------------------------------------------------
 * Dual drop-off strategy works as follows. 
 *
 * 1) Theresholding in L and U as set by tol. Any element whose size
 *    is less than some tolerance (relative to the norm of current
 *    row in u) is dropped.
 *
 * 2) Keeping only the largest lfil elements in the i-th row of L
 *    and the largest lfil elements in the i-th row of U.
 *
 * Flexibility: one can use tol=0 to get a strategy based on keeping the
 * largest elements in each row of L and U. Taking tol .ne. 0 but lfil=n
 * will give the usual threshold strategy (however, fill-in is then
 * impredictible).
 *--------------------------------------------------------------------------*/
int itsol_pc_vbilutC(ITS_VBSparMat *vbmat, ITS_VBILUSpar *lu, int lfil, double tol, ITS_BData * w, FILE * fp)
{
    int n = vbmat->n, *bsz = vbmat->bsz, ierr;
    double one = 1.0, zero = 0.0;
    int dim, szjrow, sz, len, lenu, lenl, col, jpos, jrow, upos, para;
    int max_blk_sz = ITS_MAX_BLOCK_SIZE * ITS_MAX_BLOCK_SIZE * sizeof(double);
    int nzcount, *ja, *jbuf, *iw, i, j, k, kk;
    double t, tnorm, tolnorm, *xnrm, *wn;
    ITS_VBSparMat *L, *U;
    ITS_BData *ba, *D, buf_fact, buf_ns;

    if (lfil < 0) {
        fprintf(fp, "vbilut: Illegal value for lfil.\n");
        return -1;
    }

    itsol_setupVBILU(lu, n, bsz);
    L = lu->L;
    U = lu->U;
    D = lu->D;

    iw = lu->work;
    jbuf = (int *)itsol_malloc(n * sizeof(int), "vbilut");
    buf_fact = (ITS_BData) itsol_malloc(max_blk_sz, "vbilut");
    buf_ns = (ITS_BData) itsol_malloc(max_blk_sz, "vbilut");
    xnrm = (double *)itsol_malloc(n * sizeof(double), "vbilut");
    wn = (double *)itsol_malloc(n * sizeof(double), "vbilut");

    /* set indicator array jw to -1 */
    for (i = 0; i < n; i++)
        iw[i] = -1;

    /* beginning of main loop */
    for (i = 0; i < n; i++) {
        dim = ITS_B_DIM(bsz, i);    /* number of rows of blocks in i-th row */
        nzcount = vbmat->nzcount[i];
        ja = vbmat->ja[i];
        ba = vbmat->ba[i];
        tnorm = 0;
        for (j = 0; j < nzcount; j++) {
            sz = ITS_B_DIM(bsz, ja[j]);
            t = itsol_vbnorm2(dim * sz, ba[j]);
            tnorm = its_max(t, tnorm);
        }
        if (tnorm == 0.0) {
            fprintf(fp, "vbilut:  zero row encountered.\n");
            return -2;
        }
        tolnorm = tol * tnorm;

        /* unpack L-part and U-part of row of A in arrays w */
        lenu = 1;
        lenl = 0;
        jbuf[i] = i;
        itsol_zrmC(dim, dim, w[i]);
        iw[i] = i;

        for (j = 0; j < nzcount; j++) {
            col = ja[j];
            sz = ITS_B_DIM(bsz, col);
            t = itsol_vbnorm2(dim * sz, ba[j]);
            if (t < tolnorm && col != i)
                continue;
            if (col < i) {
                jbuf[lenl] = col;
                iw[col] = lenl;
                itsol_copyBData(dim, sz, w[lenl], ba[j], 0);
                lenl++;
            }
            else if (col == i) {
                itsol_copyBData(dim, dim, w[i], ba[j], 0);
            }
            else {
                jpos = i + lenu;
                jbuf[jpos] = col;
                iw[col] = jpos;
                itsol_copyBData(dim, sz, w[jpos], ba[j], 0);
                lenu++;
            }
        }

        j = -1;
        len = 0;
        /* eliminate previous rows */
        while (++j < lenl) {
            /*----------------------------------------------------------------------
             *   in order to do the elimination in the correct order we must select
             *   the smallest column index among jw(k), k=j+1, ..., lenl.
             *--------------------------------------------------------------------*/
            jrow = jbuf[j];
            jpos = j;
            /* determine smallest column index */
            for (k = j + 1; k < lenl; k++) {
                if (jbuf[k] < jrow) {
                    jrow = jbuf[k];
                    jpos = k;
                }
            }
            szjrow = ITS_B_DIM(bsz, jrow);
            if (jpos != j) {
                col = jbuf[j];
                jbuf[j] = jbuf[jpos];
                jbuf[jpos] = col;
                iw[jrow] = j;
                iw[col] = jpos;
                sz = ITS_B_DIM(bsz, col);
                itsol_copyBData(dim, sz, buf_ns, w[j], 0);
                itsol_copyBData(dim, szjrow, w[j], w[jpos], 0);
                itsol_copyBData(dim, sz, w[jpos], buf_ns, 0);
            }
            /* get the multiplier for row to be eliminated (jrow). */
            /* fact = w(jj)*alu(jrow)                              */
            FC_FUNC(itsol_bxinv,ITSOL_BXINV)(&dim, &szjrow, D[jrow], w[j], buf_fact);

            /* zero out element in row by resetting jw(n+jrow) to -1 */
            iw[jrow] = -1;

            if (itsol_vbnorm2(dim * szjrow, buf_fact) * xnrm[jrow] <= tolnorm)
                continue;

            /* combine current row and row jrow */
            nzcount = U->nzcount[jrow];
            ja = U->ja[jrow];
            ba = U->ba[jrow];
            for (k = 0; k < nzcount; k++) {
                col = ja[k];
                sz = ITS_B_DIM(bsz, col);
                FC_FUNC(dgemm,DGEMM)("n", "n", &dim, &sz, &szjrow, &one, buf_fact, &dim,
                        ba[k], &szjrow, &zero, buf_ns, &dim);
                jpos = iw[col];

                /* if fill-in element is small then disregard: */
                if (itsol_vbnorm2(dim * sz, buf_ns) < tolnorm && jpos == -1)
                    continue;

                if (col >= i) {
                    /* dealing with upper part */
                    //          if( jpos == -1 ) {
                    if (jpos == -1 && itsol_vbnorm2(dim * sz, buf_ns) > tolnorm) {
                        /* this is a fill-in element */
                        upos = i + lenu;
                        jbuf[upos] = col;
                        iw[col] = upos;
                        itsol_copyBData(dim, sz, w[upos], buf_ns, 0);
                        lenu++;
                    }
                    else {
                        /*-------------------- this is not a fill-in element */
                        for (kk = 0; kk < dim * sz; kk++)
                            w[jpos][kk] += buf_ns[kk];

                    }
                }
                else {
                    /* dealing with lower part */
                    if (jpos == -1) {
                        /* this is a fill-in element */
                        jbuf[lenl] = col;
                        iw[col] = lenl;
                        itsol_copyBData(dim, sz, w[lenl], buf_ns, 0);
                        lenl++;
                    }
                    else {
                        /*-------------------- this is not a fill-in element */
                        for (kk = 0; kk < dim * sz; kk++)
                            w[jpos][kk] += buf_ns[kk];
                    }
                }
                }
                itsol_copyBData(dim, szjrow, w[len], buf_fact, 1);
                jbuf[len] = jrow;
                len++;
            }

            /* update l-matrix */
            lenl = len;
            len = its_min(lenl, lfil);
            for (j = 0; j < lenl; j++) {
                sz = ITS_B_DIM(bsz, jbuf[j]);
                wn[j] = itsol_vbnorm2(dim * sz, w[j]);
                iw[j] = j;
            }
            FC_FUNC(itsol_qsplit,ITSOL_QSPLIT)(wn, iw, &lenl, &len);
            L->nzcount[i] = len;
            if (len > 0) {
                L->ja[i] = (int *)itsol_malloc(len * sizeof(int), "vbilut");
                L->ba[i] = (ITS_BData *) itsol_malloc(len * sizeof(ITS_BData), "vbilut");
            }
            ja = L->ja[i];
            ba = L->ba[i];
            for (j = 0; j < len; j++) {
                jpos = iw[j];
                ja[j] = jbuf[jpos];
                sz = ITS_B_DIM(bsz, ja[j]);
                ba[j] = (ITS_BData) itsol_malloc(dim * sz * sizeof(double), "vbilut");
                itsol_copyBData(dim, sz, ba[j], w[jpos], 0);
            }
            for (j = 0; j < lenl; j++)
                iw[j] = -1;

            /* update u-matrix */
            len = its_min(lenu, lfil);
            for (j = 1; j < lenu; j++) {
                jpos = i + j;
                sz = ITS_B_DIM(bsz, jbuf[jpos]);
                wn[j - 1] = itsol_vbnorm2(dim * sz, w[jpos]);
                iw[j - 1] = jpos;
            }
            para = lenu - 1;
            FC_FUNC(itsol_qsplit,ITSOL_QSPLIT)(wn, iw, &para, &len);
            nzcount = U->nzcount[i] = len - 1;
            if (nzcount > 0) {
                U->ja[i] = (int *)itsol_malloc(nzcount * sizeof(int), "vbilut");
                U->ba[i] = (ITS_BData *) itsol_malloc(nzcount * sizeof(ITS_BData), "vbilut");
            }
            ja = U->ja[i];
            ba = U->ba[i];
            t = itsol_vbnorm2(dim * dim, w[i]);
            for (j = 0; j < nzcount; j++) {
                jpos = iw[j];
                ja[j] = jbuf[jpos];
                sz = ITS_B_DIM(bsz, ja[j]);
                ba[j] = (ITS_BData) itsol_malloc(dim * sz * sizeof(double), "vbilut");
                itsol_copyBData(dim, sz, ba[j], w[jpos], 0);
                t = its_max(t, wn[j]);
            }
            for (j = 0; j < lenu - 1; j++)
                iw[j] = -1;
            /* save norm in xnrm. Norm = average abs value. divide by len+1
             * instead of len to avoid division by zero. */
            xnrm[i] = t;

            /* store inverse of diagonal element of u */
            D[i] = (ITS_BData) itsol_malloc(dim * dim * sizeof(double), "vbilut");
            itsol_copyBData(dim, dim, D[i], w[i], 0);

            FC_FUNC(itsol_gauss,ITSOL_GAUSS)(&dim, D[i], &ierr);

            if (ierr != 0) {
                fprintf(fp, "singular block encountered.\n");
                for (j = i + 1; j < n; j++) {
                    D[j] = NULL;
                    L->ja[j] = NULL;
                    L->ba[j] = NULL;
                    U->ja[j] = NULL;
                    U->ba[j] = NULL;
                }
                return -2;
            }
            for (j = 0; j < lenu; j++) {
                iw[jbuf[i + j]] = -1;
            }
        }
        lu->DiagOpt = 1;
        free(jbuf);
        free(buf_fact);
        free(buf_ns);
        free(xnrm);
        free(wn);

        return 0;
    }

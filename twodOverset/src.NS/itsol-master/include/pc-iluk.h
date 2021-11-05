
#ifndef ITSOL_ILUK_H__
#define ITSOL_ILUK_H__

#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

int itsol_pc_lofC(int lofM, ITS_SparMat *csmat, ITS_ILUSpar *lu, FILE *fp); 
int itsol_pc_ilukC(int lofM, ITS_SparMat *csmat, ITS_ILUSpar *lu, FILE *fp);

#ifdef __cplusplus
}
#endif
#endif

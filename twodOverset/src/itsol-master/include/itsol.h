
#ifndef ITSOL_INCLUDED_PROTOS_H__
#define ITSOL_INCLUDED_PROTOS_H__

#include "solver-fgmres.h"
#include "solver-bicgstab.h"
#include "solver-bicgstabl.h"

#include "pc-arms2.h"
#include "pc-iluk.h"
#include "pc-ilutc.h"
#include "pc-ilut.h"
#include "pc-ilutpc.h"
#include "pc-pilu.h"
#include "pc-vbiluk.h"
#include "pc-vbilut.h"

#ifdef __cplusplus
extern "C" {
#endif

void itsol_solver_initialize(ITS_SOLVER *s, ITS_SOLVER_TYPE stype, ITS_PC_TYPE pctype, ITS_CooMat *A);
void itsol_solver_finalize(ITS_SOLVER *s);

int itsol_solver_assemble(ITS_SOLVER *s);

int itsol_solver_solve(ITS_SOLVER *s, double *x, double *rhs);

void itsol_pc_initialize(ITS_PC *pc, ITS_PC_TYPE pctype);
void itsol_pc_finalize(ITS_PC *pc);
int itsol_pc_assemble(ITS_SOLVER *s);

void itsol_solver_set_pars(ITS_SOLVER *s, ITS_PARS par);
void itsol_solver_init_pars(ITS_PARS *par);

#ifdef __cplusplus
}
#endif
#endif 

#ifndef SOLVERELAXED_CALLBACK
#define SOLVERELAXED_CALLBACK

#include "Decomp.h"
#include "DecompAlgo.h"
#include "MCF_Instance.h"

extern MCF_Instance* instance;
extern DecompAlgo* algo;
extern DecompApp* app;

DecompSolverStatus MyRelaxedSolver(const int whichBlock,
                                           const double *redCostX,
                                           const double target,
                                           DecompVarList &varList);

#endif

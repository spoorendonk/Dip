//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Authors: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)      //
//          Ted Ralphs, Lehigh University (ted@lehigh.edu)                   //
//          Jiadong Wang, Lehigh University (jiw408@lehigh.edu)              //
//                                                                           //
// Copyright (C) 2002-2019, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#include "Dip_C_Interface.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct instance
{
        int numCommodities;
        int demand[4];
        int source[4];
        int sink[4];
        int numNodes;
        int numArcs;
        int tail[12];
        int head[12];
        int weight[12];

} instance;

// p frac2 6 12 4
// d 0 3 1
// d 1 2 2
// d 2 1 2
// d 4 5 1
// a 0 1 0 1 1
// a 0 2 0 1 1
// a 1 2 0 1 1
// a 1 3 0 1 1
// a 1 5 0 1 1
// a 2 1 0 1 1
// a 2 3 0 1 1
// a 2 5 0 1 1
// a 3 0 0 1 1
// a 4 1 0 1 1
// a 4 2 0 1 1
// a 5 4 0 1 1

instance small = {.numCommodities = 4,
                  .demand = {1, 2, 2, 1},
                  .source = {0, 1, 2, 4},
                  .sink = {3, 2, 1, 5},
                  .numNodes = 6,
                  .numArcs = 12,
                  .tail = {0, 0, 1, 1, 1, 2, 2, 2, 3, 4, 4, 5},
                  .head = {1, 2, 2, 3, 5, 1, 3, 5, 0, 1, 2, 4},
                  .weight = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}};

int callback_SolvedRelaxed(const int whichBlock,
                           const double *redCostX,
                           const double target,
                           DecompVar ***varList,
                           int *n)
{
        printf("Crazy");
        int len = 2;
        DecompVar **vars = malloc((unsigned int)len * sizeof(DecompVar *));

        for (int i = 0; i < len; i++)
        {
                int len = 1;
                int ind[1] = {i};
                double els[1] = {i};
                double origCost = i;
                double redCost = i * 10;

                DecompVar *var = Dip_DecompVar_new(len, ind, els, origCost, redCost, 0);
                vars[i] = var;
        }

        *n = len;
        *varList = vars;
        return 0;
}

int main(int argc, char *argv[])
{
        struct UtilParameters *utilParam = Dip_UtilParameters_new();
        Dip_UtilParameters_ScanCmdLineArgs(utilParam, argc, argv);
        int doCut = Dip_UtilParameters_getSetting(utilParam, "doCut", 1, 0);
        int doPriceCut = Dip_UtilParameters_getSetting(utilParam, "doPriceCut", 1, 0);

        printf("doCut %d\n", doCut);
        printf("doPriceCut %d\n", doPriceCut);

        struct DecompApp *app = Dip_DecompApp_new(utilParam);

        double infinity = Dip_DecompApp_getInfinity(app);
        int numCols = small.numCommodities * small.numArcs;

        // objective
        double objective[numCols];
        int colIndex = 0;
        for (int k = 0; k < small.numCommodities; k++)
                for (int a = 0; a < small.numArcs; a++)
                {
                        objective[colIndex++] = small.weight[a] * small.demand[k];
                }
        Dip_DecompApp_setModelObjective(app, objective, numCols);

        // models
        DecompConstraintSet *modelCore = Dip_DecompConstraintSet_new();
        Dip_DecompConstraintSet_init(modelCore, numCols, small.numCommodities + small.numArcs);

        // capacity
        for (int a = 0; a < small.numArcs; a++)
        {
                int size = 0;
                int inds[small.numArcs];
                double elems[small.numArcs];

                for (int k = 0; k < small.numCommodities; k++)
                {
                        int colIndex = k * small.numArcs + a;
                        inds[size] = colIndex;
                        elems[size++] = 1.0 * small.demand[k];
                }

                char rowName[20] = "e(";
                char buffer[16];
                sprintf(buffer, "%d_", a);
                strcat(rowName, buffer);
                sprintf(buffer, "%d,", small.tail[a]);
                strcat(rowName, buffer);
                sprintf(buffer, "%d)", small.head[a]);
                strcat(rowName, buffer);

                Dip_DecompConstraintSet_appendRow(modelCore, size, inds, elems, -infinity, small.weight[a], rowName);
        }

        // ensure demand
        for (int k = 0; k < small.numCommodities; k++)
        {
                int source = small.source[k];

                int size = 0;
                int inds[small.numArcs];
                double elems[small.numArcs];

                for (int a = 0; a < small.numArcs; a++)
                {
                        int tail = small.tail[a];
                        colIndex = k * small.numArcs + a;

                        if (tail == source)
                        {
                                inds[size] = colIndex;
                                elems[size++] = 1.0;
                        }

                        char colName[20] = "x(";
                        char buffer2[16];
                        sprintf(buffer2, "%d_", k);
                        strcat(colName, buffer2);
                        sprintf(buffer2, "%d_", a);
                        strcat(colName, buffer2);
                        sprintf(buffer2, "%d,", small.tail[a]);
                        strcat(colName, buffer2);
                        sprintf(buffer2, "%d)", small.head[a]);
                        strcat(colName, buffer2);

                        Dip_DecompConstraintSet_pushCol(modelCore, 0.0, 1.0, colName, 0, -1);
                }
                char rowName[20] = "d(";
                char buffer[16];
                sprintf(buffer, "%d_", k);
                strcat(rowName, buffer);
                sprintf(buffer, "%d)", source);
                strcat(rowName, buffer);

                Dip_DecompConstraintSet_appendRow(modelCore, size, inds, elems, 1.0, 1.0, rowName);
        }

        Dip_DecompApp_setModelCore(app, modelCore, "modelCore");

        // relaxed models
        DecompConstraintSet **modelsRelaxed = malloc((unsigned int)small.numCommodities * sizeof(DecompConstraintSet *));
        for (int k = 0; k < small.numCommodities; k++)
        {
                DecompConstraintSet *modelRelaxed = Dip_DecompConstraintSet_new();
                Dip_DecompConstraintSet_init(modelRelaxed, small.numArcs, small.numNodes);
                Dip_DecompConstraintSet_setSparse(modelRelaxed, numCols);

                for (int i = 0; i < small.numNodes; i++)
                {
                        int size = 0;
                        int inds[small.numArcs];
                        double elems[small.numArcs];

                        for (int a = 0; a < small.numArcs; a++)
                        {
                                int tail = small.tail[a];
                                int head = small.head[a];

                                if (head == i)
                                {
                                        inds[size] = a;
                                        elems[size++] = 1.0;
                                }
                                else if (tail == i)
                                {
                                        inds[size] = a;
                                        elems[size++] = -1.0;
                                }
                        }

                        if (i == small.source[k])
                                Dip_DecompConstraintSet_appendRow(modelRelaxed, size, inds, elems, -1.0, -1.0, "out");
                        else if (i == small.sink[k])
                                Dip_DecompConstraintSet_appendRow(modelRelaxed, size, inds, elems, 1.0, 1.0, "in");
                        else
                        {
                                Dip_DecompConstraintSet_appendRow(modelRelaxed, size, inds, elems, 0.0, 0.0, "balance");
                        }
                }

                int origColIndex = k * small.numArcs;

                for (int a = 0; a < small.numArcs; a++)
                {
                        char colName[20] = "x(";
                        char buffer[16];
                        sprintf(buffer, "%d_", a);
                        strcat(colName, buffer);
                        sprintf(buffer, "%d,", small.tail[a]);
                        strcat(colName, buffer);
                        sprintf(buffer, "%d)", small.head[a]);
                        strcat(colName, buffer);

                        Dip_DecompConstraintSet_pushCol(modelRelaxed, 0.0, 1.0, colName, 0, origColIndex);
                        origColIndex++;
                }

                char string[20] = "relaxed";
                char buffer2[16];
                sprintf(buffer2, "%d", k);
                strcat(string, buffer2);

                Dip_DecompApp_setModelRelax(app, modelRelaxed, string, k);
                modelsRelaxed[k] = modelRelaxed;
        }

        DecompAlgo *algo;

        // setting up algo
        if (doCut)
        {
                algo = Dip_DecompAlgoC_new(app, utilParam);
        }

        if (doPriceCut)
        {
                algo = Dip_DecompAlgoPC_new(app, utilParam);

                Dip_DecompApp_setCallbackSolveRelaxed(app, callback_SolvedRelaxed);
        }

        if (!algo)
        {
                printf("No algo choosen");
                exit(1);
        }

        AlpsDecompModel *alps = Dip_AlpsDecompModel_new(utilParam, algo);
        int status = Dip_AlpsDecompModel_solve(alps);

        printf("status %d", status);

        Dip_AlpsDecompModel_delete(alps);
        Dip_DecompAlgo_delete(algo);
        Dip_DecompApp_delete(app);
        Dip_UtilParameters_delete(utilParam);
        free(modelCore);
        for (int k = 0; k < small.numCommodities; k++)
                free(modelsRelaxed[k]);
        free(modelsRelaxed);
}
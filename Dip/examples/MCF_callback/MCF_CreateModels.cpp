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

//===========================================================================//
#include "DecompVar.h"
#include "MCF_CreateModels.h"

//===========================================================================//
void MCF_CreateModels::initialize()
{
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "initializeApp()", m_appParam.LogLevel, 2);

   //---
   //--- read problem instance
   //---
   string instanceFile   = m_appParam.DataDir
                           + UtilDirSlash() + m_appParam.Instance;
   int rc = m_instance.readInstance(instanceFile, false);

   if (rc)
      throw UtilException("Error in readInstance",
                          "initializeApp", "MCF_DecompApp");

   //---
   //--- create models
   //---
   createModels();
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "initializeApp()", m_appParam.LogLevel, 2);
}


//===========================================================================//
void MCF_CreateModels::createModels()
{
   //---
   //--- This function does the work to create the different models
   //---  that will be used. This memory is owned by the user. It will
   //---  be passed to the application interface and used by the algorithms.
   //---
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "createModels()", m_appParam.LogLevel, 2);
   //---
   //--- (Integer) Multi-Commodity Flow Problem (MCF).
   //---
   //--- We are given:
   //---    (1) a directed graph G=(N,A),
   //---    (2) a set of commodities K, where each commodity is
   //---         a source-sink pair.
   //---
   //--- min  sum{k in K} sum{(i,j) in A} d[i,k] w[i,j] x[k,i,j]
   //--- s.t. sum{(j,i) in A} x[k,i,j] -
   //---        sum{(i,j) in A} x[k,i,j] = d[i,k],  for all i in N, k in K
   //---      sum{k in K} d[i,k] x[k,i,j] <= u[i,j],       for all (i,j) in A
   //---      sum{(i,j) in delta(s)} x[k,i,j] = 1, for all k in K
   //---      x[k,i,j] >= 0 <= 1, for all (i,j) in A
   //--- For k=(s,t) in K,
   //---    d[i,k] = -1 if i=s
   //---           =  1 if i=t
   //---           =  0, otherwise
   //---
   //---
   //--- The decomposition is formed as:
   //---
   //--- MASTER (A''):
   //---      sum{k in K} d[i,k] x[k,i,j] <= u[i,j],       for all (i,j) in A
   //---      sum{(i,j) in delta(s)} x[k,i,j] = 1, for all k in K 
   //---      x[k,i,j] >= 0 <= 1, for all (i,j) in A
   //---
   //--- SUBPROBLEM (A'): (one block for each k in K)
   //---      sum{(j,i) in A} x[k,i,j] -
   //---         sum{(i,j) in A} x[k,i,j] = d[i,k], for all i in N
   //---      x[k,i,j] >= 0 <= 1, for all (i,j) in A
   //--- For k=(s,t) in K,
   //---    d[i,k] = -1 if i=s
   //---           =  1 if i=t
   //---           =  0, otherwise
   //---
   //---
   //--- Get information about this problem instance.
   //---
   int   k, a, colIndex;
   int   numCommodities = m_instance.m_numCommodities;
   int   numArcs        = m_instance.m_numArcs;
   int numCols        = numCommodities * numArcs;
   MCF_Instance::arc* arcs = m_instance.m_arcs;
   MCF_Instance::commodity* commodities = m_instance.m_commodities;
   //---
   //--- Construct the objective function and set it
   //---    columns indexed as [k,a]= k*numArcs + a
   //---
   objective = new double[numCols];
   n = numCols;

   if (!objective) {
      throw UtilExceptionMemory("createModels", "MCF_DecompApp");
   }

   colIndex = 0;

   for (k = 0; k < numCommodities; k++)
      for (a = 0; a < numArcs; a++) {
         objective[colIndex++] = arcs[a].weight * commodities[k].demand;
      }

   //---
   //--- set the objective
   //---
   //setModelObjective(objective, numCols);
   //---
   //--- create the core/master model and set it
   //---
   modelCore = new DecompConstraintSet();
   createModelCore(modelCore);

   //---
   //--- create the relaxed/subproblem models and set them
   //---
   for (k = 0; k < numCommodities; k++) {
      modelRelax = new DecompConstraintSet();
      string                modelName  = "relax" + UtilIntToStr(k);

      if (m_appParam.UseSparse) {
         createModelRelaxSparse(modelRelax, k);
      } else {
         createModelRelax(modelRelax, k);
      }

      m_models.push_back(modelRelax);
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModels()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void MCF_CreateModels::createModelCore(DecompConstraintSet* model)
{
   //---
   //--- MASTER (A''):
   //---      sum{k in K} d[i,k] x[k,i,j] <= u[i,j],       for all (i,j) in A
   //---      sum{(i,j) in delta(s)} x[k,i,j] = 1, for all k in K 
   //---      x[k,i,j] >= 0 <= 1, for all (i,j) in A
   //---
   int   k, a, colIndex;
   int   numCommodities = m_instance.m_numCommodities;
   int   numArcs        = m_instance.m_numArcs;
   int   numCols        = numCommodities * numArcs;
   int   numRows        = 2              * numArcs;
   MCF_Instance::arc* arcs = m_instance.m_arcs;
   MCF_Instance::commodity* commodities = m_instance.m_commodities;
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "createModelCore()", m_appParam.LogLevel, 2);
   //---
   //--- create space for the model matrix (row-majored)
   //---
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);

   if (!model->M) {
      throw UtilExceptionMemory("createModelCore", "MCF_DecompApp");
   }

   model->M->setDimensions(0, numCols);
   model->reserve(numRows, numCols);
   //---
   //--- create the rows and set the col/row bounds
   //---
   UtilFillN(model->colLB, numCols,  0.0);
   UtilFillN(model->colUB, numCols,  1.0);

   for (a = 0; a < numArcs; a++) {
      CoinPackedVector row;
      double           arcLB = arcs[a].lb;
      double           arcUB = arcs[a].ub;

      for (k = 0; k < numCommodities; k++) {
         colIndex = k * numArcs + a;
         model->colLB[colIndex] = 0.0;
         model->colUB[colIndex] = 1.0;
         row.insert(colIndex, 1.0 * commodities[k].demand);
      }

      model->appendRow(row, -m_infinity, arcUB);
      string rowNameUB = "e(" +
                         UtilIntToStr(a)            + "_" +
                         UtilIntToStr(arcs[a].tail) + "," +
                         UtilIntToStr(arcs[a].head) + ")";
      model->rowNames.push_back(rowNameUB);
   }

   for (k = 0; k < numCommodities; k++) {
      int source = commodities[k].source;
      
      CoinPackedVector row;

      for (a = 0; a < numArcs; a++) {
         int tail = arcs[a].tail;
         colIndex = k * numArcs + a;

         if (tail == source)
            row.insert(colIndex, 1.0);
      }
      model->appendRow(row, 1.0, 1.0);
      string rowNameUB = "d(" +
                         UtilIntToStr(a)            + "_" +
                         UtilIntToStr(arcs[a].tail) + ")";
      model->rowNames.push_back(rowNameUB);
   }

   //---
   //--- create column names (helps with debugging)
   //---
   for (k = 0; k < numCommodities; k++) {
      for (a = 0; a < numArcs; a++) {
         string colName = "x(comm_" + UtilIntToStr(k) + "," +
                          UtilIntToStr(a)            + "_" +
                          UtilIntToStr(arcs[a].tail) + "," +
                          UtilIntToStr(arcs[a].head) + ")";
         model->colNames.push_back(colName);
         cout << colName << endl;
      }
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "Size cols x rows = " + UtilIntToStr(numCols) + " x " + UtilIntToStr(numRows), m_appParam.LogLevel, 2);
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModelCore()", m_appParam.LogLevel, 2);
}


//===========================================================================//
void MCF_CreateModels::createModelRelax(DecompConstraintSet* model,
                                     int                   commId)
{
   //--- SUBPROBLEM (A'): (one block for each k in K)
   //---      sum{(j,i) in A} x[k,i,j] -
   //---         sum{(i,j) in A} x[k,i,j] = d[i,k], for all i in N
   //---      x[k,i,j] >= 0 <= 1, for all (i,j) in A
   //--- For k=(s,t) in K,
   //---    d[i,k] = -1 if i=s
   //---           =  1 if i=t
   //---           =  0, otherwise
   int         a, i, head, tail, colIndex, source, sink;
   int         numCommodities = m_instance.m_numCommodities;
   int         numArcs        = m_instance.m_numArcs;
   int         numNodes       = m_instance.m_numNodes;
   int         numCols        = numCommodities * numArcs;
   int         numRows        = numNodes;
   MCF_Instance::arc*        arcs        = m_instance.m_arcs;
   MCF_Instance::commodity* commodities = m_instance.m_commodities;
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "createModelRelax()", m_appParam.LogLevel, 2);
   //---
   //--- create space for the model matrix (row-majored)
   //---
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);

   if (!model->M) {
      throw UtilExceptionMemory("createModelCore", "MCF_DecompApp");
   }

   model->M->setDimensions(0, numCols);
   model->reserve(numRows, numCols);
   //---
   //--- get this commodity's source and sink node
   //---
   source = commodities[commId].source;
   sink   = commodities[commId].sink;

   //---
   //--- create the rows
   //---   NOTE: this is somewhat inefficient (but simple)
   //---
   for (i = 0; i < numNodes; i++) {
      CoinPackedVector row;

      for (a = 0; a < numArcs; a++) {
         tail = arcs[a].tail;
         head = arcs[a].head;

         if (head == i) {
            colIndex = commId * numArcs + a;
            row.insert(colIndex, 1.0);
         } else if (tail == i) {
            colIndex = commId * numArcs + a;
            row.insert(colIndex, -1.0);
         }
      }

      if (i == source)
         model->appendRow(row,
                          -1,
                          -1);
      else if (i == sink)
         model->appendRow(row,
                          1,
                          1);
      else {
         model->appendRow(row, 0.0, 0.0);
      }

      string rowName = "flow(" +
                       UtilIntToStr(commId) + "_" +
                       UtilIntToStr(i)      + "_" +
                       UtilIntToStr(source) + "," +
                       UtilIntToStr(sink)   + ")";
      model->rowNames.push_back(rowName);
   }

   //---
   //--- If using callback function, then skip setup of subproblem constraints and only add bound constraints
   //---
   // for (a = 0; a < numArcs; a++) {
   //    CoinPackedVector row;
   //    tail = arcs[a].tail;
   //    head = arcs[a].head;
   //    row.insert(a, 1.0);
   //    model->appendRow(row, 0.0, m_infinity);
   // }

   //---
   //--- create a list of the "active" columns (those related
   //---   to this commmodity) all other columns are fixed to 0
   //---
   UtilFillN(model->colLB, numCols,  0.0);
   UtilFillN(model->colUB, numCols,  0.0);
   colIndex = commId * numArcs;

   for (a = 0; a < numArcs; a++) {
      double           arcLB = 0.0;
      double           arcUB = 1.0;
      model->colLB[colIndex] = arcLB;
      model->colUB[colIndex] = arcUB;
      model->activeColumns.push_back(colIndex);
      colIndex++;
   }

   //---
   //--- set the indices of the integer variables of model
   //---
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModelRelax()", m_appParam.LogLevel, 2);
}

//===========================================================================//
void MCF_CreateModels::createModelRelaxSparse(DecompConstraintSet* model,
      int                   commId)
{
   //--- SUBPROBLEM (A'): (one block for each k in K)
   //---      sum{(j,i) in A} x[k,i,j] -
   //---         sum{(i,j) in A} x[k,i,j] = d[i,k], for all i in N
   //---      x[k,i,j] >= 0 <= 1, for all (i,j) in A
   //--- For k=(s,t) in K,
   //---    d[i,k] = -1 if i=s
   //---           =  1 if i=t
   //---           =  0, otherwise
   int         a, i, head, tail, origColIndex, source, sink;
   int         numArcs        = m_instance.m_numArcs;
   int         numNodes       = m_instance.m_numNodes;
   int         numCommodities = m_instance.m_numCommodities;
   int         numCols        = numArcs;
   int         numRows        = numNodes;
   int         numColsOrig    = numArcs * numCommodities;
   MCF_Instance::arc*        arcs        = m_instance.m_arcs;
   MCF_Instance::commodity* commodities = m_instance.m_commodities;
   UtilPrintFuncBegin(m_osLog, m_classTag,
                      "createModelRelaxSparse()", m_appParam.LogLevel, 2);
   //---
   //--- create space for the model matrix (row-majored)
   //---
   model->M = new CoinPackedMatrix(false, 0.0, 0.0);

   if (!model->M) {
      throw UtilExceptionMemory("createModelCore", "MCF_DecompApp");
   }

   model->M->setDimensions(0, numCols);
   model->reserve(numCols, numCols);
   model->setSparse(numColsOrig);
   //---
   //--- get this commodity's source and sink node
   //---
   source = commodities[commId].source;
   sink   = commodities[commId].sink;

   //---
   //--- create the rows
   //---   NOTE: this is somewhat inefficient (but simple)
   //---
   for (i = 0; i < numNodes; i++) {
      CoinPackedVector row;

      for (a = 0; a < numArcs; a++) {
         tail = arcs[a].tail;
         head = arcs[a].head;

         if (head == i) {
            row.insert(a, 1.0);
         } else if (tail == i) {
            row.insert(a, -1.0);
         }
      }

      if (i == source)
         model->appendRow(row,
                          -1, //commodities[commId].demand,
                          -1); //commodities[commId].demand);
      else if (i == sink)
         model->appendRow(row,
                          1, // commodities[commId].demand,
                          1); //commodities[commId].demand);
      else {
         model->appendRow(row, 0.0, 0.0);
      }
   }

   //---
   //--- If using callback function, then skip setup of subproblem constraints and only add bound constraints
   //---
   // for (a = 0; a < numArcs; a++) {
   //    CoinPackedVector row;
   //    tail = arcs[a].tail;
   //    head = arcs[a].head;
   //    row.insert(a, 1.0);
   //    model->appendRow(row, 0.0, m_infinity);
   // }

   //---
   //--- set the colLB, colUB, integerVars and sparse mapping
   //---
   origColIndex = commId * numArcs;

   for (a = 0; a < numArcs; a++) {
      double           arcLB = 0;
      double           arcUB = 1;
      model->pushCol(arcLB, arcUB, false, origColIndex);
      origColIndex++;
   }

   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "Size cols x rows = " + UtilIntToStr(numCols) + " x " + UtilIntToStr(numRows), m_appParam.LogLevel, 2);
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "createModelRelaxSparse()", m_appParam.LogLevel, 2);
}

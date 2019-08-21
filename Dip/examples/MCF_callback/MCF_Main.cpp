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
#include "UtilParameters.h"
//===========================================================================//
#include "MCF_CreateModels.h"
//===========================================================================//
#include "AlpsDecompModel.h"
//===========================================================================//
#include "DecompAlgoC.h"
#include "DecompAlgoPC.h"
//===========================================================================//
#include "UtilTimer.h"

#include "solveRelaxed_callback.h"

class callbacks{
   public:
DecompSolverStatus MyRelaxedSolverCallback(const DecompApp* app,
   const int whichBlock,
                                   const double *redCostX,
                                   const double target,
                                   DecompVarList &varList)
{
   std::cout << "MyRelaxedSolverCallback" << std::endl;
   return DecompSolStatNoSolution;
}

int MyGenerateCutsCallback( const DecompApp *app,
    const double *x,
    DecompCutList &newCuts)
{
   std::cout << "MyGenerateCutsCallback" << std::endl;
   return 0;
}
};



class DecompCallbackDataSolStatus : public DecompCallbackData
{
private:
   const double *m_x;
   int m_numCols;
   double m_tolZero;
   bool m_isFeasible = true;

public:
   void set(const double *x,
            const int numCols,
            const double tolZero)
   {
      m_x = x;
      m_numCols = numCols;
      m_tolZero = tolZero;
      m_isFeasible = true;
   }

   bool get(){
      return m_isFeasible;
   }

   const int getInt(DecompCallbackWhere where, DecompCallbackWhat what)
   {
      if (where == MIPSOL_FEAS && what == MIPSOL_NUMCOLS)
      {
         return m_numCols;
      }

      return DecompCallbackData::getInt(where, what);
   }

   const double getDouble(DecompCallbackWhere where, DecompCallbackWhat what)
   {
      if (where == MIPSOL_FEAS && what == MIPSOL_TOLZERO)
      {
         return m_tolZero;
      }

      return DecompCallbackData::getDouble(where, what);
   }
   const double *getDoubles(DecompCallbackWhere where, DecompCallbackWhat what)
   {
      if (where == MIPSOL_FEAS && what == MIPSOL_X)
      {
         return m_x;
      }

      return DecompCallbackData::getDoubles(where, what);
   }

   void setSolutionStatus(DecompCallbackWhere where, DecompCallbackWhat what, bool isFeasible)
   {
      if(where == MIPSOL_FEAS && what == MIPSOL_STATUS){
         m_isFeasible = isFeasible;
      }
      
      DecompCallbackData::setSolutionStatus(where, what, isFeasible);
   }
};

// class DecompCallbackDataSol : public DecompCallbackData{
//    public:
//       DecompCallbackWhere where = GENERATE_VARS;
// };

// class DecompCallbackDataVars : public DecompCallbackData{
//    public:
//       DecompCallbackWhere where = GENERATE_VARS;
// };

// class DecompCallbackDataCuts : public DecompCallbackData{
//    public:
//       DecompCallbackWhere where = GENERATE_CUTS;
// };

void myCallback(DecompApp* app, DecompCallbackWhere where, DecompCallbackData& data){
   std::cout << "myCallback " << where << std::endl;
         
   if(where == DecompCallbackWhere::MIPSOL_FEAS){
      const int numCols = data.getInt(where, DecompCallbackWhat::MIPSOL_NUMCOLS);
      const double* x = data.getDoubles(where, DecompCallbackWhat::MIPSOL_X);
      const double tolZero = data.getDouble(where, DecompCallbackWhat::MIPSOL_TOLZERO);

      // dostuff
      bool isFeasible = true;

      data.setSolutionStatus(where, DecompCallbackWhat::MIPSOL_STATUS, isFeasible);
   }

   if(where == DecompCallbackWhere::MIPSOL_HEUR){
      const double* xhat = data.getDoubles(where, DecompCallbackWhat::MIPSOL_X);
      const double* origCost = data.getDoubles(where, DecompCallbackWhat::MIPSOL_XOBJ);

      // dostuff
      DecompSolution *sol;

      data.setSolution(where, DecompCallbackWhat::MIPSOL_ADD, sol);
   }

   if(where == DecompCallbackWhere::INITVARS){
      // dostuff
      DecompVar *var;

      data.addVar(where, DecompCallbackWhat::MIPVARS_ADD, var);
   }

   if(where == DecompCallbackWhere::MIPVARS){
      const int whichBlock = data.getInt(where, DecompCallbackWhat::MIPVARS_BLOCK);
      const double* redCostX = data.getDoubles(where, DecompCallbackWhat::MIPVARS_REDCOSTX);
		const double target = data.getDouble(where, DecompCallbackWhat::MIPVARS_TARGET);

      // dostuff
      DecompVar *var;
      DecompSolverStatus status = DecompSolStatNoSolution;

      data.addVar(where, DecompCallbackWhat::MIPVARS_ADD, var);
      data.setSolverStatus(where, DecompCallbackWhat::MIPVARS_ADD, status);
   }

   if(where == DecompCallbackWhere::MIPCUTS){
      const double* x = data.getDoubles(where, DecompCallbackWhat::MIPCUTS_X);

      // dostuff
      DecompCut *cut;

      data.addCut(where, DecompCallbackWhat::MIPCUTS_ADD, cut);
   }
} 

//===========================================================================//
int main(int argc, char** argv)
{
   try {
      //---
      //--- create the utility class for parsing parameters
      //---
      UtilParameters utilParam(argc, argv);
      bool doCut          = utilParam.GetSetting("doCut",          true);
      bool doPriceCut     = utilParam.GetSetting("doPriceCut",     false);
      bool doDirect       = utilParam.GetSetting("doDirect",       false);
      UtilTimer timer;
      double    timeSetupReal = 0.0;
      double    timeSetupCpu  = 0.0;
      double    timeSolveReal = 0.0;
      double    timeSolveCpu  = 0.0;
      //---
      //--- start overall timer
      //---
      timer.start();
      //---
      //--- create the user application (a DecompApp)
      //---
      MCF_CreateModels *model = new MCF_CreateModels(utilParam);
      instance = &(model->m_instance);

      app = new DecompApp(utilParam);

      app->setModelObjective(model->objective, model->n);
      app->setModelCore(model->modelCore, "core");

      int k = 0;
      for(auto x: model->m_models){
         string modelName  = "relax" + UtilIntToStr(k);
         app->setModelRelax(x, modelName, k);
         k++;     
      }

      // callbacks
      //app->setCallbackSolveRelaxed(MyRelaxedSolver);

      DecompCallbackDataSolStatus data;
      double x[] = {1,2,3};
      data.set(x,3,0.1);

      app->callback(DecompCallbackWhere::MIPSOL_FEAS, data);
      app->setCallback(myCallback);
      std::cout << data.get() << std::endl;
      app->callback(DecompCallbackWhere::MIPSOL_FEAS, data);
      std::cout << data.get() << std::endl;




      // callbacks c;
      // DecompCallbackSolveRelaxedWrong fn = std::bind(&callbacks::MyRelaxedSolverCallback, &c, 
      //       std::placeholders::_1, std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::placeholders::_5);
      // app->setCallback(fn);
      //app->setCallback(c.MyGenerateCutsCallback);
      
      // -----------------------------------------------------------------------------------------------------
      //---
      //--- create the algorithm (a DecompAlgo)
      //---
      
      assert(doCut + doPriceCut == 1);

      //---
      //--- create the CPM algorithm object
      //---
      if (doCut) {
         algo = new DecompAlgoC(app, utilParam);
      }

      //---
      //--- create the PC algorithm object
      //---
      if (doPriceCut) {
         algo = new DecompAlgoPC(app, utilParam);
      }

      if (doCut && doDirect) {
         timer.stop();
         timeSetupCpu  = timer.getCpuTime();
         timeSetupReal = timer.getRealTime();
         //---
         //--- solve
         //---
         timer.start();
         algo->solveDirect();
         timer.stop();
         timeSolveCpu  = timer.getCpuTime();
         timeSolveReal = timer.getRealTime();
      } else {
         //---
         //--- create the driver AlpsDecomp model
         //---
         int             status = 0;
         AlpsDecompModel alpsModel(utilParam, algo);
         timer.stop();
         timeSetupCpu  = timer.getCpuTime();
         timeSetupReal = timer.getRealTime();
         //---
         //--- solve
         //---
         timer.start();
         status = alpsModel.solve();
         timer.stop();
         timeSolveCpu  = timer.getCpuTime();
         timeSolveReal = timer.getRealTime();
         //---
         //--- sanity check
         //---
         cout << setiosflags(ios::fixed | ios::showpoint);
         cout << "Status= "   << status
              << " BestLB=  " << setw(10)
              << UtilDblToStr(alpsModel.getGlobalLB(), 5)
              << " BestUB= " << setw(10)
              << UtilDblToStr(alpsModel.getGlobalUB(), 5)
              << " Nodes= " << setw(6)
              << alpsModel.getNumNodesProcessed()
              << " SetupCPU= "  << timeSetupCpu
              << " SolveCPU= "  << timeSolveCpu
              << " TotalCPU= "  << timeSetupCpu + timeSolveCpu
              << " SetupReal= " << timeSetupReal
              << " SolveReal= " << timeSolveReal
              << " TotalReal= " << timeSetupReal + timeSolveReal
              << endl;
         //---
         //--- free local memory
         //---
         delete algo;
      }
   } catch (CoinError& ex) {
      cerr << "COIN Exception [ " << ex.message() << " ]"
           << " at " << ex.fileName()  << ":L" << ex.lineNumber()
           << " in " << ex.className() << "::" << ex.methodName() << endl;
      return 1;
   }

   return 0;
}


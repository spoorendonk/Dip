//===========================================================================//
// This file is part of the DIP Solver Framework.                            //
//                                                                           //
// DIP is distributed under the Eclipse Public License as part of the        //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Conceptual Design: Matthew Galati, SAS Institute Inc.                     //
//                    Ted Ralphs, Lehigh University                          //
//                                                                           //
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

//===========================================================================//
#ifndef DecompAlgoPC_h_
#define DecompAlgoPC_h_

//===========================================================================//
/**
 * \class DecompAlgoPC
 * \brief Class for DECOMP algorithm Price and Cut.
 *
 */
//===========================================================================//

//===========================================================================//
#include "DecompAlgo.h"

//===========================================================================//
class DecompAlgoPC : public DecompAlgo {
private:

   //----------------------------------------------------------------------//
   /**
    * @name Data.
    * @{
    */
   //----------------------------------------------------------------------//
   /**
    * Store the name of the class (for logging/debugging) - "who am I?"
    */
   string m_classTag;
   vector<double> m_dual;    //duals from stabilized (if bound improved)
   vector<double> m_dualRM;  //duals from restricted master
   vector<double> m_dualST;  //duals from stabilized method

   /**
    * @}
    */

   //-----------------------------------------------------------------------//
   /**
    * @name Derived from pure virtual functions of DecompAlgo.
    * @{
    */
   //-----------------------------------------------------------------------//
   /**
    * Create the master problem (all algorithms must define this function).
    */
   virtual void createMasterProblem(DecompVarList & initVars){
      DecompAlgo::createMasterProblem(initVars);
   }
   virtual int generateVarsFea(DecompVarList    & newVars, 
			       double           & mostNegReducedCost){
      return DecompAlgo::generateVarsFea(newVars, mostNegReducedCost);
   }
   virtual void phaseInit(DecompPhase & phase);


   virtual const double * getMasterDualSolution()  {
      //---
      //--- return the duals to be used in pricing step
      //---
      if(m_param.DualStab){
	 //---
	 //--- resize dual vectors
	 //---
	 int nRows = static_cast<int>(m_masterSI->getNumRows());
	 m_dual.resize(nRows);
	 m_dualRM.resize(nRows);
	 m_dualST.resize(nRows);

	 printf("cutCallsTotal=%d, priceCallsTotal=%d\n",
		m_nodeStats.cutCallsTotal,
		m_nodeStats.priceCallsTotal);
	 //---
	 //--- calculate smoothed dual
	 //---    pi_ST = alpha * pi_Bar + (1-alpha) * pi_RM
	 //--- this is dual feasible because it is taking
	 //---    a convex combination of previously dual feasible
	 //---    vectors
	 //--- need to be careful here, as the init dual is 0, which
	 //---    might not be dual feasible, therefore, in the first
	 //---    iteration, we need to skip the smoothing and enforce
	 //---    that the first dual be set to dualRM
	 //---
	 int            r;
	 const double * u      = &m_dualSolution[0];
	 double         alpha  = m_param.DualStabAlpha;
	 double         alpha1 = 1.0 - alpha; 
	 copy(u, u + nRows, m_dualRM.begin()); //copy for sake of debugging

	 //---
	 //--- for both the first PhaseI and first PhaseII calls,
	 //---   be sure to set the dual vector to dualRM as dual=0
	 //---   might not be feasible
	 //---
	 if(((m_nodeStats.cutCallsTotal + 
	      m_nodeStats.priceCallsTotal) == 1) || m_firstPhase2Call){
	    (*m_osLog) << "Init dual to dualRM" << endl;
	    copy(m_dualRM.begin(), m_dualRM.end(), m_dual.begin());
	 }
	 
	 for(r = 0; r < nRows; r++){
	    m_dualST[r] = (alpha * m_dual[r]) + (alpha1 * m_dualRM[r]);
	 }      
	 
	 const vector<string> & rowNames = m_masterSI->getRowNames();
	 if(m_param.LogDebugLevel >= 3){
	    for(r = 0; r < m_masterSI->getNumRows(); r++){
	       if(!(UtilIsZero(m_dual[r]) && 
		    UtilIsZero(m_dualRM[r]) && UtilIsZero(m_dualST[r]))){
		  if(r < static_cast<int>(rowNames.size())){
		     (*m_osLog) << "MASTER " 
				<< DecompRowTypeStr[m_masterRowType[r]]
				<< " DUAL[ " << r << "->" << rowNames[r]
				<< "] = " << m_dual[r] << " RM = " 
				<< m_dualRM[r] << " ST = " << m_dualST[r]
				<< endl;
		  }
		  else
		     (*m_osLog) << "MASTER " 
				<< DecompRowTypeStr[m_masterRowType[r]]
				<< " DUAL[ " << r
				<< "] = " << m_dual[r] << " RM = " 
				<< m_dualRM[r] << " ST = " << m_dualST[r]
				<< endl;
	       }
	    }
	 }
	 return &m_dualST[0];
      }
      else{
	 return &m_dualSolution[0];
      }
   }
   

   virtual void setObjBoundLB(const double thisBound,
			      const double thisBoundUB){
      UtilPrintFuncBegin(m_osLog, m_classTag,
			 "setObjBoundLB()", m_param.LogDebugLevel, 2);
      if(m_param.DualStab){
	 if(thisBound > m_nodeStats.objBest.first){
	    (*m_osLog) << "Bound improved " << m_nodeStats.objBest.first
		       << " to " << thisBound << " , update duals" << endl;
	    copy(m_dualST.begin(), m_dualST.end(), m_dual.begin());
	 }
      }
      DecompAlgo::setObjBoundLB(thisBound, thisBoundUB);
      UtilPrintFuncEnd(m_osLog, m_classTag,
		       "setObjBoundLB()", m_param.LogDebugLevel, 2);
   }
   
   /**
    * @}
    */
   
   //-----------------------------------------------------------------------//
   /**
    * @name Derived from virtual functions of DecompAlgo
    * @{
    */
   //-----------------------------------------------------------------------//
   //TODO
   void addCutsToPool(const double  *  x,
                      DecompCutList & newCuts,
                      int           & n_newCuts);

   //TODO
   void phaseDone();
   int  addCutsFromPool();
   void solutionUpdateAsIP();
   int  adjustColumnsEffCnt();
   int  compressColumns    ();


   
   /**
    * @}
    */

   
   //-----------------------------------------------------------------------//
   /**
    * @name Constructors and destructor.
    * @{
    */
   //-----------------------------------------------------------------------//
public:
   /**
    * Default constructors.
    */   
   DecompAlgoPC(DecompApp      * app,
                UtilParameters * utilParam,
                bool             doSetup    = true) :
      DecompAlgo(PRICE_AND_CUT, app, utilParam),
      m_classTag("D-ALGOPC") {

      //---
      //--- do any parameter overrides of the defaults here
      //---    by default turn off gomory cuts for PC
      //---
      m_param.CutCglGomory = 0;

      //---
      //--- run init setup
      //---
      if(doSetup){
         string paramSection = DecompAlgoStr[PRICE_AND_CUT];
         initSetup(utilParam, paramSection);
      }
   }
      
   DecompAlgoPC(DecompApp      * app,
                UtilParameters * utilParam,
                string         & paramSection,
                bool             doSetup    = true) :
      //is utilParam used in base class?
      DecompAlgo(PRICE_AND_CUT, app, utilParam),
      m_classTag("D-ALGOPC") {

      //---
      //--- do any parameter overrides of the defaults here
      //---    by default turn off gomory cuts for PC
      //---
      m_param.CutCglGomory = 0;

      //---
      //--- run init setup
      //---
      if(doSetup)
         initSetup(utilParam, paramSection);
   }
      
      
   /**
    * Destructor.
    */
   ~DecompAlgoPC(){}
   /**
    * @}
    */
};

#endif

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

#ifndef MCF_DECOMPAPP2_INCLUDED
#define MCF_DECOMPAPP2_INCLUDED

#include "DecompConstraintSet.h"

//===========================================================================//
#include "MCF_Instance.h"
#include "MCF_Param.h"
//===========================================================================//

//===========================================================================//
/*!
 * \class MCF_DecompApp
 * A DecompApp for solving the
 *     (Integer) Multi-Commodity Flow Problem (MCF)
 *
 * \see
 * DecompApp
 *
 */

//===========================================================================//
class MCF_CreateModels {
public:
   /** Class id tag (for log / debugging). */
   const string m_classTag;

   /** Application specific parameters. */
   MCF_Param m_appParam;

   /** MCF problem instance data */
   MCF_Instance m_instance;

public:
   /** The model objective coefficients (original space). */
   double* objective;

   /** Model constraint systems. */
   vector<DecompConstraintSet*> m_models;

   DecompConstraintSet* modelRelax;
   DecompConstraintSet* modelCore;

   std::ostream* m_osLog;
   double m_infinity = OsiClpInfinity;
   int n;

public:
   /** @name Helper functions (public). */


   /** Initialize application. */
   void initialize();

   /* Create models. */
   void createModels();
   void createModelCore(DecompConstraintSet* model);
   void createModelRelax(DecompConstraintSet* model,
                         int                   commId);
   void createModelRelaxSparse(DecompConstraintSet* model,
                               int                   commId);

public:
   /** @name Constructor and Destructor */

   /** Default constructor. Takes an instance of UtilParameters */
   MCF_CreateModels(UtilParameters& utilParam) :
      m_classTag  ("MCF-MODELS"),      
      objective   (   NULL  ),
      modelRelax  (   NULL  ),
      modelCore   (   NULL  ),
      m_osLog     (&std::cout  )
   {
      //---
      //--- get application parameters
      //---
      m_appParam.getSettings(utilParam);
      
      if (m_appParam.LogLevel >= 1) {
	 m_appParam.dumpSettings(m_osLog);
      }
      
      initialize();
   }

   virtual ~MCF_CreateModels() {
      UTIL_DELARR(objective);
      UtilDeleteVectorPtr(m_models);
      UTIL_DELPTR(modelCore);
   };
};

#endif

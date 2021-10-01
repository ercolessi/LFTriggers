
////////////////////////////////////////////////////////////////////////////////////
//                                                                                //
// Macro to run strangeness analysis vs multiplicity & effective energy on grid   //
// author: Francesca Ercolessi                                                    //
//         francesca.ercolessi@cern.ch                                            //
//                                                                                //
// To Run:                                                                        //
// testmode) gridtest=KTRUE                                                       //
// gridmode) gridtest=KFALSE                                                      //
// 1)run with runmode=full and setMergeViaJDL(kTRUE) to run on grid               //
// 2)run with runmode=terminate and SetMergeViaJDL(kTRUE) to merge                //
// 3)run with runmode=terminate and SetMergeViaJDL(kFALSE) to download in local   //
//    the merged root file                                                        //
//                                                                                //
//                                                                                //
//                                                                                //              
////////////////////////////////////////////////////////////////////////////////////

#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"
//#include "AddTaskPhysicsSelection.h"
#include "AliAnalysisTask.h"
#include "AliPhysicsSelectionTask.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliAnalysisTaskWeakDecayVertexer.h"
#include "AliAnalysisTaskKzeroshort.h"
#include "AliPPVsMultUtils.h"

#endif

void run_kzeroshort()
  
{
  // if you need to access data remotely from the GRID
  Bool_t grid = 1;  
  // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
  Bool_t gridTest = kTRUE;
  
  // Load common libraries
  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");   
   
  // include the path you need to compile and the library you need
  gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");   
  //Tell root where to look for headers
  #if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
  #else
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  #endif

  // Create the analysis manager and AOD handler
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  AliAODInputHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);

  //PhysicsSelection Configuration
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* ps =  reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C(kFALSE)"));
  
  //MultSelection
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask* ms =  reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
  ms->SetAddInfo(kTRUE);
  
  //PID Response
  //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse* pid = reinterpret_cast<AliAnalysisTaskPIDResponse*>(gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(kFALSE)"));
   
  //Strangness Task 
  gROOT->LoadMacro("AliAnalysisTaskKzeroshort.cxx++g");
  gROOT->LoadMacro("AddTaskKzeroshort.C"); 
  AliAnalysisTaskKzeroshort* st = reinterpret_cast<AliAnalysisTaskKzeroshort*>(gInterpreter->ExecuteMacro("AddTaskKzeroshort.C"));
 

  if(!mgr->InitAnalysis()) return;
  mgr->SetDebugLevel(2);
  mgr->PrintStatus();

  if (grid)
    {
      AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
      alienHandler->SetOverwriteMode();

      alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        
      // Set versions of used packages
      alienHandler->SetAPIVersion("V1.1x");
      //alienHandler->SetROOTVersion("v5-26-00b-6");
      //alienHandler->SetAliROOTVersion("v4-19-21-AN");
      //Please keep this version updated 
      alienHandler->SetAliPhysicsVersion("vAN-20210826_ROOT6-1");   
      alienHandler->SetAnalysisMacro("AnalysisLeading.C"); 

      // number of files per subjob
      //alienHandler->SetSplitMaxInputFileNumber(5);
      //set Executable
      alienHandler->SetExecutable("runLeadingAOD.sh");
      //specify how many seconds your job may take
      alienHandler->SetTTL(36000);
      //set jdl name
      alienHandler->SetJDLName("runLeadingAOD.jdl");
      
      alienHandler->SetOutputToRunNo(kTRUE);
      alienHandler->SetKeepLogs(kTRUE);
   
      alienHandler->SetTerminateFiles("event_stat.root");  //to have the output file of the Physics Selection class
      alienHandler->SetInputFormat("xml-single");
      alienHandler->SetPrice(1);      
      // Optionally modify split mode (default 'se')    
      alienHandler->SetSplitMode("se");
       
      // make sure your source files get copied to grid
      alienHandler->SetAdditionalLibs("AliAnalysisTaskKzeroshort.h AliAnalysisTaskKzeroshort.cxx");
      alienHandler->SetAnalysisSource("AliAnalysisTaskKzeroshort.cxx ");
      
      //Declare input data to be processed.
      //Method 1: Create automatically XML collections using alien 'find' command.
           
      // select the input data
      alienHandler->SetGridDataDir("/alice/data/2015/LHC15o/"); 
      alienHandler->SetRunPrefix("000");
      alienHandler->SetDataPattern("pass5_lowIR/AOD194/*/AliAOD.root");
      // runnumber
      Int_t runList[11] = {244917, 244918, 244975, 244980, 244982, 244983, 245064, 245066, 245068, 246391, 246392};
      for (Int_t i= 4;i <5; i++) alienHandler->AddRunNumber(runList[i]);	
      alienHandler->SetGridWorkingDir("LFTriggers");               

     // alienHandler->SetGridDataDir("/alice/data/2018/LHC18i/");
     // alienHandler->SetRunPrefix("000");
     // alienHandler->SetDataPattern("pass1/AOD264/AOD/*/AliAOD.root");
      // runnumber
     // Int_t runList[9] = {288861, 288862, 288863, 288864, 288868, 288902, 288903, 288908, 288909};
     // for (Int_t i= 0;i <1; i++) alienHandler->AddRunNumber(runList[i]);	
     // alienHandler->SetGridWorkingDir("LFTriggers");
      
      // define the output folder
      alienHandler->SetGridOutputDir("OutputDir");
      //alienHandler->SetDefaultOutputs();
      
      //number of times to merge (if a lot of data need a higher number)
      alienHandler->SetMaxMergeStages(5);
      alienHandler->SetSplitMaxInputFileNumber(5);
      // we can specify that we want to, later on, use Grid to also merge
      // our output. to enable this, we will set 'SetMergeViaJDL' to kTRUE
      alienHandler->SetMergeViaJDL(kTRUE);
      //When all your merging jobs have finished running,
      //there is one step to be taken still, which is downloading the output
      //of the merging jobs. This is done by changing SetMergeViaJDL(kFALSE)
      //and running one last time
       
      // Connect plug-in to the analysis manager
      mgr->SetGridHandler(alienHandler);
      
      if(gridTest) {
      	// speficy on how many files you want to run
      	alienHandler->SetNtestFiles(1);
      	// and launch the analysis
      	// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
      	alienHandler->SetRunMode("test");
      	mgr->StartAnalysis("grid");
	
      } else
    	{
    	  //full grid      
    	  alienHandler->SetRunMode("full");
    	  mgr->StartAnalysis("grid");
    	}
      
    }

}


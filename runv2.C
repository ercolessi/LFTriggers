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
#include "AliAnalysisTaskQAMultistrangev2.h"
#include "AliAnalysisTaskAO2Dconverter.h"

#endif

void runv2()
  
{
  // if you need to access data remotely from the GRID
  Bool_t grid = 1;  
  // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
  Bool_t gridTest = kTRUE;
  // if the data are MC
  Bool_t isMCdata=kFALSE;
  Bool_t isMC=kTRUE; // if you have the kinematics information
  if(isMC) isMCdata=kTRUE;
  
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

  // Create the analysis manager and ESD handler
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  AliESDInputHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  // add interface to MC if requested
  if(isMC){
    AliMCEventHandler* mcH = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcH);
  }  

  Bool_t lVarFalse = kFALSE ;
  Bool_t lVarTrue = kTRUE ;  

  //PhysicsSelection Configuration
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* ps =  reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C(kTRUE)"));
  
  //MultSelection
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask* ms =  reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
  ms->SetAddInfo(kTRUE);

  //PID Response
  //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse* pid = reinterpret_cast<AliAnalysisTaskPIDResponse*>(gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(kTRUE)"));
  
  //Weak Decay Vertexer
  //gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/Run2/macros/AddTaskWeakDecayVertexer.C");
  AliAnalysisTaskWeakDecayVertexer * taskWDV = reinterpret_cast<AliAnalysisTaskWeakDecayVertexer*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/Run2/macros/AddTaskWeakDecayVertexer.C"));
  taskWDV->SetRevertexAllEvents(kTRUE);
  taskWDV->SetUseImprovedFinding();
  taskWDV->SetupLooseVertexing();
  taskWDV->SetForceResetV0s(kTRUE);
  taskWDV->SetForceResetCascades(kTRUE);

  //AO2D converter
  AliAnalysisTaskAO2Dconverter * taskAO2DConverter = reinterpret_cast<AliAnalysisTaskAO2Dconverter*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/RUN3/AddTaskAO2Dconverter.C"));
  taskAO2DConverter->SetMCMode();
  taskAO2DConverter->SetTruncation(true);
  taskAO2DConverter->SetCompression(501);
  taskAO2DConverter->SetMaxBytes(250000000);

  //QA 
  gROOT->LoadMacro("AliAnalysisTaskQAMultistrangev2.cxx++g");
  gROOT->LoadMacro("AddTaskQAMultistrangev2.C");
  AliAnalysisTaskQAMultistrangev2 * task = reinterpret_cast<AliAnalysisTaskQAMultistrangev2*>(gInterpreter->ExecuteMacro("AddTaskQAMultistrangev2.C"));
  task -> SetIsMC(true);
  task -> SetAnalysisType("ESD");
  task -> SetQualityCutTPCrefit(false);
  task -> SetQualityCutnTPCcls(false);
  task->SetQualityCutMinnTPCcls(0);
 
    

  if(!mgr->InitAnalysis()) return;
  mgr->SetDebugLevel(2);
  mgr->PrintStatus();

  if (grid)
    {
      AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
      alienHandler->SetOverwriteMode();

      alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        
      alienHandler->SetAPIVersion("V1.1x");
      alienHandler->SetAliPhysicsVersion("vAN-20201212_ROOT6-1");   
      alienHandler->SetAnalysisMacro("AnalysisLeadingMC.C"); 

      alienHandler->SetExecutable("MCrunLeadingESD.sh");
      alienHandler->SetTTL(36000);
      alienHandler->SetJDLName("MCrunLeadingESD.jdl");
      
      alienHandler->SetOutputToRunNo(kTRUE);
      alienHandler->SetKeepLogs(kTRUE);
   
      alienHandler->SetTerminateFiles("event_stat.root");  //to have the output file of the Physics Selection class
      alienHandler->SetInputFormat("xml-single");
      alienHandler->SetPrice(1);      
      alienHandler->SetSplitMode("se");
       
      alienHandler->SetAdditionalLibs("AliAnalysisTaskQAMultistrangev2.cxx AliAnalysisTaskQAMultistrangev2.h");
      alienHandler->SetAnalysisSource("AliAnalysisTaskQAMultistrangev2.cxx ");
    
      alienHandler->SetGridDataDir("/alice/sim/2020/LHC20i2b");

      alienHandler->SetDataPattern("*ESDs.root");//("pass5_lowIR/*ESDs.root");
      alienHandler->SetRunPrefix("");   // real data: 000
    
      Long_t runnumbers[] = {288909};
      Int_t numberofruns = sizeof(runnumbers)/sizeof(Long_t);
      for(int atrun = 0;atrun<numberofruns;atrun++){
          alienHandler->AddRunNumber(runnumbers[atrun]);
      }
      
      alienHandler->SetOutputToRunNo();

      alienHandler->SetGridWorkingDir("analysisV0_LHC20i2b");
      alienHandler->SetGridOutputDir("OutputDir");
      alienHandler->SetMaxMergeStages(6);
      alienHandler->SetSplitMaxInputFileNumber(100);
      alienHandler->SetMergeViaJDL(kTRUE);
      //When all your merging jobs have finished running,
      //there is one step to be taken still, which is downloading the output
      //of the merging jobs. This is done by changing SetMergeViaJDL(kFALSE)
      //and running one last time
       
      // Connect plug-in to the analysis manager
      mgr->SetGridHandler(alienHandler);
      
      if(gridTest) {
      	// speficy on how many files you want to run
      	alienHandler->SetNtestFiles(5);
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


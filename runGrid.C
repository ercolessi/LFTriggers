{
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
    gSystem->Load("libCORRFW");
    gSystem->Load("libOADB");
    gSystem->Load("libPWGLFSTRANGENESS");
    gSystem->Load("libRUN3");
    gSystem->Load("AliAnalysisTaskQAMultistrangev2.h");
  
    // Use AliRoot includes to compile our task
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    gROOT->ProcessLine(".include AliAnalysisTaskQAMultistrangev2.h");
    
    // Create and configure the alien handler plugin
    gROOT->LoadMacro("CreateAlienHandler.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandler();
    if (!alienHandler) return;
    
    // Create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
    mgr->SetDebugLevel(99);
    
    // Connect plug-in to the analysis manager
    mgr->SetGridHandler(alienHandler);
    
    //Create input handler
    AliESDInputHandler* esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);
    
    AliMCEventHandler* mcInputHandler = new AliMCEventHandler();  //remove comment if using MC data
    mgr->SetMCtruthEventHandler(mcInputHandler);
    
    //---CDB connect---------------------------------------------------------
    //gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
    //AliTaskCDBconnect* task = AddTaskCDBconnect();
    //task->SetFallBackToRaw(kTRUE);
    //------------------------------------------------------------------------
    
    Bool_t lVarFalse = kFALSE ;
    Bool_t lVarTrue = kTRUE ;
    
    //---AliPIDResponse--------------------------------------------------------
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    Bool_t isMC=lVarTrue; // kTRUE in case of MC ; kFALSE is the default
    AddTaskPIDResponse(isMC);
    //-------------------------------------------------------------------------
    
    //---AliPhysicsSelection---------------------------------------------------
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    //Arguments: this is MC, yes; Reject Background, yes; Don't necessarily compute BG
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);
    // physSelTask->SetAnalyzeMC();
    //-------------------------------------------------------------------------
    
    //-------------------------------------------------------------------------
    // AliCentrality
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    // create centrality task, add to manager
    AliMultSelectionTask* centralityTask = AddTaskMultSelection(lVarTrue);
    centralityTask->SetAddInfo(kTRUE);
    //centralityTask->SetAddInfo(kTRUE);
    //centralityTask->SetFilterMB(kFALSE);
    //TGrid::Connect("alien://");
    //centralityTask->SetAlternateOADBforEstimators("LHC15o-SuperCalibMC-HIJING");
    //centralityTask->SetOADB("alien:///alice/cern.ch/user/d/ddobrigk/OADB-LHC18q_295586_finnerBinning.root");
     
    //centralityTask-> SetAlternateOADBFullManualBypass("OADB-LHC15o_246087.root");
    //centralityTask->SetAlternateOADBforEstimators("LHC15o-DefaultMC-HIJING");
    //-------------------------------------------------------------------------
    
    //-------------------------------------------------------------------------
    // Weak Decay Vertexer
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/Run2/macros/AddTaskWeakDecayVertexer.C");
    AliAnalysisTaskWeakDecayVertexer *taskWDV = AddTaskWeakDecayVertexer();
    taskWDV->SetRevertexAllEvents(kTRUE);
    taskWDV->SetUseImprovedFinding();
    taskWDV->SetupLooseVertexing();
    taskWDV->SetForceResetV0s(kTRUE);
    taskWDV->SetForceResetCascades(kTRUE);
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // AO2D Converter
    gROOT->LoadMacro("$ALICE_PHYSICS/RUN3/AddTaskAO2Dconverter.C");
    AliAnalysisTaskAO2Dconverter *taskAO2DConverter = AddTaskAO2Dconverter();
    taskAO2DConverter->SetMCMode();
    taskAO2DConverter->SetTruncation(true);
    taskAO2DConverter->SetCompression(501);
    taskAO2DConverter->SetMaxBytes(250000000);
    //taskAO2DConverter->SetEMCALAmplitudeThreshold(0.075);
    //-------------------------------------------------------------------------
    
    cout<<" Compiling Task..."<<endl;
    // compile our task
    
    //gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/macros/AddTaskExtractCascade.C");
    
    //gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/Run2/macros/AddTaskStrangenessVsMultiplicityRsnLikeBgSub.C");
    //
    //gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/Run2/macros/AddTaskStrEffStudy.C");
    
    //gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/Run2/macros/AddTaskWeakDecayVertexerRev2.C");
    
    Short_t       lCollidingSystems=0  ;/*0 = pp, 1 = AA*/
    const TString lMasterJobSessionFlag = "";
    
    // create our task
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    //AddTask: Should take care of everything transparently...
    gROOT->LoadMacro("AliAnalysisTaskQAMultistrangev2.cxx++g");
    gROOT->LoadMacro("AddTaskQAMultistrangev2.C");
    AliAnalysisTaskQAMultistrangev2 *task = AddTaskQAMultistrangev2();
    //gROOT->LoadMacro("/home/fercoles/sw/alice/AliPhysics/PWGLF/QATasks/AliAnalysisTaskQAMultistrange.cxx++g");
    //gROOT->LoadMacro("/home/fercoles/sw/alice/AliPhysics/PWGLF/QATasks/macros/AddTaskQAMultistrange.C");
    //AliAnalysisTaskQAMultistrange *task = AddTaskQAMultistrange(lVarTrue);

    task -> SetIsMC(true);
    task -> SetAnalysisType("ESD");
    task -> SetQualityCutTPCrefit(false);
    task -> SetQualityCutnTPCcls(false);
    task->SetQualityCutMinnTPCcls(0);
        
    if (!mgr->InitAnalysis())
        return;
    
    mgr->PrintStatus();
    // Start analysis in grid.
    mgr->StartAnalysis("grid");
};


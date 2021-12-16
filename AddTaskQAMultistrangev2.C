AliAnalysisTaskQAMultistrangev2 *AddTaskQAMultistrangev2(Bool_t isMC = kFALSE) {

   // Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskQAMultistrangev2", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskQAMultistrangev2", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Create and configure the task
   AliAnalysisTaskQAMultistrangev2 *taskcheckcascade = new AliAnalysisTaskQAMultistrangev2("TaskCheckCascade");
     taskcheckcascade->SetIsMC                       (isMC);
     taskcheckcascade->SetAnalysisType               (type);

   mgr->AddTask(taskcheckcascade);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================

   // User file name (if need be)
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outputFileName += ":PWGLFStrangeness.outputCheckCascade";

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("fListHistMultistrangev2QA",
                                                             TList::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );
   
   mgr->ConnectInput( taskcheckcascade, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskcheckcascade, 1, coutput1);
 
   return taskcheckcascade;
}   

//#include "AliAnalysisManager.h"

//#include "AliAnalysisTaskHMTFMCMultEst.h"

AliAnalysisTaskHMTFMCMultEst *AddTaskHMTFMCMultEst() {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHMTFMCMultEst", "No analysis manager to connect to.");
    return NULL;
  }

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
    "Sums",
    TList::Class(), 
    AliAnalysisManager::kOutputContainer,
    mgr->GetCommonFileName() );
  AliAnalysisTaskHMTFMCMultEst *multEstTask = new AliAnalysisTaskHMTFMCMultEst("AliAnalysisTaskHMTFMCMultEst");
  if (!multEstTask) {
      Error("CreateTasks", "Failed to make my task!");  
      return NULL;
  }
  // add estimators:
  multEstTask->AddEstimator("EtaLt05");

  mgr->AddTask(multEstTask);

  mgr->ConnectInput(multEstTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(multEstTask, 1, coutput1);

  return multEstTask;
}

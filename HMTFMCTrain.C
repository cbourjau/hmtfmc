/* 
Run locally:
$ runTrain --class=HMTFMCTrain --name=myname_lite --include=. --url="lite://$PWD/input//?pattern="galice.root"&workers=2&mc&recursive#TE" --type=ESD

For debugging
 rm -rf myname_lite && runTrain --class=HMTFMCTrain --name=myname_lite --include=. --url="local://\$PWD/input/root_archive_00001/?pattern="galice.root"&mc#TE" ; cat build.log
*/
  

#ifndef __CINT__
# include <AliAnalysisManager.h>
#else 
class AliAnalysisManager;
#endif

#include "MultiplicityEstimators.h"
#include "AliAnalysisTaskHMTFMC.h"
#include "TrainSetup.C"
#include "ParUtilities.C"


class HMTFMCTrain : public TrainSetup
{
public:
  HMTFMCTrain(const char* name="HMTFMCTrain") : TrainSetup(name)
  {
    // Set Generator here?
    // fOptions.Set("type", "AOD"); // AOD input
    fOptions.Set("type", "ESD"); // ESD input
    //fOptions.Add("parameter", "VALUE", "Help on parameter", "value");
  }
  virtual AliVEventHandler* CreateMCHandler(UShort_t /*type*/, bool mc)
  {
    if (!mc) {
      Fatal("CreateMCHandler", "MC option _must_ be set");
      return 0;
    }
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mcHandler->SetReadTR(false); 
    return mcHandler;
  }
  AliVEventHandler* CreateOutputHandler(UShort_t) { return 0; }
  void CreatePhysicsSelection(Bool_t, AliAnalysisManager*) {}
  void CreateCentralitySelection(Bool_t) {}
  void CreateTasks(AliAnalysisManager*  mgr)
  {
    std::cout << "Creating Task" << std::endl;
    fRailway->LoadLibrary("pythia6_4_25");
    fRailway->LoadAux("MultiplicityEstimators.h");
    fRailway->LoadSource("MultiplicityEstimators.cxx");
    fRailway->LoadAux("AliAnalysisTaskHMTFMC.h");
    fRailway->LoadSource("AliAnalysisTaskHMTFMC.cxx");
    AliAnalysisManager::SetCommonFileName("MC_estimators.root");
 
    // AliAnalysisTaskHMTFMC *task = new AliAnalysisTaskHMTFMC("TaskdNdeta");
    Long_t ret = gROOT->ProcessLine("new AliAnalysisTaskHMTFMC(\"TaskdNdeta\")");
    if (!ret) {
      Error("CreateTasks", "Failed to make my task!");  
      return;
    }
    AliAnalysisTaskSE* task = reinterpret_cast<AliAnalysisTaskSE*>(ret);
    // add estimators here:
    // task->AddEstimator("EtaLt05");
    gROOT->ProcessLine(Form("((AliAnalysisTaskHMTFMC)%p)->AddEstimator(\"EtaLt05\")",
			    task));
    gROOT->ProcessLine(Form("((AliAnalysisTaskHMTFMC)%p)->AddEstimator(\"EtaLt08\")",
			    task));
    gROOT->ProcessLine(Form("((AliAnalysisTaskHMTFMC)%p)->AddEstimator(\"EtaLt15\")",
			    task));
    gROOT->ProcessLine(Form("((AliAnalysisTaskHMTFMC)%p)->AddEstimator(\"Eta08_15\")",
			    task));

    mgr->AddTask(task);
        
    AliAnalysisDataContainer* sums = 
      mgr->CreateContainer("Sums", TList::Class(), 
  			  AliAnalysisManager::kOutputContainer,
  			  AliAnalysisManager::GetCommonFileName());
    AliAnalysisDataContainer* results = // Needed for output from Terminate
      mgr->CreateContainer("Results", TList::Class(), 
  			  AliAnalysisManager::kParamContainer, // Important!
  			  AliAnalysisManager::GetCommonFileName());
    
    mgr->ConnectOutput(task, 1, sums);
    mgr->ConnectOutput(task, 2, results);


    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());  
  }
  const char* ClassName() const { return "HMTFMCTrain"; }

};

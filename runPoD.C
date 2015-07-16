// #include "TList.h"
// #include "TNamed.h"
// #include "TProof.h"
// #include "TROOT.h"

// #include "AliAnalysisManager.h"
// #include "AliAODInputHandler.h"

// #include "AliAnalysisTaskHMTFMCMultEst.h"

void runPoD(
  TString dataset = "Find;"
                    "BasePath=/alice/cern.ch/user/p/pwgpp_mc/2015/17_Week/TestMultiplicity/Test2/Test_1M_events_iter1/000%/;"
                    "FileName=root_archive.zip;"
                    "Anchor=galice.root;"
                    "Tree=/TE;"
                    "Mode=remote;",  // <-- much faster dataset creation
  //Bool_t usePhysicsSelection = kTRUE,
  Int_t numEvents = 999999999,
  Int_t firstEvent = 0
) {

  // Not needed on the VAF
  //gEnv->SetValue("XSec.GSI.DelegProxy","2");
   
  TString extraLibs = "ANALYSIS:ANALYSISalice:pythia6_4_25"; // extraLibs = "ANALYSIS:OADB:ANALYSISalice:CORRFW:OADB:PWGmuon";

  TList *list = new TList(); 
  list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
  list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));  // important: creates token on every PROOF worker

  // Note the difference between CAF and VAF
  //TProof::Open("alice-caf.cern.ch");
  TProof::Open("pod://");

  // Check the dataset before running the analysis!
  gProof->ShowDataSet( dataset.Data() );
  //return;  // <-- uncomment this to test search before running the analysis!

  // Not needed on the VAF
  //gProof->EnablePackage("VO_ALICE@AliRoot::v5-04-81-AN", list);

  // A single AliRoot package for *all* AliRoot versions: new on VAF
  TString aliceVafPar = "/afs/cern.ch/alice/offline/vaf/AliceVaf.par";
  gProof->UploadPackage(aliceVafPar.Data());
  gProof->EnablePackage(aliceVafPar.Data(), list);  // this "list" is the same as always

  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train");

  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

  // DON'T use double '+' when running multiple times: it uselessly recompiles everything!
  gProof->Load("MultiplicityEstimators.cxx+");
  gProof->Load("AliAnalysisTaskHMTFMCMultEst.cxx+");  
  gROOT->LoadMacro("AddTaskHMTFMCMultEst.C");

  AliAnalysisTaskHMTFMCMultEst *simplePtTask = AddTaskHMTFMCMultEst();//usePhysicsSelection);
 
  if (!mgr->InitAnalysis()) return;

  mgr->StartAnalysis("proof", dataset, numEvents, firstEvent);
}

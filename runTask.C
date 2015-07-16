#ifndef __CINT__
#include <fstream>
#include "TROOT.h"
#include "TChain.h"

#include "AliAnalysisManager.h"
#endif

void runTask(const char * incollection = 0, const char * outfile = "MC_MultEst.root")
{
  // for running with root only
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so"); 

  // load analysis framework
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("$ALICE_ROOT/lib/libpythia6_4_25.so");
  gSystem->Load("$ALICE_ROOT/lib/libpythia6.so");

  TChain * chain = new TChain ("TE");
  if (incollection == 0) {
    // if no collection is given, check for a file in workdir
    chain->Add("galice.root");
  }
  // else if (TString(incollection).Contains("xml")){
  //   TGrid::Connect("alien://");
  //   TAlienCollection * coll = TAlienCollection::Open (incollection);
  //   while(coll->Next()){
  //     chain->Add(TString("alien://")+coll->GetLFN());
  //   }
  {
    ifstream file_collect(incollection);
    TString line;
    while (line.ReadLine(file_collect) ) {
      chain->Add(line.Data());
    }
  }
  chain->GetListOfFiles()->Print();
  //
  // for includes use either global setting in $HOME/.rootrc
  // ACLiC.IncludePaths: -I$(ALICE_ROOT)/include
  // or in each macro
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("HMTFMCMultEst");

  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  mgr->SetCommonFileName(outfile);
  
  // Enable MC event handler
  AliMCEventHandler* handler = new AliMCEventHandler;
  handler->SetReadTR(kFALSE);
  mgr->SetMCtruthEventHandler(handler);

  // Create containers for input/output before adding tasks
  //AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  // Create and add tasks
  gROOT->LoadMacro("MultiplicityEstimators.cxx+g");
  gROOT->LoadMacro("AliAnalysisTaskHMTFMCMultEst.cxx+g");  
  gROOT->LoadMacro("AddTaskHMTFMCMultEst.C");
  AliAnalysisTask *task1 = AddTaskHMTFMCMultEst();
  mgr->AddTask(task1);

  // Enable debug printouts
  mgr->SetDebugLevel(0);

  if (!mgr->InitAnalysis()) return;

  mgr->PrintStatus();

  mgr->StartAnalysis("local", chain);
}

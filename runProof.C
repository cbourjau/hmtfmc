#ifndef __CINT__
#include <fstream>
#include "TROOT.h"
#include "TProof.h"
#include "TChain.h"
#include "TObject.h"

#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#endif

void loadLibs(const TString extralibs){
  // for running with root only
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so"); 

  // load extra analysis libs
  TIter it(extralibs.Tokenize(":"));
  TObjString *lib = 0;
  while ((lib = dynamic_cast<TObjString *>(it()))){
    gSystem->Load(lib->String());
  }
}

TChain* makeChain(const char* incollection) {
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
  ifstream file_collect(incollection);
  TString line;
  while (line.ReadLine(file_collect) ) {
    chain->Add(line.Data());
  }
  return chain;
}

void loadAnalysisFiles(const TString files, TString runmode ) {
  // load files from colon (:) separated string
  // first file in string is loaded first (OBS: dependencies)
  // The function does the right thing regardless if local, lite or pod
  // Do not use ++ in the string to avoid unnecessary recompiles on nodes
  
  TIter it(files.Tokenize(":"));
  TObjString *lib = 0;
  while ((lib = dynamic_cast<TObjString *>(it()))){
    Int_t status;
    if (runmode.BeginsWith("local")){
      status = gROOT->LoadMacro(lib->String());
    }
    else {
      status = gProof->Load(lib->String());
    }
    if (status !=0){
	std::cout << "Error loading " << lib->String() << std::endl;
    }
  }
}

void runProof(const TString runmode_str  = "lite",
	      Int_t max_events = -1,
	      const char * incollection = "./input/input_files.dat",
	      const char * analysisName = "hmtf_mc_mult",
	      const char * aliceExtraLibs=("libANALYSIS:"
					   "libANALYSISalice:"
					   "libpythia6_4_25"),
	      const char * analysisFiles=("MultiplicityEstimators.cxx+g:"
					  "AliAnalysisTaskHMTFMCMultEst.cxx+g"),
	      const TString addAnalysisFiles=("AddTaskHMTFMCMultEst.C"))
{
  if(!(runmode_str.BeginsWith("local") ||
       runmode_str.BeginsWith("lite") ||
       runmode_str.BeginsWith("pod"))){
    std::cout << "invalid runmode given" << std::endl;
    return;
  }
  if (runmode_str.BeginsWith("lite")) TProof::Open("lite://");
  else if (runmode_str.BeginsWith("pod")) TProof::Open("pod://");
  
  loadLibs(aliceExtraLibs);
  TChain *chain = makeChain(incollection);

  gSystem->AddIncludePath("-I$ALICE_ROOT/include");

  // Create  and setup the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager(analysisName);
  mgr->SetCommonFileName(TString(analysisName) + TString(".root"));

  // Set Handlers
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  AliMCEventHandler* handler = new AliMCEventHandler;
  handler->SetReadTR(kFALSE);
  mgr->SetMCtruthEventHandler(handler);

  loadAnalysisFiles(analysisFiles, runmode_str);

  // Add tasks
  TIter it(addAnalysisFiles.Tokenize(":"));
  TObjString *adder = 0;
  while ((adder = dynamic_cast<TObjString *>(it()))){
    mgr->AddTask(reinterpret_cast<AliAnalysisTask*>(gROOT->Macro(adder->String())));    
  }

  // Enable debug printouts
  mgr->SetDebugLevel(0);

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if (runmode_str.BeginsWith("lite")) mgr->StartAnalysis("proof", chain, max_events);
  else if (runmode_str.BeginsWith("local")) mgr->StartAnalysis("local", chain, max_events);
}

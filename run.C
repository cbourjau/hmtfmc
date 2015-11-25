#ifndef __CINT__
#include <fstream>
#include "TROOT.h"
#include "TProof.h"
#include "TChain.h"
#include "TObject.h"
#include "TList.h"

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"
#endif

void loadLibs(const TString extralibs, const TString runmode){
  if ((runmode.BeginsWith("local"))){
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
      if (gSystem->Load(lib->String()) < 0) {
	std::cout << "Error loading " << lib->String() << std::endl;
      }
    }
  }
  else{
    TList *list = new TList();
    list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extralibs.Data()));
    TString alicePar;
    
    if ((runmode.BeginsWith("lite"))){
      alicePar  = "$ALICE_ROOT/ANALYSIS/macros/AliRootProofLite.par";
    }
    
    if (runmode.BeginsWith("pod")){
      // A single AliRoot package for *all* AliRoot versions: new on VAF
      alicePar = "/afs/cern.ch/alice/offline/vaf/AliceVaf.par";
      list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));  // important: creates token on every PROOF worker
    }
    gProof->UploadPackage(alicePar.Data());
    gProof->EnablePackage(alicePar.Data(), list);  // this "list" is the same as always
    std::cout << "enabled packages: " << std::endl;
    gProof->ShowEnabledPackages();
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

TString getDatasetString(const char* incollection) {
  TString dataset("");
  if (TString(incollection).BeginsWith("pythia")){
    TString dataset = ("Find;"
		       "BasePath=/alice/cern.ch/user/p/pwgpp_mc/2015/17_Week/TestMultiplicity/Test2/Test_1M_events_iter1/%/;"
		       "FileName=root_archive.zip;"
		       "Anchor=galice.root;"
		       "Tree=/TE;"
		       "Mode=remote;"  // <-- much faster dataset creation
		       );
  }
  else if (TString(incollection).BeginsWith("phojet")){
    TString dataset = ("Find;"
		       "BasePath=/alice/sim/2014/LHC14j3c/146805/%/;"
		       "FileName=root_archive.zip;"
		       "Anchor=galice.root;"
		       "Tree=/TE;"
		       "Mode=remote;"  // <-- much faster dataset creation
		       );
  }
  return dataset;
}

void setUpIncludes(TString runmode) {
  // Set the include directories
  // gProof is without -I :P and gSystem has to be set for all runmodes
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  if (runmode.BeginsWith("lite")){
    gProof->AddIncludePath("$ALICE_ROOT/include", kTRUE);
  }
  
}
void loadAnalysisFiles(const TString files, TString runmode ) {
  // load files from colon (:) separated string
  // first file in string is loaded first (OBS: dependencies)
  // The function does the right thing regardless if local, lite or pod
  // Do not use ++ in the string to avoid unnecessary recompiles on nodes

  // loop over analysis files given by the user
  TIter it(files.Tokenize(":"));
  TObjString *lib = 0;
  while ((lib = dynamic_cast<TObjString *>(it()))){
    Int_t status;
    if (runmode.BeginsWith("local")) status = gROOT->LoadMacro(lib->String());
    else status = gProof->Load(lib->String());
    if (status !=0){
	std::cout << "Error loading " << lib->String() << std::endl;
    }
  }
}

void run(const TString runmode_str  = "lite",
	 Int_t max_events = -1,
	 //Int_t first_event= 0,
	 Int_t debug = 0,
	 const char * incollection = "./pythia/input_files.dat",
	 const char * analysisName = "hmtf_mc_mult",
	 const char * aliceExtraLibs=("libPWGHMTF:"
				      "libpythia6_4_25:"
				      "libAliPythia6"
				      ),
	 const char * analysisFiles="",//"AliMultiplicityEstimators.cxx+:"
	                               //"AliAnalysisTaskHMTFMCMultEst.cxx+"),
	 const TString adderFiles=("$ALICE_PHYSICS/PWG/HMTF/macros/AddTaskHMTFMCMultEst.C(\"\"):"
				   "$ALICE_PHYSICS/PWG/HMTF/macros/AddTaskHMTFMCMultEst.C(\"InelGt0\"):"
				   "$ALICE_PHYSICS/PWG/HMTF/macros/AddTaskHMTFMCMultEst.C(\"SphericityGt09\")"
				   ))
{
  if(!(runmode_str.BeginsWith("local") ||
       runmode_str.BeginsWith("lite") ||
       runmode_str.BeginsWith("pod"))){
    std::cout << "invalid runmode given" << std::endl;
    return;
  }
  // start proof if necessary
  if (runmode_str.BeginsWith("lite")) TProof::Open("lite://");
  else if (runmode_str.BeginsWith("pod")) TProof::Open("pod://");

  setUpIncludes(runmode_str);
  loadLibs(aliceExtraLibs, runmode_str);

  // Create  and setup the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager(analysisName);

  // Enable debug printouts
  mgr->SetDebugLevel(debug);
  

  mgr->SetCommonFileName(TString(analysisName) + TString(".root"));

  // Set Handlers
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  AliMCEventHandler* handler = new AliMCEventHandler;
  handler->SetReadTR(kFALSE);
  mgr->SetMCtruthEventHandler(handler);

  loadAnalysisFiles(analysisFiles, runmode_str);

  // Add tasks
  TIter it(adderFiles.Tokenize(":"));
  TObjString *adder = 0;
  while ((adder = dynamic_cast<TObjString *>(it()))){
    // adder files include logic to add the task to the manager, no need to do it again here
    gROOT->Macro(adder->String());    
  }

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();

  if (!(runmode_str.BeginsWith("pod"))) {
      // Process with chain
      TChain *chain = makeChain(incollection);
      if (!chain) {
	std::cout << "Dataset is empty!" << std::endl;
	return;
      }
      if (runmode_str.BeginsWith("lite")) mgr->StartAnalysis("proof", chain, max_events);

      else if (runmode_str.BeginsWith("local")) {
	// in order to be inconvinient, aliroot does interprete max_events differently for proof and local :P
	if (max_events == -1) mgr->StartAnalysis("local", chain);
	else mgr->StartAnalysis("local", chain, max_events);
      }
    }
    else {
      // process with dataset string
      TString dataset = getDatasetString(incollection);
      // Check the dataset before running the analysis!
      gProof->ShowDataSet( dataset.Data() );
      mgr->StartAnalysis("proof", dataset, max_events, 0 /*first_event*/);
    }
}

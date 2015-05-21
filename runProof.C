#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "AliAnalysisGrid.h"
#include "TSystem.h"
#include "TROOT.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisGrid.h"
#include "AliVEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisAlien.h"
#include "AliMCEventHandler.h"
#include "AliPhysicsSelectionTask.h"
#include "TRegexp.h"
#include "AliESDInputHandler.h"

#endif
void LoadLibs();

class AliAnalysisGrid; // Needed, otherwise the interpreter freaks out with the definition of CreateAlienHandler
AliAnalysisGrid* CreateAlienHandler(const char *mode, const char * rootVersion, const char * alirootVersion, const char *aliphysicsVersion, const char *proofcluster, const char *proofdataset, const char * taskname);


//______________________________________________________________________________
void runProof(
              const char *proofdataset = " /default/ilakomov/TestDownScale1", // path to dataset on proof cluster, for proof analysis
              const char *mode = "full", // Full & Test work for proof
              const Long64_t nentries = 200000, 
              const Long64_t firstentry = 0, 
              const char *proofcluster = "alice-caf.cern.ch", // which proof cluster to use in proof mode
              const char * rootVersion = "v5-34-08",
              const char * alirootVersion = "v5-06-09",
              const char * aliphysicsVersion = "vAN-20150323"
              // const char * rootVersion = "v5-34-08",
              // const char * alirootVersion = "v5-06-16",
              // const char * aliphysicsVersion = "vAN-20150503"
	      )
{

  
  // load libraries
  LoadLibs();
  // add aliroot indlude path
  //  gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_ROOT")));
        
  // analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager("AliAnalysisTaskHMTFMC");
    
  // create the alien handler and attach it to the manager
  AliAnalysisGrid *plugin = CreateAlienHandler(mode, rootVersion, alirootVersion, aliphysicsVersion, proofcluster, proofdataset, "AliAnalysisTaskHMTFMC"); 
  mgr->SetGridHandler(plugin);
    
  AliVEventHandler* iH = new AliESDInputHandler(); // WARNING: this needs to be changed for AOD!
  //    AliAODInputHandler* iH = new AliAODInputHandler();
  //    iH->SetInactiveBranches("tracks. vertices. v0s. cascades. jets. caloClusters. fmdClusters. pmdClusters. dimuons. AliAODZDC");
  //    iH->SetInactiveBranches("*");
  //    iH->SetCheckStatistics(kTRUE);
  mgr->SetInputEventHandler(iH);
        
  // mc event handlerrunEx01.C
  AliMCEventHandler* mchandler = new AliMCEventHandler();
  // Not reading track references
  mchandler->SetReadTR(kFALSE);
  mgr->SetMCtruthEventHandler(mchandler);


  // create task
  gROOT->LoadMacro("MultiplicityEstimators.cxx+g");
  gROOT->LoadMacro("AliAnalysisTaskHMTFMC.cxx+g");
  AliAnalysisTask *task = new AliAnalysisTaskHMTFMC("TaskdNdeta");
  // Enable MC event handler
  task->Print();
        
  //
  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("coutput", TList::Class(),    AliAnalysisManager::kOutputContainer, "dndeta.root");

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);

  // enable debug printouts
  mgr->SetDebugLevel(10);
  //    mgr->SetNSysInfo(100);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
   
  // start analysis
  Printf("Starting Analysis....");
  mgr->StartAnalysis("proof",nentries,firstentry);
}

//______________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler(const char *mode, const char * rootVersion, const char * alirootVersion, const char *aliphysicsVersion, const char *proofcluster, const char *proofdataset, const char * taskname)
{
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetRunMode(mode);

  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion(rootVersion);
  plugin->SetAliROOTVersion(alirootVersion);
  plugin->SetAliPhysicsVersion(aliphysicsVersion);
  //----------------------------------------------------------
  //---      PROOF MODE SPECIFIC SETTINGS         ------------
  //---------------------------------------------------------- 
  // Proof cluster
  plugin->SetProofCluster(proofcluster);
  // Dataset to be used   
  plugin->SetProofDataSet(proofdataset);
  // May need to reset proof. Supported modes: 0-no reset, 1-soft, 2-hard
  plugin->SetProofReset(0);
  // May limit number of workers
  plugin->SetNproofWorkers(1);
  // May limit the number of workers per slave
  plugin->SetNproofWorkersPerSlave(1);   // FIXME: add a parameter for this?
  // May use a specific version of root installed in proof
  plugin->SetROOTVersion(rootVersion);
  // May set the aliroot mode. Check http://aaf.cern.ch/node/83 
  plugin->SetAliRootMode("default"); // Loads AF libs by default
  // May request ClearPackages (individual ClearPackage not supported)
  plugin->SetClearPackages(kFALSE);
  // Plugin test mode works only providing a file containing test file locations, used in "local" mode also
  plugin->SetFileForTestMode("files.txt"); // file should contain path name to a local directory containg *ESDs.root etc
  // Request connection to alien upon connection to grid
  plugin->SetProofConnectGrid(kFALSE);
  // code
  plugin->SetAnalysisSource(Form("%s.cxx %s.cxx",
				 "MultiplicityEstimators",
				 taskname));
  std::cout << "Anlysis source files: " << std::endl;
  std::cout << plugin->GetAnalysisSource() << std::endl;
 
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  // List is read from the front -> important for dependencies!
  plugin->SetAdditionalLibs(Form("%s.h %s.cxx %s.h %s.cxx",
				 "MultiplicityEstimators", "MultiplicityEstimators",
				 taskname, taskname));
  std::cout << "Additional libs: " << std::endl;
  std::cout << plugin->GetAdditionalLibs() << std::endl;
  
  // Other PROOF specific parameters
  plugin->SetProofParameter("PROOF_UseMergers","-1");
  printf("Using: PROOF_UseMergers   : %s\n", plugin->GetProofParameter("PROOF_UseMergers"));
  plugin->Print();
  return plugin;
}

void LoadLibs() {
  gSystem->Load("libCore.so");  
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libProof");
  gSystem->Load("libMatrix");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTENDER");
  gSystem->Load("libCORRFW");
  gSystem->Load("libMinuit");

}

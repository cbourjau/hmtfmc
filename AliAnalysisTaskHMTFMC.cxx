#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"

#include "AliAnalysisTaskHMTFMC.h"
#include "TGraphErrors.h"
#include "AliLog.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"

#include <iostream>
#include "AliPDG.h"
#include "AliGenDPMjetEventHeader.h"
#include "TH2F.h"


#include "MultiplicityEstimators.h"

using namespace std;

ClassImp(AliAnalysisTaskHMTFMC)

AliAnalysisTaskHMTFMC::AliAnalysisTaskHMTFMC()
: AliAnalysisTaskSE(), fMyOut(0), fEstimatorsList(0) 
{

}

//________________________________________________________________________
AliAnalysisTaskHMTFMC::AliAnalysisTaskHMTFMC(const char *name) 
  : AliAnalysisTaskSE(name), fMyOut(0), fEstimatorsList(0), fEstimatorNames(0)
{
  AliPDG::AddParticlesToPdgDataBase();

  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

void AliAnalysisTaskHMTFMC::AddEstimator(const char* n)
{
  if (!fEstimatorNames.IsNull()) fEstimatorNames.Append(",");
  fEstimatorNames.Append(n);
}
    
void AliAnalysisTaskHMTFMC::InitEstimators()
{
  fEstimatorsList = new TList;
  fEstimatorsList->SetOwner();
  fEstimatorsList->SetName("estimators");
  
  TObjArray* arr = fEstimatorNames.Tokenize(",");
  TObject*   obj = 0;
  TIter      next(arr);
  while ((obj = next())) {
    MultiplicityEstimatorBase* e = MakeEstimator(obj->GetName());
    std::cout << "Init estimtor: " << obj->GetName() << std::endl;
    fEstimatorsList->Add(e);
  }
}
//________________________________________________________________________
MultiplicityEstimatorBase*
AliAnalysisTaskHMTFMC::MakeEstimator(const TString& name)
{
  if (name.BeginsWith("EtaLt05")) return new EtaBase("EtaLt05", "| #eta| #leq 0.5", -0.5, 0.0, 0.0, 0.5);
  if (name.BeginsWith("EtaLt08")) return new EtaBase("EtaLt08", "| #eta| #leq 0.8", -0.8, 0.0, 0.0, 0.8);
  if (name.BeginsWith("EtaLt15")) return new EtaBase("EtaLt15", "| #eta| #leq 1.5", -1.5, 0.0, 0.0, 1.5);
  if (name.BeginsWith("Eta08_15")) return new EtaBase("Eta08_15", "0.8 #leq | #eta| #leq 1.5",
						      -1.5, -0.8, 0.8, 1.5);
  if (name.BeginsWith("V0A")) return new EtaBase("V0A", "2.8 #leq #eta #leq 5.1",
						 0.0, 0.0, 2.8, 5.1);
  if (name.BeginsWith("V0C")) return new EtaBase("V0C", "-3.7 #leq #eta #leq -1.7",
						 -3.7, -1.7, 0.0, 0.0);
  if (name.BeginsWith("V0M")) return new EtaBase("V0M", "-3.7 #leq #eta #leq -1.7 || 2.8 #leq #eta #leq 5.1",
						 -3.7, -1.7, 2.8, 5.1);

  return 0;
}

//________________________________________________________________________
void AliAnalysisTaskHMTFMC::UserCreateOutputObjects()
{
  fMyOut = new TList();
  fMyOut->SetOwner();

  InitEstimators();

  TIter next(fEstimatorsList);
  MultiplicityEstimatorBase* e = 0;
  while ((e = static_cast<MultiplicityEstimatorBase*>(next()))) {
    e->RegisterHistograms(fMyOut);
    // putting estimators into a vector for easier looping in UserExec.
    // it is only available on the slaves
    festimators.push_back(e);
  }
  
  // Suppress annoying printout
  AliLog::SetGlobalLogLevel(AliLog::kError);
  
  PostData(1, fMyOut);
}

//________________________________________________________________________
void AliAnalysisTaskHMTFMC::UserExec(Option_t *) 
{
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
     Printf("ERROR: Could not retrieve MC event");
     return;
  }

  for (std::vector<MultiplicityEstimatorBase*>::size_type i = 0; i < festimators.size(); i++) {//estimator loop
    festimators[i]->PreEvent(mcEvent);
  }//estimator loop

  // Track loop for establishing multiplicity
  for (Int_t iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {
    AliMCParticle *track = (AliMCParticle*)mcEvent->GetTrack(iTrack);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTrack);
      continue;
    }
    if (mcEvent->Stack()->IsPhysicalPrimary(iTrack)){
      for (std::vector<MultiplicityEstimatorBase*>::size_type i = 0; i < festimators.size(); i++) { //estimator loop
	festimators[i]->ProcessTrackForMultiplicityEstimation(track);
      }//estimator loop
    }
  }  //track loop

  // Track loop with known multiplicity in each estimator
  for (Int_t iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {
    AliMCParticle *track = (AliMCParticle*)mcEvent->GetTrack(iTrack);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTrack);
      continue;
    }
    if (mcEvent->Stack()->IsPhysicalPrimary(iTrack)){
      for (std::vector<MultiplicityEstimatorBase*>::size_type i = 0; i < festimators.size(); i++) { //estimator loop
	festimators[i]->ProcessTrackWithKnownMultiplicity(track);
      }//estimator loop
    }
  }  //track loop

  // Increment eventcounters etc.
  for (std::vector<MultiplicityEstimatorBase*>::size_type i = 0; i < festimators.size(); i++) { //estimator loop
    festimators[i]->PostEvent();
  }//estimator loop

  // Post output data.
  PostData(1, fMyOut);
}      

//________________________________________________________________________
void AliAnalysisTaskHMTFMC::Terminate(Option_t *) 
{
  // recreates the fEstimatorsList
  InitEstimators();

  // This list is associated to a read only file
  fMyOut = static_cast<TList*> (GetOutputData(1));
  if (!fMyOut) {
    Error("Terminate", "Didn't get sum container");
    return;
  }
  TList* results = new TList;
  results->SetName("terminateResults");

  TIter nextEst(fEstimatorsList);
  MultiplicityEstimatorBase* e = 0;
  while ((e = static_cast<MultiplicityEstimatorBase*>(nextEst()))) {
    std::cout << "Terminating estimator " << e->GetTitle() << std::endl;
    e->Terminate(fMyOut, results);
  }
  PostData(2, results);
}

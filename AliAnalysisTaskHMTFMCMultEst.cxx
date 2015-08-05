#include <iostream>

#include "TCanvas.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TParticle.h"


#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliESDEvent.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliHeader.h"
#include "AliLog.h"
#include "AliMCEvent.h"
//#include "AliPDG.h"
#include "AliStack.h"
#include "AliVEvent.h"

#include "AliAnalysisTaskHMTFMCMultEst.h"
#include "MultiplicityEstimators.h"

using namespace std;

ClassImp(AliAnalysisTaskHMTFMCMultEst)

Int_t pdg_2_pid_enum(Int_t pdg) {
  if (pdg == 310) return kK0S;
  else if (pdg == 321) return kKPLUS;
  else if (pdg == -321) return kKMINUS;
  return -99999;
}


Bool_t IsPi0PhysicalPrimary(Int_t index, AliStack *stack)
{
  // pi0 are considered unstable and thus not primary in AliStack.
  // If the given track index is a pi0 and if that pi0 has no ancestors, it is a primary
  // Return kFALSE if the index points to anything but a pi0

  TParticle* p = stack->Particle(index);
  Int_t ist = p->GetStatusCode();

  // Initial state particle
  //if (ist > 1) return kFALSE; // pi0 are not initial state (?)
  // pi0 are unstable and thus this returned false in the original implementation.
  //if (!stack->IsStable(pdg)) return kFALSE;

  Int_t pdg = TMath::Abs(p->GetPdgCode());

  // The function is only for pi0's so I'm out of here if its not a pi'0 we are looking at
  if (pdg != kPi0) return kFALSE;
  
  if (index < stack->GetNprimary()) {
    // Particle produced by generator
    return kTRUE;
  }
  else {
    // Particle produced during transport
    Int_t imo =  p->GetFirstMother();
    TParticle* pm  = stack->Particle(imo);
    Int_t mpdg = TMath::Abs(pm->GetPdgCode());
    // Check for Sigma0 
    if ((mpdg == 3212) &&  (imo <  stack->GetNprimary())) return kTRUE;
    
    // Check if it comes from a pi0 decay
    if ((mpdg == kPi0) && (imo < stack->GetNprimary()))   return kTRUE; 

    // Check if this is a heavy flavor decay product
    Int_t mfl  = Int_t (mpdg / TMath::Power(10, Int_t(TMath::Log10(mpdg))));

    // Light hadron
    if (mfl < 4) return kFALSE;

    // Heavy flavor hadron produced by generator
    if (imo <  stack->GetNprimary()) {
      return kTRUE;
    }
	
    // To be sure that heavy flavor has not been produced in a secondary interaction
    // Loop back to the generated mother
    while (imo >=  stack->GetNprimary()) {
      imo = pm->GetFirstMother();
      pm  =  stack->Particle(imo);
    }
    mpdg = TMath::Abs(pm->GetPdgCode());
    mfl  = Int_t (mpdg / TMath::Power(10, Int_t(TMath::Log10(mpdg))));

    if (mfl < 4) {
      return kFALSE;
    } else {
      return kTRUE;
    } 
    } // produced by generator ?
} 


AliAnalysisTaskHMTFMCMultEst::AliAnalysisTaskHMTFMCMultEst()
  : AliAnalysisTaskSE(), fMyOut(0), fEstimatorsList(0), fEstimatorNames(0),
    festimators(0), fRequireINELgt0(kTRUE), fRunconditions(0)
    fdebugk0(0),fdebugkm(0), fdebugkp(0), fdebugnch_gt_30(0), fdebugk0_vs_kch(0),
    fdebugnch_vs_k0(0), fdebugnch_vs_kch(0)
{

}

//________________________________________________________________________
AliAnalysisTaskHMTFMCMultEst::AliAnalysisTaskHMTFMCMultEst(const char *name) 
  : AliAnalysisTaskSE(name), fMyOut(0), fEstimatorsList(0), fEstimatorNames(0),
    festimators(0), fRequireINELgt0(kTRUE), fRunconditions(0)
    fdebugk0(0),fdebugkm(0), fdebugkp(0), fdebugnch_gt_30(0), fdebugk0_vs_kch(0),
    fdebugnch_vs_k0(0), fdebugnch_vs_kch(0)
{
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

void AliAnalysisTaskHMTFMCMultEst::AddEstimator(const char* n)
{
  if (!fEstimatorNames.IsNull()) fEstimatorNames.Append(",");
  fEstimatorNames.Append(n);
}
    
void AliAnalysisTaskHMTFMCMultEst::InitEstimators()
{
  fEstimatorsList = new TList;
  fEstimatorsList->SetOwner();
  fEstimatorsList->SetName("estimators");
  
  TObjArray* arr = fEstimatorNames.Tokenize(",");
  TObject*   obj = 0;
  TIter      next(arr);
  std::cout << "Init estimators... " << std::endl;
  while ((obj = next())) {
    MultiplicityEstimatorBase* e = MakeEstimator(obj->GetName());
    fEstimatorsList->Add(e);
  }
}
//________________________________________________________________________
MultiplicityEstimatorBase*
AliAnalysisTaskHMTFMCMultEst::MakeEstimator(const TString& name)
{
  if (name.BeginsWith("Total")) return new EtaBase("Total", "full #eta coverage ");
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
void AliAnalysisTaskHMTFMCMultEst::UserCreateOutputObjects()
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

  fdebugk0 = new TProfile("debug_nch_k0", "debug nch K0", 100, 0, 100);
  fdebugkp = new TProfile("debug_nch_kplus", "debug nch K+", 100, 0, 100);
  fdebugkm = new TProfile("debug_nch_kminus", "debug nch K-", 100, 0, 100);
  fdebugnch_gt_30 = new TH1F("debug_nch_gt_30", "debug nch > 30", kNPID, 0, kNPID);
  fdebugk0_vs_kch = new TH2F("debug_k0_vs_kch", "fdebugk0_vs_kch", 100, -.5, 99.5, 100, -.5, 99.5);
  fdebugnch_vs_kch = new TH2F("debug_nch_vs_kch", "fdebugnch_vs_kch", 100, -.5, 99.5, 100, -.5, 99.5);
  fdebugnch_vs_k0 = new TH2F("debug_nch_vs_k0", "fdebugnch_vs_k0", 100, -.5, 99.5, 100, -.5, 99.5);
  // fdebugk0->Sumw2(); fdebugk0->SetDirectory(0); fMyOut->Add(fdebugk0);
  // fdebugkp->Sumw2(); fdebugkp->SetDirectory(0); fMyOut->Add(fdebugkp);
  // fdebugkm->Sumw2(); fdebugkm->SetDirectory(0); fMyOut->Add(fdebugkm);
  // fdebugnch_gt_30->Sumw2(); fdebugnch_gt_30->SetDirectory(0); fMyOut->Add(fdebugnch_gt_30);
  // fdebugk0_vs_kch->Sumw2(); fdebugk0_vs_kch->SetDirectory(0); fMyOut->Add(fdebugk0_vs_kch);
  // fdebugnch_vs_kch->Sumw2(); fdebugnch_vs_kch->SetDirectory(0); fMyOut->Add(fdebugnch_vs_kch);
  // fdebugnch_vs_k0->Sumw2(); fdebugnch_vs_k0->SetDirectory(0); fMyOut->Add(fdebugnch_vs_k0);
  PostData(1, fMyOut);
}

//________________________________________________________________________
void AliAnalysisTaskHMTFMCMultEst::UserExec(Option_t *) 
{
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
     Printf("ERROR: Could not retrieve MC event");
     return;
  }

  for (std::vector<MultiplicityEstimatorBase*>::size_type i = 0; i < festimators.size(); i++) {//estimator loop
    festimators[i]->PreEvent(mcEvent);
  }//estimator loop

  // Track loop for establishing multiplicity and checking for INEL > 0
  Bool_t isINEL_gt_0(kFALSE);
  Int_t nch_05 = 0;
  Int_t nch_05_k0 = 0;
  Int_t nch_05_kplus = 0;
  Int_t nch_05_kminus = 0;
  for (Int_t iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {
    AliMCParticle *track = (AliMCParticle*)mcEvent->GetTrack(iTrack);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTrack);
      continue;
    }
    // Pass the particle on to the estimators if it is a primary. Extra check for pi0's is needed since they are unstable
    if (mcEvent->Stack()->IsPhysicalPrimary(iTrack) ||
	IsPi0PhysicalPrimary(iTrack, mcEvent->Stack())){
      if (TMath::Abs(track->Eta()) < 1) isINEL_gt_0 = kTRUE;

      // Debug for K0 ratios -----------------------
      if (TMath::Abs(track->Eta()) < .5){
	if (track->Charge() != 0) nch_05++;
	if (track->PdgCode() == 310) nch_05_k0++;
	if (track->PdgCode() == 321) nch_05_kplus++;
	if (track->PdgCode() == -321) nch_05_kminus++;
      }
      // Debug ends here ----------------------------
      
      for (std::vector<MultiplicityEstimatorBase*>::size_type i = 0; i < festimators.size(); i++) { //estimator loop
	festimators[i]->ProcessTrackForMultiplicityEstimation(track);
      }//estimator loop
    }
  }  //track loop

  // Debug for K0 ratios -----------------------
  fdebugk0->Fill(nch_05, nch_05_k0);
  fdebugkp->Fill(nch_05, nch_05_kplus);
  fdebugkm->Fill(nch_05, nch_05_kminus);
  fdebugk0_vs_kch->Fill(nch_05_k0, (nch_05_kplus+nch_05_kminus));
  fdebugnch_vs_kch->Fill(nch_05, (nch_05_kplus+nch_05_kminus));
  fdebugnch_vs_k0->Fill(nch_05, (nch_05_k0));
  // Debug ends here ----------------------------
  
  // Track loop with known multiplicity in each estimator
  // Skip this if event is not INEL>0 and it is required to be so
  if (!fRequireINELgt0 || isINEL_gt_0){
    for (Int_t iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {
      AliMCParticle *track = (AliMCParticle*)mcEvent->GetTrack(iTrack);
      if (!track) {
	Printf("ERROR: Could not receive track %d", iTrack);
	continue;
      }
      // Debug starts here ----------------------------
      if (nch_05 > 30) {
	fdebugnch_gt_30->Fill(pdg_2_pid_enum(track->PdgCode()), mcEvent->GenEventHeader()->EventWeight());
      }
      // Debug ends here ----------------------------

      if (mcEvent->Stack()->IsPhysicalPrimary(iTrack) ||
	  IsPi0PhysicalPrimary(iTrack, mcEvent->Stack())){

	for (std::vector<MultiplicityEstimatorBase*>::size_type i = 0; i < festimators.size(); i++) { //estimator loop
	  festimators[i]->ProcessTrackWithKnownMultiplicity(track);
	}//estimator loop
      }
    }  //track loop
  
    // Increment eventcounters etc.
    for (std::vector<MultiplicityEstimatorBase*>::size_type i = 0; i < festimators.size(); i++) { //estimator loop
      festimators[i]->PostEvent();
    }//estimator loop
  }
  // Post output data.
  PostData(1, fMyOut);
}      

//________________________________________________________________________
void AliAnalysisTaskHMTFMCMultEst::Terminate(Option_t *) 
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
  std::cout << "Terminating estimators..." << std::endl;
  while ((e = static_cast<MultiplicityEstimatorBase*>(nextEst()))) {
    e->Terminate(fMyOut);
  }
  PostData(1, fMyOut);

  fRunconditions = new TList;
  fRunconditions->SetName("rclist");
  fRunconditions->SetOwner(0);
  TObjString *rcinfo = new TObjString;
  if (fRequireINELgt0) rcinfo->SetString("INELgt0_true");
  else rcinfo->SetString("INELgt0_false");
  fRunconditions->Add(rcinfo);
  PostData(2, fRunconditions);
}

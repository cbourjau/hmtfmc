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
: AliAnalysisTaskSE(), fMyOut(0)
{

}

//________________________________________________________________________
AliAnalysisTaskHMTFMC::AliAnalysisTaskHMTFMC(const char *name) 
  : AliAnalysisTaskSE(name), fMyOut(0)
    
{

  AliPDG::AddParticlesToPdgDataBase();

  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskHMTFMC::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  std::cout << "creating user output object" << std::endl;
  fMyOut = new TList();

  festimators.push_back(new EtaLt05());
  for (std::vector<MultiplicityEstimatorBase*>::size_type i = 0; i < festimators.size(); i++) {
    festimators[i]->RegisterHistograms(fMyOut);
  }

  fMyOut->SetOwner();
  
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

  // Track loop
  Int_t nchMidEta = 0;
  for (Int_t iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {
    AliMCParticle *track = (AliMCParticle*)mcEvent->GetTrack(iTrack);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTrack);
      continue;
    }
    for (std::vector<MultiplicityEstimatorBase*>::size_type i = 0; i < festimators.size(); i++) { //estimator loop
      festimators[i]->ProcessTrack(track, iTrack);
    }//estimator loop
  }  //track loop

  for (std::vector<MultiplicityEstimatorBase*>::size_type i = 0; i < festimators.size(); i++) {//estimator loop
    festimators[i]->PostEvent();
  }//estimator loop


  // fHistNch->Fill(nchMidEta,evWeight);
  // fHistNchUnweighted->Fill(nchMidEta);

  // Post output data.
  PostData(1, fMyOut);
}      

//________________________________________________________________________
void AliAnalysisTaskHMTFMC::Terminate(Option_t *) 
{
  std::cout << "Terminating" << std::endl;
  // Draw result to the screen
  // Called once at the end of the query   
  //  return;
  //  fHistPt = dynamic_cast<TH1F*> (GetOutputData(1));
  //fMyOut  = dynamic_cast<TList*> (GetOutputData(1));
  
  //TH2F* fdNdeta  = (TH2F*) fMyOut->FindObject("fdNdetaMCInel_EtaLt05");
  //std::cout << "integral: " << fdNdeta << std::endl;


  // std::cout << "Processed " << fHistIev->GetBinContent(0) << " events, equivalent to "<< fHistIev->GetBinContent(1) << " generated events."  << std::endl;
  
}

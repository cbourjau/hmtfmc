#include <iostream>

#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"

// #include "AliLog.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"

#include "MultiplicityEstimators.h"

using namespace std;

ClassImp(MultiplicityEstimatorBase)

MultiplicityEstimatorBase::MultiplicityEstimatorBase()
: TNamed()
{
}
MultiplicityEstimatorBase::MultiplicityEstimatorBase(const char* name, const char* title)
: TNamed(name, title), fdNdeta(0)
{
}

void MultiplicityEstimatorBase::RegisterHistograms(TList *outputList){
  std::cout << "Registering Histogram" << std::endl;
  TString name_appendix = TString("_") + fName;
  TString title_appendix = TString(" ") + fTitle;
    fdNdeta = new TH2F("fdNdetaMCInel" + name_appendix,
		       "MC $dN/d\\eta$ Inel," + title_appendix,
		     200, -10.0, 10.0,
		     festimator_bins, 0, 100);    
    fdNdeta->GetXaxis()->SetTitle("#eta");
    fdNdeta->GetYaxis()->SetTitle("$N_{ch}$ in " + title_appendix);
    fdNdeta->GetZaxis()->SetTitle("$N_{ch}$ per $\eta$ bin");
    fdNdeta->SetMarkerStyle(kFullCircle);
    fdNdeta->Sumw2();
    fdNdeta->SetDirectory(0);
    outputList->Add(fdNdeta);

    // fdNdeta_stack = new THStack("dndeta_stack" + name_appendix,
    // 				   "$dN/d\\eta$" + title_appendix);
    // outputList->Add(fdNdeta_stack);
    
    // // Create N_{ch} distribution histograms
    // fHistNch = new TH1F ("fNch" + name_appendix,
    // 			    "Multiplicity distribution" + title_appendix,
    // 			    festimator_bins, 0.0, 100);
    // fHistNch->Sumw2();
    // outputList->Add(fHistNch);

    
    // // Setup counter histograms
    // // Use double in order to not saturate!
    // fEventCounter = new TH2D("fEventCounter" + name_appendix,
    // 				"Number of (weighted) events per cent. bin" + title_appendix,
    // 				2, -0.5, 1.5,  // 1: processed; 2: weighted events
    // 				festimator_bins, 0, 100);
    // fEventCounter->GetXaxis()->SetTitle("processed and weighted");
    // fEventCounter->GetYaxis()->SetTitle("Multiplicity");
    // outputList->Add(fEventCounter);
}


void MultiplicityEstimatorBase::ReadEventHeaders(AliMCEvent *event){
  // Here we need some generator-dependent logic, in case we need to extract useful information from the headers.
  AliGenPythiaEventHeader * headPy  = 0;
  AliGenDPMjetEventHeader * headPho = 0;
  AliGenEventHeader * htmp = event->GenEventHeader();
  if(!htmp) {
    //AliError("Cannot Get MC Header!!");
    return;
  }
  if( TString(htmp->IsA()->GetName()) == "AliGenPythiaEventHeader") {
    headPy =  (AliGenPythiaEventHeader*) htmp;
  } else if (TString(htmp->IsA()->GetName()) == "AliGenDPMjetEventHeader") {
    headPho = (AliGenDPMjetEventHeader*) htmp;
  } else {
    cout << "Unknown header" << endl;
  }
  fevent  = event;
  feventWeight = htmp->EventWeight();
  fheader = event->Header();
  fstack = fheader->Stack();
}

////////////////////////////////////////
// Eta based multiplicity estimators  //
////////////////////////////////////////

//______________________________________________________________________________________
EtaBase::EtaBase() : MultiplicityEstimatorBase() {}

EtaBase::EtaBase(const char* name, const char* title)
  : MultiplicityEstimatorBase(name, title)
{
  festimator_bins = 10;
}

void EtaBase::PreEvent(AliMCEvent *event){
  ReadEventHeaders(event);
  // Room for more header logic
}

void EtaBase::ProcessTrack(AliMCParticle *track, Int_t iTrack){
  //std::cout << track->IsPrimary() << std::endl;
  Bool_t isPrimary = fevent->Stack()->IsPhysicalPrimary(iTrack);
  if (isPrimary && track->Charge() != 0){
    if(TMath::Abs(track->Eta()) < eta_max) nch_current_event++;
  }
  eta_values_current_event.push_back(track->Eta());
}

void EtaBase::PostEvent(){
  for (vector<Float_t>::size_type track = 0; track != eta_values_current_event.size(); track++) {
    fdNdeta->Fill(eta_values_current_event[track], nch_current_event, feventWeight);
  }
  //std::cout << fdNdeta->Integral() << std::endl;
  //std::cout << nch_current_event <<  feventWeight << std::endl;
  //fHistNch->Fill(nch_current_event,feventWeight);
  //fEventCounter->Fill(0.0, nch_current_event);             // Number of processed events
  //fEventCounter->Fill(1, nch_current_event, feventWeight); // Number of weighted events (accounts for downscalling)

  // Clear counters and chaches for the next event:
  eta_values_current_event.clear();
  nch_current_event = 0;
}


EtaLt05::EtaLt05() : EtaBase("EtaLt05", "$\\eta \\lt 0.5") {
  eta_min = 0.0;
  eta_max = 0.5;
}

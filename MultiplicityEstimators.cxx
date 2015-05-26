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
  : TNamed(), fdNdeta(0), fdNdeta_stack(0), fEventCounter(0), fEventCounterUnweighted(0)
   
{
}

MultiplicityEstimatorBase::MultiplicityEstimatorBase(const char* name, const char* title)
  : TNamed(name, title), fdNdeta(0), fdNdeta_stack(0), fEventCounter(0), fEventCounterUnweighted(0)
{
}

void MultiplicityEstimatorBase::RegisterHistograms(TList *outputList){
  std::cout << "Registering Histogram" << std::endl;
  fdNdeta = new TH2F("fdNdetaMCInel" + GetNamePostfix(),
		     "MC $dN/d\\eta$ Inel," + GetTitlePostfix(),
		     200, -10.0, 10.0,
		     festimator_bins, 0, 100);    
  fdNdeta->GetXaxis()->SetTitle("#eta");
  fdNdeta->GetYaxis()->SetTitle("$N_{ch}$ in " + GetTitlePostfix());
  fdNdeta->GetZaxis()->SetTitle("$N_{ch}$ per $\eta$ bin");
  fdNdeta->SetMarkerStyle(kFullCircle);
  fdNdeta->Sumw2();
  fdNdeta->SetDirectory(0);
  outputList->Add(fdNdeta);

  // Histogram stack to make some nice plots in terminate
  fdNdeta_stack = new THStack("dndeta_stack" + GetNamePostfix(),
			      "$dN/d\\eta$" + GetTitlePostfix());
  outputList->Add(fdNdeta_stack);

  // Setup counter histograms
  // Use double in order to not saturate!
  fEventCounter = new TH1D ("fEventCounter" + GetNamePostfix(),
			    "Multiplicity distribution" + GetTitlePostfix(),
			    festimator_bins, 0.0, 100);
  fEventCounter->Sumw2();
  outputList->Add(fEventCounter);
    
  fEventCounterUnweighted = new TH1D ("fEventCounter_unweighted" + GetNamePostfix(),
				      "Multiplicity distribution (unweighted)" + GetTitlePostfix(),
				      festimator_bins, 0.0, 100);
  fEventCounterUnweighted->GetXaxis()->SetTitle("processed and weighted");
  fEventCounterUnweighted->GetYaxis()->SetTitle("Multiplicity");
  outputList->Add(fEventCounterUnweighted);
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
    if(TMath::Abs(track->Eta()) < eta_max && TMath::Abs(track->Eta()) > eta_min) nch_in_estimator_region++;
    eta_values_current_event.push_back(track->Eta());
  }
}

void EtaBase::PostEvent(){
  for (vector<Float_t>::size_type track = 0; track != eta_values_current_event.size(); track++) {
    fdNdeta->Fill(eta_values_current_event[track], nch_in_estimator_region, feventWeight);
  }
  fEventCounter->Fill(nch_in_estimator_region, feventWeight);
  fEventCounterUnweighted->Fill(nch_in_estimator_region);

  // Clear counters and chaches for the next event:
  eta_values_current_event.clear();
  nch_in_estimator_region = 0;
}

void EtaBase::Terminate(TList* outputlist,TList* results){
  // recover pointers to histograms since they might get lost in PROOF
  fEventCounter = static_cast<TH1D*>(outputlist->FindObject("fEventCounter" + GetNamePostfix()));
  fEventCounterUnweighted = static_cast<TH1D*>(outputlist->FindObject("fEventCounter" + GetNamePostfix()));
  fdNdeta = static_cast<TH2F*>(outputlist->FindObject("fdNdetaMCInel" + GetNamePostfix()));
  fdNdeta_stack = static_cast<THStack*>(outputlist->FindObject("dndeta_stack" + GetNamePostfix()));

  // loop over multiplicity bins
  Int_t mult_bin_max = fdNdeta->GetYaxis()->GetNbins();
  for (Int_t mult_bin = 1; mult_bin < mult_bin_max; mult_bin++) {
    Double_t nevents_unweighted = fEventCounterUnweighted->Integral(mult_bin, mult_bin);
    Double_t nevents_weighted = fEventCounter->Integral(mult_bin, mult_bin);
    if (nevents_weighted > 1.0) {
      fdNdeta->GetYaxis()->SetRange(mult_bin, mult_bin);
      // cannot use scale here since that ignores the range. Because ROOT hates me...
      TH1D *temp = (TH1D*) fdNdeta->ProjectionX()->Clone();
      // See histogram def in Begin for which bin is which
      // Divide by number of events in each multiplicity bin
      temp->Scale(1.0 / nevents_weighted);
      std::cout << temp->Integral() << std::endl;
      // used colors between 800 and 900
      temp->SetMarkerColor(800 + (Int_t)((100.0/mult_bin_max)*(mult_bin-1) + 1));
      char title[40];
      sprintf(title, "%.2f $ \\le N_{ch} < $ %.2f",
	      fdNdeta->GetYaxis()->GetBinLowEdge(mult_bin),
	      fdNdeta->GetYaxis()->GetBinUpEdge(mult_bin));
      temp->SetTitle(title + TString("(weighted events"));
      fdNdeta_stack->Add(temp);
    }
  }  
}


EtaLt05::EtaLt05() : EtaBase("EtaLt05", "$\\eta \\lt 0.5$") {
  eta_min = 0.0;
  eta_max = 0.5;
}

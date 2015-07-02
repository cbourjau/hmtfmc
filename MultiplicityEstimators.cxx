#include <iostream>

#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TString.h"

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

Int_t pid_enum_to_pdg(Int_t pid_enum) {
  if (pid_enum == kPROTON) return 2212;
  else if (pid_enum == kLAMBDA) return 3122;
  else if (pid_enum == kK0S) return 310;
  else if (pid_enum == kKPLUS) return 321;
  else if (pid_enum == kKMINUS) return -321;
  else if (pid_enum == kPIPLUS) return 211;
  else if (pid_enum == kPIMINUS) return -211;
  else if (pid_enum == kPI0) return 111;
  else if (pid_enum == kXI) return 3322;
  else if (pid_enum == kOMEGAMINUS) return 3334;
  else if (pid_enum == kOMEGAPLUS) return -3334;
  else return 99999;
}


MultiplicityEstimatorBase::MultiplicityEstimatorBase()
  : TNamed(), fweight_esti(0)
{
}

MultiplicityEstimatorBase::MultiplicityEstimatorBase(const char* name, const char* title)
  : TNamed(name, title), fweight_esti(0)
{
}

void MultiplicityEstimatorBase::RegisterHistograms(TList *outputList){
  std::cout << "Registering Histogram: " << GetName() << std::endl;
  // Put all histograms of one estimator in their own sub-list
  TList *curr_est = new TList();
  curr_est->SetName(GetName());
  outputList->Add(curr_est);


  for (Int_t weighted_or_not = 0; weighted_or_not <= 1; weighted_or_not++) {
    TString postfix("");
    if (weighted_or_not == kUnweighted) {postfix = TString("_unweighted");}

    fdNdeta[weighted_or_not] = new TH2F("fdNdeta" + postfix ,
		       "dN/d\\eta \\ Inel,\\ " + GetTitlePostfix(),
		       200, -10.0, 10.0,
		       festimator_bins, 0, festimator_bins);
    fdNdeta[weighted_or_not]->GetXaxis()->SetTitle("\\eta");
    fdNdeta[weighted_or_not]->GetYaxis()->SetTitle("N_{ch} \\ in \\ " + GetTitlePostfix());
    fdNdeta[weighted_or_not]->GetZaxis()->SetTitle("N_{ch} \\ per\\ \\eta \\ bin");
    fdNdeta[weighted_or_not]->SetMarkerStyle(kFullCircle);
    fdNdeta[weighted_or_not]->Sumw2();
    fdNdeta[weighted_or_not]->SetDirectory(0);
    curr_est->Add(fdNdeta[weighted_or_not]);

    // Setup counter histograms
    // Use double in order to not saturate!
    fEventcounter[weighted_or_not] = new TH1D ("fEventcounter" + postfix,
					       Form("Multiplicity\\ distribution\\ %s", GetTitlePostfix().Data()),
					       festimator_bins, 0.0, festimator_bins);
    fEventcounter[weighted_or_not]->Sumw2();
    fEventcounter[weighted_or_not]->GetXaxis()->SetTitle("Multiplicity\\ in\\ estimator");
    fEventcounter[weighted_or_not]->GetYaxis()->SetTitle("Events");
    curr_est->Add(fEventcounter[weighted_or_not]);

    // 3D: Mult, pT, PID
    festi_pT_pid[weighted_or_not] =
      new TH3F("festi_pT_pid" + postfix,
	       Form("Event\\ class\\ vs.\\ p_{T}\\ vs.\\ pid,\\ %s", GetTitlePostfix().Data()),
	       festimator_bins, 0.0, festimator_bins,
	       20, 0, 20,
	       kNPID, -.5, kNPID - 0.5);
    festi_pT_pid[weighted_or_not]->GetXaxis()->SetTitle("Multiplicity");
    festi_pT_pid[weighted_or_not]->GetYaxis()->SetTitle("p_{T}\\ [GeV]}");
    festi_pT_pid[weighted_or_not]->GetZaxis()->SetTitle("PID");
    // Name bins with pdg code:
    for (Int_t ipid = 0; ipid < kNPID; ipid++) {
      festi_pT_pid[weighted_or_not]->GetZaxis()->SetBinLabel(ipid + 1, Form("%d",pid_enum_to_pdg(ipid)));
    }
    festi_pT_pid[weighted_or_not]->Sumw2();
    festi_pT_pid[weighted_or_not]->SetDirectory(0);			    
    curr_est->Add(festi_pT_pid[weighted_or_not]);
  }

  fNchInEstimatorRegion = new TNtuple("fevent_counter", "Nch in estimator region", "nch:ev_weight");
  curr_est->Add(fNchInEstimatorRegion);

  fweight_esti = new TH2D("fweight_esti", "Distribution of weights in each mult class",
			  10000, 0, 10000,
			  festimator_bins, 0.0, festimator_bins);
  fweight_esti->SetDirectory(0);
  fweight_esti->Sumw2();
  fweight_esti->GetXaxis()->SetTitle("Event weight");
  fweight_esti->GetYaxis()->SetTitle("Multiplicity in estimator region");
  curr_est->Add(fweight_esti);
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

EtaBase::EtaBase(const char* name, const char* title,
		 Float_t feta_min_backwards, Float_t feta_max_backwards, 
		 Float_t feta_min_forwards, Float_t feta_max_forwards)
  : MultiplicityEstimatorBase(name, title),
    feta_min_backwards(feta_min_backwards), feta_max_backwards(feta_max_backwards),
    feta_min_forwards(feta_min_forwards), feta_max_forwards(feta_max_forwards)
{
  festimator_bins = 100;
}

void EtaBase::PreEvent(AliMCEvent *event){
  // Room for more header logic
  ReadEventHeaders(event);

  // Clear counters and chaches for the following event:
  fnch_in_estimator_region = 0;
  memset(fn_pid_in_event, 0, kNPID*sizeof(*fn_pid_in_event));
}

/*
  Count tracks to establish multiplicity in this estimator
*/
void EtaBase::ProcessTrackForMultiplicityEstimation(AliMCParticle *track){
  if (track->Charge() != 0){
    if((track->Eta() >= feta_min_backwards &&
       track->Eta() <= feta_max_backwards) ||
       (track->Eta() >= feta_min_forwards &&
	track->Eta() <= feta_max_forwards))
      {
	fnch_in_estimator_region++;
      }
  }
}

// loop over tracks again, now that mult is known:
void EtaBase::ProcessTrackWithKnownMultiplicity(AliMCParticle *track){
  if (track->Charge() != 0){
    fdNdeta[kUnweighted]->Fill(track->Eta(), fnch_in_estimator_region);
    fdNdeta[kWeighted]->Fill(track->Eta(), fnch_in_estimator_region, feventWeight);

    // y axis are the different particles defined as enum in the header file.
    // ipid is not the histogram bin number! The bins of the histogram are chosen to consume the enum
    // value for the pid.
    Int_t pdgCode = track->PdgCode();
    for (Int_t ipid = 0; ipid < kNPID; ipid++) {
      if (pdgCode == pid_enum_to_pdg(ipid)){
	festi_pT_pid[kUnweighted]->Fill(fnch_in_estimator_region,
					track->Pt(),
					ipid);
	festi_pT_pid[kWeighted]->Fill(fnch_in_estimator_region,
				      track->Pt(),
				      ipid,
				      feventWeight);
	break;
      }
    }
  }
}

void EtaBase::PostEvent(){
  // Fill event counters
  fEventcounter[kUnweighted]->Fill(fnch_in_estimator_region);
  fEventcounter[kWeighted]->Fill(fnch_in_estimator_region, feventWeight);
  fNchInEstimatorRegion->Fill(fnch_in_estimator_region, feventWeight);
  fweight_esti->Fill(feventWeight, fnch_in_estimator_region);
}

void EtaBase::Terminate(TList* outputlist,TList* results){
  // recover pointers to histograms since they are null on master
  std::cout << "Terminate " << fName << std::endl;
  TList *curr_est = static_cast<TList*>(outputlist->FindObject(GetName()));
  fEventcounter[kWeighted] = static_cast<TH1D*>(curr_est->FindObject("fEventcounter" ));
  fEventcounter[kUnweighted] = static_cast<TH1D*>(curr_est->FindObject("fEventcounter_unweighted" ));
  fdNdeta[kWeighted] = static_cast<TH2F*>(curr_est->FindObject("fdNdeta" ));
  fdNdeta[kUnweighted] = static_cast<TH2F*>(curr_est->FindObject("fdNdeta_unweighted" ));
}



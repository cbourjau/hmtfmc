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
  else if (pid_enum == kANTIPROTON) return -2212;
  else if (pid_enum == kLAMBDA) return 3122;
  else if (pid_enum == kANTILAMBDA) return -3122;
  else if (pid_enum == kK0S) return 310;
  else if (pid_enum == kKPLUS) return 321;
  else if (pid_enum == kKMINUS) return -321;
  else if (pid_enum == kPIPLUS) return 211;
  else if (pid_enum == kPIMINUS) return -211;
  else if (pid_enum == kPI0) return 111;
  else if (pid_enum == kXI) return 3322;
  else if (pid_enum == kANTIXI) return -3322;
  else if (pid_enum == kOMEGAMINUS) return 3334;
  else if (pid_enum == kOMEGAPLUS) return -3334;
  else return 99999;
}


MultiplicityEstimatorBase::MultiplicityEstimatorBase()
  : TNamed(), fweight_esti(0), fuseWeights(kTRUE)
{
}

MultiplicityEstimatorBase::MultiplicityEstimatorBase(const char* name, const char* title)
  : TNamed(name, title), fweight_esti(0), fuseWeights(kTRUE)
{
}

void MultiplicityEstimatorBase::RegisterHistograms(TList *outputList){
  std::cout << "Registering Histogram: " << GetName() << std::endl;
  // Put all histograms of one estimator in their own sub-list
  TList *curr_est = new TList();
  curr_est->SetName(GetName());
  outputList->Add(curr_est);

  TString postfix = fuseWeights?"":postfix = fuseWeights?"":"_unweighted";
  fdNdeta = new TH2F("fdNdeta" + postfix ,
		     "dN/d#eta Inel, " + GetTitlePostfix(),
		     200, -10.0, 10.0,
		     festimator_bins, 0, festimator_bins);
  fdNdeta->GetXaxis()->SetTitle("#eta");
  fdNdeta->GetYaxis()->SetTitle("N_{ch} in " + GetTitlePostfix());
  fdNdeta->GetZaxis()->SetTitle("N_{ch} per #eta bin");
  fdNdeta->SetMarkerStyle(kFullCircle);
  fdNdeta->Sumw2();
  fdNdeta->SetDirectory(0);
  curr_est->Add(fdNdeta);

  // Setup counter histograms
  // Use double in order to not saturate!
  fEventcounter = new TH1D ("fEventcounter" + postfix,
			    Form("Multiplicity distribution %s", GetTitlePostfix().Data()),
			    festimator_bins, 0.0, festimator_bins);
  fEventcounter->Sumw2();
  fEventcounter->GetXaxis()->SetTitle("Multiplicity in estimator");
  fEventcounter->GetYaxis()->SetTitle("Events");
  curr_est->Add(fEventcounter);

  // 3D: Mult, pT, PID
  Float_t mult_bin_edges[festimator_bins + 1];
  for (Int_t idx = 0; idx < festimator_bins+1; idx++){
    mult_bin_edges[idx] = idx;
  }
  Float_t pt_bin_edges[] =
    {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,  // 0.1 * 20
     2.2,2.4,2.6,2.8,3.0,                                                                // 0.2 * 5
     3.3,3.6,3.9,4.2,                                                                    // 0.3 * 4
     4.6,5,5.4,                                                                          // 0.4 * 3
     5.9,
     6.5,7,7.5,8,8.5,
     9.2,
     10,11,12,
     13.5,15,
     17,20
    };
    
  Float_t pid_bin_edges[kNPID + 1];
  for (Int_t idx = 0; idx < kNPID +1; idx++) {
    pid_bin_edges[idx] = Float_t(idx) - 0.5;  // shift bins to center around int values
  }

  festi_pT_pid = new TH3F("festi_pT_pid" + postfix,
			  Form("Event class vs. p_{T} vs. pid, %s", GetTitlePostfix().Data()),
			  festimator_bins, mult_bin_edges,
			  sizeof(pt_bin_edges)/sizeof(*pt_bin_edges) - 1, pt_bin_edges,
			  kNPID, pid_bin_edges);
  festi_pT_pid->GetXaxis()->SetTitle("Multiplicity");
  festi_pT_pid->GetYaxis()->SetTitle("p_{T} [GeV]}");
  festi_pT_pid->GetZaxis()->SetTitle("PID");
  // Name bins with pdg code:
  for (Int_t ipid = 0; ipid < kNPID; ipid++) {
    festi_pT_pid->GetZaxis()->SetBinLabel(ipid + 1, Form("%d",pid_enum_to_pdg(ipid)));
  }
  festi_pT_pid->Sumw2();
  festi_pT_pid->SetDirectory(0);			    
  curr_est->Add(festi_pT_pid);

  fNchInEstimatorRegion = new TNtuple("fevent_counter", "Nch in estimator region", "nch:ev_weight:nmpi");
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
    fnMPI = headPy->GetNMPI();
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
    feta_min_forwards(feta_min_forwards), feta_max_forwards(feta_max_forwards),
    fbypass_eta_selection(false)
{
  festimator_bins = 200;
}

// Constructor for bypassing the eta selection to get the full range
EtaBase::EtaBase(const char* name, const char* title)
  : MultiplicityEstimatorBase(name, title),
    feta_min_backwards(0), feta_max_backwards(0),
    feta_min_forwards(0), feta_max_forwards(0),
    fbypass_eta_selection(true)
{
  festimator_bins = 200;
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
    if(fbypass_eta_selection ||
       ((track->Eta() >= feta_min_backwards &&
       track->Eta() <= feta_max_backwards) ||
       (track->Eta() >= feta_min_forwards &&
	track->Eta() <= feta_max_forwards)))
      {
	fnch_in_estimator_region++;
      }
  }
}

// loop over tracks again, now that mult is known:
void EtaBase::ProcessTrackWithKnownMultiplicity(AliMCParticle *track){
  if (track->Charge() != 0){
    // only enforce charged tracks for dN/deta!
    fdNdeta->Fill(track->Eta(), fnch_in_estimator_region, fuseWeights?feventWeight:1);
  }
  // y axis are the different particles defined as enum in the header file.
  // ipid is not the histogram bin number! The bins of the histogram are chosen to consume the enum
  // value for the pid.
  Int_t pdgCode = track->PdgCode();
  for (Int_t ipid = 0; ipid < kNPID; ipid++) {
    // these tracks might be uncharged!
    if (pdgCode == pid_enum_to_pdg(ipid)){
      festi_pT_pid->Fill(fnch_in_estimator_region,
			 track->Pt(),
			 ipid,
			 fuseWeights?feventWeight:1);
      break;
    }
  }
}

void EtaBase::PostEvent(){
  // Fill event counters
  fEventcounter->Fill(fnch_in_estimator_region);
  fEventcounter->Fill(fnch_in_estimator_region, fuseWeights?feventWeight:1);
  fNchInEstimatorRegion->Fill(fnch_in_estimator_region, feventWeight, fnMPI);
  fweight_esti->Fill(feventWeight, fnch_in_estimator_region);
}

void EtaBase::Terminate(TList* outputlist){
  // recover pointers to histograms since they are null on master
  std::cout << "Terminate " << fName << std::endl;
  TList *curr_est = static_cast<TList*>(outputlist->FindObject(GetName()));
  fEventcounter = static_cast<TH1D*>(curr_est->FindObject("fEventcounter" ));
  fdNdeta = static_cast<TH2F*>(curr_est->FindObject("fdNdeta" ));
}



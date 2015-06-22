#ifndef MultiplicityEstimators_cxx
#define MultiplicityEstimators_cxx

//class TH1F;
//class TH1I;
//class TGraphErrors;
class AliMCEvent;
class AliMCParticle;
class AliHeader;
class AliStack;
class AliMCParticle;

#include "TNamed.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"

enum {
  kPROTON,
  kLAMBDA,
  kK0S,
  kKPLUS,
  kKMINUS,
  kPIPLUS,
  kPIMINUS,
  kPI0,
  kXI,
  kOMEGAMINUS,
  kOMEGAPLUS,
  kNPID
};

class MultiplicityEstimatorBase : public TNamed {
 public:
  MultiplicityEstimatorBase();
  MultiplicityEstimatorBase(const char* name, const char* title);
  virtual ~MultiplicityEstimatorBase() {}
  // Available estimators
  enum {kEtaLt05,
	kEtaLt08};
  void RegisterHistograms(TList* outputList);
  //get the ending of the name common to all histograms from this estimator:
  TString GetNamePostfix() {return TString("_") + fName;};
  //get the ending of the title common to all histograms from this estimator:
  TString GetTitlePostfix() {return TString(" ") + fTitle;};
  virtual void PreEvent(AliMCEvent* event) = 0;
  virtual void ProcessTrackForMultiplicityEstimation(AliMCParticle* track) = 0;
  virtual void ProcessTrackWithKnownMultiplicity(AliMCParticle* track) = 0;
  virtual void PostEvent() = 0;
  virtual void Terminate(TList* sum, TList* results) = 0;
  
 protected:
  /*
    Header dependent logic is done here. This function is supposed to be called in the
    beginning of PreEvent().
   */
  void ReadEventHeaders(AliMCEvent* event);
  Int_t festimator_bins;
  // histograms for weighted [0] and unweighted [1] are created where appropriate
  enum {kWeighted,
	kUnweighted};
  TH2F  *fdNdeta[2];          // dNdEta distributions; multiplicity is on the y-axis
  TH3F  *festi_pT_pid[2];  // multiplicity class; pT; pid
  TH1D  *fEventcounter[2];
  TH2D  *fweight_esti;  // Distribution of weights in each multiplicity class
  AliHeader *fheader;       // Event header
  AliMCEvent *fevent;       // current event
  AliStack  *fstack;
  Float_t feventWeight;        // weight of the event read from the header

 private:
  MultiplicityEstimatorBase(const MultiplicityEstimatorBase&);           // not implemented
  MultiplicityEstimatorBase& operator=(const MultiplicityEstimatorBase&);// not implemented


  ClassDef(MultiplicityEstimatorBase, 1); // example of analysis
};


class EtaBase : public MultiplicityEstimatorBase {
 public:
  EtaBase();
  EtaBase(const char* name, const char* title, Float_t eta_min, Float_t eta_max);
 protected:
  void PreEvent(AliMCEvent* event);
  void ProcessTrackForMultiplicityEstimation(AliMCParticle* track);
  void ProcessTrackWithKnownMultiplicity(AliMCParticle *track);
  void PostEvent();
  void Terminate(TList* sum, TList* results);
  Int_t fnch_in_estimator_region;   // counter for charged particles in current event
  Int_t n_pid_in_event[kNPID];     // counter for PID'ed particles in this event
  Float_t eta_min, eta_max;  // range in eta for mult. estimation
  ClassDef(EtaBase, 1)
};

#endif

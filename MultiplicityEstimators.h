#ifndef MultiplicityEstimators_cxx
#define MultiplicityEstimators_cxx

class TH1F;
class TH1I;
class TGraphErrors;

#include "TNamed.h"
#include "THStack.h"
#include "TH2F.h"


class MultiplicityEstimatorBase : public TNamed {
 public:
  MultiplicityEstimatorBase();
  MultiplicityEstimatorBase(const char* name, const char* title);
  //virtual ~MultiplicityEstimatorBase() {}
  // Available estimators
  enum {kEtaLt05,
	kEtaLt08};
  void RegisterHistograms(TList* outputList);
  virtual void PreEvent(AliMCEvent* event) {};
  virtual void ProcessTrack(AliMCParticle* track, Int_t iTrack) {};
  virtual void PostEvent() {};
  //MultiplicityEstimatorBase(const MultiplicityEstimatorBase&);           // not implemented
  //MultiplicityEstimatorBase& operator=(const MultiplicityEstimatorBase&);// not implemented
  
 protected:
  /*
    Header dependent logic is done here. This function is supposed to be called in the
    beginning of PreEvent().
   */
  void ReadEventHeaders(AliMCEvent* event);
  Int_t festimator_bins;
  TH2F  *fdNdeta;          // dNdEta distributions; multiplicity is on the y-axis
  /* THStack *fdNdeta_stack;  // 1D dNdeta histograms scaled to number of events */
  /* TH1F  *fHistNch;         // Multiplicity distribution  */
  /* TH1F  *fHistNchUnweighted;           // Multiplicity distribution  */
  /* TH2D  *fEventCounter;    // Event counter, xaxis: #processed/weighted; yaxis:cent. bins */
  AliHeader *fheader;       // Event header
  AliMCEvent *fevent;       // current event
  AliStack  *fstack;
  Float_t feventWeight;        // weight of the event read from the header

  ClassDef(MultiplicityEstimatorBase, 1); // example of analysis
};


class EtaBase : public MultiplicityEstimatorBase {
 public:
  EtaBase();
  EtaBase(const char* name, const char* title);
 protected:
  void PreEvent(AliMCEvent* event);
  void ProcessTrack(AliMCParticle* track, Int_t itrack);
  void PostEvent();
  Int_t nch_current_event;   // counter for charged particles in current event
  std::vector<Float_t> eta_values_current_event;
  Float_t eta_min, eta_max;  // range in eta for mult. estimation
  ClassDef(EtaBase, 1)
};

class EtaLt05 : public EtaBase {
 public:
  EtaLt05();
  ClassDef(EtaLt05, 1);
};
#endif

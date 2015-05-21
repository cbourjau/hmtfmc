#ifndef AliAnalysisTaskHMTFMC_cxx
#define AliAnalysisTaskHMTFMC_cxx


class TH1F;
class TH1I;
class TGraphErrors;

enum {kHistINEL,kHistNSD,kHistND,kHistSiD,kNHist};
#include "AliAnalysisTaskSE.h"
#include "MultiplicityEstimators.h"
#include <iostream>
#include <stdio.h>


class AliAnalysisTaskHMTFMC : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHMTFMC();
  AliAnalysisTaskHMTFMC(const char *name );
  virtual ~AliAnalysisTaskHMTFMC() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  // AliAnalysisTaskHMTFMC(const AliAnalysisTaskHMTFMC&); // not implemented
  // AliAnalysisTaskHMTFMC& operator=(const AliAnalysisTaskHMTFMC&); // not implemented

 private:
  TList *fMyOut;             // Output list
  MultiplicityEstimatorBase* festi;
  std::vector<MultiplicityEstimatorBase*> festimators;
  
  ClassDef(AliAnalysisTaskHMTFMC, 1); // example of analysis
};

#endif

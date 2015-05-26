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
  virtual ~AliAnalysisTaskHMTFMC() {};

  void AddEstimator(const char* n);
  void InitEstimators();
  MultiplicityEstimatorBase* MakeEstimator(const TString& name);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

 private:
  TList *fMyOut;             // Output list
  TList *fEstimatorsList;   // List to get the estimators out in terminate
  TString fEstimatorNames;
  // MultiplicityEstimatorBase* festi;

  std::vector<MultiplicityEstimatorBase*> festimators;

  // Declaring these shuts up warnings from Weffc++
  AliAnalysisTaskHMTFMC(const AliAnalysisTaskHMTFMC&); // not implemented
  AliAnalysisTaskHMTFMC& operator=(const AliAnalysisTaskHMTFMC&); // not implemented

  ClassDef(AliAnalysisTaskHMTFMC, 1); // example of analysis
};

#endif

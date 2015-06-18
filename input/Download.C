#include "TString.h"

void Download(Bool_t unpack=true)
{
  TString base = "alien:///alice/cern.ch/user/p/pwgpp_mc/2015/17_Week/TestMultiplicity/Test2/Test_1M_events_iter1/";
  TString runs = "";
  Int_t start_file = 301;
  Int_t end_file = 500;
  for (Int_t i=start_file; i <= end_file; i++) {
    if (i<10) runs.Append(Form("0000%d ", i));
    else if (i<100) runs.Append(Form("000%d ", i));
    else if (i<1000) runs.Append(Form("00%d ", i));
    else if (i<10000) runs.Append(Form("0%d ", i));
  }
  std::cout << runs << std::endl;
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/trains/GridDownload.C");

  GridDownload(base, runs, unpack);
}
// EOF


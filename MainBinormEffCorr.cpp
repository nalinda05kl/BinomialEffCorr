#include "BinormEffCorr.h"

#include <iostream>
#include <string>

#include "TFile.h"
#include "TString.h"

using namespace std;

void MainProgram()
{
  cout << "\n\t\t\t\t\t~~~-<><>- Main program starts here -<><>-~~~ \n" << endl;

  //_____________________________________________________
  Int_t *BIN_ENTRIES_REFMULT;
  BIN_ENTRIES_REFMULT = new Int_t[1000];
  Emptying1DArrayInt(BIN_ENTRIES_REFMULT, 1000);

  //_____________________________________________________
  Double_t **CUMULANTS_BINORM_EFF_CORR_REFMULT;
  CUMULANTS_BINORM_EFF_CORR_REFMULT = new Double_t*[1000];
  for(Int_t i = 0; i < 1000; i++) CUMULANTS_BINORM_EFF_CORR_REFMULT[i] = new Double_t[3];
  Emptying2DArray(CUMULANTS_BINORM_EFF_CORR_REFMULT, 1000, 3);

  CumulantCalc_from_BinomEff_correction(BIN_ENTRIES_REFMULT, CUMULANTS_BINORM_EFF_CORR_REFMULT);

  //_____________________________________________________
  Double_t **CBWC_CUMULANTS_BINORM_EFF_CORR;
  CBWC_CUMULANTS_BINORM_EFF_CORR = new Double_t*[9];
  for(Int_t i = 0; i < 9; i++) CBWC_CUMULANTS_BINORM_EFF_CORR[i] = new Double_t[3];
  Emptying2DArray(CBWC_CUMULANTS_BINORM_EFF_CORR, 9, 3);

  CBWC_CumulantCalc_from_BinomEff_correction(CUMULANTS_BINORM_EFF_CORR_REFMULT, CBWC_CUMULANTS_BINORM_EFF_CORR, BIN_ENTRIES_REFMULT, 9, 0, 1000, 1);

  //_____________________________________________________
  DeleteDynamic1DArrayInt(BIN_ENTRIES_REFMULT, "BIN_ENTRIES_REFMULT");
  DeleteDynamic2DArray(CUMULANTS_BINORM_EFF_CORR_REFMULT, 1000, "CUMULANTS_BINORM_EFF_CORR_REFMULT");
  DeleteDynamic2DArray(CBWC_CUMULANTS_BINORM_EFF_CORR, 9, "CBWC_CUMULANTS_BINORM_EFF_CORR");

  cout << "\n\t\t\t\t\t________________________________________" << endl;
  cout << "\t\t\t\t\t-<><><>- Main program ends here -<><><>- \n\n" << endl;
}

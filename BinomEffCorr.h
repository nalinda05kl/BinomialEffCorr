#ifndef BINORMEFFCORR_H_
#define BINORMEFFCORR_H_

#include "TCanvas.h"
#include "TString.h"

void Emptying1DArrayInt(Int_t *, Int_t);
void Emptying1DArray(Double_t *, Int_t);
void Emptying2DArray(Double_t **, Int_t, Int_t);
void DeleteDynamic1DArrayInt(Int_t *, TString);
void DeleteDynamic1DArray(Double_t *, TString);
void DeleteDynamic2DArray(Double_t **, Int_t, TString);

void CumulantCalc_from_BinomEff_correction(Int_t *, Double_t **);
void CBWC_CumulantCalc_from_BinomEff_correction(Double_t **, Double_t **, Int_t *, Int_t, Int_t, Int_t, Int_t);

TString InRootFile = "Merged_file_200GeV_FULL_BinormEffCorr_FD15.root";

TString AddName     = ""; // modify in the systematic error calcualtion. Eg: "_DaDb1". For default, leave this blank (ie. like "").
TString hh_q11      = "h_q11" + AddName;
TString hh_q12      = "h_q12" + AddName;
TString hh_q13      = "h_q13" + AddName;
TString hh_q21      = "h_q21" + AddName;
TString hh_q22      = "h_q22" + AddName;
TString hh_q23      = "h_q23" + AddName;
TString hh_q31      = "h_q31" + AddName;
TString hh_q32      = "h_q32" + AddName;
TString hh_q33      = "h_q33" + AddName;
TString hh_q11_2    = "h_q11_2" + AddName;
TString hh_q11_3    = "h_q11_3" + AddName;
TString hh_q11_q21  = "h_q11_q21" + AddName;
TString hh_q11_q22  = "h_q11_q22" + AddName;

Double_t RefMult3[10]           = {16, 34, 67, 120, 196, 301, 440, 618, 725, 1000};
Double_t RefMult3BinCntre[9]    = {24.0074, 48.8787, 91.3678, 155.454, 245.675, 367.351, 526.062, 670.213, 897.869};

#endif

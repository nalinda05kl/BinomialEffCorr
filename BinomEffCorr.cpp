#include "BinormEffCorr.h"

// headers
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TString.h"

#include "Math/IFunction.h"

// C headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

void Emptying1DArrayInt(Int_t *array, Int_t raws)
{
	for(Int_t i = 0; i < raws; i++) array[i] = 0;
}

void Emptying1DArray(Double_t *array, Int_t raws)
{
	for(Int_t i = 0; i < raws; i++) array[i] = 0;
}

void Emptying2DArray(Double_t **array, Int_t raws, Int_t columns)
{
 	for(Int_t i = 0; i < raws; i++)
 		for(Int_t j = 0; j < columns; j++) array[i][j] = 0.;
}

void DeleteDynamic1DArrayInt(Int_t *Array, TString ArrayName)
{
  delete[] Array; cout << "\n[ The pointer to 1-d array named, " << ArrayName << " is being freed ]" << endl;
}

void DeleteDynamic1DArray(Double_t *Array, TString ArrayName)
{
  delete[] Array; cout << "\n[ The pointer to 1-d array named, " << ArrayName << " is being freed ]" << endl;
}


void DeleteDynamic2DArray(Double_t **Array, Int_t rows, TString ArrayName)
{
 	for(Int_t i = 0; i < rows; i++) delete[] Array[i];
 	delete[] Array;
 	cout << "\n[ The pointer to 2-d array named, " << ArrayName << " is being freed ]" << endl;
}

void CumulantCalc_from_BinomEff_correction(Int_t *BinEntriesRefMult, Double_t **CumulantBinorm)
{

  TFile *LambdaMom_file = new TFile(InRootFile, "read");

  TProfile *Histo_q11     = (TProfile*) LambdaMom_file->Get(hh_q11);
  TProfile *Histo_q12     = (TProfile*) LambdaMom_file->Get(hh_q12);
  TProfile *Histo_q13     = (TProfile*) LambdaMom_file->Get(hh_q13);
  TProfile *Histo_q21     = (TProfile*) LambdaMom_file->Get(hh_q21);
  TProfile *Histo_q22     = (TProfile*) LambdaMom_file->Get(hh_q22);
  TProfile *Histo_q23     = (TProfile*) LambdaMom_file->Get(hh_q23);
  TProfile *Histo_q31     = (TProfile*) LambdaMom_file->Get(hh_q31);
  TProfile *Histo_q32     = (TProfile*) LambdaMom_file->Get(hh_q32);
  TProfile *Histo_q33     = (TProfile*) LambdaMom_file->Get(hh_q33);
  TProfile *Histo_q11_2   = (TProfile*) LambdaMom_file->Get(hh_q11_2);
  TProfile *Histo_q11_3   = (TProfile*) LambdaMom_file->Get(hh_q11_3);
  TProfile *Histo_q11_q21 = (TProfile*) LambdaMom_file->Get(hh_q11_q21);
  TProfile *Histo_q11_q22 = (TProfile*) LambdaMom_file->Get(hh_q11_q22);

  TH1D *RefMultBinEntries = new TH1D("RefMultBinEntries", "RefMultBinEntries", 1000,0,1000);

  Double_t q11 = 0.0, q12 = 0.0, q13 = 0.0, q21 = 0.0, q22 = 0.0, q23 = 0.0, q31 = 0.0, q32 = 0.0, q33 = 0.0, q11_2 = 0.0, q11_3 = 0.0, q11_q21 = 0.0, q11_q22 =0.0;
  Double_t C1 = 0.0, C2 = 0.0, C3 = 0.0;

  Int_t BinEntries = 0;

  for(Int_t i = 0; i < 1000; i++)
  {
  	BinEntriesRefMult[i] = 0;
  	for(Int_t j = 0; j < 3; j++)
  	{
  		CumulantBinorm[i][j] = 0.;
  	}
  }

  //cout << " RefMult3 " << "   \tC1 " << "   \tC2 " << "   \tEntries " << endl;
  for(Int_t RefMultBin = 0; RefMultBin < 1000; RefMultBin++)
  {

  	q11        = Histo_q11 -> GetBinContent(RefMultBin);
  	BinEntries = Histo_q11 -> GetBinEntries(RefMultBin);
  	BinEntriesRefMult[RefMultBin] = BinEntries;

    	RefMultBinEntries   -> SetBinContent(RefMultBin, BinEntries);

  	q12       = Histo_q12   -> GetBinContent(RefMultBin);
    	q13       = Histo_q13   -> GetBinContent(RefMultBin);
    	q21       = Histo_q21   -> GetBinContent(RefMultBin);
    	q22       = Histo_q22   -> GetBinContent(RefMultBin);
    	q23       = Histo_q23   -> GetBinContent(RefMultBin);
    	q31       = Histo_q31   -> GetBinContent(RefMultBin);
    	q32       = Histo_q32   -> GetBinContent(RefMultBin);
    	q33       = Histo_q33   -> GetBinContent(RefMultBin);
    	q11_2     = Histo_q11_2 -> GetBinContent(RefMultBin);
    	q11_3     = Histo_q11_3 -> GetBinContent(RefMultBin);
    	q11_q21   = Histo_q11_q21 -> GetBinContent(RefMultBin);
    	q11_q22   = Histo_q11_q22 -> GetBinContent(RefMultBin);

    	// Calculation of cumulants.
    	C1 = q11;
    	C2 = q11_2 + q21 - q22;
    	C3 = q11_3 + 3*q11_q21 - 3*q11_q22 + q31 - 3*q32 + 2*q33;

    	//cout << " " << RefMultBin << " \t" << C1 << "     \t" << C2 << "     \t" << BinEntries << endl;

    	CumulantBinorm[RefMultBin][0] = C1;
    	CumulantBinorm[RefMultBin][1] = C2;
    	CumulantBinorm[RefMultBin][2] = C3;
  }
  LambdaMom_file -> Close();

}

void CBWC_CumulantCalc_from_BinomEff_correction(Double_t **mu, Double_t **cbwc_mu, Int_t *events, Int_t CentBins, Int_t RefMultLow, Int_t RefMultHigh, Int_t Type)
{
	Int_t RefMultBinWidth = 0;
	if(Type == 0)
	{
		RefMultBinWidth = (Int_t) ((RefMultHigh - RefMultLow)/CentBins);
		cout << CentBins << " " << RefMultLow << " " << RefMultHigh << " " << RefMultBinWidth << endl;
	}

	Int_t binLow  = 0;
	Int_t binHigh = 0;

	if(Type == 0){binLow  = RefMultLow; binHigh = RefMultLow + RefMultBinWidth;}
	else{binLow  = RefMult3[0]; binHigh = RefMult3[1];}

  	Double_t Sum_MU_1 = 0., Sum_MU_2 = 0., Sum_MU_3 = 0.;
  	Double_t Cum_1 = 0., Cum_2 = 0., Cum_3 = 0., Cum_21 = 0., Cum_32 = 0.;
	Int_t eventsSum = 0;

	for(Int_t i = 0; i < CentBins; i++)
	{
		cout << "\n\t______ [" << binLow << " , "<< binHigh << "] ______" << endl;

		for(Int_t j = binLow; j < binHigh; j++)
		{
			Sum_MU_1 += mu[j][0]*events[j];
			Sum_MU_2 += mu[j][1]*events[j];
      			Sum_MU_3 += mu[j][2]*events[j];
			eventsSum += events[j];
		}
		cbwc_mu[i][0] =  Sum_MU_1/eventsSum;
		cbwc_mu[i][1] =  Sum_MU_2/eventsSum;
    		cbwc_mu[i][2] =  Sum_MU_3/eventsSum;

    		Sum_MU_1 = 0.; Sum_MU_2 = 0.; Sum_MU_3 = 0.;
    		Cum_1 = 0.; Cum_2 = 0.; Cum_3 = 0.; Cum_21 = 0.; Cum_32 = 0.;
		eventsSum = 0;
    		Cum_1   = cbwc_mu[i][0];
    		Cum_2   = (cbwc_mu[i][1] - pow(cbwc_mu[i][0],2));
    		Cum_3   = (2*pow(cbwc_mu[i][0],3) - 3*cbwc_mu[i][0]*cbwc_mu[i][1] + cbwc_mu[i][2]);
    		Cum_21  = Cum_2/Cum_1;
    		Cum_32  = Cum_3/Cum_2;

		cout << "\n\tC_1 \t\t= " << Cum_1 << endl;
		cout << "\tC_2 \t\t= " << Cum_2 << endl;
    		cout << "\tC_3 \t\t= " << Cum_3 << endl;
    		cout << "\tC_2/C_1 \t= " << Cum_21 << endl;
    		cout << "\tC_3/C_2 \t= " << Cum_32 << endl;


		if(Type == 0){binLow  += RefMultBinWidth; binHigh += RefMultBinWidth;}
		else{binLow  = RefMult3[i+1]; binHigh = RefMult3[i+2];}
	}
  	cout << "\n    ~~~~~~~~~~~~ CBWC END ~~~~~~~~~~~~" << endl;
}

string base="/unix/dune/tof/";
string date="2018Jul16";
string BARS[]={ "D1714"};
string POSITIONS[]={"16cm", "32cm", "48cm", "64cm", "80cm", "96cm", "112cm"};
string THRESHOLDS[]={"10mV", "15mV", "20mV"};
int colors []      ={kBlue+1, kCyan+1, kGreen+1, kOrange+1, kRed, kViolet};
string CHANNELS[]={"PMTA", "PMTB"};
string divisions="20mVdiv";

#include "functions.h"

void buildPositionEstimator(){


  gStyle->SetOptStat(0);

  TH2D *hEstimator = new TH2D ("hEstimator", "", 100, -2e-8, 2e-8, 100, -1, 1);
  TH2D *hDenominator = new TH2D ("hDenominator", "", 100, -2e-8, 2e-8, 100, -1, 1);
  TH2D *hEstimator2 = new TH2D ("hEstimator2", "", 100, -2e-8, 2e-8, 100, -1, 1);
  TH2D *hDenominator2 = new TH2D ("hDenominator2", "", 100, -2e-8, 2e-8, 100, -1, 1);

  int npos = sizeof(POSITIONS)/sizeof(POSITIONS[0]);
  int nthresh = sizeof(THRESHOLDS)/sizeof(THRESHOLDS[0]);
  int nch = sizeof(CHANNELS)/sizeof(CHANNELS[0]);
  int nbars = sizeof(BARS)/sizeof(BARS[0]);

  for (int ibar=0; ibar<nbars; ibar++){
  for (int ich=0; ich<2; ich++){

    for (int ithresh=0; ithresh<nthresh; ithresh++){
   
      for (int ipos=0; ipos<npos; ipos++){
      
	string barname = BARS[ibar];

	string basename   =   base + date + "/Bar" + barname + "/SourceAt" + POSITIONS[ipos] + "/Trig" + CHANNELS[ich] + "_thres" + THRESHOLDS[ithresh] + ".root";
      
	TFile *fin = new TFile(basename.c_str(), "read");

	TH2D *htemp = (TH2D*)fin->Get("hChargeAsymmetryVSdeltaT");
	TH2D *htemp2 = (TH2D*)fin->Get("hChargeAsymmetryVSdeltaT2");

	double pos = 16.*(ipos+1);

	hEstimator->Add(htemp, pos);
	hDenominator->Add(htemp, +1);

	hEstimator2->Add(htemp2, pos);
	hDenominator2->Add(htemp2, +1);
      


	fin->Close();

      }
      }   
    }
  }

  TCanvas *c = new TCanvas("c", "", 1200, 800);
  c->Divide(2,2);
  c->cd(1);
  hEstimator->Divide(hDenominator);
  hEstimator->Smooth();
  hEstimator->Draw("colz");
  c->cd(3);
  TF2 *f1 = new TF2("f1", "[0]+[1]*( x[0] + x[1]*[2] )", -1e-10, 1e-10, -1, 1);
  f1->FixParameter(0, 64.);
  f1->SetParLimits(1, 5e9, 6e9);
  // f1->SetParLimits(2, TMath::Pi()/4., TMath::Pi()*2.);
  hEstimator->Fit(f1, "", "r", -1e10, 1e10);

  f1->Draw("surf2");
  c->cd(4);

  hEstimator->GetXaxis()->SetRangeUser(-1e10, +1e10);
  hEstimator->Draw("surf2");


  cout << f1->Eval(0,0) << " " << f1->Eval(6.5e-9, -0.55) << " " << f1->Eval(-6.5e-9, +0.55) << endl;


  TCanvas *c2 = new TCanvas("c2", "", 1200, 500);
  c2->Divide(2,1);
  c2->cd(1);
  hEstimator2->SetTitle(";#Delta t [s];Charge asymmetry");
  hEstimator2->Divide(hDenominator2);
  hEstimator2->Smooth();
  hEstimator2->Draw("colz");

  TF2 *f2 = new TF2("f2", "[0]+[1]*( x[0] + x[1]*[2] )", -1e-10, 1e-10, -1, 1);
  f2->FixParameter(0, 64.);
  f2->SetParLimits(1, 5e9, 6e9);
  // f2->SetParLimits(2, TMath::Pi()/4., TMath::Pi()*2.);
  hEstimator2->Fit(f2, "", "r", -1e10, 1e10);



  //  f2->Draw("surf2");
  c2->cd(2);
  hEstimator2->GetXaxis()->SetRangeUser(-1e10, +1e10);
  hEstimator2->Draw("surf2");

  cout << f2->Eval(0,0) << " " << f2->Eval(6.5e-9, -0.55) << " " << f2->Eval(-6.5e-9, +0.55) << endl;

}

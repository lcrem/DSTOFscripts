#include <sys/stat.h>
#include <time.h>
#include <stdio.h>
#include "FFTtools.h"
#include "functions.h"

string trigChan[2] = {"Ch3", "Ch4"};
string chanNice[2] = {"PMTA", "PMTB"};
string thresholds[3] = {"10mV", "15mV", "20mV"};

void findCalibrationNumbers(string base, string date, string barname, string div){

  string sourcePos = "64cm";
  

  string whichPMT[2];
  double gains[2];

  int foundPMT = findWhichPMT(barname, whichPMT);
  if (foundPMT==0) return;
  int foundGains = findGains(whichPMT, gains);
  if (foundGains==0) return;

  //  example file path /unix/dune/tof/2018Jul06/BarD1719/SourceAt64cm/TrigCh3_thres50mV_50mVdiv.ch3.traces.root
  // string base      = "/unix/dune/tof/";
  // string date      = "2018Jul06/";
  // string barname   = "D1719";
  // string sourcePos = "64cm";
  // string trigChan  = "Ch3";
  // string thresh    = "50mV";
  // string div       = "50mVdiv";

  TH1D *hChargeAsymmetry = new TH1D ("hChargeAsymmetry", "", 100, -1, 1);
  TH2D *hChargeAsymmetryVSdeltaT = new TH2D ("hChargeAsymmetryVSdeltaT", "", 100, -2e-8, 2e-8, 100, -1, 1);
  TH1D *hDeltaT = new TH1D ("hDeltaT", "", 50, -2e-8, 2e-8);
  TH1D *hDeltaTsimple = new TH1D("hDeltaTsimple", "", 50, -2e-8, 2e-8);

  for (int ichan=0; ichan<2; ichan++){

    for (int ithresh=0; ithresh<3; ithresh++){


      string inputname1 = base + date + "/Bar" + barname + "/SourceAt" + sourcePos + "/Trig" + trigChan[ichan] + "_thres" + thresholds[ithresh] + "_" + div + ".ch3.traces";
      string inputname2 = base + date + "/Bar" + barname + "/SourceAt" + sourcePos + "/Trig" + trigChan[ichan] + "_thres" + thresholds[ithresh] + "_" + div + ".ch4.traces";

      cout << inputname1 << endl;
      cout << inputname2 << endl;



      int ngraph = 10000;
      double R = 50; // Ohm
      double electron = 1.60217661E-19; // Coulomb
      double numPhotoEle1, numPhotoEle2;
      double min1, min2;
      double mean1, mean2;
      double integral1, integral2;
      double integralNoise1, integralNoise2;
      double deltat, deltatold;
      double maxTime1, maxTime2;
      double chargeAsymmetry;

      TFile *f1signal  = new TFile ((inputname1+".root").c_str(), "read");
      TFile *f2signal  = new TFile ((inputname2+".root").c_str(), "read");

      for (int igraph=1; igraph<=ngraph; igraph++){
  
	TGraph *g1 = (TGraph*)f1signal->Get(Form("graph%d", igraph));
	if (!g1)  continue;
	zeroGraph(g1);

	min1 = TMath::MinElement(g1->GetN(), g1->GetY());

	TGraph *g2 = (TGraph*)f2signal->Get(Form("graph%d", igraph));
	if (!g2)  continue;
	zeroGraph(g2);

	min2 = TMath::MinElement(g2->GetN(), g2->GetY());

	integral1 = getIntegralFromHisto(g1);
	integral2 = getIntegralFromHisto(g2);

	numPhotoEle1 = integral1/(gains[0]*R*electron);
	numPhotoEle2 = integral2/(gains[1]*R*electron);
     
  
	maxTime1 = g1->GetX()[TMath::LocMin(g1->GetN(), g1->GetY())] ;
	maxTime2 = g2->GetX()[TMath::LocMin(g2->GetN(), g2->GetY())] ;
	deltatold = maxTime1 - maxTime2;
 
     
	TGraph *gCor = FFTtools::getCorrelationGraph(g1, g2);
	double *x = gCor->GetX();
	double *y = gCor->GetY();
	int ntot = gCor->GetN();
	deltat = x[TMath::LocMax(ntot, y)];
     
	//     deltat = x[TMath::LocMax(ntot, y)];
	//     cout << igraph << " " << deltat << " " << x[TMath::LocMax(ntot, y)] << " " << deltat - x[TMath::LocMax(ntot, y)] << endl;
	hDeltaT->Fill(deltat);

	chargeAsymmetry = (numPhotoEle1-numPhotoEle2)/(numPhotoEle1+numPhotoEle2);
	hChargeAsymmetry->Fill(chargeAsymmetry);
	hChargeAsymmetryVSdeltaT->Fill(deltat, chargeAsymmetry);
	hDeltaTsimple->Fill(deltatold);

	delete g1;
	delete g2;
	delete gCor;
      }
      
      f1signal->Close();
      f2signal->Close();

    }
  }

  string titleStr = barname + " ^{90}Sr at " + sourcePos; 

  hDeltaT->SetTitle(Form("%s;#Delta t [s];Entries", titleStr.c_str()));
  hDeltaT->SetLineWidth(2);
  
  gStyle->SetOptFit(1);

  TCanvas *c2 = new TCanvas("c2", "", 1200, 500);
  c2->Divide(2,2);
  c2->cd(1);
  hDeltaT->Draw("e");
  TF1 *func = new TF1("func", "gaus");
  hDeltaT->Fit(func);

  hDeltaTsimple->SetLineColor(kRed);
  hDeltaTsimple->Draw("e same");

  c2->cd(2);
  hChargeAsymmetry->SetTitle(Form("%s;Charge Asymmetry;Entries", titleStr.c_str()));
  hChargeAsymmetry->SetLineWidth(2);
  hChargeAsymmetry->Draw("e");
  TF1 *func2 = new TF1("func2", "gaus");
  hChargeAsymmetry->Fit(func2);

  c2->cd(3);

  hChargeAsymmetryVSdeltaT->Draw("colz");

  // delete hDeltaT;
  // delete hDeltaTsimple;
  // delete hChargeAsymmetry;
  // delete hChargeAsymmetryVSdeltaT;

}


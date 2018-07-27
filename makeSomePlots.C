#include <sys/stat.h>
#include <time.h>
#include <stdio.h>
#include "FFTtools.h"
#include "functions.h"


string trigChan[2] = {"Ch3", "Ch4"};
//string trigChan[2] = {"Ch2", "Ch3"};
string chanNice[2] = {"PMTA", "PMTB"};

void makeSomePlots(string base, string date, string barname, string sourcePos, int ichan, string thresh, string div){

  TF2 *fPosEstimator = new TF2("fPosEstimator", "[0]+[1]*( x[0] + x[1]*[2] )", -1e-10, 1e-10, -1, 1);
  fPosEstimator->FixParameter(0, 64.);
  fPosEstimator->FixParameter(1, 5e9);
  fPosEstimator->FixParameter(2, -4.3e-9);

  TF2 *fPosEstimator2 = new TF2("fPosEstimator2", "[0]+[1]*( x[0] + x[1]*[2] )", -1e-10, 1e-10, -1, 1);
  fPosEstimator2->FixParameter(0, 64.);
  fPosEstimator2->FixParameter(1, 5e9);
  fPosEstimator2->FixParameter(2, -5.21e-9);

  string whichPMT[2];
  double gains[2];
  double calvar[2];
  int foundPMT = findWhichPMT(barname, whichPMT);
  if (foundPMT==0) return;
  int foundCal = findCalibration(barname, calvar);
  if (foundCal==0) return;
  int foundGains = findGains(whichPMT, gains);
  if (foundGains==0) return;


  double calDeltaT = calvar[0];
  double calAsymmetry = calvar[1];

  //  example file path /unix/dune/tof/2018Jul06/BarD1719/SourceAt64cm/TrigCh3_thres50mV_50mVdiv.ch3.traces.root
  // string base      = "/unix/dune/tof/";
  // string date      = "2018Jul06/";
  // string barname   = "D1719";
  // string sourcePos = "64cm";
  // string trigChan  = "Ch3";
  // string thresh    = "50mV";
  // string div       = "50mVdiv";

  string inputname1 = base + date + "/Bar" + barname + "/SourceAt" + sourcePos + "/Trig" + trigChan[ichan] + "_thres" + thresh + "_" + div + ".ch3.traces";
  string inputname2 = base + date + "/Bar" + barname + "/SourceAt" + sourcePos + "/Trig" + trigChan[ichan] + "_thres" + thresh + "_" + div + ".ch4.traces";

  string inputnoise1 = base + date + "/Bar" + barname + "/NoSource/Trig" + trigChan[ichan] + "_thres" + thresh + "_" + div + ".ch3.traces";
  string inputnoise2 = base + date + "/Bar" + barname + "/NoSource/Trig" + trigChan[ichan] + "_thres" + thresh + "_" + div + ".ch4.traces";

  string basename   =   base + date + "/Bar" + barname + "/SourceAt" + sourcePos + "/Trig" + chanNice[ichan] + "_thres" + thresh ;
  string outputname =   basename + ".root";

  string titleStr = barname + " ^{90}Sr at " + sourcePos + ", trig " + chanNice[ichan] + " " + thresh + " thresh"; 

  cout << inputname1 << endl;
  cout << inputname2 << endl;
  cout << inputnoise1 << endl;
  cout << inputnoise2 << endl;
  cout << outputname << endl;

  TH1D *hPeak1 = new TH1D ("hPeak1", "", 50, 0, 0.2);
  TH1D *hPeak2 = new TH1D ("hPeak2", "", 50, 0, 0.2);
  TH1D *hIntegral1 = new TH1D ("hIntegral1", "", 100, 0, 500);
  TH1D *hIntegral2 = new TH1D ("hIntegral2", "", 100, 0, 500);
  TH1D *hDifference1 = new TH1D ("hDifference1", "", 100, 0, 500);
  TH1D *hDifference2 = new TH1D ("hDifference2", "", 100, 0, 500);
  TH1D *hPeakNoise1 = new TH1D ("hPeakNoise1", "", 50, 0, 0.2);
  TH1D *hPeakNoise2 = new TH1D ("hPeakNoise2", "", 50, 0, 0.2);
  TH1D *hIntegralNoise1 = new TH1D ("hIntegralNoise1", "", 100, 0, 500);
  TH1D *hIntegralNoise2 = new TH1D ("hIntegralNoise2", "", 100, 0, 500);
  TH1D *hChargeAsymmetry = new TH1D ("hChargeAsymmetry", "", 100, -1, 1);
  TH2D *hChargeAsymmetryVSdeltaT = new TH2D ("hChargeAsymmetryVSdeltaT", "", 100, -2e-8, 2e-8, 100, -1, 1);
  TH1D *hChargeAsymmetry2 = new TH1D ("hChargeAsymmetry2", "", 100, -1, 1);
  TH2D *hChargeAsymmetryVSdeltaT2 = new TH2D ("hChargeAsymmetryVSdeltaT2", "", 100, -2e-8, 2e-8, 100, -1, 1);
  TH1D *hDeltaT = new TH1D ("hDeltaT", "", 50, -2e-8, 2e-8);
  TH1D *hDeltaTsimple = new TH1D("hDeltaTsimple", "", 50, -2e-8, 2e-8);
  TH1D *hPosResFromDeltaT = new TH1D("hPosResFromDeltaT", "", 100, -40, 40);
  TH1D *hPosResFromAsymm  = new TH1D("hPosResFromAsymm", "", 100, -40, 40);
  TH1D *hPosResFromAsymm2  = new TH1D("hPosResFromAsymm2", "", 100, -40, 40);
  TH1D *hPosRes = new TH1D("hPosRes", "", 100, -40, 40);
  TH1D *hPosRes2 = new TH1D("hPosRes2", "", 100, -40, 40);


  TH2D *hPeakAvsPeakB = new TH2D("hPeakAvsPeakB", "", 50, 0, 500, 50, 0, 500);
  TH2D *hPeakAvsTimeA = new TH2D("hPeakAvsTimeA", "", 50, 0, 500, 50, -20e-9, 20e-9);
  TH2D *hPeakBvsTimeB = new TH2D("hPeakBvsTimeB", "", 50, 0, 500, 50, -20e-9, 20e-9);
  TH2D *hPeakAvsDeltaT = new TH2D("hPeakAvsDeltaT", "", 50, 0, 500, 50, -20e-9, 20e-9);
  TH2D *hPeakBvsDeltaT = new TH2D("hPeakBvsDeltaT", "", 50, 0, 500, 50, -20e-9, 20e-9);

  TH2D *hExpAvsDeltaT = new TH2D("hExpAvsDeltaT", "", 50, -100, 100, 50, -20e-9, 20e-9);
  TH2D *hExpBvsDeltaT = new TH2D("hExpBvsDeltaT", "", 50, -100, 100, 50, -20e-9, 20e-9);


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
  double chargeAsymmetry, chargeAsymmetry2;

  double posFromDeltaT, posFromAsymm, posFromAsymm2, posEstim, posEstim2;

  double expectedPos =  stod(sourcePos.substr (0, sourcePos.size()-2));

  TFile *f1signal  = new TFile ((inputname1+".root").c_str(), "read");
  TFile *f2signal  = new TFile ((inputname2+".root").c_str(), "read");
  TFile *f1noise = new TFile ((inputnoise1+".root").c_str(), "read");
  TFile *f2noise = new TFile ((inputnoise2+".root").c_str(), "read");

   for (int igraph=1; igraph<=ngraph; igraph++){
  
     TGraph *g1 = (TGraph*)f1signal->Get(Form("graph%d", igraph));
     if (!g1)  continue;
     zeroGraph(g1);

     min1 = TMath::MinElement(g1->GetN(), g1->GetY());

     TGraph *g2 = (TGraph*)f2signal->Get(Form("graph%d", igraph));
     if (!g2)  continue;
     zeroGraph(g2);

     min2 = TMath::MinElement(g2->GetN(), g2->GetY());

     hPeak1->Fill(TMath::Abs(min1));
     hPeak2->Fill(TMath::Abs(min2));

     integral1 = getIntegralFromHisto(g1);
     integral2 = getIntegralFromHisto(g2);

     numPhotoEle1 = integral1/(gains[0]*R*electron);
     numPhotoEle2 = integral2/(gains[1]*R*electron);
     
  
     maxTime1 = g1->GetX()[TMath::LocMin(g1->GetN(), g1->GetY())] ;
     maxTime2 = g2->GetX()[TMath::LocMin(g2->GetN(), g2->GetY())] ;
     deltatold = maxTime1 - maxTime2 - calDeltaT;

  
     
     TGraph *gCor = FFTtools::getCorrelationGraph(g1, g2);
     double *x = gCor->GetX();
     double *y = gCor->GetY();
     int ntot = gCor->GetN();
     deltat = x[TMath::LocMax(ntot, y)] - calDeltaT;
     chargeAsymmetry = (numPhotoEle1-numPhotoEle2)/(numPhotoEle1+numPhotoEle2) - calAsymmetry;
     chargeAsymmetry2 = (integral1-integral2)/(integral1+integral2);
     
     posFromDeltaT = getPositionFromDeltaT(deltat);
     posFromAsymm = getPositionFromAsymmetry(chargeAsymmetry);
     posFromAsymm2 = getPositionFromAsymmetry(chargeAsymmetry2);
     posEstim = fPosEstimator->Eval(deltat, chargeAsymmetry);
     posEstim2 = fPosEstimator2->Eval(deltat, chargeAsymmetry2);

     if (qualityCut(ichan, sourcePos, numPhotoEle1, numPhotoEle2, deltat)) continue;

     //     deltat = x[TMath::LocMax(ntot, y)];
     //     cout << igraph << " " << deltat << " " << x[TMath::LocMax(ntot, y)] << " " << deltat - x[TMath::LocMax(ntot, y)] << endl;
     hDeltaT->Fill(deltat);
  

     hChargeAsymmetry->Fill(chargeAsymmetry);
     hChargeAsymmetryVSdeltaT->Fill(deltat, chargeAsymmetry);
     hChargeAsymmetry2->Fill(chargeAsymmetry2);
     hChargeAsymmetryVSdeltaT2->Fill(deltat, chargeAsymmetry2);

     hPosResFromDeltaT->Fill(posFromDeltaT-expectedPos);
     hPosResFromAsymm->Fill(posFromAsymm-expectedPos);
     hPosResFromAsymm2->Fill(posFromAsymm2-expectedPos);
     hPosRes->Fill(posEstim - expectedPos);
     hPosRes2->Fill(posEstim2 - expectedPos);
     //     cout << expectedPos << " " << posFromAsymm << " " << posFromDeltaT << " " << pos << endl;


     hIntegral1->Fill(numPhotoEle1);
     hIntegral2->Fill(numPhotoEle2);
     hDifference1->Fill(numPhotoEle1);
     hDifference2->Fill(numPhotoEle2);

     hDeltaTsimple->Fill(deltatold);

     hPeakAvsPeakB->Fill(numPhotoEle1, numPhotoEle2);
     hPeakAvsTimeA->Fill(numPhotoEle1, maxTime1);
     hPeakBvsTimeB->Fill(numPhotoEle2, maxTime2);

     hPeakAvsDeltaT->Fill(numPhotoEle1, deltat);
     hPeakBvsDeltaT->Fill(numPhotoEle2, deltat);

     hExpAvsDeltaT->Fill(numPhotoEle1/maxTime1/1e9, deltat);
     hExpBvsDeltaT->Fill(numPhotoEle2/maxTime2/1e9, deltat);
     

     // if (deltat>5e-9){
     // TCanvas *c = new TCanvas("c");
     // c->Divide(2,2);
     // c->cd(1);
     // g1->Draw("Al");
     // c->cd(2);
     // g2->Draw("Al");
     // c->cd(3);
     // gCor->Draw("Al");
     // }

     delete g1;
     delete g2;
     delete gCor;
   }


   // for (int igraph=1; igraph<=ngraph; igraph++){

   //  TGraph *g1n = (TGraph*)f1noise->Get(Form("graph%d", igraph));
   //  if (!g1n){
   //    continue;
   //  }
   //  min1 = TMath::MinElement(g1n->GetN(), g1n->GetY());

   //  TGraph *g2n = (TGraph*)f2noise->Get(Form("graph%d", igraph));
   //  if (!g2n){
   //    continue;
   //  }
   //  min2 = TMath::MinElement(g2n->GetN(), g2n->GetY());

   //  integralNoise1 = getIntegralFromHisto(g1n);
   //  integralNoise2 = getIntegralFromHisto(g2n);
   //  numPhotoEle1 = integralNoise1/(gains[0]*R*electron);
   //  numPhotoEle2 = integralNoise2/(gains[1]*R*electron);
    
   //  TGraph *gCor = FFTtools::getCorrelationGraph(g1n, g2n);
   //  double *x = gCor->GetX();
   //  double *y = gCor->GetY();
   //  int ntot = gCor->GetN();
   //  deltat = x[TMath::LocMax(ntot, y)];
    
   //  if (qualityCut(ichan, sourcePos, numPhotoEle1, numPhotoEle2, deltat)) continue;
    
    
   //  hPeakNoise1->Fill(TMath::Abs(min1));
   //  hPeakNoise2->Fill(TMath::Abs(min2));

   //  hIntegralNoise1->Fill(numPhotoEle1);
   //  hIntegralNoise2->Fill(numPhotoEle2);

   //  delete g1n;
   //  delete g2n;
   //  delete gCor;

   // }

  f1signal->Close();
  f2signal->Close();
  f1noise->Close();
  f2noise->Close();


  double deltaTsource = getDeltaT(inputname1);

  double deltaTnoise = getDeltaT(inputnoise1);
  
  cout << "Delta T source and noise : " << deltaTsource << " " << deltaTnoise << endl;

  double scale = 0.; // deltaTsource/deltaTnoise;

  hPeakNoise1->Scale(scale);
  hPeakNoise2->Scale(scale);
  hIntegralNoise1->Scale(scale);
  hIntegralNoise2->Scale(scale);
 
  hDifference1->Add(hIntegralNoise1, -1);
  hDifference2->Add(hIntegralNoise2, -1);


  hIntegral1->SetTitle(Form("Channel 2: %s;p.e.;Entries", titleStr.c_str()));
  hIntegral2->SetTitle(Form("Channel 3: %s;p.e.;Entries", titleStr.c_str()));

  hIntegral1->SetLineWidth(2);
  hIntegral1->SetFillStyle(3003);
  hIntegral1->SetFillColor(kBlue);
  hIntegral1->SetLineColor(kBlue);
  hIntegralNoise1->SetLineWidth(2);
  hIntegralNoise1->SetFillStyle(3003);
  hIntegralNoise1->SetFillColor(kRed);
  hIntegralNoise1->SetLineColor(kRed);
  hIntegral2->SetLineWidth(2);
  hIntegral2->SetFillStyle(3003);
  hIntegral2->SetFillColor(kBlue);
  hIntegral2->SetLineColor(kBlue);
  hIntegralNoise2->SetLineWidth(2);
  hIntegralNoise2->SetFillStyle(3003);
  hIntegralNoise2->SetFillColor(kRed);
  hIntegralNoise2->SetLineColor(kRed);

  hDifference1->SetLineWidth(2);
  hDifference1->SetFillStyle(3003);
  hDifference1->SetFillColor(kBlue);
  hDifference1->SetLineColor(kBlue);
  hDifference2->SetLineWidth(2);
  hDifference2->SetFillStyle(3003);
  hDifference2->SetFillColor(kBlue);
  hDifference2->SetLineColor(kBlue);


  hDeltaT->SetTitle(Form("%s;#Delta t [s];Entries", titleStr.c_str()));
  hDeltaT->SetLineWidth(2);
  
  gStyle->SetOptFit(1);

  TFile *fout = new TFile(outputname.c_str(), "recreate");
  hPeak1->Write("hPeakA");
  hPeak2->Write("hPeakB");
  hIntegral1->Write("hIntegralA");
  hIntegral2->Write("hIntegralB");
  hPeakNoise1->Write("hPeakNoiseA");
  hPeakNoise2->Write("hPeakNoiseB");
  hIntegralNoise1->Write("hIntegralNoiseA");
  hIntegralNoise2->Write("hIntegralNoiseB");
  hDeltaT->Write("hDeltaT");
  hChargeAsymmetry->Write("hChargeAsymmetry");
  hChargeAsymmetryVSdeltaT->Write("hChargeAsymmetryVSdeltaT");
  hChargeAsymmetry2->Write("hChargeAsymmetry2");
  hChargeAsymmetryVSdeltaT2->Write("hChargeAsymmetryVSdeltaT2");

  hDifference1->Write("hDifferenceA");
  hDifference2->Write("hDifferenceB");

  hPeakAvsPeakB->SetTitle(";Photo-electrons PMT A;Photo-electrons PMT B");
  hPeakAvsTimeA->SetTitle(";Photo-electrons PMT A;Peak time PMT A");
  hPeakBvsTimeB->SetTitle(";Photo-electrons PMT B;Peak time PMT B");
  hPeakAvsDeltaT->SetTitle(";Photo-electrons PMT A;Time difference PMT A - PMT B");
  hPeakBvsDeltaT->SetTitle(";Photo-electrons PMT B;Time difference PMT A - PMT B");


  hPeakAvsPeakB->Write("hPeakAvsPeakB");
  hPeakAvsTimeA->Write("hPeakAvsTimeA");
  hPeakBvsTimeB->Write("hPeakBvsTimeB");
  hPeakAvsDeltaT->Write("hPeakAvsDeltaT");
  hPeakBvsDeltaT->Write("hPeakBvsDeltaT");
  hExpAvsDeltaT->Write("hExpAvsDeltaT");
  hExpBvsDeltaT->Write("hExpBvsDeltaT");

  hPosResFromDeltaT->Write("hPosResFromDeltaT");
  hPosResFromAsymm->Write("hPosResFromAsymm");
  hPosResFromAsymm2->Write("hPosResFromAsymm2");
  hPosRes->Write("hPosRes");
  hPosRes2->Write("hPosRes2");

  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
  c1->Divide(2,2);
  c1->cd(1);
  hIntegral1->Draw("histo");
  hIntegralNoise1->Draw("histo same");

  c1->cd(3);
  hDifference1->Draw("histo");

  c1->cd(2);
  hIntegral2->Draw("histo");
  hIntegralNoise2->Draw("histo same");

  c1->cd(4);
  hDifference2->Draw("histo");

  string outputc1 = basename+"_c1";
  c1->Write("c1");
  c1->Print(Form("%s.png", outputc1.c_str()));
  c1->Print(Form("%s.pdf", outputc1.c_str()));


  int imax = hDeltaT->GetMaximumBin();
  double timeXcenter = hDeltaT->GetXaxis()->GetBinCenter(imax);

  TCanvas *c2 = new TCanvas("c2", "", 1200, 500);
  c2->Divide(2,1);
  c2->cd(1);
  hDeltaT->Draw("e");
  TF1 *func = new TF1("func", "gaus");
  hDeltaT->Fit(func);

  hDeltaTsimple->SetLineColor(kRed);
  hDeltaTsimple->Draw("e same");

  string outputc2 = basename+"_c2";
  
  c2->cd(2);
  hChargeAsymmetry->SetTitle(Form("%s;Charge Asymmetry;Entries", titleStr.c_str()));
  hChargeAsymmetry->SetLineWidth(2);
  hChargeAsymmetry->Draw("e");
  TF1 *func2 = new TF1("func2", "gaus");
  hChargeAsymmetry->Fit(func2);

  hChargeAsymmetry2->SetLineWidth(2);
  hChargeAsymmetry2->SetLineColor(kRed);
  hChargeAsymmetry2->Draw("e same");

  c2->Write("c2");

  c2->Print(Form("%s.png", outputc2.c_str()));
  c2->Print(Form("%s.pdf", outputc2.c_str()));

  TCanvas *c3 = new TCanvas("c3", "", 1200, 800);
  c3->Divide(2,2);
  c3->cd(1);
  hPeakAvsTimeA->Draw("colz");
  c3->cd(2);
  hPeakBvsTimeB->Draw("colz");
  c3->cd(3);
  hPeakAvsDeltaT->Draw("colz");
  c3->cd(4);
  hPeakBvsDeltaT->Draw("colz");

  string outputc3 = basename+"_c3";

  c3->Print(Form("%s.png", outputc3.c_str()));
  c3->Print(Form("%s.pdf", outputc3.c_str()));


  TCanvas *c4 = new TCanvas("c4", "", 1200, 800);
  c4->Divide(2,2);
  c4->cd(1);
  hExpAvsDeltaT->Draw("colz");
  c4->cd(2);
  hExpBvsDeltaT->Draw("colz");
  c4->cd(3);
  hPeakAvsPeakB->Draw("colz");
  c4->cd(4);
  hChargeAsymmetryVSdeltaT->Draw("colz");

  string outputc4 = basename+"_c4";

  c4->Print(Form("%s.png", outputc4.c_str()));
  c4->Print(Form("%s.pdf", outputc4.c_str()));


  TCanvas *c5 = new TCanvas ("c5", "", 1200, 800);
  c5->Divide(2,1);
  c5->cd(1);
  hPosRes->SetTitle("Position resolution;Position resolution [cm];Entries");
  hPosResFromDeltaT->SetTitle("Position resolution;Position resolution [cm];Entries");
  hPosResFromDeltaT->SetLineColor(kBlack);
  hPosResFromAsymm->SetLineColor(kBlue);
  hPosResFromAsymm2->SetLineColor(kRed);

  hPosResFromDeltaT->SetLineWidth(2);
  hPosResFromAsymm->SetLineWidth(2);
  hPosRes->SetLineWidth(2);
  hPosResFromAsymm2->SetLineWidth(2);
  hPosRes2->SetLineWidth(2);


  double ymax1 = hPosRes->GetMaximum();
  double ymax2 = hPosResFromDeltaT->GetMaximum();
  double ymax3 = hPosResFromAsymm->GetMaximum();

  double allymax = ymax1 > ymax2 ? ymax1 : ymax2;
  allymax        = ymax3 > allymax ? ymax3 : allymax;

  hPosRes->SetMaximum(allymax);
  hPosResFromDeltaT->SetMaximum(allymax);

  hPosResFromDeltaT->Draw("histo");
  hPosResFromAsymm->Draw("histo same");
  hPosResFromAsymm2->Draw("histo same");
  

  TLegend *legres = new TLegend(0.5, 0.7, 0.89, 0.89);
  legres->AddEntry(hPosResFromDeltaT, "From time difference", "l");
  legres->AddEntry(hPosResFromAsymm,  "From pe asymmetry", "l");
  legres->AddEntry(hPosResFromAsymm2, "From integral asymmetry", "l");
  legres->Draw();

  c5->cd(2);
  hPosRes2->SetLineColor(kRed);

  gStyle->SetOptFit(1);
  hPosRes->Draw("histo");
  hPosRes->Fit("gaus", "sames");
  
  hPosRes2->Draw("histo sames");
  hPosRes2->Fit("gaus", "sames");

  // TPaveStats *st = (TPaveStats*)hPosRes2->FindObject("fitstats");
  // st->SetY1NDC(0.3);
  // st->SetY2NDC(0.6);
  // st->SetLineColor(kRed);

  c5->Update();

  TLegend *legres2 = new TLegend(0.11, 0.7, 0.5, 0.89);
  legres2->AddEntry(hPosRes, "Estimator (pe, dt)", "l");
  legres2->AddEntry(hPosRes2, "Estimator (integral, dt)", "l");
  legres2->Draw();


  string outputc5 = basename+"_c5";
  c5->Print(Form("%s.png", outputc5.c_str()));
  c5->Print(Form("%s.pdf", outputc5.c_str()));

  c5->Write("c5");

  fout->Close();
  cout << outputname << endl;

  ofstream myfile;
  myfile.open (Form("%s.txt", basename.c_str()));
  myfile << "Summary output for bar " << barname << "\n";
  myfile << "PMT A : " << whichPMT[0] << " with gain " << gains[0] << "\n";
  myfile << "PMT B : " << whichPMT[1] << " with gain " << gains[1] << "\n";
  myfile << "Triggering on " << chanNice[ichan] << " with threshold -" << thresh << "\n";
  myfile << "Acquisition time source and cosmics: \n";
  myfile << deltaTsource << " " << deltaTnoise << "\n";
  myfile << "Timing resolution mean and RMS: \n";
  myfile << func->GetParameter(1) << " " << func->GetParError(1) << " " << func->GetParameter(2) << " " << func->GetParError(2) << "\n";
  myfile << "Photo-electrons mean and RMS for PMT A and B: \n";
  myfile << hDifference1->GetMean() << " " << hDifference1->GetRMS() << " " << hDifference2->GetMean() << " " << hDifference2->GetRMS() << "\n";
  myfile << "Charge Asymmetry Parameters: \n";
  myfile << func2->GetParameter(1) << " " << func2->GetParError(1) << " " << func2->GetParameter(2) << " " << func2->GetParError(2) << "\n";
  myfile.close();

  delete hPeak1;
  delete hPeak2;
  delete hPeakNoise1;
  delete hPeakNoise2;
  delete hDeltaT;
  delete hDeltaTsimple;
  delete hChargeAsymmetry;
  delete hChargeAsymmetryVSdeltaT;
  delete hChargeAsymmetry2;
  delete hChargeAsymmetryVSdeltaT2;
  delete hIntegral1;
  delete hIntegral2;
  delete hIntegralNoise1;
  delete hIntegralNoise2;
  delete hDifference1;
  delete hDifference2;
  delete hPosResFromDeltaT;
  delete hPosResFromAsymm;
  delete hPosResFromAsymm2;
  delete hPosRes;
  delete hPosRes2;


  delete hPeakAvsPeakB ;
  delete hPeakAvsTimeA ;
  delete hPeakBvsTimeB ;
  delete hPeakAvsDeltaT;
  delete hPeakBvsDeltaT;
  delete hExpAvsDeltaT ;
  delete hExpBvsDeltaT ;


}


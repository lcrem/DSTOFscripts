#include <sys/stat.h>
#include <time.h>
#include <stdio.h>
#include "FFTtools.h"

double getDeltaT(string inputname1);
double getIntegralFromHisto(TGraph *g);
int findWhichPMT(string barname, string whichPMT[2]);
int findGains(string whichPMT[2], double gains[2]);
void zeroGraph(TGraph *g);
bool qualityCut(int ichan, string sourcePos, double peakA, double peakB, double deltat);


string trigChan[2] = {"Ch3", "Ch4"};
string chanNice[2] = {"PMTA", "PMTB"};

void makeSomePlots(string base, string date, string barname, string sourcePos, int ichan, string thresh, string div){

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
  TH1D *hIntegral1 = new TH1D ("hIntegral1", "", 100, 0, 400);
  TH1D *hIntegral2 = new TH1D ("hIntegral2", "", 100, 0, 400);
  TH1D *hDifference1 = new TH1D ("hDifference1", "", 100, 0, 400);
  TH1D *hDifference2 = new TH1D ("hDifference2", "", 100, 0, 400);
  TH1D *hPeakNoise1 = new TH1D ("hPeakNoise1", "", 50, 0, 0.2);
  TH1D *hPeakNoise2 = new TH1D ("hPeakNoise2", "", 50, 0, 0.2);
  TH1D *hIntegralNoise1 = new TH1D ("hIntegralNoise1", "", 100, 0, 400);
  TH1D *hIntegralNoise2 = new TH1D ("hIntegralNoise2", "", 100, 0, 400);
  TH1D *hDeltaT = new TH1D ("hDeltaT", "", 100, -2e-8, 2e-8);
  TH1D *hDeltaTsimple = new TH1D("hDeltaTsimple", "", 100, -2e-8, 2e-8);


  TH2D *hPeakAvsPeakB = new TH2D("hPeakAvsPeakB", "", 50, 0, 400, 50, 0, 400);
  TH2D *hPeakAvsTimeA = new TH2D("hPeakAvsTimeA", "", 50, 0, 400, 50, -20e-9, 20e-9);
  TH2D *hPeakBvsTimeB = new TH2D("hPeakBvsTimeB", "", 50, 0, 400, 50, -20e-9, 20e-9);
  TH2D *hPeakAvsDeltaT = new TH2D("hPeakAvsDeltaT", "", 50, 0, 400, 50, -20e-9, 20e-9);
  TH2D *hPeakBvsDeltaT = new TH2D("hPeakBvsDeltaT", "", 50, 0, 400, 50, -20e-9, 20e-9);

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
     deltatold = maxTime1 - maxTime2;

  
     
     TGraph *gCor = FFTtools::getInterpolatedCorrelationGraph(g1, g2, (1/2.6)*1e-9);
     double *x = gCor->GetX();
     double *y = gCor->GetY();
     int ntot = gCor->GetN();
     deltat = x[TMath::LocMax(ntot, y)];
     
     if (qualityCut(ichan, sourcePos, numPhotoEle1, numPhotoEle2, deltat)) continue;

     //     deltat = x[TMath::LocMax(ntot, y)];
     //     cout << igraph << " " << deltat << " " << x[TMath::LocMax(ntot, y)] << " " << deltat - x[TMath::LocMax(ntot, y)] << endl;
     hDeltaT->Fill(deltat);
  
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


   for (int igraph=1; igraph<=ngraph; igraph++){

    TGraph *g1n = (TGraph*)f1noise->Get(Form("graph%d", igraph));
    if (!g1n){
      continue;
    }
    min1 = TMath::MinElement(g1n->GetN(), g1n->GetY());

    TGraph *g2n = (TGraph*)f2noise->Get(Form("graph%d", igraph));
    if (!g2n){
      continue;
    }
    min2 = TMath::MinElement(g2n->GetN(), g2n->GetY());

    integralNoise1 = getIntegralFromHisto(g1n);
    integralNoise2 = getIntegralFromHisto(g2n);
    numPhotoEle1 = integralNoise1/(gains[0]*R*electron);
    numPhotoEle2 = integralNoise2/(gains[1]*R*electron);
    
    TGraph *gCor = FFTtools::getInterpolatedCorrelationGraph(g1n, g2n, (1/2.6)*1e-9);
    double *x = gCor->GetX();
    double *y = gCor->GetY();
    int ntot = gCor->GetN();
    deltat = x[TMath::LocMax(ntot, y)];
    
    if (qualityCut(ichan, sourcePos, numPhotoEle1, numPhotoEle2, deltat)) continue;
    
    
    hPeakNoise1->Fill(TMath::Abs(min1));
    hPeakNoise2->Fill(TMath::Abs(min2));

    hIntegralNoise1->Fill(numPhotoEle1);
    hIntegralNoise2->Fill(numPhotoEle2);

    delete g1n;
    delete g2n;
    delete gCor;

   }

  f1signal->Close();
  f2signal->Close();
  f1noise->Close();
  f2noise->Close();


  double deltaTsource = getDeltaT(inputname1);

  double deltaTnoise = getDeltaT(inputnoise1);
  
  cout << "Delta T source and noise : " << deltaTsource << " " << deltaTnoise << endl;

  double scale = deltaTsource/deltaTnoise;

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

  TCanvas *c2 = new TCanvas("c2");
  hDeltaT->Draw("e");
  TF1 *func = new TF1("func", "gaus");
  hDeltaT->Fit(func, "r", "", timeXcenter-4e-9, timeXcenter+4e-9);

  hDeltaTsimple->SetLineColor(kRed);
  hDeltaTsimple->Draw("e same");

  string outputc2 = basename+"_c2";
  
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

  string outputc4 = basename+"_c4";

  c4->Print(Form("%s.png", outputc4.c_str()));
  c4->Print(Form("%s.pdf", outputc4.c_str()));

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
  myfile.close();

  delete hPeak1;
  delete hPeak2;
  delete hPeakNoise1;
  delete hPeakNoise2;
  delete hDeltaT;
  delete hDeltaTsimple;
  delete hIntegral1;
  delete hIntegral2;
  delete hIntegralNoise1;
  delete hIntegralNoise2;
  delete hDifference1;
  delete hDifference2;

  delete hPeakAvsPeakB ;
  delete hPeakAvsTimeA ;
  delete hPeakBvsTimeB ;
  delete hPeakAvsDeltaT;
  delete hPeakBvsDeltaT;
  delete hExpAvsDeltaT ;
  delete hExpBvsDeltaT ;


}


bool qualityCut(int ichan, string sourcePos, double peakA, double peakB, double deltat){

  int length = stoi(sourcePos.substr (0, sourcePos.size()-2))+6;
  int barlength = 140;
  double c = 3.e8/1.6;
  // Add 2 ns to deltat expected to account for difference in cable lengths
  double deltaTexpected = (length - (barlength - length))*0.01/c - 4.5E-9;
  double mindeltat = deltaTexpected - 2e-9;
  double maxdeltat = deltaTexpected + 2e-9;
  double expectedPeakFraction = (barlength-length)*1.0/length;
  double minPeak = expectedPeakFraction*0.75;
  double maxPeak = expectedPeakFraction*1.25;
  
  bool cut = false;

  
  //  cout << deltaTexpected << endl;

  if ( ( peakA/peakB<minPeak  || peakA/peakB>maxPeak) &&
       (deltat < mindeltat || deltat>maxdeltat) )  cut = true;
  
  

  return cut;
}

double getDeltaT(string inputname1){
  
  ifstream f0 (Form("%s.times", inputname1.c_str()));
  string line;
  double acqTot = 0.;
  double acq0;
  while ( getline(f0, line)){
     std::size_t pos = line.find(")");
     std::string add0 = line.substr (pos+1, line.size()-1);
     acq0 = stod(add0);
     acqTot += acq0;
     //     cout << add0 << " " << line << endl;
  }
  f0.close();

  return acqTot;

}


double getIntegralFromHisto(TGraph *g){

  double *x = g->GetX();
  double *y = g->GetY();
  int N = g->GetN();
  TH1D *h = new TH1D("h", "", 1000, x[0], x[N-1]);
  double integral=0;
  for (int i=N*0; i<N; i++){

    //  integral+=y[i];
    h->Fill(x[i], -y[i]);
  }

  integral = h->Integral("width");
  
  delete h;

  return integral;
}


int findWhichPMT(string barname, string whichPMT[2]){

  ifstream f0 ("barToPMTmap.txt");
  string line;
  string thisbar;
  while ( getline(f0, line)){
    stringstream ss(line);
    ss >> thisbar;
    if (thisbar==barname){
      ss >> whichPMT[0];
      ss >> whichPMT[1];
      cout << "Processing bar " << barname << endl;
      cout << "PMT A is " << whichPMT[0] << endl;
      cout << "PMT B is " << whichPMT[1] << endl;
      f0.close();
      return 1;
    }

  }

  cout << "I can't recognise bar " << barname << endl;

  return 0;
}


int findGains(string whichPMT[2], double gains[2]){

  ifstream f0 ("PMTtoGainMap.txt");
  string line;
  string linePMT[2];
  string temp;
  string thispmt;
  bool found[2];
  found[0] = false;
  found[1] = false;

  getline(f0, line);

  while ( getline(f0, line)){
    stringstream ss(line);
    ss >> thispmt;
    //    cout << thispmt << " " << whichPMT[0] << " " << whichPMT[1] << endl;
    for (int i=0; i<2; i++){
      if (thispmt==whichPMT[i]){
	linePMT[i] = line;
	ss >> temp;
	gains[i] = stod(temp);
      }
    }
    
  }
  f0.close();

  double voltage1, voltage2;
  double gainp0, gainerr0, gainp1, gainerr1;
  for (int i=0; i<2; i++){
    stringstream ss(linePMT[i]);
    ss >> thispmt >> voltage1 >> gains[i] >> voltage2;
    cout << thispmt << " " << voltage1 << " " << voltage2 << endl;
    if (voltage1!=voltage2){
      if (voltage2==0){
	cout << "You need to write the voltages used for PMT " << thispmt << endl;
	return 0;
      }
      
      ifstream f1("PMTtoGainParameterization.txt");      
      getline(f1, line);
      while ( getline(f1, line)){
	stringstream ss1(line);
	ss1 >> thispmt;
	if (thispmt==whichPMT[i]){
	  ss1 >> gainp0 >> gainerr0 >> gainp1 >> gainerr1;
	  gains[i] = gainp0*TMath::Power(voltage2, gainp1);
	  found[i] = true;	  
	}	
      }
      
      f1.close();
    } else {
      found[i] = true;
    }
  }

  if (found[0]==found[1]==true){
    cout << "PMT A gain: " << gains[0] << endl;
    cout << "PMT B gain: " << gains[1] << endl;
    return 1;
  }

  cout << "I can't recognise one of these PMTs " << whichPMT[0] << " or " << whichPMT[1] << endl;

  return 0;

}


void zeroGraph(TGraph *g){

  double *y = g->GetY();
  double mean = 0;
  int ntot = 200;

  for (int i=0; i<ntot; i++) mean+= y[i]/ntot;

  for (int i=0; i<g->GetN(); i++) g->GetY()[i]-=mean;			       

}

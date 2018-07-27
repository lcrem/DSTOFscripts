string base="/unix/dune/tof/";
string date="2018Jul10";
string barname="D1721";
string POSITIONS[]={"16cm", "32cm", "48cm", "64cm", "80cm", "96cm", "112cm"};
string THRESHOLDS[]={"10mV", "15mV", "20mV"};
int colors []      ={kBlue+1, kCyan+1, kGreen+1, kOrange+1, kRed, kViolet};
string CHANNELS[]={"PMTA", "PMTB"};
string divisions="20mVdiv";

void makeSummaryPlotsCh(int ich);

void makeSummaryPlots(){

  makeSummaryPlotsCh(0);
  makeSummaryPlotsCh(1);

}

void makeSummaryPlotsCh(int ich){

  int npos = sizeof(POSITIONS)/sizeof(POSITIONS[0]);
  int nthresh = sizeof(THRESHOLDS)/sizeof(THRESHOLDS[0]);
  int nch = sizeof(CHANNELS)/sizeof(CHANNELS[0]);
  

  TGraph *gAcqTime[10];
  TGraphErrors *gTimeDifference[10];
  TGraphErrors *gTimeResolution[10];
  TGraphErrors *gAsymmetryPeak[10];
  TGraphErrors *gAsymmetryWidth[10];
  TGraphErrors *gPEdistPMTA[10];
  TGraphErrors *gPEdistPMTB[10];



  double x[10], xerr[10], acqTime[10];
  double timesource, timecosmics, maxtime;
  maxtime=0;
  double resmean[10], resmeanerr[10], resrms[10], resrmserr[10];
  double asymean[10], asymeanerr[10], asyrms[10], asyrmserr[10];
  double pemeanA[10], permsA[10], pemeanB[10], permsB[10];
  double maxpeA, maxpeB;
  maxpeA=maxpeB=0;

  string baseoutputname =   base + date + "/Bar" + barname + "/Trig" + CHANNELS[ich]  ;

  TLegend *leg = new TLegend(0.8, 0.75, 0.99, 0.99, "Trigger thresholds");

  for (int ithresh=0; ithresh<nthresh; ithresh++){
   
    for (int ipos=0; ipos<npos; ipos++){
      
      timesource=timecosmics=0;
      string basename   =   base + date + "/Bar" + barname + "/SourceAt" + POSITIONS[ipos] + "/Trig" + CHANNELS[ich] + "_thres" + THRESHOLDS[ithresh] ;
      ifstream f0 (Form("%s.txt", basename.c_str()));
      string line;
      getline(f0, line);
      getline(f0, line);
      getline(f0, line);
      getline(f0, line);
      getline(f0, line);
      getline(f0, line);
      stringstream ss(line);
      ss >> timesource >> timecosmics;
      
      x[ipos] = 16*(ipos+1);
      xerr[ipos] = 0;
      acqTime[ipos] = timesource;
      if (acqTime[ipos]>maxtime) maxtime=acqTime[ipos];

      getline(f0, line);
      getline(f0, line);
      stringstream ss1(line);
      ss1 >> resmean[ipos] >> resmeanerr[ipos] >> resrms[ipos] >> resrmserr[ipos];

      getline(f0, line);
      getline(f0, line);
      stringstream ss2(line);
      ss2 >> pemeanA[ipos] >> permsA[ipos] >> pemeanB[ipos] >> permsB[ipos];
      if (pemeanA[ipos]>maxpeA) maxpeA = pemeanA[ipos];
      if (pemeanB[ipos]>maxpeB) maxpeB = pemeanB[ipos];

      getline(f0, line);
      getline(f0, line);
      stringstream ss3(line);
      ss3 >> asymean[ipos] >> asymeanerr[ipos] >> asyrms[ipos] >> asyrmserr[ipos];

      f0.close();


    }   

    gAcqTime[ithresh] = new TGraph(npos, x, acqTime);
    gAcqTime[ithresh]->SetLineColor(colors[ithresh]);
    gAcqTime[ithresh]->SetLineWidth(2);
    leg->AddEntry(gAcqTime[ithresh], THRESHOLDS[ithresh].c_str(), "l");

    gTimeDifference[ithresh] = new TGraphErrors(npos, x, resmean, xerr, resmeanerr);
    gTimeDifference[ithresh]->SetLineColor(colors[ithresh]);
    gTimeDifference[ithresh]->SetLineWidth(2);

    gTimeResolution[ithresh] = new TGraphErrors(npos, x, resrms, xerr, resrmserr);
    gTimeResolution[ithresh]->SetLineColor(colors[ithresh]);
    gTimeResolution[ithresh]->SetLineWidth(2);

    gAsymmetryPeak[ithresh] = new TGraphErrors(npos, x, asymean, xerr, asymeanerr);
    gAsymmetryPeak[ithresh]->SetLineColor(colors[ithresh]);
    gAsymmetryPeak[ithresh]->SetLineWidth(2);

    gAsymmetryWidth[ithresh] = new TGraphErrors(npos, x, asyrms, xerr, asyrmserr);
    gAsymmetryWidth[ithresh]->SetLineColor(colors[ithresh]);
    gAsymmetryWidth[ithresh]->SetLineWidth(2);

    gPEdistPMTA[ithresh] = new TGraphErrors(npos, x, pemeanA, xerr, permsA);
    gPEdistPMTA[ithresh]->SetLineColor(colors[ithresh]);
    gPEdistPMTA[ithresh]->SetLineWidth(2);

    gPEdistPMTB[ithresh] = new TGraphErrors(npos, x, pemeanB, xerr, permsB);
    gPEdistPMTB[ithresh]->SetLineColor(colors[ithresh]);
    gPEdistPMTB[ithresh]->SetLineWidth(2);

  }


  TCanvas *cAcqTime = new TCanvas("cAcqTime");
  cAcqTime->SetLogy();
  
  for (int ithresh=0; ithresh<nthresh; ithresh++){
    if (ithresh==0){
      gAcqTime[ithresh]->SetMaximum(maxtime*2);
      gAcqTime[ithresh]->SetTitle(Form("Bar %s, triggering on %s;Source distance [cm];Acquisition time [s]", barname.c_str(), CHANNELS[ich].c_str()));
      gAcqTime[ithresh]->Draw("Al");
    } else {
      gAcqTime[ithresh]->Draw("l");
    }
  }
  leg->Draw();
  cAcqTime->Print((baseoutputname+"_acqTime.png").c_str());
  cAcqTime->Print((baseoutputname+"_acqTime.pdf").c_str());


  TCanvas *cTimeDifference = new TCanvas("cTimeDifference");
 
  for (int ithresh=0; ithresh<nthresh; ithresh++){
    if (ithresh==0){
      gTimeDifference[ithresh]->SetMinimum(-1.5e-8);
      gTimeDifference[ithresh]->SetMaximum(+1.5e-8);
      gTimeDifference[ithresh]->SetTitle(Form("Bar %s, triggering on %s;Source distance [cm];Time Difference [s]", barname.c_str(), CHANNELS[ich].c_str()));
      gTimeDifference[ithresh]->Draw("Al");
    } else {
      gTimeDifference[ithresh]->Draw("l");
    }
    gTimeDifference[ithresh]->Fit("pol1");
  }
  leg->Draw();
  cTimeDifference->Print((baseoutputname+"_timeDifference.png").c_str());
  cTimeDifference->Print((baseoutputname+"_timeDifference.pdf").c_str());


  TCanvas *cTimeResolution = new TCanvas("cTimeResolution");
 
  for (int ithresh=0; ithresh<nthresh; ithresh++){
    if (ithresh==0){
      gTimeResolution[ithresh]->SetMinimum(0);
      gTimeResolution[ithresh]->SetMaximum(+2e-9);
      gTimeResolution[ithresh]->SetTitle(Form("Bar %s, triggering on %s;Source distance [cm];Time Resolution [s]", barname.c_str(), CHANNELS[ich].c_str()));
      gTimeResolution[ithresh]->Draw("Al");
    } else {
      gTimeResolution[ithresh]->Draw("l");
    }
  }
  leg->Draw();
  cTimeResolution->Print((baseoutputname+"_timeResolution.png").c_str());
  cTimeResolution->Print((baseoutputname+"_timeResolution.pdf").c_str());

 TCanvas *cAsymmetryPeak = new TCanvas("cAsymmetryPeak");
 
  for (int ithresh=0; ithresh<nthresh; ithresh++){
    if (ithresh==0){
      gAsymmetryPeak[ithresh]->SetMinimum(-1);
      gAsymmetryPeak[ithresh]->SetMaximum(+1);
      gAsymmetryPeak[ithresh]->SetTitle(Form("Bar %s, triggering on %s;Source distance [cm];Charge asymmetry peak", barname.c_str(), CHANNELS[ich].c_str()));
      gAsymmetryPeak[ithresh]->Draw("Al");
    } else {
      gAsymmetryPeak[ithresh]->Draw("l");
    }
    gAsymmetryPeak[ithresh]->Fit("pol1");
  }
  leg->Draw();
  cAsymmetryPeak->Print((baseoutputname+"_asymmetryPeak.png").c_str());
  cAsymmetryPeak->Print((baseoutputname+"_asymmetryPeak.pdf").c_str());


 TCanvas *cAsymmetryWidth = new TCanvas("cAsymmetryWidth");
 
  for (int ithresh=0; ithresh<nthresh; ithresh++){
    if (ithresh==0){
      gAsymmetryWidth[ithresh]->SetMinimum(0);
      gAsymmetryWidth[ithresh]->SetMaximum(0.5);
      gAsymmetryWidth[ithresh]->SetTitle(Form("Bar %s, triggering on %s;Source distance [cm];Charge asymmetry width", barname.c_str(), CHANNELS[ich].c_str()));
      gAsymmetryWidth[ithresh]->Draw("Al");
    } else {
      gAsymmetryWidth[ithresh]->Draw("l");
    }
  }
  leg->Draw();
  cAsymmetryWidth->Print((baseoutputname+"_asymmetryWidth.png").c_str());
  cAsymmetryWidth->Print((baseoutputname+"_asymmetryWidth.pdf").c_str());


  TCanvas *cPEpmtA = new TCanvas("cPEpmtA");
 
  for (int ithresh=0; ithresh<nthresh; ithresh++){
    if (ithresh==0){
      gPEdistPMTA[ithresh]->SetMaximum(maxpeA*1.5);
      gPEdistPMTA[ithresh]->SetTitle(Form("Bar %s, triggering on %s, PMT A;Source distance [cm];Photo-electrons", barname.c_str(), CHANNELS[ich].c_str()));
      gPEdistPMTA[ithresh]->Draw("Al");
    } else {
      gPEdistPMTA[ithresh]->Draw("l");
    }
  }
  leg->Draw();
  cPEpmtA->Print((baseoutputname+"_photoelectronsPMTA.png").c_str());
  cPEpmtA->Print((baseoutputname+"_photoelectronsPMTA.pdf").c_str());

  TCanvas *cPEpmtB = new TCanvas("cPEpmtB");
 
  for (int ithresh=0; ithresh<nthresh; ithresh++){
    if (ithresh==0){
      gPEdistPMTB[ithresh]->SetMaximum(maxpeB*1.5);
      gPEdistPMTB[ithresh]->SetTitle(Form("Bar %s, triggering on %s, PMT B;Source distance [cm];Photo-electrons", barname.c_str(), CHANNELS[ich].c_str()));
      gPEdistPMTB[ithresh]->Draw("Al");
    } else {
      gPEdistPMTB[ithresh]->Draw("l");
    }
  }
  leg->Draw();
  cPEpmtB->Print((baseoutputname+"_photoelectronsPMTB.png").c_str());
  cPEpmtB->Print((baseoutputname+"_photoelectronsPMTB.pdf").c_str());



}


double barlength = 140.;
double c = 3.e8/2.;
double interpolationDeltaT = 0.01e-9;

double getDeltaT(string inputname1);
double getIntegralFromHisto(TGraph *g);
int findWhichPMT(string barname, string whichPMT[2]);
int findCalibration(string barname, double calvar[2]);
int findGains(string whichPMT[2], double gains[2]);
void zeroGraph(TGraph *g);
bool qualityCut(int ichan, string sourcePos, double peakA, double peakB, double deltat);
double getPositionFromDeltaT(double deltat);
double getPositionFromasymmetry(double deltat);



bool qualityCut(int ichan, string sourcePos, double peakA, double peakB, double deltat){

  double length = stoi(sourcePos.substr (0, sourcePos.size()-2))+6.;
  // Add 2 ns to deltat expected to account for difference in cable lengths
  double deltaTexpected = (length - (barlength - length))*0.01/c;
  double mindeltat = deltaTexpected - 2e-9;
  double maxdeltat = deltaTexpected + 2e-9;
  double expectedPeakFraction = (barlength-length)*1.0/length;
  double minPeak = expectedPeakFraction*0.75;
  double maxPeak = expectedPeakFraction*1.25;
  
  bool cut = false;

  
  //cout << deltaTexpected << endl;

  if ( ( peakA/peakB<minPeak  || peakA/peakB>maxPeak) &&
       (deltat < mindeltat || deltat>maxdeltat) )  cut = true;
  
  

  return cut;
}


double getPositionFromDeltaT(double deltat){


  /* double position = deltat*c/0.01;  */

  /* position+=64.;  */

  //  double position = 0.5 * (deltat*c/0.01 + barlength) - 6.;

  double position = (deltat + 8.67271e-09)/(1.36627e-10);

  return position;
  
}

double getPositionFromAsymmetry(double asymm){

  /* double position = -asymm*barlength*1.; */

  /* position+= barlength*0.5; */

  double position = (asymm - 0.700352)/(-0.0113595);

  return position;

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


int findCalibration(string barname, double calvar[2]){

  ifstream f0 ("barToPMTmap.txt");
  string line;
  string thisbar;
  string temp;
  while ( getline(f0, line)){
    stringstream ss(line);
    ss >> thisbar;
    if (thisbar==barname){
      ss >> temp;
      ss >> temp;
      ss >> calvar[0];
      ss >> calvar[1];
      f0.close();
      if (calvar[0]==calvar[1]){
	cout << "I can't find calibration values for " << barname << endl;
	return 0;
      } else      
      return 1;
    }

  }

  cout << "I can't find calibration values for " << barname << endl;

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



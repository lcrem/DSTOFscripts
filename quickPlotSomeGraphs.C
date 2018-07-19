void quickPlotSomeGraphs(string name, int n=5){

  TGraph *g[100];

  TFile *f = new TFile(name.c_str(), "read");
  int color = 52;

  TCanvas *c1 = new TCanvas("c1");

  for (int i=0; i<n; i++){
    g[i] = (TGraph*)f->Get(Form("graph%d", i+1));

    color+=5;
    g[i]->SetLineColor(color);

    if (i==0) g[i]->Draw("Al");
    else g[i]->Draw("l");

  }
  
  c1->Print("Temp.png");
  
}

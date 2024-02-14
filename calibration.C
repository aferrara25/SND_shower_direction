void calibration() {

  const int NENERGIES = 5;
  const double k_fix = 0.063194; //0.0617134; //0.0579034
  const double alpha_fix = 0.0130885; //0.0111385; //0.010421
  double mean_k = 0;
  double mean_a = 0;
  std::array<int, NENERGIES> energies = {100,140,180,240,300};
  std::array<int, NENERGIES> runs = {100631,100673,100672,100647,100645}; // test: {100631,100673,100672,100647,100645} calib: {100677,100633,100671,100648,100639};
  std::array<double, NENERGIES*3> xmax = {5000,6000,7000, 8000,8000,9000, 10000,11000,12000, 12000,12000,14000, 13500,15000,15000}; //{4000,6000,7000, 8000,8000,9000, 10000,11000,12000, 12000,12000,14000, 14000,14000,15000};
  std::array<double, NENERGIES*3> xmin = {500,1000,2500, 1000,1000,3000, 1000,1000,4000, 1000,1000,4000, 1500,3000,6000}; // {500,500,500, 1000,1000,1500, 1000,1000,1500, 1000,1000,2500, 2000,1000,2500}; 
  std::array<double, NENERGIES*3> amax = {27,25,22, 21,21,20, 21,21,20, 21,21,25, 21,21,25}; //{25,25,20, 25,25,20, 25,25,20, 20,20,23, 25,25,25};
  std::array<double, NENERGIES*3> amin = {3,3,5, 3,3,10, 3,3,10, 3,3,12, 3,3,15}; // *0.001 {-3,-3,1, -5,-5,6, -5,-5,8, -5,0,8, -5,0,15};
  TFile* indata[NENERGIES];
  TGraphErrors* graph[3];
  TGraphErrors* gAlpha[3];
  TGraphErrors* gEn[3];
  TGraphErrors* gEnrel[3];
  TGraphErrors* gErr[3];
  TGraphErrors* gK[3];
  for (int k{0}; k<3; ++k) {
    graph[k] = new TGraphErrors();
    graph[k]->SetTitle(Form("start_in_SciFi%i",k+2));
    gAlpha[k] = new TGraphErrors();
    gAlpha[k]->SetTitle(Form("start_in_SciFi%i",k+2));
    gK[k] = new TGraphErrors();
    gK[k]->SetTitle(Form("start_in_SciFi%i",k+2));
    gEn[k] = new TGraphErrors();
    gEn[k]->SetTitle(Form("start_in_SciFi%i",k+2));
    gEnrel[k] = new TGraphErrors();
    gEnrel[k]->SetTitle(Form("start_in_SciFi%i",k+2));
    gErr[k] = new TGraphErrors();
    gErr[k]->SetTitle(Form("start_in_SciFi%i",k+2));
  }
  TGraphErrors* graph2[NENERGIES];
  for (int k{0}; k<NENERGIES; ++k) {
    graph2[k] = new TGraphErrors();
    graph2[k]->SetTitle(Form("%i GeV",energies[k]));
  }
  double en[NENERGIES*3];
  double SciFIQDC[NENERGIES*3];
  double err[NENERGIES*3];
  EColor colors[5] = {kRed, kBlue, kGreen, kCyan, kYellow};
  TCanvas* c[NENERGIES];
  TCanvas* cAlpha[NENERGIES];
  TCanvas* cEn[NENERGIES];
  TCanvas* cMarco = new TCanvas("c","c", 1920,1080);

  TF1 *fitFunc = new TF1("fitFunc", "[0]*x + [1]");
  fitFunc->SetParameters(-1/4, 4000);
  TF1 *gaus = new TF1("gaus", "gaus");
  gaus->SetParameters(1800, 12*0.001, 4*0.001);
  TF1 *line = new TF1("line", "[0]*x + [1]", 50, 350);
  line->SetParameters(0, 0);

  TF1 *gausM = new TF1("gaus", "gaus");
  gausM->SetParameters(12000, 180, 30);


  TFile* outf = new TFile("calibrations.root","recreate");
  for (int i{0}; i<NENERGIES; ++i) {
    std::cout<<i<<"\n";
    indata[i] = new TFile(Form("calib/TB_outputRun_%i.root",runs[i]),"read");
    c[i] = new TCanvas(Form("US QDC vs SciFi QDC %i GeV",energies[i]),Form("US QDC vs SciFi QDC %i GeV",energies[i]), 1920,1080);
    c[i]->Divide(2,2);
    cAlpha[i] = new TCanvas(Form("Alpha %i GeV",energies[i]),Form("Alpha %i GeV",energies[i]), 1920,1080);
    cAlpha[i]->Divide(3,1);
    cEn[i] = new TCanvas(Form("Reconstructed Energy %i GeV",energies[i]),Form("Reconstructed Energy %i GeV",energies[i]), 1920,1080);
    cEn[i]->Divide(3,1);
    TH1D* hAlpha[3];
    TH1D* hEn[3];
    TH1D* hMarco = new TH1D(Form("Reco_En_%iGeV", energies[i]),Form("Reco_En_%iGeV; Reconstructed energy (GeV); entries", energies[i]), 100, 50, 300);
    c[i]->cd(1);
    ((TH2D*)indata[i]->Get("GuilCut_QDCUS_vs_QDCScifi"))->DrawClone("colz");
    gPad->Modified();
    gPad->Update();
    for (int j{0}; j<3; ++j) {
      SciFIQDC[3*i+j] = ((TH1F*)indata[i]->Get(Form("GuilCut_Shower_SciFi_QDC_shStart%i",j+2)))->GetMean();
      err[3*i+j] = ((TH1F*)indata[i]->Get(Form("GuilCut_Shower_SciFi_QDC_shStart%i",j+2)))->GetStdDev();
      en[3*i+j] = energies[i] + 5*j;
      graph[j]->SetPoint(i,en[3*i+j],SciFIQDC[3*i+j]);
      graph[j]->SetPointError(i,0,err[3*i+j]);
      graph2[i]->SetPoint(j,j+2+0.1*i,SciFIQDC[3*i+j]);
      graph2[i]->SetPointError(j,0,err[3*i+j]);

      cAlpha[i]->cd(j+1);
      TH2D* h = (TH2D*)indata[i]->Get(Form("GuilCut_QDCUS_vs_QDCScifi_ShStart_st%i",j+2));
      int nx = h->GetXaxis()->FindBin(xmax[3*i+j]);
      int ny = h->GetYaxis()->GetNbins();
      int nxmax = h->GetXaxis()->GetNbins();
      hAlpha[j] = new TH1D(Form("Alpha_%iGeV_st%i", energies[i], j+2),Form("Alpha_%iGeV_st%i; Alpha (GeV/QDC); entries", energies[i], j+2), 70, -0.02, 0.04);
      hEn[j] = new TH1D(Form("Reco_En_%iGeV_st%i", energies[i], j+2),Form("Reco_En_%iGeV_st%i; Reconstructed energy (GeV); entries", energies[i], j+2), 100, 0, 500);
      for (int binX = h->GetXaxis()->FindBin(xmin[3*i+j]); binX <= nx; ++binX) {
        for (int binY = 1; binY <= ny; ++binY) {
            // Access the bin content
            double binContent = h->GetBinContent(binX, binY);
            if (binContent>0) {
              double a = (energies[i] - k_fix*(h->GetYaxis()->GetBinCenter(binY)))/(h->GetXaxis()->GetBinCenter(binX));
              hAlpha[j]->Fill(a,binContent);             
            }
        }
      }
      for (int binX = 1; binX <= nxmax; ++binX) {
        for (int binY = 1; binY <= ny; ++binY) {
            // Access the bin content
            double binContent = h->GetBinContent(binX, binY);
            if (binContent>0) {
              double recEn = k_fix*(h->GetYaxis()->GetBinCenter(binY)) + alpha_fix*(h->GetXaxis()->GetBinCenter(binX));
              hEn[j]->Fill(recEn,binContent);
              if (j<2) {hMarco->Fill(recEn,binContent);}             
            }
        }
      }
      gStyle->SetOptFit(1011);
      gaus->SetRange(amin[3*i+j]*0.001, amax[3*i+j]*0.001);
      //gaus->FixParameter(0,hAlpha[j]->GetMaximum());
      hAlpha[j]->Fit(gaus,"RQ");
      hAlpha[j]->Draw("hist colz");

      gPad->Modified();
      gPad->Update();

      gAlpha[j]->SetPoint(i,en[3*i+j],gaus->GetParameter(1));
      gAlpha[j]->SetPointError(i,0,gaus->GetParameter(2));

      cEn[i]->cd(j+1);
      hEn[j]->Draw("hist colz");
      gEn[j]->SetPoint(i,en[3*i+j],hEn[j]->GetMean() - energies[i]);
      gEn[j]->SetPointError(i,0,hEn[j]->GetStdDev());

      gEnrel[j]->SetPoint(i,en[3*i+j],(hEn[j]->GetMean() - energies[i])/energies[i]*100);
      gEnrel[j]->SetPointError(i,0,hEn[j]->GetStdDev()/energies[i]*100);

      gErr[j]->SetPoint(i,en[3*i+j],hEn[j]->GetStdDev()/energies[i]*100);

      if (i==2) {
        cMarco->cd();
        gStyle->SetOptFit(1011);
        gStyle->SetOptStat(0);
        gausM->SetRange(110, 250);
        //gaus->FixParameter(0,hAlpha[j]->GetMaximum());
        hMarco->Fit(gausM,"RQ");
        hMarco->Draw("hist colz");

        gPad->Modified();
        gPad->Update();
      }

      c[i]->cd(j+2);
      h->SetStats(0);
      gStyle->SetOptFit(1011);
      fitFunc->SetRange(xmin[3*i+j], xmax[3*i+j]);
      h->Fit(fitFunc,"RQ"); 
      h->DrawClone("colz");
      //fitFunc->Draw("same");
      gPad->Modified();
      gPad->Update();
      double alpha = -energies[i]*fitFunc->GetParameter(0)/fitFunc->GetParameter(1);
      double erralpha = alpha*(std::abs(fitFunc->GetParError(0)/fitFunc->GetParameter(0))+std::abs(fitFunc->GetParError(1)/fitFunc->GetParameter(1)));
      double k = energies[i]/fitFunc->GetParameter(1);
      double errk = k*fitFunc->GetParError(1)/fitFunc->GetParameter(1);
      // gAlpha[j]->SetPoint(i,en[3*i+j],alpha);
      // gAlpha[j]->SetPointError(i,0,erralpha);
      gK[j]->SetPoint(i,en[3*i+j],k);
      gK[j]->SetPointError(i,0,errk);
      if (j!=2) {
        mean_k += k/(NENERGIES*2);
        mean_a += gaus->GetParameter(1)/(NENERGIES*2);
      }
    }

    outf->cd();
    c[i]->Write();
    cAlpha[i]->Write();
    cEn[i]->Write();
    cMarco->Write();
    indata[i]->Close();
  }
  
  TCanvas* c1 = new TCanvas("SciFi QDC vs Energy","SciFi QDC vs Energy", 1920,1080);
  TCanvas* c2 = new TCanvas("SciFi QDC vs start","SciFi QDC vs start", 1920,1080);
  TCanvas* c3 = new TCanvas("Alpha vs Energy","Alpha vs Energy", 1920,1080);
  TCanvas* c4 = new TCanvas("K vs Energy","K vs Energy", 1920,1080);
  TCanvas* c5 = new TCanvas("E_reco - E","E_reco - E", 1920,1080);
  TCanvas* c6 = new TCanvas("Relative energy resolution","Relative energy resolution", 1920,1080);
  TCanvas* c7 = new TCanvas("(E_reco - E)/E","(E_reco - E)/E", 1920,1080);
  TMultiGraph *mg = new TMultiGraph();
  TMultiGraph *mg2 = new TMultiGraph();
  TMultiGraph *mg3 = new TMultiGraph();
  TMultiGraph *mg4 = new TMultiGraph();
  TMultiGraph *mg5 = new TMultiGraph();
  TMultiGraph *mg6 = new TMultiGraph();
  TMultiGraph *mg7 = new TMultiGraph();

  for (int k{0}; k<3; ++k) {
    graph[k]->SetMarkerColor(colors[k]);
    graph[k]->SetMarkerStyle(21);
    mg->Add(graph[k]);
    gAlpha[k]->SetMarkerColor(colors[k]);
    gAlpha[k]->SetMarkerStyle(21);
    mg3->Add(gAlpha[k]);
    gK[k]->SetMarkerColor(colors[k]);
    gK[k]->SetMarkerStyle(21);
    mg4->Add(gK[k]);
    gEn[k]->SetMarkerColor(colors[k]);
    gEn[k]->SetMarkerStyle(21);
    mg5->Add(gEn[k]);
    gErr[k]->SetMarkerColor(colors[k]);
    gErr[k]->SetMarkerStyle(21);
    mg6->Add(gErr[k]);
    gEnrel[k]->SetMarkerColor(colors[k]);
    gEnrel[k]->SetMarkerStyle(21);
    mg7->Add(gEnrel[k]);
  }
    for (int k{0}; k<NENERGIES; ++k) {
    graph2[k]->SetMarkerColor(colors[k]);
    graph2[k]->SetMarkerStyle(21);
    mg2->Add(graph2[k]);
  }
  c1->cd();
  mg->GetXaxis()->SetTitle("pion energy (GeV)");
  mg->GetYaxis()->SetTitle("mean Scifi QDC (a.u.)");
  mg->GetXaxis()->SetLimits(50,350);
  mg->SetMinimum(0.);
  mg->SetMaximum(5000.);
  mg->Draw("ape");
  c1->BuildLegend();
  c2->cd();
  mg2->GetXaxis()->SetTitle("shower start");
  mg2->GetYaxis()->SetTitle("mean Scifi QDC (a.u.)");
  mg2->GetXaxis()->SetLimits(1.5,4.5);
  mg2->SetMinimum(0.);
  mg2->SetMaximum(5000.);
  mg2->Draw("ape");
  c2->BuildLegend();
  c3->cd();
  mg3->GetXaxis()->SetTitle("pion energy (GeV)");
  mg3->GetYaxis()->SetTitle("Alpha (GeV/QDC)");
  mg3->GetXaxis()->SetLimits(50,350);
  mg3->SetMinimum(-0.002);
  mg3->SetMaximum(0.025);
  mg3->Draw("ape");
  c3->BuildLegend();
  c4->cd();
  mg4->GetXaxis()->SetTitle("pion energy (GeV)");
  mg4->GetYaxis()->SetTitle("K (GeV/QDC)");
  mg4->GetXaxis()->SetLimits(50,350);
  mg4->SetMinimum(0.03);
  mg4->SetMaximum(0.15);
  mg4->Draw("ape");
  c4->BuildLegend();
  c5->cd();
  mg5->GetXaxis()->SetTitle("pion energy (GeV)");
  mg5->GetYaxis()->SetTitle("E_reco - E (GeV)");
  mg5->GetXaxis()->SetLimits(50,350);
  mg5->SetMinimum(-100);
  mg5->SetMaximum(100);
  mg5->Draw("ape");
  c5->BuildLegend();
  line->Draw("same");
  c6->cd();
  mg6->GetXaxis()->SetTitle("pion energy (GeV)");
  mg6->GetYaxis()->SetTitle("Resolution %");
  mg6->GetXaxis()->SetLimits(50,350);
  mg6->SetMinimum(0);
  mg6->SetMaximum(40);
  mg6->Draw("ape");
  c6->BuildLegend();
  c7->cd();
  mg7->GetXaxis()->SetTitle("pion energy (GeV)");
  mg7->GetYaxis()->SetTitle("(E_reco - E)/E %");
  mg7->GetXaxis()->SetLimits(50,350);
  mg7->SetMinimum(-100);
  mg7->SetMaximum(100);
  mg7->Draw("ape");
  c7->BuildLegend();
  line->Draw("same");
  outf->cd();
  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();
  c5->Write();
  c6->Write();
  c7->Write();
  outf->Close();
  std::cout<<"k:\t"<<mean_k<<"\na:\t"<<mean_a<<"\n";
}
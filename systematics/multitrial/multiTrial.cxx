#include "../../help.h"
#include "../../effHelp.h"

void DivideAndComputeRogerBarlow(TH1F *h1, TH1F *h2);
double PassRogerBarlowCriterion(int nsigmas, Double_t dev, Double_t RBsigma);
TH1F *makeSystPlots(int num = 1);

void multiTrial()
{
  gROOT->SetBatch(kTRUE);
  const int trials = 1000;
  TH1F *h[trials];
  for (int i = 0; i < trials; i++){
    h[i] = makeSystPlots(i+1);
  }

  const int bins = h[0]->GetNbinsX();
  TH1F* hPtDev[bins];
  TF1* fgaus[bins];

  Int_t avBinning = 30;
  Int_t nBins[18] = {avBinning-15, avBinning-15, avBinning-5, avBinning, avBinning, 
                     avBinning, avBinning, avBinning, avBinning, avBinning, 
                     avBinning, avBinning, avBinning, avBinning, avBinning, 
                     avBinning-5, avBinning-5, avBinning-15};

  for (int i = 0; i < bins; i++){
    //some pt bins have different binning for the deviation histo
    hPtDev[i] = new TH1F(Form("hPtDev%i", i), Form("%.2f < #it{p}_{T} < %.2f GeV/#it{c};1 - #it{Y}_{sys}/#it{Y}_{def};Counts", h[0]->GetBinLowEdge(i+1), h[0]->GetBinLowEdge(i+2)), nBins[i], -0.15, +0.15);
    hPtDev[i]->SetStats(0);
    if (i == 18 || i == 17) {
      fgaus[i] = new TF1(Form("fgaus%i", i), "gaus", -0.15, +0.15);
    } else {
      fgaus[i] = new TF1(Form("fgaus%i", i), "gaus", -0.15, +0.15);
    }
  }

  for (int i = 0; i < trials; i++){
    for (int j = 0; j < bins; j++){
      if(h[i]->GetBinContent(j+1) != 0) {
        hPtDev[j]->Fill(h[i]->GetBinContent(j+1), 1);
      }
    }
  }

  for (int i = 0; i < bins; i++){
    fgaus[i]->SetParameter(0, hPtDev[i]->GetMaximum());
    fgaus[i]->SetParameter(1, 0);
    fgaus[i]->SetParLimits(0, 0.8*hPtDev[i]->GetMaximum(), 1.5*hPtDev[i]->GetMaximum());
    std::cout << "maximum value in bin: " << hPtDev[i]->GetMaximum() << std::endl;
    hPtDev[i]->Fit(fgaus[i], "LR");
  }

  TH1F *hSystMultiTrial = (TH1F *)h[0]->Clone("hSystMultiTrial");
  hSystMultiTrial->Reset();
  for (int i = 0; i < bins; i++) {
    hSystMultiTrial->SetBinContent(i + 1, fgaus[i]->GetParameter(2));
    hSystMultiTrial->SetBinError(i + 1, 0);
  }

  TFile* Write = new TFile("../sourcesOfSyst/multiTrialSyst-Xi.root", "RECREATE");
  for (int i = 0; i < bins; i++){
    hPtDev[i]->Write();
    fgaus[i]->Write();
  }

  hSystMultiTrial->SetTitle("");
  hSystMultiTrial->GetYaxis()->SetRangeUser(0, 0.2);
  hSystMultiTrial->GetYaxis()->SetTitle("Syst. uncertainty");
  hSystMultiTrial->Write();

  TLatex *ltx = new TLatex();
  ltx->SetTextSize(0.05);
  ltx->SetTextColor(kRed);
  gROOT->SetBatch(kFALSE);

  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1200);
  c2->SetName("MultiTrialSyst_1");
  c2->Divide(3, 4);
  StyleCanvas(c2, 0.15, 0.05, 0.05, 0.15);
  for (Int_t i = 1; i <= 12; i++) {
    c2->cd(i)->SetMargin(0.15, 0.05, 0.15, 0.1);
  }

  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 1200);
  c3->SetName("MultiTrialSyst_2");
  c3->Divide(3, 4);
  StyleCanvas(c3, 0.15, 0.05, 0.05, 0.15);
  for (Int_t i = 1; i <= 12; i++) {
    c3->cd(i)->SetMargin(0.15, 0.05, 0.15, 0.1);
  }

  for (int i = 0; i < bins; i++){
    if(i > 11){
      c3->cd(i-12+1);
      //c3->cd(i-12+1)->SetBottomMargin(0.15);
    }
    else{
      c2->cd(i+1);
      //c2->cd(i+1)->SetBottomMargin(0.15);
    }
    StyleHistoLight(hPtDev[i]);
    hPtDev[i]->GetYaxis()->SetRangeUser(0., hPtDev[i]->GetMaximum() * 1.3);
    hPtDev[i]->GetXaxis()->SetTitleSize(0.06);
    // hPtDev[i]->GetXaxis()->SetTitleOffset(1.);
    hPtDev[i]->SetMarkerStyle(kFullCircle);
    hPtDev[i]->Draw("EP SAME");
    fgaus[i]->Draw("same");

    ltx->DrawLatexNDC(0.25, 0.8, Form("#mu = %.3f", fgaus[i]->GetParameter(1)));
    ltx->DrawLatexNDC(0.25, 0.7, Form("#sigma = %.3f", fgaus[i]->GetParameter(2)));
    ltx->DrawLatexNDC(0.25, 0.6, Form("#chi^{2}/ndf = %.1f", fgaus[i]->GetChisquare() / fgaus[i]->GetNDF()));
  }
  c2->Draw();
  c3->Draw();
  c2->Write();
  c3->Write();
}

//-------------------------------------------------------------------------------
//------------------------- DEFINE FUNCTIONS ------------------------------------
//-------------------------------------------------------------------------------

//---------------------------------------------------------------
TH1F* makeSystPlots(int num = 1)
{
  TFile *fdef = TFile::Open("data/yieldsOutEffCorr/yield_XiPm_MB_inel0_DEF.root");
  TFile *fvaried = TFile::Open(Form("data/yieldsOutEffCorr/yield_XiPm_MB_inel0_%i.root", num));

  TH1F *hVariedCut = (TH1F *)fvaried->Get("effCorrYield/Yield_XiPm");
  hVariedCut->SetName(Form("hVarCutSet%i", num));
  TH1F *hDefault = (TH1F *)fdef->Get("effCorrYield/Yield_XiPm");
  hDefault->SetName(Form("hDefault"));

  TH1F *hDev = (TH1F *)hDefault->Clone("hDev");
  hDev->Reset();

  DivideAndComputeRogerBarlow(hVariedCut, hDefault);

  for (int i = 1; i <= hDev->GetNbinsX(); i++)
  {
    double dev = 1 - hVariedCut->GetBinContent(i);
    double err = hVariedCut->GetBinError(i);
    std::cout << "dev: " << dev << "err: " << err << " weight: " << PassRogerBarlowCriterion(1, dev, err) << std::endl;
    hDev->SetBinContent(i, PassRogerBarlowCriterion(1, dev, err));
  }

  hDev->GetYaxis()->SetRangeUser(-0.5, 0.5);
  //Return Max Dev Histo
  //fvaried->Close();
  //fdef->Close();
  return hDev;
}

//----------------------------------------------------------------------------------------------------
void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 ){
  //Use Roger Barlow "sigma_{delta}" as errors for ratios
  Double_t lh1NBins = h1->GetNbinsX();
  Double_t lh2NBins = h2->GetNbinsX();

  if( lh1NBins != lh2NBins ){
    cout<<"Problem! Number of bins doesn't match! "<<endl;
    return;
  }

  Double_t lSigmaDelta[100];
  for( Int_t i=1; i<h1->GetNbinsX()+1; i++){
    //Computation of roger barlow sigma_{delta}
    lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(h1->GetBinError(i),2) - TMath::Power(h2->GetBinError(i),2) ) );
    //Computation of relationship to h2 for plotting in ratio plot
    if ( h2->GetBinContent(i) > 1e-12 ){
      lSigmaDelta[i] /= h2->GetBinContent(i);
    }else{
      lSigmaDelta[i] = 0;
    }
  }
  //Regular Division
  h1->Divide(h2);
  //Replace Erorrs
  for(Int_t i=1; i<h1->GetNbinsX()+1; i++){
    h1->SetBinError(i, lSigmaDelta[i]);
  }
}

//----------------------------------------------------------------------------------------------------
double PassRogerBarlowCriterion(int nsigmas, Double_t dev, Double_t RBsigma){
  if (TMath::Abs(dev)>(nsigmas*RBsigma)) {return dev;}
  else {return 0.;}
}
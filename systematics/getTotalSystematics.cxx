#include "../help.h"
#include "../effHelp.h"

void prephisto(TH1D* h, int col);

void getTotalSystematics(const Int_t nParticle = 2) // 0-2 : xi, 3-5 : omega
{
  // Setup pt binning
  Int_t numPt = 0;
  Int_t partType;
  if (nParticle <= 2) {
    numPt = numPtXi; // Xi
    partType = 0;
  } else {
    numPt = numPtOmega; // Omega
    partType = 1;
  }

  Double_t* binpt = new Double_t[numPt + 1];

  if (nParticle <= 2) {
    for (int i = 0; i < (numPt + 1); ++i) {
        binpt[i] = binptXi[i]; // Xi
      }
  } else {
    for (int i = 0; i < (numPt + 1); ++i) {
        binpt[i] = binptOmega[i]; // Omega
      }
  }

  TH1D *hPt = new TH1D("hPt", "hPt", numPt, binpt);

  // Cut variation (multi-trial)
  TFile *fMultiTrial = TFile::Open("./sourcesOfSyst/multiTrialSyst-Xi.root");
  TH1D *hSystMultiTrial = (TH1D *)fMultiTrial->Get("hSystMultiTrial");
  prephisto(hSystMultiTrial, kGreen+1);

  // Material budget
  TFile *matbudg = TFile::Open("./sourcesOfSyst/matBudgetSystRun2-Xi.root");
  TF1 *fmatbudg = (TF1 *)matbudg->Get("funcFit");
  TH1D *hSystMatBudg = (TH1D *)hPt->Clone("hSystMatBudg");
  hSystMatBudg->Reset();
  for (int i = 1; i <= hPt->GetNbinsX(); i++)
  {
    Double_t MaterialBudget = TMath::Abs((1 - fmatbudg->Eval(hPt->GetBinCenter(i)))) / 2;
    hSystMatBudg->SetBinContent(i, MaterialBudget);
  }
  prephisto(hSystMatBudg, kBlue);

  // Efficiency multiplicity dependence
  TH1D *hSystEffMultDep = (TH1D *)hPt->Clone("hSystEffMultDep");
  for (int i = 1; i <= hSystEffMultDep->GetNbinsX(); i++)
  {
    hSystEffMultDep->SetBinContent(i, 0.02);
  }
  prephisto(hSystEffMultDep, kRed);

  // // OOB pileup
  // TFile *fOOB = TFile::Open("OOBSystOmega.root");
  // TH1D *hSystOOB = (TH1D *)fOOB->Get("hSystOOB");
  // prephisto(hSystOOB, kMagenta);

  // Bkg variation
  TFile *fBkg = TFile::Open("./sourcesOfSyst/SystBKG-Xi.root");
  TH1D *hBkgSyst = (TH1D *)fBkg->Get("hSystBKG");
  hBkgSyst->SetLineStyle(1);
  prephisto(hBkgSyst, kOrange+1);

  //Total
  TH1D* hSystTotal = (TH1D*)hPt->Clone("hSystTotal");
  hSystTotal->SetTitle("Total systematics");
  for (int i = 1; i <= hSystTotal->GetNbinsX(); i++)
  {
      double syst = 0;
      syst += hSystMatBudg->GetBinContent(hSystMatBudg->FindBin(hSystTotal->GetBinCenter(i))) * hSystMatBudg->GetBinContent(hSystMatBudg->FindBin(hSystTotal->GetBinCenter(i)));
      syst += hSystEffMultDep->GetBinContent(hSystEffMultDep->FindBin(hSystTotal->GetBinCenter(i))) * hSystEffMultDep->GetBinContent(hSystEffMultDep->FindBin(hSystTotal->GetBinCenter(i)));
      syst += hSystMultiTrial->GetBinContent(hSystMultiTrial->FindBin(hSystTotal->GetBinCenter(i))) * hSystMultiTrial->GetBinContent(hSystMultiTrial->FindBin(hSystTotal->GetBinCenter(i)));
      //syst += hSystOOB->GetBinContent(hSystOOB->FindBin(hSystTotal->GetBinCenter(i))) * hSystOOB->GetBinContent(hSystOOB->FindBin(hSystTotal->GetBinCenter(i)));
      syst += hBkgSyst->GetBinContent(hBkgSyst->FindBin(hSystTotal->GetBinCenter(i))) * hBkgSyst->GetBinContent(hBkgSyst->FindBin(hSystTotal->GetBinCenter(i)));

      syst = TMath::Sqrt(syst);
      hSystTotal->SetBinContent(i, syst);
  }
  prephisto(hSystTotal, kBlack);

  //Canvas
  // TCanvas *cTotalSyst = new TCanvas("cTotalSyst", " ", 1400, 1000);
  // cTotalSyst->SetRightMargin(0.15);
  // cTotalSyst->SetLeftMargin(0.15);
  // cTotalSyst->SetBottomMargin(0.15);

  TCanvas *cTotalSyst = new TCanvas("cTotalSyst","cTotalSyst", 800, 600);
  cTotalSyst->cd();
  StyleCanvas(cTotalSyst, 0.14, 0.05, 0.11, 0.15);
  cTotalSyst->SetGridy();
  cTotalSyst->SetGridx();
  StyleHisto(hSystTotal, 0., 0.075, color[0], MarkerMult[0], sPt, "Systematic uncertainties", "", 0, 0, 0, 1.0, 1.25, SizeMult[0], 0.04, 0.04);
  hSystTotal->SetStats(kFALSE);
  hSystTotal->SetLineWidth(2);

  // hSystTotal->SetYTitle("systematic uncertainty");
  // hSystTotal->SetXTitle("p_{T} (GeV/c)");
  // hSystTotal->SetTitle("Systematic uncertainties");
  // hSystTotal->GetYaxis()->SetTitleSize(0.05);
  // hSystTotal->GetYaxis()->SetTitleOffset(1.4);
  // hSystTotal->GetXaxis()->SetTitleSize(0.05);
  // hSystTotal->GetXaxis()->SetTitleOffset(1.);
  // hSystTotal->GetYaxis()->SetRangeUser(0., 0.15);

  hSystTotal->Draw();
  hSystMatBudg->Draw("same");
  hSystEffMultDep->Draw("same");
  hSystMultiTrial->Draw("same");
  //hSystOOB->Draw("same");
  hBkgSyst->Draw("same");

  TLegend *leg = new TLegend(0.47, 0.57, 0.95, 0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hSystTotal, "Total", "l");
  leg->AddEntry(hSystMatBudg, "Material budget", "l");
  leg->AddEntry(hSystEffMultDep, "Efficiency mult. dep.", "l");
  //leg->AddEntry(hSystOOB, "OOB pileup", "l");
  leg->AddEntry(hSystMultiTrial, "Cut variation (multi-trial)", "l");
  leg->AddEntry(hBkgSyst, "Background", "l");

  leg->Draw("same");

  TString outputfilename = Form("./sourcesOfSyst/totalSystematics.root");
  TFile *CreateFile = new TFile(outputfilename, "recreate");
  hSystTotal->Write();
  hSystMatBudg->Write();
  hSystEffMultDep->Write();
  //hSystOOB->Write();
  hSystMultiTrial->Write();
  hBkgSyst->Write();
  CreateFile->Close();
}

void prephisto(TH1D* h, int col){
  h->SetLineColor(col);
  h->SetMarkerColor(col);
  h->SetLineWidth(2);
  h->SetStats(0);
}
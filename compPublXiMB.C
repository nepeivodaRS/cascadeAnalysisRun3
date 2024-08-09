#include "help.h"

void StyleHisto(TH1F *histo, Double_t Low, Double_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, 
                Bool_t XRange, Double_t XLow, Double_t XUp, Double_t xOffset, Double_t yOffset, Double_t mSize, 
                Double_t xTitleSize, Double_t yTitleSize, Double_t xLabelSize,  Double_t yLabelSize, Double_t yLabelOffset)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)
    histo->GetXaxis()->SetRangeUser(XLow, XUp);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetLabelSize(xLabelSize);
  histo->GetXaxis()->SetTitleSize(xTitleSize);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(yTitleSize);
  histo->GetYaxis()->SetLabelSize(yLabelSize);
  histo->GetYaxis()->SetTitleOffset(yOffset);
  histo->SetTitle(title);
  histo->GetYaxis()->SetLabelOffset(yLabelOffset);
}

void compPublXiMB(const Int_t nParticle = 2, // 0-2 : xi, 3-5 : omega
                 const Int_t inel = 0, // inel > N (0/1)
                 const TString workingDir = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM",
                 const TString postFix = "")
{
  // Start of Code
  std::cout << "\x1B[1;33m"; // Set text color to yellow
  std::cout << "\n************ Starting Making Comparison with Published Results ************\n";
  std::cout << "\x1B[0m"; // Reset text color

  gStyle->SetOptStat(0);

  TString outputfilePath = workingDir + "/comparisons/" + "mbXi" + ".root";
  TFile *outputfile = new TFile(outputfilePath, "RECREATE");

  // HEP Data File 
  TString inputHEPDataPath = workingDir + "/data/published/xi.root";
  TFile* fileHEPDataIn = TFile::Open(inputHEPDataPath);
  if (!fileHEPDataIn || fileHEPDataIn->IsZombie()) {
      std::cerr << "Error opening HEP data file!" << std::endl;
      return;
  }
  
  TDirectory* yieldHEPDir = fileHEPDataIn->GetDirectory("Table 3");
  if (!yieldHEPDir)
  {
    std::cerr << "`Table 3` directory in HEP data file is not found!" << std::endl;
    return;
  }

  TH1F* hHEPYield = (TH1F *)yieldHEPDir->Get("Hist1D_y11");
  TH1F* hHEPYieldSyst = (TH1F *)hHEPYield->Clone("Hist1D_y11_Syst");
  TH1F* hHEPYielde1 = (TH1F *)yieldHEPDir->Get("Hist1D_y11_e1");
  TH1F* hHEPYielde2 = (TH1F *)yieldHEPDir->Get("Hist1D_y11_e2");
  TH1F* hHEPYielde3 = (TH1F *)yieldHEPDir->Get("Hist1D_y11_e3");
  if (!hHEPYield || !hHEPYielde1 || !hHEPYielde2)
  {
    std::cerr << "Histogram `Hist1D_y11` is not found!" << std::endl;
    return;
  }

  Double_t e1 = 0;
  Double_t e2 = 0;
  Double_t e3 = 0;

  for (Int_t iBin = 1; iBin <= hHEPYield->GetNbinsX(); iBin++) {
    e1 = hHEPYielde1->GetBinContent(iBin);
    e2 = hHEPYielde2->GetBinContent(iBin);
    e3 = hHEPYielde3->GetBinContent(iBin); // uncorrelated
    //hHEPYield->SetBinError(iBin, sqrt(pow(e1, 2) + pow(e2, 2)) + e3);
    hHEPYield->SetBinError(iBin, e1); // only statistical uncert.
    hHEPYieldSyst->SetBinError(iBin, e2); // only syst. uncert.
  }

  // Files with yields in all mult. classes + MB
  TFile* fileDataIn[numMult + 1];
  // MB
  fileDataIn[0] = TFile::Open(workingDir + "/yieldsOutEffCorr" +  "/yield_" + particleNames[nParticle] + "_MB_inel" + inel + postFix + ".root");
  if (!fileDataIn[0] || fileDataIn[0]->IsZombie()) {
    std::cerr << "Error opening input data file for MB!" << std::endl;
    return;
  } else {
    cout << "file for MB yield is opened"<< std::endl;
  }
  // in mult. classes
  for (Int_t iFile = 1; iFile < numMult + 1; iFile++) {
    TString fileInPath = workingDir + "/yieldsOutEffCorr" + "/yield_" + particleNames[nParticle] + "_" + multiplicityPerc[iFile - 1] + "-" + multiplicityPerc[iFile] + "_inel" + inel + postFix  + ".root";
    fileDataIn[iFile] = TFile::Open(fileInPath);
    if (!fileDataIn[iFile] || fileDataIn[iFile]->IsZombie()) {
      std::cerr << "Error opening input data file for mult. class: " << multiplicityPerc[iFile - 1] <<  " - " << multiplicityPerc[iFile] << std::endl;
      return;
    } else {
      cout << "file for mult. class: " << multiplicityPerc[iFile - 1] <<  " - " << multiplicityPerc[iFile] << " is opened"<< std::endl;
    }
  }

  TString fileInPathSyst = workingDir + "/systematics/sourcesOfSyst/totalSystematics.root";
  TFile* fileSystin = TFile::Open(fileInPathSyst);
  if (!fileSystin || fileSystin->IsZombie()) {
      std::cerr << "Error opening `fileSystin` data file!" << std::endl;
      return;
  }
  TH1F* hSystUncert = (TH1F *)fileSystin->Get("hSystTotal");

  TH1F* hYield[numMult + 1];
  TH1F* hYieldSyst[numMult + 1];
  TDirectory* yieldDir[numMult + 1];

  // Get yields
  for (Int_t iFile = 0; iFile < numMult + 1; iFile++) {
    yieldDir[iFile] = fileDataIn[iFile]->GetDirectory("effCorrYield");
    if (!yieldDir[iFile])
    {
      std::cerr << "`effCorrYield` directory is not found!" << std::endl;
      return;
    }

    hYield[iFile] = (TH1F *)yieldDir[iFile]->Get("Yield_" + particleNames[nParticle]);
    // Syst. histo
    hYieldSyst[iFile] = (TH1F*)hYield[iFile]->Clone(Form("YieldSys_%i_",iFile) + particleNames[nParticle]); 
    for (Int_t bin = 1; bin <= hYield[iFile]->GetNbinsX(); bin++) {
      hYieldSyst[iFile]->SetBinError(bin, hSystUncert->GetBinContent(hSystUncert->FindBin(hYield[iFile]->GetBinCenter(bin)))*hYield[iFile]->GetBinContent(bin));
    }
  }

  // Scale x 2 to comapre to + and - separately
  if(nParticle != 2) {
    hYield[0]->Scale(2.);
    hYieldSyst[0]->Scale(2.);
  }

  TH1F* hYieldsDenom = (TH1F *)hHEPYield->Clone("YieldHEPDataClone");
  TH1F* hYieldsDenomSyst = (TH1F *)hHEPYieldSyst->Clone("YieldHEPDataSystClone");
  TH1F* hYieldsRatio = (TH1F *)hYield[0]->Clone("YieldClone");
  TH1F* hYieldsRatioSyst = (TH1F *)hYieldSyst[0]->Clone("YieldSystClone");
  hYieldsRatio->Divide(hYieldsDenom);
  ErrRatioCorr(hYield[0], hYieldsDenom, hYieldsRatio, 0);
  hYieldsRatioSyst->Divide(hYieldsDenomSyst);
  ErrRatioCorr(hYieldSyst[0], hYieldsDenomSyst, hYieldsRatioSyst, 0);

  Double_t yieldYLow[2] = {0.1*1e-5, 0.2*1e-8};
  Double_t yieldYUp[2] = {0.04, 3*1e4};

  Int_t partType = 0;

  TCanvas* canvasComparison = new TCanvas("canvasComparison_Xi", "canvasComparison_Xi", 0, 70, 620, 850);
  StyleCanvas(canvasComparison, 0.15, 0.05, 0.05, 0.15);
  // Yield Pads
  TPad *padYieldUp = new TPad("padYieldUp_" + particleNames[2], "padYieldUp_" + particleNames[2], 0, 0.36, 1, 1);
  TPad *padYieldLow = new TPad("padYieldLow_" + particleNames[2], "padYieldLow_" + particleNames[2], 0, 0, 1, 0.36);
  StylePad(padYieldUp, 0.15, 0.05, 0.05, 0.);
  StylePad(padYieldLow, 0.15, 0.05, 0.02, 0.2);

  canvasComparison->cd();
  padYieldUp->Draw();
  padYieldLow->Draw();

  // Yield Legend
  TLegend *legYield = new TLegend(0.186, 0.0625, 0.893, 0.3125);
  legYield->SetHeader("Comparison with Published MB Run2 Data");
  //legYield->SetNColumns(3);
  TLegendEntry *legYieldHeader = (TLegendEntry *)legYield->GetListOfPrimitives()->First();
  legYieldHeader->SetTextSize(0.04);
  //legYield->SetTextAlign(13); // adjust text wrt markers
  StyleLegend(legYield, 0.0, 0.0);
  legYield->SetTextSize(0.035);

  TLegend *legendTitle = new TLegend(0.609, 0.693913, 0.908, 0.893);
  StyleLegend(legendTitle, 0.0, 0.0);
  legendTitle->SetTextAlign(33);
  legendTitle->SetTextSize(0.05);
  legendTitle->SetTextFont(42);
  legendTitle->SetLineColorAlpha(0.,0.);
  legendTitle->SetFillColorAlpha(0.,0.);
  legendTitle->AddEntry("", "#bf{ALICE Work In Progress}", "");
  if(nParticle != 2) {
    legendTitle->AddEntry("", particleSymnbols[nParticle] + "x2" + ", |y| < 0.5", "");
  } else {
    legendTitle->AddEntry("", particleSymnbols[2] + ", |y| < 0.5", "");
  }
  //legendTitle->AddEntry("", "Statistical uncert. only", "");

  Int_t compStyle = 3;
  StyleHisto(hHEPYield, yieldYLow[partType], yieldYUp[partType], color[compStyle], MarkerMult[compStyle], "", hHEPYield->GetYaxis()->GetTitle(), "", 0, 0, 0, 1.5, 1.0, SizeMult[compStyle], 0.0, 0.05, 0.0, 0.035, 0.005);
  StyleHisto(hHEPYieldSyst, yieldYLow[partType], yieldYUp[partType], color[compStyle], MarkerMult[compStyle], "", hHEPYieldSyst->GetYaxis()->GetTitle(), "", 0, 0, 0, 1.5, 1.0, SizeMult[compStyle], 0.0, 0.05, 0.0, 0.035, 0.005);

  // Dummy histograms to extend x-axis
  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 8.5);
  // Up
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, 1e-12);
  canvasComparison->cd();
  padYieldUp->cd();
  padYieldUp->SetLogy();
  StyleHisto(hDummy, yieldYLow[partType], yieldYUp[partType], 1, 1, sPt, sdNdPtdY, "", 0, 0, 0, 1.5, 1.0, 0, 0.0, 0.05, 0.0, 0.035, 0.005);
  SetTickLength(hDummy, 0.025, 0.03);
  TAxis *axisYieldDummy = hDummy->GetYaxis();
  axisYieldDummy->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
  hDummy->GetXaxis()->SetRangeUser(0, hYield[0]->GetXaxis()->GetBinUpEdge(hYield[0]->GetNbinsX()) + 0.5);
  hDummy->Draw("same");
  // Low
  TH1F *hDummyLow = new TH1F("hDummyLow", "hDummyLow", 10000, 0, 8.5);
  for (Int_t i = 1; i <= hDummyLow->GetNbinsX(); i++)
    hDummyLow->SetBinContent(i, 1);
  padYieldLow->cd();
  StyleHisto(hDummyLow, 0.4 , 1.6, 1, 1, sPt, "Ratio to Published", "", 0, 0, 0, 1.0, 0.7, 0, 0.08, 0.08, 0.08, 0.07, 0.01);
  SetTickLength(hDummyLow, 0.025, 0.03);
  hDummyLow->GetXaxis()->SetRangeUser(0, hYieldsRatio->GetXaxis()->GetBinUpEdge(hYieldsRatio->GetNbinsX()) + 0.5);
  hDummyLow->GetYaxis()->CenterTitle();
  TAxis *axisDummyLow = hDummyLow->GetYaxis();
  axisDummyLow->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
  axisDummyLow->ChangeLabel(-1, -1, -1, -1, -1, -1, " ");
  hDummyLow->Draw("same");

  // Plot yields
  padYieldUp->cd();
  hYield[0]->Draw("same ex0");
  hYieldSyst[0]->SetFillStyle(0);
  hYieldSyst[0]->Draw("same e2");
  legYield->AddEntry(hYield[0], "pp, #sqrt{#it{s}} = 13.6 TeV", "pef");
  hHEPYield->Draw("same ex0");
  hHEPYieldSyst->SetFillStyle(0);
  hHEPYieldSyst->Draw("same e2");
  legYield->AddEntry(hHEPYield, "Eur.Phys.J.C.80 (2020) 167, pp, #sqrt{#it{s}} = 13 TeV", "pef");

  // Plot ratios
  padYieldLow->cd();
  StyleHisto(hYieldsRatio, 0.5 , 1.1, color[0], MarkerMult[0], sPt, "Ratio to Published", "", 0, 0, 0, 1.0, 0.7, SizeMult[0], 0.08, 0.08, 0.08, 0.08, 0.005);
  hYieldsRatio->Draw("same ex0");
  hYieldsRatioSyst->SetFillStyle(0);
  hYieldsRatioSyst->Draw("same e2");

  padYieldUp->cd();
  legendTitle->Draw();
  legYield->Draw();
  gPad->Update();
  padYieldLow->cd();
  DrawHorLine(hYieldsRatio->GetXaxis()->GetBinUpEdge(hYieldsRatio->GetNbinsX()) + 0.5, 1.0);
  canvasComparison->Update();

  outputfile->cd();
  canvasComparison->Write();

  // Compare Systematics //
  TLegend *legendTitleSystComp = new TLegend(0.589, 0.625, 0.888, 0.824);
  StyleLegend(legendTitleSystComp, 0.0, 0.0);
  legendTitleSystComp->SetTextAlign(33);
  legendTitleSystComp->SetTextSize(0.05);
  legendTitleSystComp->SetTextFont(42);
  legendTitleSystComp->SetLineColorAlpha(0.,0.);
  legendTitleSystComp->SetFillColorAlpha(0.,0.);
  legendTitleSystComp->AddEntry("", "#bf{ALICE Work In Progress}", "");
  legendTitleSystComp->AddEntry("", particleSymnbols[2] + ", |y| < 0.5", "");

  TLegend *legendsystComp = new TLegend(0.182, 0.180, 0.383, 0.358);
  legendsystComp->SetFillStyle(0);
  legendsystComp->SetTextSize(0.03);

  TCanvas *canvasSystComp = new TCanvas("canvasSystComp","canvasSystComp", 0,66,857,708);
  canvasSystComp->cd();
  canvasSystComp->Draw();
  StyleCanvas(canvasSystComp, 0.14, 0.05, 0.11, 0.15);
  TH1F* hRun3Syst = (TH1F* )hYield[0]->Clone("hRun3Syst");
  TH1F* hRun2Syst = (TH1F* )hHEPYieldSyst->Clone("hRun2Syst");
  for(Int_t bin = 1; bin <= hYield[0]->GetNbinsX(); bin++) {
    hRun3Syst->SetBinContent(bin, hYieldSyst[0]->GetBinError(bin) / hYield[0]->GetBinContent(bin) * 100);
    hRun2Syst->SetBinContent(bin, hHEPYieldSyst->GetBinError(bin) / hHEPYield->GetBinContent(bin) * 100);
    hRun3Syst->SetBinError(bin, 0.);
    hRun2Syst->SetBinError(bin, 0.);
  }
  legendsystComp->AddEntry(hRun3Syst, "pp, #sqrt{#it{s}} = 13.6 TeV", "l");
  legendsystComp->AddEntry(hRun2Syst, "Eur.Phys.J.C.80 (2020) 167, pp, #sqrt{#it{s}} = 13 TeV", "l");
  StyleHisto(hRun3Syst, 0., 9, color[0], MarkerMult[0], sPt, "Systematic uncertainty (%)", "", 0, 0, 0, 1.0, 1.25, SizeMult[0], 0.04, 0.04, 0.04, 0.04, 0.01);
  StyleHisto(hRun2Syst, 0., 9, color[compStyle], MarkerMult[compStyle], sPt, "Systematic uncertainty (%)", "", 0, 0, 0, 1.0, 1.25, SizeMult[compStyle], 0.04, 0.04, 0.04, 0.04, 0.01);
  hRun3Syst->GetYaxis()->SetMaxDigits(1);
  //hRun3Syst->GetYaxis()->SetDecimals(kTRUE);
  //gPad->SetLogy();
  hRun3Syst->Draw("same");
  hRun3Syst->SetLineWidth(3);
  hRun2Syst->Draw("same");
  hRun2Syst->SetLineWidth(3);
  legendTitleSystComp->Draw("same");
  legendsystComp->Draw("same");
  outputfile->cd();
  canvasSystComp->Write();

  // End of Code
  std::cout << "\x1B[1;32m"; // Set text color to green
  std::cout << "\n************* Comparison Finished Successfully! *************\n";
  std::cout << "\x1B[0m"; // Reset text color

}
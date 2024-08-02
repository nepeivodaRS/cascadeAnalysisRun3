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

void sigmaInMult(const Int_t nParticle = 2, // 0-2 : xi, 3-5 : omega
                 const Int_t inel = 0, // inel > N (0/1)
                 const TString workingDir = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM")
{
  //gROOT->SetBatch(kTRUE);
  // Start of Code
  std::cout << "\x1B[1;33m"; // Set text color to yellow
  std::cout << "\n************ Starting Comparison of Sigma in Multiplicity Classes Plot ************\n";
  std::cout << "\x1B[0m"; // Reset text color

  gStyle->SetOptStat(0);

  TString outputfilePath = workingDir + "/sigmaInMult/" + "sigmaInMult_" + particleNames[nParticle] + "_inel" + inel + ".root";
  TFile *outputfile = new TFile(outputfilePath, "RECREATE");

  // Files with means in all mult. classes + MB
  TFile* fileDataIn[numMult + 1];
  // MB
  fileDataIn[0] = TFile::Open(workingDir + "/yieldsOut" +  "/yield_" + particleNames[nParticle] + "_MB_inel" + inel + ".root");
  if (!fileDataIn[0] || fileDataIn[0]->IsZombie()) {
    std::cerr << "Error opening input data file for MB!" << std::endl;
    return;
  } else {
    cout << "file for MB mean is opened"<< std::endl;
  }
  // in mult. classes
  for (Int_t iFile = 1; iFile < numMult + 1; iFile++) {
    TString fileInPath = workingDir + "/yieldsOut" + "/yield_" + particleNames[nParticle] + "_" + multiplicityPerc[iFile - 1] + "-" + multiplicityPerc[iFile] + "_inel" + inel + ".root";
    fileDataIn[iFile] = TFile::Open(fileInPath);
    if (!fileDataIn[iFile] || fileDataIn[iFile]->IsZombie()) {
      std::cerr << "Error opening input data file for mult. class: " << multiplicityPerc[iFile - 1] <<  " - " << multiplicityPerc[iFile] << std::endl;
      return;
    } else {
      cout << "file for mult. class: " << multiplicityPerc[iFile - 1] <<  " - " << multiplicityPerc[iFile] << " is opened"<< std::endl;
    }
  }

  TH1F* hSigma[numMult + 1];
  TH1F* hSigmaRatio[numMult + 1];
  TH1F* hSigmaDenom;
  TDirectory* sigmaDir[numMult + 1];

  // Get sigmas
  for (Int_t iFile = 0; iFile < numMult + 1; iFile++) {
    sigmaDir[iFile] = fileDataIn[iFile]->GetDirectory("fitParams");
    if (!sigmaDir[iFile])
    {
      std::cerr << "`fitParams` directory is not found!" << std::endl;
      return;
    }

    hSigma[iFile] = (TH1F *)sigmaDir[iFile]->Get("Sigma_" + particleNames[nParticle]);
    hSigmaRatio[iFile] = (TH1F *)hSigma[iFile]->Clone("SigmaClone_" + particleNames[nParticle]);
    if (iFile == 0) {
      hSigmaDenom = (TH1F *)hSigma[iFile]->Clone("SigmaCloneForDenom_" + particleNames[nParticle]);
    }

    hSigmaRatio[iFile]->Divide(hSigmaDenom);
    if (iFile != 0) {
      ErrRatioCorr(hSigma[iFile], hSigmaDenom, hSigmaRatio[iFile], 0);
    }
  }

  // Sigma canvas
  TCanvas* canvasSigma = new TCanvas("sigmaSummary_" + particleNames[nParticle], "sigmaSummary_" + particleNames[nParticle], 0, 70, 620, 850);
  StyleCanvas(canvasSigma, 0.15, 0.05, 0.05, 0.15);
  // Sigma Pads
  TPad *padSigmaUp = new TPad("padSigmaUp_" + particleNames[nParticle], "padSigmaUp_" + particleNames[nParticle], 0, 0.36, 1, 1);
  TPad *padSigmaLow = new TPad("padSigmaLow_" + particleNames[nParticle], "padSigmaLow_" + particleNames[nParticle], 0, 0, 1, 0.36);
  StylePad(padSigmaUp, 0.15, 0.05, 0.05, 0.);
  StylePad(padSigmaLow, 0.15, 0.05, 0.02, 0.2);
  canvasSigma->cd();
  padSigmaUp->Draw();
  padSigmaLow->Draw();

  // Sigma Legend
  TLegend *legSigma = new TLegend(0.186, 0.0625, 0.893, 0.3125);
  legSigma->SetHeader("FT0M Multiplicity Percentile");
  legSigma->SetNColumns(3);
  TLegendEntry *legSigmaHeader = (TLegendEntry *)legSigma->GetListOfPrimitives()->First();
  legSigmaHeader->SetTextSize(0.04);
  //legSigma->SetTextAlign(13); // adjust text wrt markers
  StyleLegend(legSigma, 0.0, 0.0);

  TLegend *legendTitle = new TLegend(0.609, 0.693913, 0.908, 0.893);
  StyleLegend(legendTitle, 0.0, 0.0);
  legendTitle->SetTextAlign(33);
  legendTitle->SetTextSize(0.05);
  legendTitle->SetTextFont(42);
  legendTitle->SetLineColorAlpha(0.,0.);
  legendTitle->SetFillColorAlpha(0.,0.);
  legendTitle->AddEntry("", "#bf{ALICE Work In Progress}", "");
  legendTitle->AddEntry("", "pp, #sqrt{#it{s}} = 13.6 TeV", "");
  legendTitle->AddEntry("", particleSymnbols[nParticle] + ", |y| < 0.5", "");
  
  TString SmoltBis[numMult + 1];
  for (Int_t m = 0; m < numMult + 1; m++) {
    if (m == 0) { 
      SmoltBis[m] += Form("%.1f#minus%.1f%s", 0.0, 100.0, "%");
    } else {
      SmoltBis[m] += Form("%.1f#minus%.1f%s", multiplicityPerc[m-1], multiplicityPerc[m], "%");
    }
  }

  Double_t sigmaYLow[2] = {0.0005, 0.0005};
  Double_t sigmaYUp[2] = {0.004, 0.01};

  Int_t partType;
  if (nParticle <= 2) {
    partType = 0; // Xi
  } else {
    partType = 1; // Omega
  }

  // Dummy histograms to extend x-axis
  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 8.5);
  // Up
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, 1e-12);
  canvasSigma->cd();
  padSigmaUp->cd();
  StyleHisto(hDummy, sigmaYLow[partType], sigmaYUp[partType], 1, 1, "#it{p}_{T} (GeV/#it{c})", "Sigma (GeV/#it{c}^{2})", "", 0, 0, 0, 1.5, 1.0, 0, 0.0, 0.05, 0.0, 0.035, 0.005);
  SetTickLength(hDummy, 0.025, 0.03);
  TAxis *axisMeanDummy = hDummy->GetYaxis();
  axisMeanDummy->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
  hDummy->GetXaxis()->SetRangeUser(0, hSigma[0]->GetXaxis()->GetBinUpEdge(hSigma[0]->GetNbinsX()) + 0.5);
  hDummy->GetYaxis()->SetMaxDigits(4);
  hDummy->GetYaxis()->SetDecimals(kTRUE);
  hDummy->Draw("same");
  // Low
  TH1F *hDummyLow = new TH1F("hDummyLow", "hDummyLow", 10000, 0, 8.5);
  for (Int_t i = 1; i <= hDummyLow->GetNbinsX(); i++)
    hDummyLow->SetBinContent(i, 1);
  padSigmaLow->cd();
  StyleHisto(hDummyLow, 0.7, 1.3, 1, 1, "#it{p}_{T} (GeV/#it{c})", "Ratio to 0-100%", "", 0, 0, 0, 1.0, 0.7, 0, 0.08, 0.08, 0.08, 0.07, 0.01);
  SetTickLength(hDummyLow, 0.025, 0.03);
  hDummyLow->GetYaxis()->SetMaxDigits(3);
  hDummyLow->GetYaxis()->SetDecimals(kTRUE);
  hDummyLow->GetXaxis()->SetRangeUser(0, hSigmaRatio[0]->GetXaxis()->GetBinUpEdge(hSigmaRatio[0]->GetNbinsX()) + 0.5);
  hDummyLow->GetYaxis()->CenterTitle();
  TAxis *axisMeanDummyLow = hDummyLow->GetYaxis();
  axisMeanDummyLow->ChangeLabel(-1, -1, -1, -1, -1, -1, " "); // last 
  axisMeanDummyLow->ChangeLabel(1, -1, -1, -1, -1, -1, " "); // first
  hDummyLow->Draw("same");

  // Define colors
  for(int iClass = 0; iClass < numMult + 1; iClass++) { 
    if(iClass==0) {
      color[iClass] = kBlack;
    } else {
      color[iClass] = numMult + FI - iClass;
    }
  }

  // Plotting
  for (Int_t iFile = 0; iFile < numMult + 1; iFile++) {
    // Up
    padSigmaUp->cd();
    StyleHisto(hSigma[iFile], sigmaYLow[partType], sigmaYUp[partType], color[iFile], MarkerMult[iFile], "", hSigma[iFile]->GetYaxis()->GetTitle(), "", 0, 0, 0, 1.5, 1.0, SizeMult[iFile], 0.0, 0.05, 0.0, 0.035, 0.005);
    hSigma[iFile]->GetYaxis()->SetMaxDigits(4);
    hSigma[iFile]->GetYaxis()->SetDecimals(kTRUE);
    legSigma->AddEntry(hSigma[iFile], SmoltBis[iFile], "pef");
    hSigma[iFile]->Draw("same ex0");
    hSigma[iFile]->SetFillStyle(0);
    //hSigma[iFile]->Draw("same e2");

    // Low
    padSigmaLow->cd();
    StyleHisto(hSigmaRatio[iFile], 0.1 , 11, color[iFile], MarkerMult[iFile], "#it{p}_{T} (GeV/#it{c})", "Ratio to 0-100%", "", 0, 0, 0, 1.0, 0.7, SizeMult[iFile], 0.08, 0.08, 0.08, 0.08, 0.005);
    if (iFile != 0) {
      hSigmaRatio[iFile]->Draw("same ex0");
      hSigmaRatio[iFile]->SetFillStyle(0);
      //hSigmaRatio[iFile]->Draw("same e2");
    }
  }

  padSigmaUp->cd();
  canvasSigma->Draw();
  legendTitle->Draw();
  legSigma->Draw();
  gPad->Update();
  canvasSigma->Update();

  padSigmaLow->cd();
  DrawHorLine(hSigmaRatio[0]->GetXaxis()->GetBinUpEdge(hSigmaRatio[0]->GetNbinsX()) + 0.5, 1.0);
  outputfile->cd();
  canvasSigma->Write();

  // End of Code
  std::cout << "\x1B[1;32m"; // Set text color to green
  std::cout << "\n************* Sigma Comparison in Multiplicity Classes Plot is Created Successfully! *************\n";
  std::cout << "\x1B[0m"; // Reset text color
}
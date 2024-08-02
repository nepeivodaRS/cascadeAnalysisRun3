#include "help.h"

void pmCompare(const Int_t partType = 0, // 0: Xi, 1: Omega
               const Int_t inel = 0, // inel > N (0/1)
               const TString workingDir = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM",
               const TString postFix = "")
{
  //gROOT->SetBatch(kTRUE);
  // Start of Code
  std::cout << "\x1B[1;33m"; // Set text color to yellow
  std::cout << "\n************ Starting Comaprison Between + and - ************\n";
  std::cout << "\x1B[0m"; // Reset text color

  gStyle->SetOptStat(0);

  // Output File //
  TString outputfilePath;
  if(partType == 0) {
    outputfilePath = workingDir + "/comparisons/pmComaprison_Xi.root";
  } else {
    outputfilePath = workingDir + "/comparisons/pmComaprison_Omega.root";
  }
  TFile *outputfile = new TFile(outputfilePath, "RECREATE");

  // Files with yields in all mult. classes + MB
  TFile* fileDataIn_Plus[numMult + 1];
  TFile* fileDataIn_Minus[numMult + 1];

  // Minus //
  cout << "Processing Minus Sign Cascades:"<< std::endl;
  Int_t nParticle = -1;
  if(partType == 0) {
    nParticle = 0;
  } else {
    nParticle = 3;
  }

  // MB
  fileDataIn_Minus[0] = TFile::Open(workingDir + "/yieldsOutEffCorr" +  "/yield_" + particleNames[nParticle] + "_MB_inel" + inel + postFix + ".root");
  if (!fileDataIn_Minus[0] || fileDataIn_Minus[0]->IsZombie()) {
    std::cerr << "Error opening input data file for MB!" << std::endl;
    return;
  } else {
    cout << "file for MB yield is opened"<< std::endl;
  }
  // in mult. classes
  for (Int_t iFile = 1; iFile < numMult + 1; iFile++) {
    TString fileInPath = workingDir + "/yieldsOutEffCorr" + "/yield_" + particleNames[nParticle] + "_" + multiplicityPerc[iFile - 1] + "-" + multiplicityPerc[iFile] + "_inel" + inel + postFix  + ".root";
    fileDataIn_Minus[iFile] = TFile::Open(fileInPath);
    if (!fileDataIn_Minus[iFile] || fileDataIn_Minus[iFile]->IsZombie()) {
      std::cerr << "Error opening input data file for mult. class: " << multiplicityPerc[iFile - 1] <<  " - " << multiplicityPerc[iFile] << std::endl;
      return;
    } else {
      cout << "file for mult. class: " << multiplicityPerc[iFile - 1] <<  " - " << multiplicityPerc[iFile] << " is opened"<< std::endl;
    }
  }

  // Plus //
  cout << "Processing Plus Sign Cascades:"<< std::endl;
  if(partType == 0) {
    nParticle = 1;
  } else {
    nParticle = 4;
  }

  // MB
  fileDataIn_Plus[0] = TFile::Open(workingDir + "/yieldsOutEffCorr" +  "/yield_" + particleNames[nParticle] + "_MB_inel" + inel + postFix + ".root");
  if (!fileDataIn_Plus[0] || fileDataIn_Plus[0]->IsZombie()) {
    std::cerr << "Error opening input data file for MB!" << std::endl;
    return;
  } else {
    cout << "file for MB yield is opened"<< std::endl;
  }
  // in mult. classes
  for (Int_t iFile = 1; iFile < numMult + 1; iFile++) {
    TString fileInPath = workingDir + "/yieldsOutEffCorr" + "/yield_" + particleNames[nParticle] + "_" + multiplicityPerc[iFile - 1] + "-" + multiplicityPerc[iFile] + "_inel" + inel + postFix  + ".root";
    fileDataIn_Plus[iFile] = TFile::Open(fileInPath);
    if (!fileDataIn_Plus[iFile] || fileDataIn_Plus[iFile]->IsZombie()) {
      std::cerr << "Error opening input data file for mult. class: " << multiplicityPerc[iFile - 1] <<  " - " << multiplicityPerc[iFile] << std::endl;
      return;
    } else {
      cout << "file for mult. class: " << multiplicityPerc[iFile - 1] <<  " - " << multiplicityPerc[iFile] << " is opened"<< std::endl;
    }
  }

  // Compute Ratios //
  TH1F* hYield_Minus[numMult + 1];
  TH1F* hYield_Plus[numMult + 1];
  TH1F* hYieldsRatio[numMult + 1];
  for (Int_t iFile = 0; iFile < numMult + 1; iFile++) {
    // Get yields
    if(partType == 0) {
      hYield_Minus[iFile] = (TH1F *)fileDataIn_Minus[iFile]->Get("effCorrYield/Yield_" + particleNames[0]);
      hYield_Plus[iFile] = (TH1F *)fileDataIn_Plus[iFile]->Get("effCorrYield/Yield_" + particleNames[1]);
    } else {
      hYield_Minus[iFile] = (TH1F *)fileDataIn_Minus[iFile]->Get("effCorrYield/Yield_" + particleNames[3]);
      hYield_Plus[iFile] = (TH1F *)fileDataIn_Plus[iFile]->Get("effCorrYield/Yield_" + particleNames[4]);
    }
    // Ratio
    hYieldsRatio[iFile] = (TH1F *)hYield_Minus[iFile]->Clone(Form("YieldMinus_%i_", iFile));
    hYieldsRatio[iFile]->Divide(hYield_Plus[iFile]);
    ErrRatioCorr(hYield_Minus[iFile], hYield_Plus[iFile], hYieldsRatio[iFile], 0);
  }

  // Scale yields for a better separation on the plot
  TString sScaleFactor[numMult + 1];
  hYield_Minus[0]->Scale(ScaleFactorMB);
  hYield_Plus[0]->Scale(ScaleFactorMB);
  sScaleFactor[0] = Form(" (x2^{%i})", int(log2(ScaleFactorMB)));
  for (Int_t iFile = 1; iFile < numMult + 1; iFile++) {
    hYield_Minus[iFile]->Scale(ScaleFactor[iFile-1]);
    hYield_Plus[iFile]->Scale(ScaleFactor[iFile-1]);
    if (int(log2(ScaleFactor[iFile-1])) == 0) {
      sScaleFactor[iFile] = "";
    } else if (int(log2(ScaleFactor[iFile-1])) == 1) {
      sScaleFactor[iFile] = " (x2)";
    } else {
      sScaleFactor[iFile] = Form(" (x2^{%i})", int(log2(ScaleFactor[iFile-1])));
    }
  }

  TString SmoltBis[numMult + 1];
  for (Int_t m = 0; m < numMult + 1; m++) {
    if (m == 0) { 
      SmoltBis[m] += Form("%.1f#minus%.1f%s", 0.0, 100.0, "%") + sScaleFactor[m];
    } else {
      SmoltBis[m] += Form("%.1f#minus%.1f%s", multiplicityPerc[m-1], multiplicityPerc[m], "%") + sScaleFactor[m];
    }
  }

  // Define colors
  for(int iClass = 0; iClass < numMult + 1; iClass++) { 
    if(iClass==0) {
      color[iClass] = kBlack;
    } else {
      color[iClass] = numMult + FI - iClass;
    }
  }

  Double_t yieldYLow[2] = {5*1e-11, 0.2*1e-8};
  Double_t yieldYUp[2] = {9*1e5, 9*1e4};

  // Setup pt binning
  Int_t numPtBins = 0;
  if (nParticle <= 2) {
    numPtBins = numPtXi; // Xi
  } else {
    numPtBins = numPtOmega; // Omega
  }

  Double_t* binpt = new Double_t[numPtBins + 1];

  if (nParticle <= 2) {
    for (int i = 0; i < (numPtBins + 1); ++i) {
        binpt[i] = binptXi[i]; // Xi
      }
  } else {
    for (int i = 0; i < (numPtBins + 1); ++i) {
        binpt[i] = binptOmega[i]; // Omega
      }
  }

  outputfile->cd();
  {
    TCanvas* canvasCompPlot = new TCanvas("canvasCompPlot_comparePM", "canvasCompPlot_comparePM", 0, 70, 620, 850);
    StyleCanvas(canvasCompPlot, 0.15, 0.05, 0.05, 0.15);
    TPad *padUp = new TPad("padUp", "padUp", 0, 0.36, 1, 1);
    TPad *padLow = new TPad("padLow", "padLow", 0, 0, 1, 0.36);
    StylePad(padUp, 0.15, 0.05, 0.05, 0.);
    StylePad(padLow, 0.15, 0.05, 0.02, 0.2);
    canvasCompPlot->cd();
    padUp->Draw();
    padLow->Draw();

    // Legend
    TLegend *legClasses = new TLegend(0.186, 0.0625, 0.893, 0.3125);
    legClasses->SetHeader("FT0M % (Open Circles: +; Full ones: -)");
    legClasses->SetNColumns(3);
    TLegendEntry *legEffAccInClassesHeader = (TLegendEntry *)legClasses->GetListOfPrimitives()->First();
    legEffAccInClassesHeader->SetTextSize(0.04);
    StyleLegend(legClasses, 0.0, 0.0);

    TLegend *legendTitle = new TLegend(0.608, 0.636, 0.907, 0.893);
    StyleLegend(legendTitle, 0.0, 0.0);
    legendTitle->SetTextAlign(33);
    legendTitle->SetTextSize(0.05);
    legendTitle->SetTextFont(42);
    legendTitle->SetLineColorAlpha(0.,0.);
    legendTitle->SetFillColorAlpha(0.,0.);
    legendTitle->AddEntry("", "#bf{ALICE Work In Progress}", "");
    legendTitle->AddEntry("", "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    if(partType == 0) {
      legendTitle->AddEntry("", particleSymnbols[0] + "; " + particleSymnbols[1] + ", |y| < 0.5", "");
    } else {
      legendTitle->AddEntry("", particleSymnbols[3] + "; " + particleSymnbols[4] + ", |y| < 0.5", "");
    }
    legendTitle->AddEntry("", "LHC22o_pass6", "");

    // Dummy histograms
    TH1F *hDummyUp = new TH1F("hDummyUp_comparePM", "hDummyUp_comparePM", 10000, 0, binpt[numPtBins] + 0.5);
    // Up
    for (Int_t i = 1; i <= hDummyUp->GetNbinsX(); i++)
      hDummyUp->SetBinContent(i, 1e-12);
    canvasCompPlot->cd();
    padUp->cd();
    padUp->SetLogy();
    StyleHisto(hDummyUp, yieldYLow[partType], yieldYUp[partType], 1, 1, sPt, sdNdPtdY, "", 0, 0, 0, 1.5, 1.0, 0, 0.0, 0.05, 0.0, 0.035, hDummyUp->GetXaxis()->GetLabelOffset(), 0.005);
    SetTickLength(hDummyUp, 0.025, 0.03);
    TAxis *axisDummyUpY = hDummyUp->GetYaxis();
    axisDummyUpY->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
    axisDummyUpY->ChangeLabel(2, -1, -1, -1, -1, -1, " ");
    hDummyUp->GetXaxis()->SetRangeUser(0, binpt[numPtBins] + 0.5);
    hDummyUp->Draw("same");
    // Low
    TH1F *hDummyLow = new TH1F("hDummyLow_comparePM", "hDummyLow_comparePM", 10000, 0, binpt[numPtBins] + 0.5);
    for (Int_t i = 1; i <= hDummyLow->GetNbinsX(); i++)
      hDummyLow->SetBinContent(i, 1);
    padLow->cd();
    StyleHisto(hDummyLow, 0.7, 1.7, 1, 1, sPt, "- / +", "", 0, 0, 0, 1.0, 0.7, 0, 0.08, 0.08, 0.08, 0.07, hDummyLow->GetXaxis()->GetLabelOffset(), 0.01);
    SetTickLength(hDummyLow, 0.025, 0.03);
    TAxis *axisDummyLowY = hDummyLow->GetYaxis();
    //padLow->SetGridy(); // Only grid on x-axis
    axisDummyLowY->ChangeLabel(-1, -1, -1, -1, -1, -1, " "); // last 
    axisDummyLowY->ChangeLabel(1, -1, -1, -1, -1, -1, " "); // first
    hDummyLow->GetXaxis()->SetRangeUser(0, binpt[numPtBins] + 0.5);
    hDummyLow->GetYaxis()->CenterTitle();
    hDummyLow->Draw("same");

    //Plotting
    for (Int_t i = 0; i < numMult + 1; i++) {
      padUp->cd();
      StyleHisto(hYield_Minus[i], yieldYLow[partType], yieldYUp[partType], color[i], 8, "", hYield_Minus[i]->GetYaxis()->GetTitle(), "", 0, 0, 0, 0, 0, SizeMult[i], 0, 0, 0, 0, 0, 0);
      hYield_Minus[i]->Draw("same");

      StyleHisto(hYield_Plus[i], yieldYLow[partType], yieldYUp[partType], color[i], 4, "", hYield_Plus[i]->GetYaxis()->GetTitle(), "", 0, 0, 0, 0, 0, SizeMult[i], 0, 0, 0, 0, 0, 0);
      hYield_Plus[i]->Draw("same");

      legClasses->AddEntry(hYield_Minus[i], SmoltBis[i], "pef");

      padLow->cd();
      StyleHisto(hYieldsRatio[i], 0.1 , 1.9, color[i], 8, sPt, "rec./gen.", "", 0, 0, 0, 0, 0, SizeMult[i], 0, 0, 0, 0, 0, 0);
      hYieldsRatio[i]->Draw("same");
    }

    padUp->cd();
    legendTitle->Draw("same");
    legClasses->Draw("same");
    canvasCompPlot->Write();
  }

  // End of Code
  std::cout << "\x1B[1;32m"; // Set text color to green
  std::cout << "\n************* Comaprison Between + and - is Performed Successfully! *************\n";
  std::cout << "\x1B[0m"; // Reset text color
}
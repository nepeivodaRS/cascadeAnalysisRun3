#include "help.h"

void yieldInMult(const Int_t nParticle = 2, // 0-2 : xi, 3-5 : omega
                 const Int_t inel = 0, // inel > N (0/1)
                 const TString workingDir = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM",
                 const TString fileData = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisPostSQM/23jul-lhc22o-pass6-test/AnalysisResults.root",
                 const TString postFix = "",
                 const TString effPostFix = "")
{
  //gROOT->SetBatch(kTRUE);
  // Start of Code
  std::cout << "\x1B[1;33m"; // Set text color to yellow
  std::cout << "\n************ Starting Making Spectra in Multiplicity Classes Plot ************\n";
  std::cout << "\x1B[0m"; // Reset text color

  gStyle->SetOptStat(0);

  TString outputfilePath = workingDir + "/yieldsInMult/" + "yieldInMult_" + particleNames[nParticle] + "_inel" + inel + postFix + ".root";
  TFile *outputfile = new TFile(outputfilePath, "RECREATE");

  // Get Mult. % event distribution
  TFile* fileDataInStat = TFile::Open(fileData);
  if (!fileDataInStat || fileDataInStat->IsZombie()) {
      std::cerr << "Error opening input data file!" << std::endl;
      return;
  }
  TDirectory* cascqaanalysisDir = fileDataInStat->GetDirectory("lf-cascqaanalysis");
  if (!cascqaanalysisDir)
  {
    std::cerr << "`lf-cascqaanalysis` directory is not found!" << std::endl;
    return;
  }
  cascqaanalysisDir->cd();
  TH2F* hCentFT0M_rec_data = (TH2F *)cascqaanalysisDir->Get("hCentFT0M_rec");
  if (!hCentFT0M_rec_data)
  {
    std::cerr << "Histogram `hCentFT0M_rec` is not found!" << std::endl;
    return;
  }

  hCentFT0M_rec_data->GetYaxis()->SetRange(2 + inel, 3); // INEL>0

  TH1F* hCentFT0M_rec_data_1D = (TH1F*)hCentFT0M_rec_data->ProjectionX();
  hCentFT0M_rec_data_1D = (TH1F*)hCentFT0M_rec_data_1D->Rebin(numMult, "hCentFT0M_rec_data_1D_rebinned", multiplicityPerc);

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

  TString fileInPathEff = workingDir + "/efficiencies" + "/efficiency_" + particleNames[nParticle] + effPostFix + postFix + ".root";
  TFile* fileEffin = TFile::Open(fileInPathEff);
  if (!fileEffin || fileEffin->IsZombie()) {
      std::cerr << "Error opening `fileEffin` data file!" << std::endl;
      return;
  }
  TH1F* cumulativeEvents = (TH1F *)fileEffin->Get("effEvent/eventCorrCumulative");

  TString fileInPathSyst = workingDir + "/systematics/sourcesOfSyst/totalSystematics.root";
  TFile* fileSystin = TFile::Open(fileInPathSyst);
  if (!fileSystin || fileSystin->IsZombie()) {
      std::cerr << "Error opening `fileSystin` data file!" << std::endl;
      return;
  }
  TH1F* hSystUncert = (TH1F *)fileSystin->Get("hSystTotal");

  TH1F* hYield[numMult + 1];
  TH1F* hYieldSyst[numMult + 1];
  TH1F* hYieldsRatio[numMult + 1];
  TH1F* hYieldsRatioSyst[numMult + 1];
  TH1F* hYieldsDenom;
  TH1F* hYieldsDenomSyst;
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
      hYieldSyst[iFile]->SetBinError(bin, hSystUncert->GetBinContent(bin)*hYield[iFile]->GetBinContent(bin));
    }

    hYieldsRatio[iFile] = (TH1F *)hYield[iFile]->Clone(Form("YieldClone_%i_", iFile) + particleNames[nParticle]);
    // Syst. histo
    hYieldsRatioSyst[iFile] = (TH1F *)hYieldSyst[iFile]->Clone(Form("YieldCloneSyst_%i_", iFile) + particleNames[nParticle]);
    if (iFile == 0) {
      hYieldsDenom = (TH1F *)hYield[0]->Clone("YieldCloneForDenom_" + particleNames[nParticle]);
      // Syst. histo
      hYieldsDenomSyst = (TH1F *)hYieldSyst[0]->Clone("YieldCloneForDenomSyst_" + particleNames[nParticle]);
    }

    hYieldsRatio[iFile]->Divide(hYieldsDenom);
    // Syst. histo
    hYieldsRatioSyst[iFile]->Divide(hYieldsDenomSyst);
    if (iFile != 0) {
      ErrRatioCorr(hYield[iFile], hYieldsDenom, hYieldsRatio[iFile], 0);
      // Syst. histo
      ErrRatioCorr(hYieldSyst[iFile], hYieldsDenomSyst, hYieldsRatioSyst[iFile], 0);
    }
  }

  TH1F* hYieldMBsummed;
  TH1F* hYieldMBsummedSyst;
  TH1F* hYieldsRatioOfMBsummed;
  TH1F* hYieldsRatioOfMBsummedSyst;
  hYieldMBsummed = (TH1F*)hYield[0]->Clone("hYieldMBsummed");
  hYieldMBsummed->Reset();
  // Syst. histo
  hYieldMBsummedSyst = (TH1F*)hYieldMBsummed->Clone("hYieldMBsummedSyst");
  Double_t sumWeights = 0;
  for (Int_t iFile = 1; iFile < numMult + 1; iFile++) {
    TH1F* hYieldPerPerc = (TH1F*)hYield[iFile]->Clone("hYieldPerPerc");
    // Syst. histo
    TH1F* hYieldPerPercSyst = (TH1F*)hYieldSyst[iFile]->Clone("hYieldPerPercSyst");
    Double_t weight = 0;
    if (iFile > 1.5) {
      weight = cumulativeEvents->GetBinContent(iFile) - cumulativeEvents->GetBinContent(iFile - 1);
      std::cout << "1: " << cumulativeEvents->GetBinContent(iFile) << std::endl;
      std::cout << "2: " << cumulativeEvents->GetBinContent(iFile - 1) << std::endl;
    } else {
      weight = cumulativeEvents->GetBinContent(iFile);
    }
    sumWeights += weight * hCentFT0M_rec_data_1D->GetBinContent(iFile);
    hYieldPerPerc->Scale(weight * hCentFT0M_rec_data_1D->GetBinContent(iFile));
    // Syst. histo
    hYieldPerPercSyst->Scale(weight * hCentFT0M_rec_data_1D->GetBinContent(iFile));

    std::cout << "Fraction: " << sumWeights/2.61149e+08*100 << std::endl;
    std::cout << "weight: " << weight * hCentFT0M_rec_data_1D->GetBinContent(iFile) << " class: " << multiplicityPerc[iFile - 1] << " - " <<  multiplicityPerc[iFile] << std::endl;
    hYieldMBsummed->Add(hYieldPerPerc);
    // Syst. histo
    hYieldMBsummedSyst->Add(hYieldPerPercSyst);
  }
  std::cout << "sum weights: " << sumWeights << std::endl;
  hYieldMBsummed->Scale(1. / sumWeights);
  // Syst. histo
  hYieldMBsummedSyst->Scale(1. / sumWeights);
  hYieldsRatioOfMBsummed = (TH1F*)hYieldMBsummed->Clone("hYieldsRatioToMBsummed");
  hYieldsRatioOfMBsummed->Divide(hYield[0]);
  ErrRatioCorr(hYieldMBsummed, hYield[0], hYieldsRatioOfMBsummed, 0);
  // Syst. histo
  hYieldsRatioOfMBsummedSyst = (TH1F*)hYieldMBsummedSyst->Clone("hYieldsRatioToMBsummedSyst");
  hYieldsRatioOfMBsummedSyst->Divide(hYieldSyst[0]);
  ErrRatioCorr(hYieldMBsummedSyst, hYieldSyst[0], hYieldsRatioOfMBsummedSyst, 0);

  // Scale yields for a better separation on the plot
  TString sScaleFactor[numMult + 1];
  hYield[0]->Scale(ScaleFactorMB);
  hYieldSyst[0]->Scale(ScaleFactorMB);
  hYieldMBsummed->Scale(ScaleFactorMB);
  hYieldMBsummedSyst->Scale(ScaleFactorMB);
  sScaleFactor[0] = Form(" (x2^{%i})", int(log2(ScaleFactorMB)));
  for (Int_t iFile = 1; iFile < numMult + 1; iFile++) {
    hYield[iFile]->Scale(ScaleFactor[iFile-1]);
    hYieldSyst[iFile]->Scale(ScaleFactor[iFile-1]);
    if (int(log2(ScaleFactor[iFile-1])) == 0) {
      sScaleFactor[iFile] = "";
    } else if (int(log2(ScaleFactor[iFile-1])) == 1) {
      sScaleFactor[iFile] = " (x2)";
    } else {
      sScaleFactor[iFile] = Form(" (x2^{%i})", int(log2(ScaleFactor[iFile-1])));
    }
  }

  // Yield canvas
  TCanvas* canvasYield = new TCanvas("yieldSummary_" + particleNames[nParticle], "yieldSummary_" + particleNames[nParticle], 0, 70, 620, 850);
  StyleCanvas(canvasYield, 0.15, 0.05, 0.05, 0.15);
  // Yield Pads
  TPad *padYieldUp = new TPad("padYieldUp_" + particleNames[nParticle], "padYieldUp_" + particleNames[nParticle], 0, 0.36, 1, 1);
  TPad *padYieldLow = new TPad("padYieldLow_" + particleNames[nParticle], "padYieldLow_" + particleNames[nParticle], 0, 0, 1, 0.36);
  StylePad(padYieldUp, 0.15, 0.05, 0.05, 0.);
  StylePad(padYieldLow, 0.15, 0.05, 0.02, 0.2);
  canvasYield->cd();
  padYieldUp->Draw();
  padYieldLow->Draw();

  // Yield Legend
  TLegend *legYield = new TLegend(0.186, 0.0625, 0.893, 0.3125);
  legYield->SetHeader("FT0M Multiplicity Percentile");
  legYield->SetNColumns(3);
  TLegendEntry *legYieldHeader = (TLegendEntry *)legYield->GetListOfPrimitives()->First();
  legYieldHeader->SetTextSize(0.04);
  //legYield->SetTextAlign(13); // adjust text wrt markers
  StyleLegend(legYield, 0.0, 0.0);

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
      SmoltBis[m] += Form("%.1f#minus%.1f%s", 0.0, 100.0, "%") + sScaleFactor[m];
    } else {
      SmoltBis[m] += Form("%.1f#minus%.1f%s", multiplicityPerc[m-1], multiplicityPerc[m], "%") + sScaleFactor[m];
    }
  }

  Double_t yieldYLow[2] = {1e-9, 0.2*1e-8};
  Double_t yieldYUp[2] = {9*1e5, 9*1e4};

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
  canvasYield->cd();
  padYieldUp->cd();
  padYieldUp->SetLogy();
  StyleHisto(hDummy, yieldYLow[partType], yieldYUp[partType], 1, 1, sPt, sdNdPtdY, "", 0, 0, 0, 1.5, 1.0, 0, 0.0, 0.05, 0.0, 0.035, 0.005, 0.005);
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
  padYieldLow->SetLogy();
  StyleHisto(hDummyLow, 0.09, 11, 1, 1, sPt, "Ratio to 0-100%", "", 0, 0, 0, 1.0, 0.7, 0, 0.08, 0.08, 0.08, 0.07, 0.01, 0.01);

  hDummyLow->GetXaxis()->SetRangeUser(0, hYieldsRatio[0]->GetXaxis()->GetBinUpEdge(hYieldsRatio[0]->GetNbinsX()) + 0.5);
  hDummyLow->GetYaxis()->CenterTitle();
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
    padYieldUp->cd();
    StyleHisto(hYield[iFile], yieldYLow[partType], yieldYUp[partType], color[iFile], MarkerMult[iFile], "", hYield[iFile]->GetYaxis()->GetTitle(), "", 0, 0, 0, 1.5, 1.0, SizeMult[iFile], 0.0, 0.05, 0.0, 0.035, 0.005, 0.005);
    StyleHisto(hYieldSyst[iFile], yieldYLow[partType], yieldYUp[partType], color[iFile], MarkerMult[iFile], "", hYield[iFile]->GetYaxis()->GetTitle(), "", 0, 0, 0, 1.5, 1.0, SizeMult[iFile], 0.0, 0.05, 0.0, 0.035, 0.005, 0.005);
    hYield[iFile]->GetYaxis()->SetMaxDigits(4);
    hYield[iFile]->GetYaxis()->SetDecimals(kTRUE);
    legYield->AddEntry(hYield[iFile], SmoltBis[iFile], "pef");
    hYield[iFile]->Draw("same ex0");
    hYield[iFile]->SetFillStyle(0);
    hYieldSyst[iFile]->SetFillStyle(0);
    hYieldSyst[iFile]->Draw("same e2"); // e2 - syst

    // Low
    padYieldLow->cd();
    StyleHisto(hYieldsRatio[iFile], 0.1 , 11, color[iFile], MarkerMult[iFile], sPt, "Ratio to 0-100%", "", 0, 0, 0, 1.0, 0.7, SizeMult[iFile], 0.08, 0.08, 0.08, 0.07, 0.005, 0.01);
    StyleHisto(hYieldsRatioSyst[iFile], 0.1 , 11, color[iFile], MarkerMult[iFile], sPt, "Ratio to 0-100%", "", 0, 0, 0, 1.0, 0.7, SizeMult[iFile], 0.08, 0.08, 0.08, 0.07, 0.005, 0.01);
    if (iFile != 0) {
      hYieldsRatio[iFile]->Draw("same ex0");
      hYieldsRatio[iFile]->SetFillStyle(0);
      hYieldsRatioSyst[iFile]->SetFillStyle(0);
      hYieldsRatioSyst[iFile]->Draw("same e2");
    }
  }

  // Plot summed MB //
  // Up
  padYieldUp->cd();
  StyleHisto(hYieldMBsummed, yieldYLow[partType], yieldYUp[partType], kGray, MarkerMult[numMult + 2], "", "", "", 0, 0, 0, 1.5, 1.0, 0., 0.0, 0.05, 0.0, 0.035, 0.005, 0.005);
  StyleHisto(hYieldMBsummedSyst, yieldYLow[partType], yieldYUp[partType], kGray, MarkerMult[numMult + 2], "", "", "", 0, 0, 0, 1.5, 1.0, 0., 0.0, 0.05, 0.0, 0.035, 0.005, 0.005);
  legYield->AddEntry(hYieldMBsummed, "MB (Sum)", "pef");
  hYieldMBsummed->Draw("same ex0");
  hYieldMBsummed->SetFillStyle(0);
  hYieldMBsummedSyst->SetFillStyle(0);
  hYieldMBsummedSyst->Draw("same e2"); // e2 - syst
  // Low
  padYieldLow->cd();
  StyleHisto(hYieldsRatioOfMBsummed, 0.1 , 11, kGray, MarkerMult[numMult + 2], sPt, "Ratio to 0-100%", "", 0, 0, 0, 1.0, 0.7, 0.7, 0.08, 0.08, 0.08, 0.07, 0.005, 0.01);
  StyleHisto(hYieldsRatioOfMBsummedSyst, 0.1 , 11, kGray, MarkerMult[numMult + 2], sPt, "Ratio to 0-100%", "", 0, 0, 0, 1.0, 0.7, 0.7, 0.08, 0.08, 0.08, 0.07, 0.005, 0.01);
  hYieldsRatioOfMBsummed->Draw("same ex0");
  hYieldsRatioOfMBsummed->SetFillStyle(0);
  hYieldsRatioOfMBsummedSyst->SetFillStyle(0);
  hYieldsRatioOfMBsummedSyst->Draw("same e2"); // e2 - syst

  padYieldUp->cd();
  legendTitle->Draw("same");
  legYield->Draw("same");
  gPad->Update();
  canvasYield->Update();

  padYieldLow->cd();
  DrawHorLine(hYieldsRatio[0]->GetXaxis()->GetBinUpEdge(hYieldsRatio[0]->GetNbinsX()) + 0.5, 1.0);
  outputfile->cd();
  canvasYield->Write();

  // End of Code
  std::cout << "\x1B[1;32m"; // Set text color to green
  std::cout << "\n************* Spectra in Multiplicity Classes Plot is Created Successfully! *************\n";
  std::cout << "\x1B[0m"; // Reset text color
}
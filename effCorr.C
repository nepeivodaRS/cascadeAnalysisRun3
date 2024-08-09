#include "help.h"

void effCorr(const  Int_t nParticle = 2, // 0-2 : xi, 3-5 : omega
            const Int_t inel = 0, // inel > N (0/1)
            const TString fileMC = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisPostSQM/24jul-lhc24b1b/AnalysisResults.root",
            const TString fileMCPP = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisPostSQM/run3_13tev/xi/lhc24b1b/2024-07-24/AnalysisResults.root",
            const TString fileData = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisPostSQM/23jul-lhc22o-pass6-test/AnalysisResults.root",
            const TString outputDir = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisPostSQM",
            const TString yieldPostFix = "", // postfix of yields (needed to compute MB signal loss)
            const TString postFix = "_LHC24b1b") // postfix to the output file
{
  // Start of Code
  std::cout << "\x1B[1;33m"; // Set text color to yellow
  std::cout << "\n************ Starting Computing Eff. x Acc.,  Signal Loss and Event Factor ************\n";
  std::cout << "\x1B[0m"; // Reset text color

  // Set Style
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat("me");
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);

  // Define colors
  for(int iClass = 0; iClass < numMultEff + 1; iClass++) { 
    if(iClass==0) {
      color[iClass] = kBlack;
    } else {
      color[iClass] = numMultEff + FIEff - iClass;
    }
  }

  Double_t binPrecision = 1e-6;

  // File with calibration histograms (cascqaanalysis output of MC)
  TFile* fileMCin = TFile::Open(fileMC);
  if (!fileMCin || fileMCin->IsZombie()) {
      std::cerr << "Error opening `fileMC` data file!" << std::endl;
      return;
  }

  // File with main histograms (postprocessing output of MC)
  TFile* fileMCPPIn = TFile::Open(fileMCPP);
  if (!fileMCPPIn || fileMCPPIn->IsZombie()) {
      std::cerr << "Error opening input `fileMCPP` file!" << std::endl;
      return;
  }

  // Open yield files
  TFile *inputFile[numMult + 1];
  TString inputFilePath;
  for (Int_t iFile = 0; iFile < numMult + 1; iFile++) {
    if (iFile == 0) {
      inputFilePath = outputDir + "/yieldsOut/" + "yield_" + particleNames[nParticle] + "_MB_inel" + inel + yieldPostFix + ".root";
      inputFile[0] = TFile::Open(inputFilePath);
      if (!inputFile[0] || inputFile[0]->IsZombie()) {
        std::cerr << "Error opening input data file for MB!" << std::endl;
        return;
      } else {
        cout << "file for MB yield is opened"<< std::endl;
      }
    } else {
      inputFilePath = outputDir + "/yieldsOut/" + "yield_" + particleNames[nParticle] + "_" + multiplicityPerc[iFile - 1] + "-" + multiplicityPerc[iFile] + "_inel" + inel + yieldPostFix + ".root";
      inputFile[iFile] = TFile::Open(inputFilePath);
      if (!inputFile[iFile] || inputFile[iFile]->IsZombie()) {
        std::cerr << "Error opening input data file for mult. class: " << multiplicityPerc[iFile - 1] <<  " - " << multiplicityPerc[iFile] << std::endl;
      } else {
        cout << "file for mult. class: " << multiplicityPerc[iFile - 1] <<  " - " << multiplicityPerc[iFile] << " is opened"<< std::endl;
      }
    }
  }

  TH1F* hYield[numMult + 1];
  TH1F* hYieldCorrected[numMult + 1];
  TDirectory* yieldDir[numMult + 1];

  // Get yields
  for (Int_t iFile = 0; iFile < numMult + 1; iFile++) {
    yieldDir[iFile] = inputFile[iFile]->GetDirectory("fitParams");
    if (!yieldDir[iFile])
    {
      std::cerr << "`fitParams` directory is not found!" << std::endl;
      return;
    }
    hYield[iFile] = (TH1F *)yieldDir[iFile]->Get("Yield_" + particleNames[nParticle]);
    if (!hYield[iFile])
    {
      std::cerr << "Yield histogram is not found!" << std::endl;
      return;
    }
    hYieldCorrected[iFile] = (TH1F*)hYield[iFile]->Clone(Form("hYields_%d", iFile));
    hYieldCorrected[iFile]->SetName("Yield_" + particleNames[nParticle]);
  }

  // Get Mult. % event distribution
  TFile* fileDataIn = TFile::Open(fileData);
  if (!fileDataIn || fileDataIn->IsZombie()) {
      std::cerr << "Error opening input data file!" << std::endl;
      return;
  }
  TDirectory* cascqaanalysisDir = fileDataIn->GetDirectory("lf-cascqaanalysis");
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
  hCentFT0M_rec_data->GetXaxis()->SetRangeUser(0., 100.);
  TH1F* hCentFT0M_rec_data_1D = (TH1F*)hCentFT0M_rec_data->ProjectionX();

  TH1F* hCentFT0M_rec_data_1D_finest = (TH1F*)hCentFT0M_rec_data_1D->Clone("hCentFT0M_rec_data_1D_finest");
  TH1F* hCentFT0M_rec_data_1D_anal = (TH1F*)hCentFT0M_rec_data_1D->Clone("hCentFT0M_rec_data_1D_anal");
  hCentFT0M_rec_data_1D_anal = (TH1F*)hCentFT0M_rec_data_1D_anal->Rebin(numMult, "hCentFT0M_rec_data_1D_anal", multiplicityPerc);
  hCentFT0M_rec_data_1D_finest = (TH1F*)hCentFT0M_rec_data_1D_finest->Rebin(numMultFinest, "hCentFT0M_rec_data_1D_finest", multiplicityFinest);

  // Output File
  TString outputfilePath = outputDir + "/efficiencies/" + "efficiency_" + particleNames[nParticle] + postFix + ".root";
  TFile *outputfile = new TFile(outputfilePath, "RECREATE");

  // Setup pt binning
  Int_t numPtBins = 0;
  Int_t partType;
  if (nParticle <= 2) {
    numPtBins = numPtXi; // Xi
    partType = 0;
  } else {
    numPtBins = numPtOmega; // Omega
    partType = 1;
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

  // Histos of generated cascades from generated events with accepted z vrtx + chosen event type (evSelFlag)
  TH3F* hCascadePlusTrue;
  TH3F* hCascadeMinusTrue;
  TH3F* hCascadeSumTrue;

  // Histos of generated cascades from generated events with accepted z vrtx + chosen event type (evSelFlag) + associated to the accepted reconstructed event of the same type
  TH3F* hCascadePlusTrueAssoc;
  TH3F* hCascadeMinusTrueAssoc;
  TH3F* hCascadeSumTrueAssoc;

  if(partType == 1){
    hCascadePlusTrueAssoc = (TH3F*)fileMCPPIn->Get("lf-cascpostprocessing/hPtOmegaPlusTrueAssocWithSelColl");
    hCascadeMinusTrueAssoc = (TH3F*)fileMCPPIn->Get("lf-cascpostprocessing/hPtOmegaMinusTrueAssocWithSelColl");
    hCascadePlusTrue = (TH3F*)fileMCPPIn->Get("lf-cascpostprocessing/hPtOmegaPlusTrue");
    hCascadeMinusTrue = (TH3F*)fileMCPPIn->Get("lf-cascpostprocessing/hPtOmegaMinusTrue");
  } else {
    hCascadePlusTrueAssoc = (TH3F*)fileMCPPIn->Get("lf-cascpostprocessing/hPtXiPlusTrueAssocWithSelColl");
    hCascadeMinusTrueAssoc = (TH3F*)fileMCPPIn->Get("lf-cascpostprocessing/hPtXiMinusTrueAssocWithSelColl");
    hCascadePlusTrue = (TH3F*)fileMCPPIn->Get("lf-cascpostprocessing/hPtXiPlusTrue");
    hCascadeMinusTrue = (TH3F*)fileMCPPIn->Get("lf-cascpostprocessing/hPtXiMinusTrue");
  }

  if (!hCascadePlusTrueAssoc || !hCascadeMinusTrueAssoc || !hCascadePlusTrue || !hCascadeMinusTrue)
  {
    std::cerr << "At least one of the histograms is not found in `fileMCPPIn`!" << std::endl;
    return;
  }

  // Reconstructed cascades (only true ones)
  TH3F* hCascadePlusTrueRec = (TH3F*)fileMCPPIn->Get("lf-cascpostprocessing/hPtCascPlusTrueRec");
  TH3F* hCascadeMinusTrueRec = (TH3F*)fileMCPPIn->Get("lf-cascpostprocessing/hPtCascMinusTrueRec");
  TH3F* hCascadeSumTrueRec = (TH3F*)hCascadePlusTrueRec->Clone("SumRec");
  hCascadeSumTrueRec->Add(hCascadeMinusTrueRec);
  hCascadeSumTrueRec->SetName("hPtRecSum");

  hCascadeSumTrue = (TH3F*)hCascadePlusTrue->Clone("SumTrue");
  hCascadeSumTrue->Add(hCascadeMinusTrue);
  hCascadeSumTrue->SetName("hPtSumTrue");

  hCascadeSumTrueAssoc = (TH3F*)hCascadePlusTrueAssoc->Clone("SumTrueAssoc");
  hCascadeSumTrueAssoc->Add(hCascadeMinusTrueAssoc);
  hCascadeSumTrueAssoc->SetName("hPtSumTrueAssoc");

  TDirectory* dirSupport = outputfile->mkdir("supportHistos");
  TDirectory* dirEffAcc = outputfile->mkdir("effAcc");
  TDirectory* dirSignalLoss = outputfile->mkdir("effSignalLoss");

  Double_t leftRap = -0.5;
  Double_t rightRap = 0.5;

  TH1D* hEffAccInClasses[numMultEff + 1];
  TH1D* hEffAccInClassesRatio[numMultEff + 1];
  TH1D* hSignalLossInClasses[numMultEff + 1];
  TH1D* hSignalLossInClassesRatio[numMultEff + 1];

  for (Int_t i = 0; i < numMultEff + 1; i++) {
    dirSupport->cd();
    // True lvl Not Assoc //
    TH3F* hCascadeTrue;
    if (nParticle == 2 || nParticle == 5) {
      hCascadeTrue = (TH3F*)hCascadeSumTrue->Clone(Form("cascadeTrue_%i", i));
    } else if (nParticle == 0 || nParticle == 3) {
      hCascadeTrue = (TH3F*)hCascadeMinusTrue->Clone(Form("cascadeTrue_%i", i));
    } else if (nParticle == 1 || nParticle == 4) {
      hCascadeTrue = (TH3F*)hCascadePlusTrue->Clone(Form("cascadeTrue_%i", i));
    }
    hCascadeTrue->GetXaxis()->SetRange(hCascadeTrue->GetXaxis()->FindBin(binpt[0] + binPrecision), hCascadeTrue->GetXaxis()->FindBin(binpt[numPtBins] - binPrecision)); // Set pt interval
    hCascadeTrue->GetYaxis()->SetRange(hCascadeTrue->GetYaxis()->FindBin(leftRap + binPrecision), hCascadeTrue->GetYaxis()->FindBin(rightRap - binPrecision)); // Set rapidity interval
    if(i!=0){
      hCascadeTrue->GetZaxis()->SetRange(hCascadeTrue->GetZaxis()->FindBin(multiplicityPercEff[i - 1] + binPrecision), hCascadeTrue->GetZaxis()->FindBin(multiplicityPercEff[i] - binPrecision)); // Set mult. % interval
    }
    else {
      hCascadeTrue->GetZaxis()->SetRange(hCascadeTrue->GetZaxis()->FindBin(multiplicityPercEff[0] + binPrecision), hCascadeTrue->GetZaxis()->FindBin(multiplicityPercEff[numMultEff] - binPrecision)); // Set mult. % interval
    }
    hCascadeTrue->Sumw2();
    hCascadeTrue->Write();
    TH1D* hCascadeTrue_1D = static_cast<TH1D*>(hCascadeTrue->Project3D("x"));

    // True lvl Assoc //
    TH3F* hCascadeTrueAssoc;
    if (nParticle == 2 || nParticle == 5) {
      hCascadeTrueAssoc = (TH3F*)hCascadeSumTrueAssoc->Clone(Form("cascadeTrueAssoc_%i", i));
    } else if (nParticle == 0 || nParticle == 3) {
      hCascadeTrueAssoc = (TH3F*)hCascadeMinusTrueAssoc->Clone(Form("cascadeTrueAssoc_%i", i));
    } else if (nParticle == 1 || nParticle == 4) {
      hCascadeTrueAssoc = (TH3F*)hCascadePlusTrueAssoc->Clone(Form("cascadeTrueAssoc_%i", i));
    }
    hCascadeTrueAssoc->GetXaxis()->SetRange(hCascadeTrueAssoc->GetXaxis()->FindBin(binpt[0] + binPrecision), hCascadeTrueAssoc->GetXaxis()->FindBin(binpt[numPtBins] - binPrecision)); // Set pt interval
    hCascadeTrueAssoc->GetYaxis()->SetRange(hCascadeTrueAssoc->GetYaxis()->FindBin(leftRap + binPrecision), hCascadeTrueAssoc->GetYaxis()->FindBin(rightRap - binPrecision)); // Set rapidity interval
    if(i!=0){
      hCascadeTrueAssoc->GetZaxis()->SetRange(hCascadeTrueAssoc->GetZaxis()->FindBin(multiplicityPercEff[i - 1] + binPrecision), hCascadeTrueAssoc->GetZaxis()->FindBin(multiplicityPercEff[i] - binPrecision)); // Set mult. % interval
    }
    else {
      hCascadeTrueAssoc->GetZaxis()->SetRange(hCascadeTrueAssoc->GetZaxis()->FindBin(multiplicityPercEff[0] + binPrecision), hCascadeTrueAssoc->GetZaxis()->FindBin(multiplicityPercEff[numMultEff] - binPrecision)); // Set mult. % interval
    }
    hCascadeTrueAssoc->Sumw2();
    hCascadeTrueAssoc->Write();
    TH1D* hCascadeTrueAssoc_1D = static_cast<TH1D*>(hCascadeTrueAssoc->Project3D("x"));

    // True Rec //
    TH3F* hCascadeTrueRec;
    if (nParticle == 2 || nParticle == 5) {
      hCascadeTrueRec = (TH3F*)hCascadeSumTrueRec->Clone(Form("cascadeTrueRec_%i", i));
    } else if (nParticle == 0 || nParticle == 3) {
      hCascadeTrueRec = (TH3F*)hCascadeMinusTrueRec->Clone(Form("cascadeTrueRec_%i", i));
    } else if (nParticle == 1 || nParticle == 4) {
      hCascadeTrueRec = (TH3F*)hCascadePlusTrueRec->Clone(Form("cascadeTrueRec_%i", i));
    }
    hCascadeTrueRec->GetXaxis()->SetRange(hCascadeTrueRec->GetXaxis()->FindBin(binpt[0] + binPrecision), hCascadeTrueRec->GetXaxis()->FindBin(binpt[numPtBins] - binPrecision)); // Set pt interval
    hCascadeTrueRec->GetYaxis()->SetRange(hCascadeTrueRec->GetYaxis()->FindBin(leftRap + binPrecision), hCascadeTrueRec->GetYaxis()->FindBin(rightRap - binPrecision)); // Set rapidity interval
    if(i!=0){
      hCascadeTrueRec->GetZaxis()->SetRange(hCascadeTrueRec->GetZaxis()->FindBin(multiplicityPercEff[i - 1] + binPrecision), hCascadeTrueRec->GetZaxis()->FindBin(multiplicityPercEff[i] - binPrecision)); // Set mult. % interval
    }
    else {
      hCascadeTrueRec->GetZaxis()->SetRange(hCascadeTrueRec->GetZaxis()->FindBin(multiplicityPercEff[0] + binPrecision), hCascadeTrueRec->GetZaxis()->FindBin(multiplicityPercEff[numMultEff] - binPrecision)); // Set mult. % interval
    }
    hCascadeTrueRec->Sumw2();
    hCascadeTrueRec->Write();
    TH1D* hCascadeTrueRec_1D = static_cast<TH1D*>(hCascadeTrueRec->Project3D("x"));

    //Rebin to match analysis binning
    hCascadeTrue_1D = (TH1D*)hCascadeTrue_1D->Rebin(numPtBins, Form("hCascadeTrue_1D_rebinned_%i", i), binpt);
    hCascadeTrueAssoc_1D = (TH1D*)hCascadeTrueAssoc_1D->Rebin(numPtBins, Form("hCascadeTrueAssoc_1D_rebinned_%i", i), binpt);
    hCascadeTrueRec_1D = (TH1D*)hCascadeTrueRec_1D->Rebin(numPtBins, Form("hCascadeTrueRec_1D_rebinned_%i", i), binpt);
    hCascadeTrue_1D->Write();
    hCascadeTrueAssoc_1D->Write();
    hCascadeTrueRec_1D->Write();

    dirEffAcc->cd();

    // Efficiency x Acceptance //
    TCanvas* canvasEffAcc = new TCanvas(Form("canvasEffAcc_%i", i), Form("canvasEffAcc_%i", i), 0,100,1200,850);
    hEffAccInClasses[i] = (TH1D*)hCascadeTrueRec_1D->Clone(Form("effAcc_%i", i));
    hEffAccInClasses[i]->Divide(hCascadeTrueAssoc_1D);
    SetErrorDivide(hEffAccInClasses[i], hCascadeTrueRec_1D, hCascadeTrueAssoc_1D);
    hEffAccInClasses[i]->SetStats(0);
    StyleCanvas(canvasEffAcc, 0.15, 0.05, 0.10, 0.15);
    StyleHisto(hEffAccInClasses[i], 0, 1.2 * hEffAccInClasses[i]->GetBinContent(hEffAccInClasses[i]->GetMaximumBin()), kBlack, 20, sPt, "eff. x acc.", "", 0, 0, 0, 1.0, 1.25, 1,  0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
    hEffAccInClasses[i]->Draw();
    canvasEffAcc->Write();

    // Signal Loss //
    TCanvas* canvasSignalLoss = new TCanvas(Form("canvasSignalLoss_%i", i), Form("canvasSignalLoss_%i", i), 0,100,1200,850);
    hSignalLossInClasses[i] = (TH1D*)hCascadeTrueAssoc_1D->Clone(Form("effAcc_%i", i));
    hSignalLossInClasses[i]->Divide(hCascadeTrue_1D);
    SetErrorDivide(hSignalLossInClasses[i], hCascadeTrueAssoc_1D, hCascadeTrue_1D);
    hSignalLossInClasses[i]->SetStats(0);
    StyleCanvas(canvasSignalLoss, 0.15, 0.05, 0.10, 0.15);
    StyleHisto(hSignalLossInClasses[i], 0, 1.2 * hSignalLossInClasses[i]->GetBinContent(hSignalLossInClasses[i]->GetMaximumBin()), kBlack, 20, sPt, "Signal Loss", "", 0, 0, 0, 1.0, 1.25, 1,  0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
    hSignalLossInClasses[i]->Draw();
    canvasSignalLoss->Write();
  }

  // Efficiency x Acceptance //
  TString SmoltBis[numMultEff + 1];
  for (Int_t m = 0; m < numMultEff + 1; m++) {
    if (m == 0) { 
      SmoltBis[m] += Form("%.1f#minus%.1f%s", 0.0, 100.0, "%");
    // } else if(m == 1) {
    //   SmoltBis[m] += Form("%.3f#minus%.3f%s", multiplicityPercEff[m-1], multiplicityPercEff[m], "%");
    } else {
      SmoltBis[m] += Form("%.1f#minus%.1f%s", multiplicityPercEff[m-1], multiplicityPercEff[m], "%");
    }
  }
  gROOT->SetBatch(kFALSE);
  {
    TCanvas* canvasCompPlot = new TCanvas("canvasCompPlot_EffAcc", "canvasCompPlot_EffAcc", 0, 70, 620, 850);
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
    legClasses->SetHeader("FT0M Multiplicity Percentile");
    legClasses->SetNColumns(3);
    TLegendEntry *legEffAccInClassesHeader = (TLegendEntry *)legClasses->GetListOfPrimitives()->First();
    legEffAccInClassesHeader->SetTextSize(0.04);
    StyleLegend(legClasses, 0.0, 0.0);

    TLegend *legendTitleMC = new TLegend(0.608, 0.636, 0.907, 0.893);
    StyleLegend(legendTitleMC, 0.0, 0.0);
    legendTitleMC->SetTextAlign(33);
    legendTitleMC->SetTextSize(0.05);
    legendTitleMC->SetTextFont(42);
    legendTitleMC->SetLineColorAlpha(0.,0.);
    legendTitleMC->SetFillColorAlpha(0.,0.);
    legendTitleMC->AddEntry("", "#bf{ALICE Work In Progress}", "");
    legendTitleMC->AddEntry("", "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    legendTitleMC->AddEntry("", particleSymnbols[nParticle] + ", |y| < 0.5", "");
    legendTitleMC->AddEntry("", "LHC24b1b", "");

    // Dummy histograms
    TH1F *hDummyUp = new TH1F("hDummyUp_EffAcc", "hDummyUp_EffAcc", 10000, 0, binpt[numPtBins] + 0.5);
    // Up
    for (Int_t i = 1; i <= hDummyUp->GetNbinsX(); i++)
      hDummyUp->SetBinContent(i, 1000);
    canvasCompPlot->cd();
    padUp->cd();
    StyleHisto(hDummyUp, -0.1, 0.22, 1, 1, sPt, "eff. x. acc.", "", 0, 0, 0, 1.5, 1.0, 0, 0.0, 0.05, 0.0, 0.035, hDummyUp->GetXaxis()->GetLabelOffset(), 0.005);
    SetTickLength(hDummyUp, 0.025, 0.03);
    TAxis *axisDummyUpY = hDummyUp->GetYaxis();
    axisDummyUpY->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
    axisDummyUpY->ChangeLabel(2, -1, -1, -1, -1, -1, " ");
    hDummyUp->GetXaxis()->SetRangeUser(0, binpt[numPtBins] + 0.5);
    //padUp->SetLogy();
    hDummyUp->Draw("same");
    // Low
    TH1F *hDummyLow = new TH1F("hDummyLow_EffAcc", "hDummyLow_EffAcc", 10000, 0, binpt[numPtBins] + 0.5);
    for (Int_t i = 1; i <= hDummyLow->GetNbinsX(); i++)
      hDummyLow->SetBinContent(i, 1);
    padLow->cd();
    StyleHisto(hDummyLow, 0.2, 1.8, 1, 1, sPt, "Ratio to 0-100%", "", 0, 0, 0, 1.0, 0.7, 0, 0.08, 0.08, 0.08, 0.07, hDummyLow->GetXaxis()->GetLabelOffset(), 0.01);
    SetTickLength(hDummyLow, 0.025, 0.03);
    TAxis *axisDummyLowY = hDummyLow->GetYaxis();
    //padLow->SetGridy(); // Only grid on x-axis
    axisDummyLowY->ChangeLabel(-1, -1, -1, -1, -1, -1, " "); // last 
    axisDummyLowY->ChangeLabel(1, -1, -1, -1, -1, -1, " "); // first
    hDummyLow->GetXaxis()->SetRangeUser(0, binpt[numPtBins] + 0.5);
    hDummyLow->GetYaxis()->CenterTitle();
    hDummyLow->Draw("same");

    //Plotting
    for (Int_t i = 0; i < numMultEff + 1; i++) {
      padUp->cd();
      StyleHisto(hEffAccInClasses[i], 0, hEffAccInClasses[i]->GetMaximum()*1e2, color[i], MarkerMult[i], "", hEffAccInClasses[i]->GetYaxis()->GetTitle(), "", 0, 0, 0, 0, 0, SizeMult[i], 0, 0, 0, 0, 0, 0);
      hEffAccInClasses[i]->Draw("same");
      legClasses->AddEntry(hEffAccInClasses[i], SmoltBis[i], "pef");

      padLow->cd();
      if(i!=0){
        hEffAccInClassesRatio[i] = (TH1D*)hEffAccInClasses[i]->Clone(Form("hEffAccInClassesRatio_%i", i));
        hEffAccInClassesRatio[i]->Divide(hEffAccInClasses[0]);
        StyleHisto(hEffAccInClassesRatio[i], 0.1 , 1.9, color[i], MarkerMult[i], sPt, "Ratio to 0-100%", "", 0, 0, 0, 0, 0, SizeMult[i], 0, 0, 0, 0, 0, 0);
        hEffAccInClassesRatio[i]->Draw("same");
      }
    }

    padUp->cd();
    legendTitleMC->Draw("same");
    legClasses->Draw("same");
    canvasCompPlot->Write();
  }
  // Signal Loss // (global so we can add MB later)
  TCanvas* canvasCompPlot_SignalLoss = new TCanvas("canvasCompPlot_SignalLoss", "canvasCompPlot_SignalLoss", 0, 70, 620, 850);
  StyleCanvas(canvasCompPlot_SignalLoss, 0.15, 0.05, 0.05, 0.15);
  TPad *padUp_SignalLoss = new TPad("padUp_SignalLoss", "padUp_SignalLoss", 0, 0.36, 1, 1);
  TPad *padLow_SignalLoss = new TPad("padLow_SignalLoss", "padLow_SignalLoss", 0, 0, 1, 0.36);
  StylePad(padUp_SignalLoss, 0.15, 0.05, 0.05, 0.);
  StylePad(padLow_SignalLoss, 0.15, 0.05, 0.02, 0.2);
  canvasCompPlot_SignalLoss->cd();
  padUp_SignalLoss->Draw();
  padLow_SignalLoss->Draw();

  // Legend
  TLegend *legClasses_SignalLoss = new TLegend(0.186, 0.0625, 0.893, 0.3125);
  legClasses_SignalLoss->SetHeader("FT0M Multiplicity Percentile");
  legClasses_SignalLoss->SetNColumns(3);
  TLegendEntry *legEffAccInClassesHeader = (TLegendEntry *)legClasses_SignalLoss->GetListOfPrimitives()->First();
  legEffAccInClassesHeader->SetTextSize(0.04);
  StyleLegend(legClasses_SignalLoss, 0.0, 0.0);

  TLegend *legendTitleMC_SignalLoss = new TLegend(0.608, 0.636, 0.907, 0.893);
  StyleLegend(legendTitleMC_SignalLoss, 0.0, 0.0);
  legendTitleMC_SignalLoss->SetTextAlign(33);
  legendTitleMC_SignalLoss->SetTextSize(0.05);
  legendTitleMC_SignalLoss->SetTextFont(42);
  legendTitleMC_SignalLoss->SetLineColorAlpha(0.,0.);
  legendTitleMC_SignalLoss->SetFillColorAlpha(0.,0.);
  legendTitleMC_SignalLoss->AddEntry("", "#bf{ALICE Work In Progress}", "");
  legendTitleMC_SignalLoss->AddEntry("", "pp, #sqrt{#it{s}} = 13.6 TeV", "");
  legendTitleMC_SignalLoss->AddEntry("", particleSymnbols[nParticle] + ", |y| < 0.5", "");
  legendTitleMC_SignalLoss->AddEntry("", "LHC24b1b", "");

  // Dummy histograms
  TH1F *hDummyUp_SignalLoss = new TH1F("hDummyUp_SignalLoss", "hDummyUp_SignalLoss", 10000, 0, binpt[numPtBins] + 0.5);
  // Up
  for (Int_t i = 1; i <= hDummyUp_SignalLoss->GetNbinsX(); i++)
    hDummyUp_SignalLoss->SetBinContent(i, 1000);
  canvasCompPlot_SignalLoss->cd();
  padUp_SignalLoss->cd();
  StyleHisto(hDummyUp_SignalLoss, -0.8, 2.0, 1, 1, sPt, "Signal Loss", "", 0, 0, 0, 1.5, 1.0, 0, 0.0, 0.05, 0.0, 0.035, hDummyUp_SignalLoss->GetXaxis()->GetLabelOffset(), 0.005);
  SetTickLength(hDummyUp_SignalLoss, 0.025, 0.03);
  TAxis *axisDummyUpY_SignalLoss = hDummyUp_SignalLoss->GetYaxis();
  axisDummyUpY_SignalLoss->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
  //axisDummyUpY_SignalLoss->ChangeLabel(2, -1, -1, -1, -1, -1, " ");
  hDummyUp_SignalLoss->GetXaxis()->SetRangeUser(0, binpt[numPtBins] + 0.5);
  //padUp_SignalLoss->SetLogy();
  hDummyUp_SignalLoss->Draw("same");
  // Low
  TH1F *hDummyLow_SignalLoss = new TH1F("hDummyLow_SignalLoss", "hDummyLow_SignalLoss", 10000, 0, binpt[numPtBins] + 0.5);
  for (Int_t i = 1; i <= hDummyLow_SignalLoss->GetNbinsX(); i++)
    hDummyLow_SignalLoss->SetBinContent(i, 1);
  padLow_SignalLoss->cd();
  StyleHisto(hDummyLow_SignalLoss, 0.1, 1.2, 1, 1, sPt, "Ratio to 0-100%", "", 0, 0, 0, 1.0, 0.7, 0, 0.08, 0.08, 0.08, 0.07, hDummyLow_SignalLoss->GetXaxis()->GetLabelOffset(), 0.01);
  SetTickLength(hDummyLow_SignalLoss, 0.025, 0.03);
  TAxis *axisDummyLowY_SignalLoss = hDummyLow_SignalLoss->GetYaxis();
  //padLow_SignalLoss->SetGridy(); // Only grid on x-axis
  axisDummyLowY_SignalLoss->ChangeLabel(-1, -1, -1, -1, -1, -1, " "); // last 
  axisDummyLowY_SignalLoss->ChangeLabel(1, -1, -1, -1, -1, -1, " "); // first
  hDummyLow_SignalLoss->GetXaxis()->SetRangeUser(0, binpt[numPtBins] + 0.5);
  hDummyLow_SignalLoss->GetYaxis()->CenterTitle();
  hDummyLow_SignalLoss->Draw("same");

  //Plotting
  for (Int_t i = 1; i < numMultEff + 1; i++) {
    padUp_SignalLoss->cd();
    StyleHisto(hSignalLossInClasses[i], 0, hSignalLossInClasses[i]->GetMaximum()*1e2, color[i], MarkerMult[i], "", hSignalLossInClasses[i]->GetYaxis()->GetTitle(), "", 0, 0, 0, 0, 0, SizeMult[i], 0, 0, 0, 0, 0, 0);
    hSignalLossInClasses[i]->Draw("same");
    legClasses_SignalLoss->AddEntry(hSignalLossInClasses[i], SmoltBis[i], "pef");
  }

  padUp_SignalLoss->cd();
  legendTitleMC_SignalLoss->Draw("same");
  legClasses_SignalLoss->Draw("same");

  // Event Correction //
  TDirectory* dirEventCorr = outputfile->mkdir("effEvent");
  dirEventCorr->cd();

  TH2F* hCentFT0M_rec = (TH2F*)fileMCin->Get("lf-cascqaanalysis/hCentFT0M_rec");
  TH2F* hCentFT0M_genMC = (TH2F*)fileMCin->Get("lf-cascqaanalysis/hCentFT0M_genMC");
  hCentFT0M_rec->GetYaxis()->SetRange(2 + inel, 3); // INEL>0
  hCentFT0M_genMC->GetYaxis()->SetRange(2 + inel, 3); // INEL>0
  hCentFT0M_rec->GetXaxis()->SetRangeUser(0., 100.);
  hCentFT0M_genMC->GetXaxis()->SetRange(0., 100.);

  TH2F* hCentFT0M_rec_to_rebin = (TH2F*)hCentFT0M_rec->Clone("hCentFT0M_rec_to_rebin");
  TH2F* hCentFT0M_genMC_to_rebin = (TH2F*)hCentFT0M_genMC->Clone("hCentFT0M_genMC_to_rebin");

  TH2F* hCentFT0M_rec_to_finest = (TH2F*)hCentFT0M_rec->Clone("hCentFT0M_rec_to_finest");
  TH2F* hCentFT0M_genMC_to_finest = (TH2F*)hCentFT0M_genMC->Clone("hCentFT0M_genMC_to_finest");

  TH1D* eventCorr;
  TH1D* eventCorr_finest;

  TH1D* hCentFT0M_rec_1D_rebinned = (hCentFT0M_rec_to_rebin->ProjectionX());
  TH1D* hCentFT0M_genMC_1D_rebinned = (hCentFT0M_genMC_to_rebin->ProjectionX());

  TH1D* hCentFT0M_rec_1D_finest = (hCentFT0M_rec_to_finest->ProjectionX());
  TH1D* hCentFT0M_genMC_1D_finest = (hCentFT0M_genMC_to_finest->ProjectionX());

  hCentFT0M_rec_1D_rebinned = (TH1D*)hCentFT0M_rec_1D_rebinned->Rebin(numMult, "hCentFT0M_rec_1D_rebinned", multiplicityPerc);
  hCentFT0M_genMC_1D_rebinned = (TH1D*)hCentFT0M_genMC_1D_rebinned->Rebin(numMult, "hCentFT0M_genMC_1D_rebinned", multiplicityPerc);

  hCentFT0M_rec_1D_finest = (TH1D*)hCentFT0M_rec_1D_finest->Rebin(numMultFinest, "hCentFT0M_rec_1D_finest", multiplicityFinest);
  hCentFT0M_genMC_1D_finest = (TH1D*)hCentFT0M_genMC_1D_finest->Rebin(numMultFinest, "hCentFT0M_genMC_1D_finest", multiplicityFinest);

  // TCanvas *canvasTest = new TCanvas("canvasTest","canvasTest", 800, 600);
  // canvasTest->cd();
  // StyleCanvas(canvasTest, 0.14, 0.05, 0.11, 0.15);
  // hCentFT0M_rec_1D_rebinned->Draw();
  // canvasTest->Draw();

  {
    TLegend *LegendTitleEventCorr = new TLegend(0.497, 0.626, 0.918, 0.833);
    LegendTitleEventCorr->SetFillStyle(0);
    LegendTitleEventCorr->SetTextAlign(33);
    LegendTitleEventCorr->SetTextSize(0.04);
    LegendTitleEventCorr->SetTextFont(42);
    LegendTitleEventCorr->SetLineColorAlpha(0.,0.);
    LegendTitleEventCorr->SetFillColorAlpha(0.,0.);
    LegendTitleEventCorr->SetBorderSize(0.);
    LegendTitleEventCorr->AddEntry("", "#bf{ALICE Work In Progress}", "");
    LegendTitleEventCorr->AddEntry("", "Event Correction", "");
    LegendTitleEventCorr->AddEntry("", "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    LegendTitleEventCorr->AddEntry("", "LHC24b1b, rec. true INEL > 0", "");

    TCanvas *canvasEventCorr = new TCanvas("canvasEventCorr","canvasEventCorr", 800, 600);
    canvasEventCorr->cd();
    StyleCanvas(canvasEventCorr, 0.14, 0.05, 0.11, 0.15);

    Double_t multiplicityPercEventFacotr[numMult + 1 + 1] = {0}; 
    for (int i = 0; i < numMult + 1; ++i) {
        multiplicityPercEventFacotr[i] = multiplicityPerc[i];
    }
    multiplicityPercEventFacotr[numMult + 1] = 150;
    eventCorr = new TH1D("eventCorr", "eventCorr", numMult + 1, multiplicityPercEventFacotr);
    std::cout << "\n************* Event Corrections in Classes *************\n";
    long int totalRec = 0;
    long int totalGen = 0;
    for (int i = 1; i <= eventCorr->GetNbinsX(); ++i) {
      Double_t k = 0;
      Double_t n = 0;
      if (i == eventCorr->GetNbinsX()) { // MB
        k = totalRec;
        n = totalGen;
      } else {
        k = hCentFT0M_rec_1D_rebinned->GetBinContent(i);
        n = hCentFT0M_genMC_1D_rebinned->GetBinContent(i);
        totalRec += k;
        totalGen += n;
      }
      double ratio = (n != 0) ? (k / n) : 0.0;
      Double_t errorInBin = sqrt(((Double_t)k+1)*((Double_t)k+2)/(n+2)/(n+3) - pow((Double_t)(k+1),2)/pow(n+2,2));
      std::cout << "mult. %: " << multiplicityPercEventFacotr[i-1] << " -- " << multiplicityPercEventFacotr[i] <<  " k = " << k << " n = " << n << " ratio: " << ratio << " error: " << errorInBin << std::endl;
      eventCorr->SetBinContent(i, ratio);
      eventCorr->SetBinError(i, errorInBin);
    }
    // MB efficiency is a weighted by the class width
    Double_t totWeights = 0;
    Double_t MBeventCorr = 0;
    Double_t errorMBeventCorr = 0;
    for (int j = 1; j <= eventCorr->GetNbinsX() - 1; ++j) {
      MBeventCorr += eventCorr->GetBinContent(j)*hCentFT0M_rec_data_1D_anal->GetBinContent(j);
      totWeights += hCentFT0M_rec_data_1D_anal->GetBinContent(j);
      errorMBeventCorr += pow(eventCorr->GetBinError(j)*hCentFT0M_rec_data_1D_anal->GetBinContent(j), 2);
      std::cout << "number of events in class " << j << " is: " << hCentFT0M_rec_data_1D_anal->GetBinContent(j) << std::endl;
      std::cout << "error is: " << errorMBeventCorr << std::endl;
    }
    MBeventCorr = MBeventCorr / totWeights;
    std::cout << "total weights: " << totWeights << std::endl;
    eventCorr->SetBinContent(eventCorr->GetNbinsX(), MBeventCorr);
    eventCorr->SetBinError(eventCorr->GetNbinsX(), sqrt(errorMBeventCorr) / totWeights);

    // MB efficiency finest multiplicity bins
    eventCorr_finest = new TH1D("eventCorr_finest", "eventCorr_finest", numMultFinest, multiplicityFinest);
    for (int i = 1; i <= hCentFT0M_rec_1D_finest->GetNbinsX(); ++i) {
      Double_t k = hCentFT0M_rec_1D_finest->GetBinContent(i);
      Double_t n = hCentFT0M_genMC_1D_finest->GetBinContent(i);
      double ratio = (n != 0) ? (k / n) : 0.0;
      Double_t errorInBin = sqrt(((Double_t)k+1)*((Double_t)k+2)/(n+2)/(n+3) - pow((Double_t)(k+1),2)/pow(n+2,2));
      std::cout << "Finest bin number: " << i << " k = " << k << " n = " << n << " ratio: " << ratio << " error: " << errorInBin << std::endl;
      eventCorr_finest->SetBinContent(i, ratio);
      eventCorr_finest->SetBinError(i, errorInBin);
    }
    totWeights = 0;
    MBeventCorr = 0;
    errorMBeventCorr = 0;
    for (int j = 1; j <= eventCorr_finest->GetNbinsX(); j++) {
      // MBeventCorr += eventCorr_finest->GetBinContent(j)*hCentFT0M_rec_data_1D_finest->GetBinContent(j);
      // totWeights += hCentFT0M_rec_data_1D_finest->GetBinContent(j);
      MBeventCorr += hCentFT0M_rec_data_1D_finest->GetBinContent(j);
      totWeights += hCentFT0M_rec_data_1D_finest->GetBinContent(j) / eventCorr_finest->GetBinContent(j);
      //errorMBeventCorr += pow(eventCorr_finest->GetBinError(j)*hCentFT0M_rec_data_1D_finest->GetBinContent(j), 2);
      errorMBeventCorr += pow(eventCorr_finest->GetBinError(j) / eventCorr_finest->GetBinContent(j), 2);
      std::cout << "number of events in finest class " << j << " is: " << hCentFT0M_rec_data_1D_finest->GetBinContent(j) << std::endl;
      std::cout << "error is: " << errorMBeventCorr << std::endl;
    }
    MBeventCorr = MBeventCorr / totWeights;
    std::cout << "total weights (finest): " << totWeights << std::endl;
    eventCorr->SetBinContent(eventCorr->GetNbinsX(), MBeventCorr);
    //eventCorr->SetBinError(eventCorr->GetNbinsX(), sqrt(errorMBeventCorr) / totWeights);
    eventCorr->SetBinError(eventCorr->GetNbinsX(), sqrt(errorMBeventCorr) * MBeventCorr);

    eventCorr->Draw();

    TAxis *axisEventCorr = eventCorr->GetXaxis();
    axisEventCorr->ChangeLabel(-1, -1, -1, -1, -1, -1, " "); // last
    axisEventCorr->ChangeLabel(-2, -1, -1, -1, -1, -1, " "); // pre-last
    StyleHisto(eventCorr, 0., 2. * eventCorr->GetBinContent(eventCorr->GetMaximumBin()), color[0], MarkerMult[0], "FT0M Multiplicity percentile", "Event Factor", "", 0, 0, 0, 1.0, 1.25, 1.1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
    TH1F* eventCorrClone = (TH1F*)eventCorr->Clone("eventCorrClone");
    StyleHisto(eventCorrClone, 0., 1.4 * eventCorrClone->GetBinContent(eventCorr->GetMaximumBin()), color[0], 0, "FT0M Multiplicity percentile", "Event Factor", "", 0, 0, 0, 1.0, 1.25, SizeMult[0], 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
    eventCorrClone->GetXaxis()->SetRangeUser(40, 150);
    eventCorr->Draw("same HIST");
    eventCorrClone->Draw("same TEXT90");
    TLatex labelMB;
    labelMB.SetTextSize(0.05);
    labelMB.SetTextAlign(12);
    labelMB.DrawLatex(104.9,.087,"Minimum Bias");

    LegendTitleEventCorr->Draw("same");
    canvasEventCorr->Write();
  }

  // Weights rescaling //
  {
    TCanvas *canvasEventCorrInverse = new TCanvas("canvasEventCorrInverse","canvasEventCorrInverse", 800, 600);
    canvasEventCorrInverse->cd();
    StyleCanvas(canvasEventCorrInverse, 0.14, 0.05, 0.11, 0.15);
    TH1D* eventCorrInverse = new TH1D("eventCorrInverse", "eventCorrInverse", numMult, multiplicityPerc);
    Double_t totEventCorrInverse = 0;
    for (int j = 1; j <= eventCorr->GetNbinsX() - 1; ++j) {
      eventCorrInverse->SetBinContent(j, 1. / eventCorr->GetBinContent(j));
      totEventCorrInverse += eventCorrInverse->GetBinContent(j);
    }
    StyleHisto(eventCorrInverse, 0., 1.2 * eventCorrInverse->GetBinContent(eventCorrInverse->GetMaximumBin()), color[0], MarkerMult[0], "FT0M Multiplicity percentile", "1. / Event Factor", "", 0, 0, 0, 1.0, 1.25, 1.1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
    eventCorrInverse->Draw();

    TCanvas *canvasEventCorrCumulative = new TCanvas("canvasEventCorrCumulative","canvasEventCorrCumulative", 800, 600);
    canvasEventCorrCumulative->cd();
    StyleCanvas(canvasEventCorrCumulative, 0.14, 0.05, 0.11, 0.15);
    TH1D* eventCorrCumulative = new TH1D("eventCorrCumulative", "eventCorrCumulative", numMult, multiplicityPerc);
    Double_t cumulativeSum = 0;
    for (int j = 1; j <= eventCorrCumulative->GetNbinsX(); ++j) {
      cumulativeSum += eventCorrInverse->GetBinContent(j);
      eventCorrCumulative->SetBinContent(j, cumulativeSum / totEventCorrInverse);
    }
    StyleHisto(eventCorrCumulative, 0., 1.1, color[0], MarkerMult[0], "FT0M Multiplicity percentile", "Cumulative", "", 0, 0, 0, 1.0, 1.25, 1.1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
    eventCorrCumulative->Draw("same");
    eventCorrCumulative->Write();
  }

  // Correct spectra in classes for Efficiency x Acceptance, Signal Loss and Event Factor //
  // Create output files
  TFile *outputfileCorrYield[numMult + 1];
  for (Int_t iFile = 0; iFile < numMult + 1; iFile++) {
    if (iFile == 0) {
      outputfilePath = outputDir + "/yieldsOutEffCorr/" + "yield_" + particleNames[nParticle] + "_MB_inel" + inel + yieldPostFix + ".root";
    } else {
      outputfilePath = outputDir + "/yieldsOutEffCorr/" + "yield_" + particleNames[nParticle] + "_" + multiplicityPerc[iFile - 1] + "-" + multiplicityPerc[iFile] + "_inel" + inel + yieldPostFix + ".root";
    }
    outputfileCorrYield[iFile] = new TFile(outputfilePath, "RECREATE");
  }

  TDirectory* baseHistsOutDir[numMult + 1];
  TDirectory* corrYieldOutDir[numMult + 1];
  Double_t totYieldCorrected[numPtBins][numMult + 1]; // 0 - MB
  for (Int_t iClass = 0; iClass < numMult + 1; iClass++) {
    for (Int_t b = 0; b < numPtBins; b++) {
      totYieldCorrected[b][iClass] = 0;
    }
  }

  TH1F* hSignalLoss;
  TH1F* hEffAcc;

  for (Int_t iFile = 1; iFile < numMult + 1; iFile++) {
    //cout << "iFile = " << iFile << std::endl;
    cout << "Correction in mult. %: " << multiplicityPerc[iFile - 1] << " -- " << multiplicityPerc[iFile] << endl; 
    // Choose Eff. x Acc. and Signal Loss histogram
    for(Int_t i = 1; i < numMultEff + 1; i++){
      if(multiplicityPercEff[i] > multiplicityPerc[iFile] - binPrecision){
        hSignalLoss = (TH1F *)hSignalLossInClasses[i]->Clone(Form("SignalLoss_%i",iFile));
        hEffAcc = (TH1F *)hEffAccInClasses[0]->Clone(Form("EffAcc_%i",iFile)); // 0 (MB) or i
        StyleHisto(hSignalLoss, 0, 1.2 * hSignalLoss->GetBinContent(hSignalLoss->GetMaximumBin()), kBlack, 20, sPt, "Signal Loss", "", 0, 0, 0, 1.0, 1.25, 1,  0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
        StyleHisto(hEffAcc, 0, 1.2 * hEffAcc->GetBinContent(hEffAcc->GetMaximumBin()), kBlack, 20, sPt, "Eff. x Acc.", "", 0, 0, 0, 1.0, 1.25, 1,  0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
        cout << "i = " << i << std::endl;
        break;
      }
    }

    TF1 *fitExp = new TF1("fitExp", "[0]*exp([1]*x-[2]) + [3]", binpt[0], binpt[numPtBins]);
    //fitExp->SetParameters(-0.4, -5.7, -0.1, 1);
    fitExp->SetParameters(-0.4, -0.1, -0.3, 1);
    TF1 *fitPol1 = new TF1("fitPol1", "[0]*x + [1]", binpt[0], binpt[numPtBins]);
    fitPol1->SetParameters(0.01, 0.9);
    TF1 *fitInverse = new TF1("fitInverse", "[0]/(x + [1]) + [2]", binpt[0], binpt[numPtBins]);
    fitInverse->SetParameters(0.2, 0.3, 0.9);

    if (multiplicityPerc[iFile] > 1e-3) { // perform exp fit if mult.% > x - 1e-3
      hSignalLoss->Fit(fitExp, "QR0"); // Fit with exponent
      padUp_SignalLoss->cd();
      fitExp->Draw("same");
      cout << "Chi2/NDF " << fitExp->GetChisquare() / fitExp->GetNDF() << endl;

      // hSignalLoss->Fit(fitInverse, "R0"); // Fit with inverse
      // padUp_SignalLoss->cd();
      // fitInverse->Draw("same");
      // cout << "Chi2/NDF " << fitInverse->GetChisquare() / fitInverse->GetNDF() << endl;
    } else {
      hSignalLoss->Fit(fitPol1, "QR0"); // Fit wil pol1
      padUp_SignalLoss->cd();
      fitPol1->Draw("same");
      //cout << "Chi2/NDF " << fitPol1->GetChisquare() / fitPol1->GetNDF() << endl;      
    }

    for (Int_t b = 1; b <= numPtBins; b++) {
      Double_t relErr = 0;
      Double_t relErrEffAcc = 0;
      // Eff. x Acc.
      relErr = hYieldCorrected[iFile]->GetBinError(b) / hYieldCorrected[iFile]->GetBinContent(b);
      relErrEffAcc = hEffAcc->GetBinError(hEffAcc->FindBin(hYieldCorrected[iFile]->GetBinCenter(b))) / hEffAcc->GetBinContent(hEffAcc->FindBin(hYieldCorrected[iFile]->GetBinCenter(b)));
      hYieldCorrected[iFile]->SetBinContent(b, hYieldCorrected[iFile]->GetBinContent(b) / hEffAcc->GetBinContent(hEffAcc->FindBin(hYieldCorrected[iFile]->GetBinCenter(b))));
      hYieldCorrected[iFile]->SetBinError(b, sqrt(pow(relErr, 2) + pow(relErrEffAcc, 2)) * hYieldCorrected[iFile]->GetBinContent(b));
      //cout << "eff. x acc.: " << hEffAcc->GetBinContent(hEffAcc->FindBin(hYield[iFile]->GetBinCenter(b))) << " relErrEffAcc: " << relErrEffAcc << endl;
      //std::cout << "yield after eff. x acc.: " << hYieldCorrected[iFile]->GetBinContent(b) << std::endl;
      // Signal Loss
      relErr = hYieldCorrected[iFile]->GetBinError(b) / hYieldCorrected[iFile]->GetBinContent(b);
      Double_t relErrSigLoss = 0;
      Double_t sigLoss = 0;

      Double_t xValue = hSignalLoss->GetXaxis()->GetBinCenter(b);

      if (multiplicityPerc[iFile] > 1e-3) {
        // Take values from histogram
        sigLoss = hSignalLoss->GetBinContent(b);
        relErrSigLoss = hSignalLoss->GetBinError(hSignalLoss->FindBin(hYieldCorrected[iFile]->GetXaxis()->GetBinCenter(b))) / sigLoss;
        // Take values from fit
        // Exp
        // sigLoss = fitExp->Eval(xValue);

        // Double_t withoutConst = sigLoss - fitExp->GetParameter(3);
        // Double_t sig0 = (withoutConst / fitExp->GetParameter(0)) * fitExp->GetParError(0);
        // Double_t sig1 = withoutConst * xValue * fitExp->GetParError(1);
        // Double_t sig2 = withoutConst * fitExp->GetParError(2);
        // Double_t sig3 = fitExp->GetParError(3);

        // relErrSigLoss = sqrt(pow(sig0,2) + pow(sig1,2) + pow(sig2,2) + pow(sig3,2)) / sigLoss;
        // Invertse
        // sigLoss = fitInverse->Eval(xValue);

        // Double_t withoutA = (sigLoss - fitInverse->GetParameter(2)) / fitInverse->GetParameter(0);
        // Double_t sig0 = withoutA * fitInverse->GetParError(0);
        // Double_t sig1 = withoutA * withoutA * fitInverse->GetParameter(0) * fitInverse->GetParError(1);
        // Double_t sig2 = fitInverse->GetParError(2);

        // relErrSigLoss = sqrt(pow(sig0,2) + pow(sig1,2) + pow(sig2,2)) / sigLoss;
      } else {
        // Take values from histogram
        sigLoss = hSignalLoss->GetBinContent(b);
        relErrSigLoss = hSignalLoss->GetBinError(hSignalLoss->FindBin(hYieldCorrected[iFile]->GetXaxis()->GetBinCenter(b))) / sigLoss;
        // Take values from fit
        // sigLoss = fitPol1->Eval(xValue);

        // Double_t sig0 = xValue * fitPol1->GetParError(0);
        // Double_t sig1 = fitPol1->GetParError(1);

        // relErrSigLoss = sqrt(pow(sig0,2) + pow(sig1,2)) / sigLoss;
      }

      hYieldCorrected[iFile]->SetBinContent(b, hYieldCorrected[iFile]->GetBinContent(b) / sigLoss);
      //std::cout << "yield after signal loss: " << hYieldCorrected[iFile]->GetBinContent(b) << std::endl;
      hYieldCorrected[iFile]->SetBinError(b, sqrt(pow(relErr, 2) + pow(relErrSigLoss, 2)) * hYieldCorrected[iFile]->GetBinContent(b));
      //cout << "sigLoss: " << sigLoss << " relErrSigLoss: " << relErrSigLoss<< endl;

      // Event splitting/loss correction
      relErr = hYieldCorrected[iFile]->GetBinError(b) / hYieldCorrected[iFile]->GetBinContent(b);
      Double_t relErrEventFactor = eventCorr->GetBinError(iFile)/eventCorr->GetBinContent(iFile);
      hYieldCorrected[iFile]->SetBinContent(b, hYieldCorrected[iFile]->GetBinContent(b) * eventCorr->GetBinContent(iFile));
      //std::cout << "yield after event corr. in bin " << b << " is: " << hYieldCorrected[iFile]->GetBinContent(b) << std::endl;
      hYieldCorrected[iFile]->SetBinError(b, sqrt(pow(relErr, 2) + pow(relErrEventFactor, 2)) * hYieldCorrected[iFile]->GetBinContent(b));
      //cout << "event correction: " << eventCorr[iFile] << " relErrEventFactor: " << relErrEventFactor<< endl;

      // Compute "integrated yiled"
      totYieldCorrected[b-1][iFile] += hYieldCorrected[iFile]->GetBinContent(b);
    }
    corrYieldOutDir[iFile] = outputfileCorrYield[iFile]->mkdir("effCorrYield");
    corrYieldOutDir[iFile]->cd();
    hYieldCorrected[iFile]->Write();

    baseHistsOutDir[iFile] = outputfileCorrYield[iFile]->mkdir("baseHists");
    baseHistsOutDir[iFile]->cd();
    hYield[iFile]->Write();
    eventCorr->Write();
    hSignalLoss->Write();
    //canvasSignalLossFit->Write();
    hEffAcc->Write();
  }

  // Compute MB Signal Loss //
  Double_t cascadesInClass[numPtBins][numMultEff]; // cascades in class per percentile of mult. class of efficiencies
  for (Int_t iClass = 0; iClass < numMultEff; iClass++) {
    for (Int_t b = 0; b < numPtBins; b++) {
      cascadesInClass[b][iClass] = 0;
    }
  }
  for (Int_t iFile = 1; iFile < numMult + 1; iFile++) {
    Int_t multBin = -1;
    for (Int_t b = 1; b <= numPtBins; b++) {
      for(Int_t i = 1; i < numMultEff + 1; i++){
        if(multiplicityPercEff[i] > multiplicityPerc[iFile] - binPrecision){
          multBin = i - 1;
          break;
        }
      }
      cascadesInClass[b-1][multBin] += totYieldCorrected[b-1][iFile] * hCentFT0M_rec_data_1D_anal->GetBinContent(iFile) / eventCorr->GetBinContent(iFile) * hYieldCorrected[iFile]->GetBinWidth(b);
      //std::cout << "totYieldCorrected in bin " << b-1 << " is: " << totYieldCorrected[b-1][iFile] << " " << multiplicityPerc[iFile-1] << " -- " << multiplicityPerc[iFile] <<  " number of cascades: " << totYieldCorrected[b-1][iFile] * hCentFT0M_rec_data_1D_anal->GetBinContent(iFile) / eventCorr->GetBinContent(iFile) * hYieldCorrected[iFile]->GetBinWidth(b) << " pt width: " << hYieldCorrected[iFile]->GetBinWidth(b) << std::endl;
    }
    //std::cout << "multBin: " <<  multBin << std::endl;
  }

  // Compute MB Signal Loss //
  hSignalLossInClasses[0]->Reset();
  for (Int_t b = 1; b <= numPtBins; b++) {
    Double_t signalLossInBin = 0; 
    Double_t sumWeights = 0;
    Double_t errorInBin = 0;
    for (Int_t i = 0; i < numMultEff; i++) {
      signalLossInBin += hSignalLossInClasses[i+1]->GetBinContent(b)*cascadesInClass[b-1][i];
      sumWeights += cascadesInClass[b-1][i];
      errorInBin += pow(hSignalLossInClasses[i+1]->GetBinError(b)*cascadesInClass[b-1][i], 2);
      //std::cout << "Total number of cascades in pt bin " << b-1 << " is: " <<  cascadesInClass[b-1][i] << " in class: " << i << std::endl;
    }
    //std::cout << "Total number of cascades in pt bin " << b-1 << " is: " <<  sumWeights << std::endl;
    hSignalLossInClasses[0]->SetBinContent(b, signalLossInBin / sumWeights);
    hSignalLossInClasses[0]->SetBinError(b, sqrt(errorInBin) / sumWeights);
  }

  // Plot MB Signal loss efficiency together with the ones in classes and compute ratios to MB
  padUp_SignalLoss->cd();
  StyleHisto(hSignalLossInClasses[0], 0, hSignalLossInClasses[0]->GetMaximum()*1e2, color[0], MarkerMult[0], "", hSignalLossInClasses[0]->GetYaxis()->GetTitle(), "", 0, 0, 0, 0, 0, SizeMult[0], 0, 0, 0, 0, 0, 0);
  hSignalLossInClasses[0]->Draw("same");
  legClasses_SignalLoss->AddEntry(hSignalLossInClasses[0], SmoltBis[0], "pef");

  for (Int_t i = 1; i < numMultEff + 1; i++) {
    padLow_SignalLoss->cd();
    hSignalLossInClassesRatio[i] = (TH1D*)hSignalLossInClasses[i]->Clone(Form("hEffAccInClassesRatio_%i", i));
    hSignalLossInClassesRatio[i]->Divide(hSignalLossInClasses[0]);
    StyleHisto(hSignalLossInClassesRatio[i], 0.1 , 1.9, color[i], MarkerMult[i], sPt, "Ratio to 0-100%", "", 0, 0, 0, 0, 0, SizeMult[i], 0, 0, 0, 0, 0, 0);
    hSignalLossInClassesRatio[i]->Draw("same");
  }

  // Corrected yield in MB //
  {
    hSignalLoss = (TH1F*)hSignalLossInClasses[0]->Clone("hSignalLossMB");
    StyleHisto(hSignalLoss, 0, 1.2 * hSignalLoss->GetBinContent(hSignalLoss->GetMaximumBin()), kBlack, 20, sPt, "Signal Loss", "", 0, 0, 0, 1.0, 1.25, 1,  0.05, 0.05, 0.05, 0.05, 0.007, 0.007);

    // TF1 *fitPol1 = new TF1("fitPol1", "[0]*x + [1]", binpt[0], binpt[numPtBins]);
    // fitPol1->SetParameters(0.01, 0.9);
    // padUp_SignalLoss->cd();
    // hSignalLoss->Fit(fitPol1, "QR0"); // Fit wil pol1
    // fitPol1->Draw("same");
    //cout << "Chi2/NDF " << fitPol1->GetChisquare() / fitPol1->GetNDF() << endl;

    TF1 *fitExp = new TF1("fitExp", "[0]*exp([1]*x-[2]) + [3]", binpt[0], binpt[numPtBins]);
    fitExp->SetParameters(-0.4, -0.1, -0.3, 1);
    padUp_SignalLoss->cd();
    hSignalLoss->Fit(fitExp, "QR0"); // Fit low-mult. with exponent
    fitExp->Draw("same");
    //cout << "Chi2/NDF " << fitExp->GetChisquare() / fitExp->GetNDF() << endl;

    // TF1 *fitInverse = new TF1("fitInverse", "[0]/(x + [1]) + [2]", binpt[0], binpt[numPtBins]);
    // fitInverse->SetParameters(0.2, 0.3, 0.9);
    // padUp_SignalLoss->cd();
    // hSignalLoss->Fit(fitInverse, "R0"); // Fit low-mult. with exponent
    // fitInverse->Draw("same");
    //cout << "Chi2/NDF " << fitInverse->GetChisquare() / fitInverse->GetNDF() << endl;

    // Compute it like a separate class //
    for (Int_t b = 1; b <= hYield[0]->GetNbinsX(); b++) {
      Double_t relErr = 0;
      Double_t relErrEffAcc = 0;
      // Eff. x Acc.
      relErr = hYieldCorrected[0]->GetBinError(b) / hYieldCorrected[0]->GetBinContent(b);
      relErrEffAcc = hEffAcc->GetBinError(hEffAcc->FindBin(hYieldCorrected[0]->GetBinCenter(b))) / hEffAcc->GetBinContent(hEffAcc->FindBin(hYieldCorrected[0]->GetBinCenter(b)));
      hYieldCorrected[0]->SetBinContent(b, hYieldCorrected[0]->GetBinContent(b) / hEffAcc->GetBinContent(hEffAcc->FindBin(hYieldCorrected[0]->GetBinCenter(b))));
      hYieldCorrected[0]->SetBinError(b, sqrt(pow(relErr, 2) + pow(relErrEffAcc, 2)) * hYieldCorrected[0]->GetBinContent(b));
      //cout << "eff. x acc.: " << hEffAcc->GetBinContent(hEffAcc->FindBin(hYield[0]->GetBinCenter(b))) << " relErrEffAcc: " << relErrEffAcc << endl;
      //std::cout << "yield after eff. x acc.: " << hYieldCorrected[0]->GetBinContent(b) << std::endl;
      // Signal Loss
      relErr = hYieldCorrected[0]->GetBinError(b) / hYieldCorrected[0]->GetBinContent(b);
      Double_t relErrSigLoss = 0;
      Double_t sigLoss = 0;

      Double_t xValue = hSignalLoss->GetXaxis()->GetBinCenter(b);
      // Take values from histogram
      sigLoss = hSignalLoss->GetBinContent(b);
      relErrSigLoss = hSignalLoss->GetBinError(hSignalLoss->FindBin(hYieldCorrected[0]->GetXaxis()->GetBinCenter(b))) / sigLoss;
      // Take values from fit
      // Pol1
      // sigLoss = fitPol1->Eval(xValue);

      // Double_t sig0 = xValue * fitPol1->GetParError(0);
      // Double_t sig1 = fitPol1->GetParError(1);

      // relErrSigLoss = sqrt(pow(sig0,2) + pow(sig1,2)) / sigLoss;
      // EXP
      // sigLoss = fitExp->Eval(xValue);

      // Double_t withoutConst = sigLoss - fitExp->GetParameter(3);
      // Double_t sig0 = (withoutConst / fitExp->GetParameter(0)) * fitExp->GetParError(0);
      // Double_t sig1 = withoutConst * xValue * fitExp->GetParError(1);
      // Double_t sig2 = withoutConst * fitExp->GetParError(2);
      // Double_t sig3 = fitExp->GetParError(3);

      // relErrSigLoss = sqrt(pow(sig0,2) + pow(sig1,2) + pow(sig2,2) + pow(sig3,2)) / sigLoss;
      // Invertse
      // sigLoss = fitInverse->Eval(xValue);

      // Double_t withoutA = (sigLoss - fitInverse->GetParameter(2)) / fitInverse->GetParameter(0);
      // Double_t sig0 = withoutA * fitInverse->GetParError(0);
      // Double_t sig1 = withoutA * withoutA * fitInverse->GetParameter(0) * fitInverse->GetParError(1);
      // Double_t sig2 = fitInverse->GetParError(2);

      // relErrSigLoss = sqrt(pow(sig0,2) + pow(sig1,2) + pow(sig2,2)) / sigLoss;

      hYieldCorrected[0]->SetBinContent(b, hYieldCorrected[0]->GetBinContent(b) / sigLoss);
      //std::cout << "yield after signal loss: " << hYieldCorrected[0]->GetBinContent(b) << std::endl;
      hYieldCorrected[0]->SetBinError(b, sqrt(pow(relErr, 2) + pow(relErrSigLoss, 2)) * hYieldCorrected[0]->GetBinContent(b));
      //cout << "sigLoss: " << sigLoss << " relErrSigLoss: " << relErrSigLoss<< endl;

      // Event splitting/loss correction
      relErr = hYieldCorrected[0]->GetBinError(b) / hYieldCorrected[0]->GetBinContent(b);
      Double_t relErrEventFactor = eventCorr->GetBinError(eventCorr->GetNbinsX())/eventCorr->GetBinContent(eventCorr->GetNbinsX());
      hYieldCorrected[0]->SetBinContent(b, hYieldCorrected[0]->GetBinContent(b) * eventCorr->GetBinContent(eventCorr->GetNbinsX()));
      //std::cout << "yield after event corr.: " << hYieldCorrected[0]->GetBinContent(b) << std::endl;
      hYieldCorrected[0]->SetBinError(b, sqrt(pow(relErr, 2) + pow(relErrEventFactor, 2)) * hYieldCorrected[0]->GetBinContent(b));
      //cout << "event correction: " << eventCorr->GetBinContent(eventCorr->GetNbinsX()) << " relErrEventFactor: " << relErrEventFactor<< endl;
    }
    corrYieldOutDir[0] = outputfileCorrYield[0]->mkdir("effCorrYield");
    corrYieldOutDir[0]->cd();
    hYieldCorrected[0]->Write();

    baseHistsOutDir[0] = outputfileCorrYield[0]->mkdir("baseHists");
    baseHistsOutDir[0]->cd();
    hYield[0]->Write();
    eventCorr->Write();
    hSignalLoss->Write();
    hEffAcc->Write();
  }

  dirSignalLoss->cd();
  canvasCompPlot_SignalLoss->Write();

  delete[] binpt;

  // End of Code
  std::cout << "\x1B[1;32m"; // Set text color to green
  std::cout << "\n************* Eff. x Acc.,  Signal Loss and Event Factor are Obtained Successfully! *************\n";
  std::cout << "\x1B[0m"; // Reset text color
}

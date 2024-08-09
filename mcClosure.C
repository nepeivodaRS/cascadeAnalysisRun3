#include "help.h"

void mcClosure(const Int_t nParticle = 2, // 0-2 : xi, 3-5 : omega
               const Int_t inel = 0, // inel > N (0/1)
               const TString workingDir = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisPostSQM",
               const TString fileMC = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisPostSQM/24jul-lhc24b1b/AnalysisResults.root",
               const TString fileMCPP = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisPostSQM/run3_13tev/xi/lhc24b1b/2024-07-24/AnalysisResults.root",
               const TString postFix = "")
{
  // Start of Code
  std::cout << "\x1B[1;33m"; // Set text color to yellow
  std::cout << "\n************ Starting Making MC Closure ************\n";
  std::cout << "\x1B[0m"; // Reset text color

  TString outputfilePath = workingDir + "/mc-closures/" + particleNames[nParticle] + ".root";
  TFile *outputfile = new TFile(outputfilePath, "RECREATE");

  // Files with yields in all mult. classes + MB
  TFile* fileDataIn[numMult + 1];
  // MB
  fileDataIn[0] = TFile::Open(workingDir + "/yieldsOutEffCorr" +  "/yield_" + particleNames[nParticle] + "_MB_inel" + inel + postFix  + ".root");
  if (!fileDataIn[0] || fileDataIn[0]->IsZombie()) {
    std::cerr << "Error opening input data file for MB!" << std::endl;
    return;
  } else {
    cout << "file for MB yield is opened"<< std::endl;
  }
  // in mult. classes
  for (Int_t iFile = 1; iFile < numMult + 1; iFile++) {
    TString fileInPath = workingDir + "/yieldsOutEffCorr" + "/yield_" + particleNames[nParticle] + "_" + multiplicityPerc[iFile - 1] + "-" + multiplicityPerc[iFile] + "_inel" + inel + postFix + ".root";
    fileDataIn[iFile] = TFile::Open(fileInPath);
    if (!fileDataIn[iFile] || fileDataIn[iFile]->IsZombie()) {
      std::cerr << "Error opening input data file for mult. class: " << multiplicityPerc[iFile - 1] <<  " - " << multiplicityPerc[iFile] << std::endl;
      return;
    } else {
      cout << "file for mult. class: " << multiplicityPerc[iFile - 1] <<  " - " << multiplicityPerc[iFile] << " is opened"<< std::endl;
    }
  }

  TH1F* hYield[numMult + 1];
  TH1F* hYieldsRatio[numMult + 1];
  TH1F* hYieldsDenom;
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
    hYieldsRatio[iFile] = (TH1F *)hYield[iFile]->Clone(Form("YieldClone_%i_", iFile) + particleNames[nParticle]);

    if (iFile == 0) {
      hYieldsDenom = (TH1F *)hYield[0]->Clone("YieldCloneForDenom_" + particleNames[nParticle]);
    }
  }

  // Get Mult. % event distribution
  TFile* fileMCIn = TFile::Open(fileMC);
  if (!fileMCIn || fileMCIn->IsZombie()) {
      std::cerr << "Error opening input mc file!" << std::endl;
      return;
  }

  TH2F* hCentFT0M_genMC = (TH2F*)fileMCIn->Get("lf-cascqaanalysis/hCentFT0M_genMC");
  hCentFT0M_genMC->GetYaxis()->SetRange(2 + inel, 3);
  hCentFT0M_genMC->GetXaxis()->SetRangeUser(0., 100.);
  TH1D* hCentFT0M_genMC_1D = (hCentFT0M_genMC->ProjectionX());
  hCentFT0M_genMC_1D = (TH1D*)hCentFT0M_genMC_1D->Rebin(numMult, "hCentFT0M_genMC_1D_rebinned", multiplicityPerc);

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

  // Get gen. distributions //
  Double_t binPrecision = 1e-6;

  // File with gen. histograms (postprocessing output of MC)
  TFile* fileMCPPIn = TFile::Open(fileMCPP);
  if (!fileMCPPIn || fileMCPPIn->IsZombie()) {
      std::cerr << "Error opening input `fileMCPP` file!" << std::endl;
      return;
  }

  // Histos of generated cascades from generated events with accepted z vrtx + chosen event type (evSelFlag)
  TH3F* hCascadePlusTrue;
  TH3F* hCascadeMinusTrue;
  TH3F* hCascadeSumTrue;

  if(partType == 1){
    hCascadePlusTrue = (TH3F*)fileMCPPIn->Get("lf-cascpostprocessing/hPtOmegaPlusTrue");
    hCascadeMinusTrue = (TH3F*)fileMCPPIn->Get("lf-cascpostprocessing/hPtOmegaMinusTrue");
  } else {
    hCascadePlusTrue = (TH3F*)fileMCPPIn->Get("lf-cascpostprocessing/hPtXiPlusTrue");
    hCascadeMinusTrue = (TH3F*)fileMCPPIn->Get("lf-cascpostprocessing/hPtXiMinusTrue");
  }

  hCascadeSumTrue = (TH3F*)hCascadePlusTrue->Clone("SumTrue");
  hCascadeSumTrue->Add(hCascadeMinusTrue);
  hCascadeSumTrue->SetName("hPtSumTrue");

  TH1F* hGenYield[numMult + 1];

  TDirectory* dirSupport = outputfile->mkdir("supportHistos");
  Double_t leftRap = -0.5;
  Double_t rightRap = 0.5;


  for (Int_t i = 0; i < numMult + 1; i++) {
    std::cout << "Processing class number " << i << std::endl;
    dirSupport->cd();
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

    hGenYield[i] = static_cast<TH1F*>(hCascadeTrue->Project3D("x"));
    //Rebin to match analysis binning
    hGenYield[i] = (TH1F*)hGenYield[i]->Rebin(numPtBins, Form("hCascadeTrue_1D_rebinned_%i", i), binpt);
    long int nEvents = 0;
    if(i!=0){
      nEvents = hCentFT0M_genMC_1D->GetBinContent(i);
    } else {
      for (Int_t j = 1; j <= hCentFT0M_genMC_1D->GetNbinsX(); j++) {
        nEvents += hCentFT0M_genMC_1D->GetBinContent(j); // ne rabotaet verno
      }
    }

    // normalize MC gen.
    for (Int_t b = 1; b <= hGenYield[i]->GetNbinsX(); b++) {
      Int_t k = hGenYield[i]->GetBinContent(b);
      Int_t n = nEvents;
      Double_t errorInBin = sqrt(((Double_t)k+1)*((Double_t)k+2)/(n+2)/(n+3) - pow((Double_t)(k+1),2)/pow(n+2,2));
      hGenYield[i]->SetBinContent(b, (Double_t)k/(Double_t)n * 1./hGenYield[i]->GetBinWidth(b));
      hGenYield[i]->SetBinError(b, errorInBin * 1./hGenYield[i]->GetBinWidth(b));
    }

    hYieldsRatio[i]->Divide(hGenYield[i]);
    ErrRatioCorr(hYield[i], hGenYield[i], hYieldsRatio[i], 0);
  }

  // Scale yields for a better separation on the plot
  TString sScaleFactor[numMult + 1];
  hYield[0]->Scale(ScaleFactorMB);
  hGenYield[0]->Scale(ScaleFactorMB);
  sScaleFactor[0] = Form(" (x2^{%i})", int(log2(ScaleFactorMB)));
  for (Int_t iFile = 1; iFile < numMult + 1; iFile++) {
    hYield[iFile]->Scale(ScaleFactor[iFile-1]);
    hGenYield[iFile]->Scale(ScaleFactor[iFile-1]);
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

  Double_t yieldYLow[2] = {1e-9, 0.2*1e-8};
  Double_t yieldYUp[2] = {9*1e5, 9*1e4};

  // Check for negative/0 values
  for (Int_t i = 0; i < numMult + 1; i++) {
    for (Int_t b = 1; b <= hGenYield[i]->GetNbinsX(); b++) {
      if(hGenYield[i]->GetBinContent(b) <= 0 || hYield[i]->GetBinContent(b) <= 0){
        std::cout << "Bin content gen: " << hGenYield[i]->GetBinContent(b) << std::endl;
        std::cout << "Bin content rec: " << hYield[i]->GetBinContent(b) << std::endl;
      }
    }
  }

  Int_t numMultPlotted = 0; // 0 - only MB
  std::cout << "Number of plotted classes + MB: " << numMultPlotted + 1 << std::endl;

  {
    TCanvas* canvasCompPlot = new TCanvas("canvasCompPlot_MCclosure", "canvasCompPlot_MCclosure", 0, 70, 620, 850);
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
    legClasses->SetHeader("FT0M % (Open Circles: Gen.; Full ones: Rec.)");
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
    TH1F *hDummyUp = new TH1F("hDummyUp_MCclosure", "hDummyUp_MCclosure", 10000, 0, binpt[numPtBins] + 0.5);
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
    TH1F *hDummyLow = new TH1F("hDummyLow_MCclosure", "hDummyLow_MCclosure", 10000, 0, binpt[numPtBins] + 0.5);
    for (Int_t i = 1; i <= hDummyLow->GetNbinsX(); i++)
      hDummyLow->SetBinContent(i, 1);
    padLow->cd();
    StyleHisto(hDummyLow, 0.01, 1.5, 1, 1, sPt, "rec./gen.", "", 0, 0, 0, 1.0, 0.7, 0, 0.08, 0.08, 0.08, 0.07, hDummyLow->GetXaxis()->GetLabelOffset(), 0.01);
    SetTickLength(hDummyLow, 0.025, 0.03);
    TAxis *axisDummyLowY = hDummyLow->GetYaxis();
    //padLow->SetGridy(); // Only grid on x-axis
    axisDummyLowY->ChangeLabel(-1, -1, -1, -1, -1, -1, " "); // last 
    axisDummyLowY->ChangeLabel(1, -1, -1, -1, -1, -1, " "); // first
    hDummyLow->GetXaxis()->SetRangeUser(0, binpt[numPtBins] + 0.5);
    hDummyLow->GetYaxis()->CenterTitle();
    hDummyLow->Draw("same");

    //Plotting
    for (Int_t i = 0; i < numMultPlotted + 1; i++) {
      padUp->cd();
      StyleHisto(hYield[i], yieldYLow[partType], yieldYUp[partType], color[i], 8, "", hYield[i]->GetYaxis()->GetTitle(), "", 0, 0, 0, 0, 0, SizeMult[i], 0, 0, 0, 0, 0, 0);
      hYield[i]->Draw("same");

      StyleHisto(hGenYield[i], yieldYLow[partType], yieldYUp[partType], color[i], 4, "", hGenYield[i]->GetYaxis()->GetTitle(), "", 0, 0, 0, 0, 0, SizeMult[i], 0, 0, 0, 0, 0, 0);
      hGenYield[i]->Draw("same");

      legClasses->AddEntry(hYield[i], SmoltBis[i], "pef");

      padLow->cd();
      StyleHisto(hYieldsRatio[i], 0.1 , 1.9, color[i], 8, sPt, "rec./gen.", "", 0, 0, 0, 0, 0, SizeMult[i], 0, 0, 0, 0, 0, 0);
      hYieldsRatio[i]->Draw("same");
    }

    padUp->cd();
    legendTitleMC->Draw("same");
    legClasses->Draw("same");
    canvasCompPlot->Write();
  }

  // End of Code
  std::cout << "\x1B[1;32m"; // Set text color to green
  std::cout << "\n************* MC Closure is Finished! *************\n";
  std::cout << "\x1B[0m"; // Reset text color
}
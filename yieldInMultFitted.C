#include "help.h"
#include "fitting.h"

void yieldInMultFitted(const Int_t nParticle = 2, // 0-2 : xi, 3-5 : omega
                       const Int_t inel = 0, // inel > N (0/1)
                       const Int_t typefit = 3, // 0 - mT scaling, 1 - Boltzmann, 2 - Fermi-Dir, 3 - Levi, 4 - Bose-Einstein
                       const TString workingDir = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM",
                       const TString postFix = "")
{
  // Start of Code
  std::cout << "\x1B[1;33m"; // Set text color to yellow
  std::cout << "\n************ Starting Fitting the Spectra in Multiplicity Classes ************\n";
  std::cout << "\x1B[0m"; // Reset text color

  gStyle->SetOptStat(0);

  TString outputfilePath = workingDir + "/yieldInMultFitted/" + "yieldInMultFitted_" + particleNames[nParticle] + "_inel" + inel + postFix + ".root";
  TFile *outputfile = new TFile(outputfilePath, "RECREATE");

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

  // Run 3 systematics
  TString fileInPathSyst = workingDir + "/systematics/sourcesOfSyst/totalSystematics.root";
  TFile* fileSystin = TFile::Open(fileInPathSyst);
  if (!fileSystin || fileSystin->IsZombie()) {
      std::cerr << "Error opening `fileSystin` data file!" << std::endl;
      return;
  }
  TH1F* hSystUncert = (TH1F *)fileSystin->Get("hSystTotal");

  // HEP Data File (MB integrated yield)
  TString inputHEPDataPath = workingDir + "/data/published/integrated-xi.root";
  TFile* fileHEPDataIn = TFile::Open(inputHEPDataPath);
  if (!fileHEPDataIn || fileHEPDataIn->IsZombie()) {
      std::cerr << "Error opening HEP data file!" << std::endl;
      return;
  }
  
  TDirectory* yieldHEPDir = fileHEPDataIn->GetDirectory("Table 11a");
  if (!yieldHEPDir)
  {
    std::cerr << "`Table 11a` directory in HEP data file is not found!" << std::endl;
    return;
  }

  TH1F* hHEPYield = (TH1F *)yieldHEPDir->Get("Hist1D_y3");
  TH1F* hHEPYieldSyst = (TH1F *)hHEPYield->Clone("Hist1D_y3_Syst");
  TH1F* hHEPYielde1 = (TH1F *)yieldHEPDir->Get("Hist1D_y3_e1");
  TH1F* hHEPYielde2 = (TH1F *)yieldHEPDir->Get("Hist1D_y3_e2");
  if (!hHEPYield || !hHEPYielde1 || !hHEPYielde2)
  {
    std::cerr << "Histogram `Hist1D_y3` is not found!" << std::endl;
    return;
  }

  Double_t e1 = 0;
  Double_t e2 = 0;

  for (Int_t iBin = 1; iBin <= hHEPYield->GetNbinsX(); iBin++) {
    e1 = hHEPYielde1->GetBinContent(iBin);
    e2 = hHEPYielde2->GetBinContent(iBin);
    //hHEPYield->SetBinError(iBin, sqrt(pow(e1, 2) + pow(e2, 2)) + e3);
    hHEPYield->SetBinError(iBin, e1); // only statistical uncert.
    hHEPYieldSyst->SetBinError(iBin, e2); // only syst. uncert.
  }

  // HEP Data File (integrated yield in V0M classes)
  TString inputHEPDataInClassesPath = workingDir + "/data/published/integratedInClasses-xi.root";
  TFile* fileHEPDataInClassesIn = TFile::Open(inputHEPDataInClassesPath);
  if (!fileHEPDataInClassesIn || fileHEPDataInClassesIn->IsZombie()) {
      std::cerr << "Error opening HEP data file (integrated in classes)!" << std::endl;
      return;
  }
  
  TDirectory* yieldHEP2Dir = fileHEPDataInClassesIn->GetDirectory("Table 8c");
  if (!yieldHEP2Dir)
  {
    std::cerr << "`Table 8c` directory in HEP data file is not found!" << std::endl;
    return;
  }

  TH1F* hHEPYieldInClasses = (TH1F *)yieldHEP2Dir->Get("Hist1D_y3");
  TH1F* hHEPYieldInClassese1 = (TH1F *)yieldHEP2Dir->Get("Hist1D_y3_e1");
  TH1F* hHEPYieldInClassese2 = (TH1F *)yieldHEP2Dir->Get("Hist1D_y3_e2");
  TH1F* hHEPYieldInClassese3 = (TH1F *)yieldHEP2Dir->Get("Hist1D_y3_e3");

  // HEP Data File (integrated yield in tracklet classes )
  TString inputHEPDataInTrackletClassesPath = workingDir + "/data/published/integratedInClassesTracklets-xi.root";
  TFile* fileHEPDataInTrackletClassesIn = TFile::Open(inputHEPDataInTrackletClassesPath);
  if (!fileHEPDataInTrackletClassesIn || fileHEPDataInTrackletClassesIn->IsZombie()) {
      std::cerr << "Error opening HEP data file (integrated in classes)!" << std::endl;
      return;
  }
  
  TDirectory* yieldHEP3Dir = fileHEPDataInTrackletClassesIn->GetDirectory("Table 8a");
  if (!yieldHEP3Dir)
  {
    std::cerr << "`Table 8c` directory in HEP data file is not found!" << std::endl;
    return;
  }

  TH1F* hHEPYieldInTrackletClasses = (TH1F *)yieldHEP3Dir->Get("Hist1D_y3");
  TH1F* hHEPYieldInTrackletClassese1 = (TH1F *)yieldHEP3Dir->Get("Hist1D_y3_e1");
  TH1F* hHEPYieldInTrackletClassese2 = (TH1F *)yieldHEP3Dir->Get("Hist1D_y3_e2");
  TH1F* hHEPYieldInTrackletClassese3 = (TH1F *)yieldHEP3Dir->Get("Hist1D_y3_e3");

  TH1F* hYield[numMult + 1];
  TH1F* hYieldSyst[numMult + 1];
  TH1F* hYieldTotal[numMult + 1];
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
    hYieldTotal[iFile] = (TH1F*)hYield[iFile]->Clone(Form("YieldTotal_%i_",iFile) + particleNames[nParticle]); // eff. corrected yield with total uncertainty -> to be used for fitting
    for (Int_t bin = 1; bin <= hYield[iFile]->GetNbinsX(); bin++) {
      hYieldSyst[iFile]->SetBinError(bin, hSystUncert->GetBinContent(bin)*hYield[iFile]->GetBinContent(bin));
      hYieldTotal[iFile]->SetBinError(bin, sqrt(pow(hYield[iFile]->GetBinError(bin),2) + pow(hYieldSyst[iFile]->GetBinError(bin),2)));
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
  
  // Scale yields for a better separation on the plot
  TString sScaleFactor[numMult + 1];
  sScaleFactor[0] = Form(" (x2^{%i})", int(log2(ScaleFactorMB)));
  for (Int_t iFile = 1; iFile < numMult + 1; iFile++) {
    if (int(log2(ScaleFactor[iFile-1])) == 0) {
      sScaleFactor[iFile] = "";
    } else if (int(log2(ScaleFactor[iFile-1])) == 1) {
      sScaleFactor[iFile] = " (x2)";
    } else {
      sScaleFactor[iFile] = Form(" (x2^{%i})", int(log2(ScaleFactor[iFile-1])));
    } 
  }

  TString SmoltBis[numMult + 1];
  TString SmoltBisNoFactor[numMult + 1];
  for (Int_t m = 0; m < numMult + 1; m++) {
    if (m == 0) { 
      SmoltBis[m] += Form("%.1f#minus%.1f%s", 0.0, 100.0, "%") + sScaleFactor[m];
      SmoltBisNoFactor[m] += Form("%.1f#minus%.1f%s", 0.0, 100.0, "%");
    } else {
      SmoltBis[m] += Form("%.1f#minus%.1f%s", multiplicityPerc[m-1], multiplicityPerc[m], "%") + sScaleFactor[m];
      SmoltBisNoFactor[m] += Form("%.1f#minus%.1f%s", multiplicityPerc[m-1], multiplicityPerc[m], "%");
    }
  }

  Double_t yieldYLow[2] = {0.2*1e-9, 0.2*1e-8};
  Double_t yieldYUp[2] = {3*1e5, 3*1e4};

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
  StyleHisto(hDummy, yieldYLow[partType], yieldYUp[partType], 1, 1, "#it{p}_{T} (GeV/#it{c})", sdNdPtdY, "", 0, 0, 0, 1.5, 1.0, 0, 0.0, 0.05, 0.0, 0.035, hDummy->GetXaxis()->GetLabelOffset(), 0.005);
  SetTickLength(hDummy, 0.025, 0.03);
  TAxis *axisYieldDummy = hDummy->GetYaxis();
  axisYieldDummy->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
  hDummy->GetXaxis()->SetRangeUser(0, hYield[0]->GetXaxis()->GetBinUpEdge(hYield[0]->GetNbinsX()) + 0.5);
  hDummy->GetYaxis()->SetMaxDigits(4);
  hDummy->GetYaxis()->SetDecimals(kTRUE);
  hDummy->Draw("same");
  // Low
  TH1F *hDummyLow = new TH1F("hDummyLow", "hDummyLow", 10000, 0, 8.5);
  for (Int_t i = 1; i <= hDummyLow->GetNbinsX(); i++)
    hDummyLow->SetBinContent(i, 1);
  padYieldLow->cd();
  StyleHisto(hDummyLow, 0.1, 1.9, 1, 1, "#it{p}_{T} (GeV/#it{c})", "Data/Fit", "", 0, 0, 0, 1.0, 0.7, 0, 0.08, 0.08, 0.08, 0.07, hDummyLow->GetXaxis()->GetLabelOffset(), 0.01);
  SetTickLength(hDummyLow, 0.025, 0.03);
  hDummyLow->GetXaxis()->SetRangeUser(0, hYieldsRatio[0]->GetXaxis()->GetBinUpEdge(hYieldsRatio[0]->GetNbinsX()) + 0.5);
  hDummyLow->GetYaxis()->CenterTitle();
  TAxis *axisDummyLow = hDummyLow->GetYaxis();
  axisDummyLow->ChangeLabel(-1, -1, -1, -1, -1, -1, " "); // last 
  axisDummyLow->ChangeLabel(1, -1, -1, -1, -1, -1, " "); // first
  hDummyLow->Draw("same");

  // Define colors
  for(int iClass = 0; iClass < numMult + 1; iClass++) { 
    if(iClass==0) {
      color[iClass] = kBlack;
    } else {
      color[iClass] = numMult + FI - iClass;
    }
  }

  // Setup Fitting
  TString namepwgfunc[numMult + 1];
  // Setup variables/arrays
  TFitResultPtr fFitResultPtr[numMult + 1];
  TF1 *fitFunction[numMult + 1];
  TF1 *fitFunctionScaled[numMult + 1];
  TH1F* hFitsRatio[numMult + 1];
  TH1F* hYieldForStats[numMult + 1] = {0};
  TH1F* hYieldForStatsSyst[numMult + 1] = {0};

  Double_t chi2NDF[numMult + 1] = {0};

  // Fit limits
  Double_t fitLow[numFitFuncs] = {0.6, 0.6, 0.6, 0.6, 0.6, 0.6};
  Double_t fitUp[numFitFuncs] = {2.5, 2.5, 2.2, 8.0, 2.5, 2.5};

  gROOT->SetBatch(kFALSE);
  // Plotting yeilds with systematic and fitting //
  for (Int_t iFile = 0; iFile < numMult + 1; iFile++) {
    // Up
    StyleHisto(hYield[iFile], yieldYLow[partType], yieldYUp[partType], color[iFile], MarkerMult[iFile], "", hYield[iFile]->GetYaxis()->GetTitle(), "", 0, 0, 0, 1.5, 1.0, SizeMult[iFile], 0.0, 0.05, 0.0, 0.035, hYield[iFile]->GetXaxis()->GetLabelOffset(), 0.005);
    StyleHisto(hYieldSyst[iFile], yieldYLow[partType], yieldYUp[partType], color[iFile], MarkerMult[iFile], "", hYield[iFile]->GetYaxis()->GetTitle(), "", 0, 0, 0, 1.5, 1.0, SizeMult[iFile], 0.0, 0.05, 0.0, 0.035, hYieldSyst[iFile]->GetXaxis()->GetLabelOffset(), 0.005);
    legYield->AddEntry(hYield[iFile], SmoltBis[iFile], "pef");

    // Fitting
    namepwgfunc[iFile] = Form("fitpwgfunc_m%i_fit%i", iFile, typefit);
    std::cout << "\n********* Fitting pt spectra with: " << nameFit[typefit] << " *********\n";
    if (iFile > 0) {
      cout << "Mult. class: " << multiplicityPerc[iFile-1] <<  " - " << multiplicityPerc[iFile] << std::endl;
    } else {
      cout << "Minimum Bias" << std::endl;
    }
    std::cout << "\n****** Name of the fitting function: " << namepwgfunc[iFile] << " ******\n";
  
    // Save Yield histogram to make to-fit ratio after fitting
    hFitsRatio[iFile] = (TH1F *)hYield[iFile]->Clone(Form("hFitsRatio_m%i_typefit%i", iFile, typefit));
    // Save not scaled yield histogram to pass it then to YieldMean function
    hYieldForStats[iFile] = (TH1F *)hYield[iFile]->Clone(Form("hYieldForStats_m%i_fit%i", iFile, typefit));
    hYieldForStatsSyst[iFile] = (TH1F *)hYieldSyst[iFile]->Clone(Form("hYieldForStatsSyst_m%i_fit%i", iFile, typefit));
    // Setup fit function
    Double_t specNorm = hYield[iFile]->GetBinContent(hYield[iFile]->GetMaximumBin());
    switch (typefit) {
      case 0:
        fitFunction[iFile] = GetMTExpdNdptTimesPt(pdgMass[partType], 0.1, specNorm, namepwgfunc[iFile]); // mass, T, norm, name
        break;
      case 1:
        fitFunction[iFile] = GetBoltzmanndNdptTimesPt(pdgMass[partType], 0.1, specNorm, namepwgfunc[iFile]);
        break;
      case 2:
        fitFunction[iFile] = GetFermiDiracdNdptTimesPt(pdgMass[partType], 0.1, specNorm, namepwgfunc[iFile]);
        break;
      case 3:
        fitFunction[iFile] = GetLevidNdptTimesPt(pdgMass[partType], 0.1, 4, specNorm, namepwgfunc[iFile]); // m, T, n, norm
        fitFunction[iFile]->SetParLimits(0, 0, specNorm * 10); // norm
        fitFunction[iFile]->SetParLimits(1, 2, 30); // n
        fitFunction[iFile]->SetParLimits(2, 0.01, 10); // T (GeV)
        fitFunction[iFile]->SetParameter(2, 0.5); // T
        break;
      case 4:
        fitFunction[iFile] = GetBoseEinsteindNdptTimesPt(pdgMass[partType], 0.1, specNorm, namepwgfunc[iFile]); // mass, T, norm, name
        break;
      // case 5:
      //   fitFunction[iFile] = GetBGBWdNdptTimesPt(pdgMass[partType], 0.5, 0.1, 0.03, 0.04 * factor, namepwgfunc[iFile]); // doesn't work
      //   break;
      default:
        std::cerr << "Fit type is not implemented" << std::endl;
        break;
    }

    // Fit spectrum in range
    fitFunction[iFile]->SetRange(fitLow[typefit], fitUp[typefit]); // Fit Range
    fFitResultPtr[iFile] = hYieldTotal[iFile]->Fit(fitFunction[iFile], "SR0IQ"); // fit with the total uncert., do not draw
    switch (typefit) {
      case 1:
      case 2:
      case 4:
        std::cout << "norm: " << fitFunction[iFile]->GetParameter(0) << std::endl;
        std::cout << "T: " << fitFunction[iFile]->GetParameter(1) << std::endl;
        break;
      case 3:
        std::cout << "norm: " << fitFunction[iFile]->GetParameter(0) << std::endl;
        std::cout << "n: " << fitFunction[iFile]->GetParameter(1) << std::endl;
        std::cout << "T: " << fitFunction[iFile]->GetParameter(2) << std::endl;
        std::cout << "mass: " << fitFunction[iFile]->GetParameter(3) << std::endl;
        break;
      default:
        std::cerr << "Fit type is not implemented" << std::endl;
        break;
    }

    chi2NDF[iFile] = fitFunction[iFile]->GetChisquare() / fitFunction[iFile]->GetNDF();
    std::cout << "chi/ndf: " << chi2NDF[iFile] << std::endl;

    // Scaling
    fitFunctionScaled[iFile] = (TF1 *)fitFunction[iFile]->Clone(namepwgfunc[iFile] + Form("_Scaled_m%i", iFile));
    fitFunctionScaled[iFile]->SetLineColor(color[iFile]);
    fitFunctionScaled[iFile]->SetLineWidth(2);
    fitFunctionScaled[iFile]->SetLineStyle(7);
    fitFunctionScaled[iFile]->SetRange(0., fitUp[typefit]); // Display range
    if (iFile == 0) { // MB
      hYield[0]->Scale(ScaleFactorMB);
      hYieldSyst[0]->Scale(ScaleFactorMB);
      fitFunctionScaled[0]->SetParameter(0, fitFunction[0]->GetParameter(0) * ScaleFactorMB); // make sure param #0 is just normalization!
    } else { // in classes
      hYield[iFile]->Scale(ScaleFactor[iFile-1]);
      hYieldSyst[iFile]->Scale(ScaleFactor[iFile-1]);
      fitFunctionScaled[iFile]->SetParameter(0, fitFunction[iFile]->GetParameter(0) * ScaleFactor[iFile-1]); // make sure param #0 is just normalization!
    }

    padYieldUp->cd();
    hYield[iFile]->Draw("same ex0");
    hYieldSyst[iFile]->SetFillStyle(0);
    hYieldSyst[iFile]->Draw("same e2");
    fitFunctionScaled[iFile]->Draw("same");

    // Low
    // data yield/fit yield (ratio of integrated yield bin by bin)
    for (Int_t b = 1; b <= hFitsRatio[iFile]->GetNbinsX(); b++)
    {
      Double_t integralFit = fitFunction[iFile]->Integral(hFitsRatio[iFile]->GetXaxis()->GetBinLowEdge(b), hFitsRatio[iFile]->GetXaxis()->GetBinUpEdge(b));
      Double_t errorYield = hYieldTotal[iFile]->GetBinError(b) / hYieldTotal[iFile]->GetBinContent(b); // stat. + syst.
      Double_t errorFit = fitFunction[iFile]->IntegralError(hFitsRatio[iFile]->GetXaxis()->GetBinLowEdge(b), hFitsRatio[iFile]->GetXaxis()->GetBinUpEdge(b)) / integralFit;
      hFitsRatio[iFile]->SetBinContent(b, hYieldTotal[iFile]->GetBinContent(b) * hYieldTotal[iFile]->GetBinWidth(b) / integralFit);
      hFitsRatio[iFile]->SetBinError(b, sqrt(pow(errorYield, 2) + pow(errorFit, 2))*hFitsRatio[iFile]->GetBinContent(b));
    }
    StyleHisto(hFitsRatio[iFile], 0.1 , 1.9, color[iFile], MarkerMult[iFile], "#it{p}_{T} (GeV/#it{c})", "Ratio to 0-100%", "", 0, 0, 0, 1.0, 0.7, SizeMult[iFile], 0.08, 0.08, 0.08, 0.08, hFitsRatio[iFile]->GetXaxis()->GetLabelOffset(), 0.005);
    padYieldLow->cd();
    hFitsRatio[iFile]->Draw("same e0x0");
  }

  padYieldUp->cd();
  canvasYield->Draw();
  legendTitle->Draw();

  TLegend *legendFit = new TLegend(0.19, 0.297, 0.383, 0.398);
  legendFit->SetFillStyle(0);
  legendFit->SetTextSize(0.04);
  legendFit->AddEntry(fitFunctionScaled[0], nameFit[typefit] + " fit", "l");
  legendFit->Draw();

  legYield->Draw();
  gPad->Update();
  canvasYield->Update();

  padYieldLow->cd();
  DrawHorLine(hFitsRatio[0]->GetXaxis()->GetBinUpEdge(hFitsRatio[0]->GetNbinsX()) + 0.5, 1.0);
  outputfile->cd();
  canvasYield->Write();

  // Analysis Of Integrated Yields
  TH1 *hhout[numMult + 1];
  // TString Titlehhout[9] = {"kYield",
  //                          "kYieldStat",
  //                          "kYieldSysHi",
  //                          "kYieldSysLo",
  //                          "kMean",
  //                          "kMeanStat",
  //                          "kMeanSysHi",
  //                          "kMeanSysLo",
  //                          "kExtra"};

  // Double_t YieldExtr[numMult + 1] = {0};
  // Double_t YieldErrStat[numMult + 1] = {0};
  // Double_t YieldErrSistHi[numMult + 1] = {0};
  // Double_t YieldErrSistLow[numMult + 1] = {0};

  // Double_t Mean[numMult + 1] = {0};
  // Double_t MeanErrStat[numMult + 1] = {0};
  // Double_t MeanErrSistHi[numMult + 1] = {0};
  // Double_t MeanErrSistLow[numMult + 1] = {0};

  // Double_t Temp[numMult + 1] = {0};
  // Double_t TempError[numMult + 1] = {0};

  gROOT->SetBatch(kTRUE);
  Float_t integratedYield[numMult + 1] = {0};
  Float_t integratedYieldStat[numMult + 1] = {0};
  Float_t integratedYieldSystHigh[numMult + 1] = {0};
  Float_t integratedYieldSystLow[numMult + 1] = {0};

  TH1F* hIntegratedYield = new TH1F("hIntegratedYield", "hIntegratedYield", numMult + 1, 0, numMult + 1);
  TH1F* hIntegratedYieldSyst = new TH1F("hIntegratedYieldSyst", "hIntegratedYieldSyst", numMult + 1, 0, numMult + 1);

  TH1F* hIntegratedYieldRun2 = new TH1F("hIntegratedYieldRun2", "hIntegratedYieldRun2", numMult + 1, 0, numMult + 1);
  TH1F* hIntegratedYieldRun2Syst = new TH1F("hIntegratedYieldRun2Syst", "hIntegratedYieldRun2Syst", numMult + 1, 0, numMult + 1);

  for (Int_t iFile = 0; iFile < numMult + 1; iFile++) {
    cout << "\n*** Calling YieldMean macro ***\n";
    hhout[iFile] = YieldMean(hYieldForStats[iFile], hYieldForStatsSyst[iFile], fitFunction[iFile], 0, 20, 0.01, 0.1, "0qI", "logs/yieldMean-log.root", fitLow[typefit], fitUp[typefit]);
    cout << "\n*** End of call ***\n";

    // hhout[iFile]->SetLineColor(ColorFit[typefit]);
    // hhout[iFile]->SetName("hhout_" + nameFit[typefit] + Form("_m%i", iFile));
    // hhout[iFile]->GetYaxis()->SetRangeUser(0, 2);
    // for (Int_t b = 1; b <= hhout[iFile]->GetNbinsX(); b++)
    // {
    //   hhout[iFile]->GetXaxis()->SetBinLabel(b, Titlehhout[b - 1]);
    // }

    // cout << "iFile " << iFile << " typefit " << typefit << endl;
    integratedYield[iFile] = hhout[iFile]->GetBinContent(1);           // yield (spectra + extrapolated one)
    // YieldExtr[iFile] = hhout[iFile]->GetBinContent(9);       // extrapolated yield
    integratedYieldStat[iFile] = hhout[iFile]->GetBinContent(2);    // stat error
    integratedYieldSystHigh[iFile] = hhout[iFile]->GetBinContent(3);  // syst error
    integratedYieldSystLow[iFile] = hhout[iFile]->GetBinContent(4); // syst error
    // Mean[iFile] = hhout[iFile]->GetBinContent(5);            // mean pt
    // MeanErrStat[iFile] = hhout[iFile]->GetBinContent(6);     // stat error
    // MeanErrSistHi[iFile] = hhout[iFile]->GetBinContent(7);   // syst error
    // MeanErrSistLow[iFile] = hhout[iFile]->GetBinContent(8);  // syst error
    // cout << "************************************" << endl;
    // if (iFile > 0) {
    //   cout << "mult.%: " << multiplicityPerc[iFile - 1] <<  " - " << multiplicityPerc[iFile] << std::endl;
    // } else {
    //   cout << " MB " << std::endl;
    // }
    // cout << "integratedYield: " << integratedYield[iFile] << " +- " << integratedYieldStat[iFile] << " (stat.) + " << integratedYieldSystHigh[iFile] << " - " << integratedYieldSystLow[iFile] << " (syst.) " << endl;
    // cout << "************************************" << endl;
    // cout << "Mean: " << Mean[iFile] << " +- " << MeanErrStat[iFile] << endl;

    // if (typefit == 3)
    // {
    //   Temp[iFile] = fitFunction[iFile]->GetParameter(2);
    //   TempError[iFile] = fitFunction[iFile]->GetParError(2);
    // }
    // else
    // {
    //   Temp[iFile] = fitFunction[iFile]->GetParameter(1);
    //   TempError[iFile] = fitFunction[iFile]->GetParError(1);
    // }
  }

  TLegend *legendTitleIntegrated = new TLegend(0.576, 0.629, 0.871, 0.829);
  StyleLegend(legendTitleIntegrated, 0.0, 0.0);
  legendTitleIntegrated->SetTextAlign(33);
  legendTitleIntegrated->SetTextSize(0.05);
  legendTitleIntegrated->SetTextFont(42);
  legendTitleIntegrated->SetLineColorAlpha(0.,0.);
  legendTitleIntegrated->SetFillColorAlpha(0.,0.);
  legendTitleIntegrated->AddEntry("", "#bf{ALICE Work In Progress}", "");
  legendTitleIntegrated->AddEntry("", "pp, #sqrt{#it{s}} = 13.6 TeV", "");
  legendTitleIntegrated->AddEntry("", particleSymnbols[nParticle] + ", |y| < 0.5", "");

  TLegend *legendIntegrated = new TLegend(0.236, 0.212, 0.437, 0.311);
  legendIntegrated->SetFillStyle(0);
  legendIntegrated->SetTextSize(0.04);
  legendIntegrated->Draw();

  for (Int_t m = numMult; m >= 0; m--) {
    hIntegratedYield->GetXaxis()->SetBinLabel(m + 1, SmoltBisNoFactor[m]);
    hIntegratedYield->SetBinContent(m + 1, integratedYield[m]);
    hIntegratedYield->SetBinError(m + 1, integratedYieldStat[m]);
    hIntegratedYieldSyst->GetXaxis()->SetBinLabel(m + 1, SmoltBisNoFactor[m]);
    hIntegratedYieldSyst->SetBinContent(m + 1, integratedYield[m]);
    hIntegratedYieldSyst->SetBinError(m + 1, 1. / 2 * (integratedYieldSystHigh[m] + integratedYieldSystLow[m])); // proxy
    if(m == 0) { // take only MB from Run 2
      hIntegratedYieldRun2->SetBinContent(m + 1, hHEPYield->GetBinContent(1));
      hIntegratedYieldRun2->SetBinError(m + 1, hHEPYield->GetBinError(1));
      hIntegratedYieldRun2Syst->SetBinContent(m + 1, hHEPYieldSyst->GetBinContent(1));
      hIntegratedYieldRun2Syst->SetBinError(m + 1, hHEPYieldSyst->GetBinError(1));
    } else {
      hIntegratedYieldRun2->SetBinContent(m + 1, 0);
      hIntegratedYieldRun2->SetBinError(m + 1, 0);
      hIntegratedYieldRun2Syst->SetBinContent(m + 1, 0);
      hIntegratedYieldRun2Syst->SetBinError(m + 1, 0);
    }
    // hChi2->GetXaxis()->SetBinLabel(m + 1, SmoltBis[m] + "%");
    // hChi2->SetBinContent(m + 1, Chi2NDF[m]);
    // hChi2->SetBinError(m + 1, 0);
    // hTemp->GetXaxis()->SetBinLabel(m + 1, SmoltBis[m] + "%");
    // hTemp->SetBinContent(m + 1, Temp[m]);
    // hTemp->SetBinError(m + 1, TempError[m]);
    // hFracExtrYield->GetXaxis()->SetBinLabel(m + 1, SmoltBis[m] + "%");
    // hFracExtrYield->SetBinContent(m + 1, YieldExtr[m] / Yield[m]);
    // hFracExtrYield->SetBinError(m + 1, 0);
    cout << "\n\n************************************" << endl;
    cout << "Multiplicity class: " << SmoltBis[m] << endl;
    cout << "integratedYield: " << integratedYield[m] << " +- " << integratedYieldStat[m] << " (stat.) + " << integratedYieldSystHigh[m] << " - " << integratedYieldSystLow[m] << " (syst.) " << endl;
    // cout << "Mean: " << Mean[m] << " +- " << MeanErrStat[m] << endl;
  }
  gROOT->SetBatch(kFALSE);

  TCanvas *canvasIntegratedYields = new TCanvas("canvasIntegratedYields","canvasIntegratedYields", 800, 600);
  canvasIntegratedYields->cd();
  canvasIntegratedYields->Draw();
  StyleCanvas(canvasIntegratedYields, 0.14, 0.10, 0.11, 0.15);
  StyleHisto(hIntegratedYield, 0, 1.2 * hIntegratedYield->GetBinContent(hIntegratedYield->GetMaximumBin()), color[0], MarkerMult[0], "FT0M Multiplicity percentile", sdNdY, "", 0, 0, 0, 1.35, 1.25, 1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  StyleHisto(hIntegratedYieldSyst, 0, 1.2 * hIntegratedYieldSyst->GetBinContent(hIntegratedYieldSyst->GetMaximumBin()), color[0], MarkerMult[0], "FT0M Multiplicity percentile", sdNdY, "", 0, 0, 0, 1.35, 1.25, 1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  StyleHisto(hIntegratedYieldRun2, 0, 1.2 * hIntegratedYieldRun2->GetBinContent(hIntegratedYieldRun2->GetMaximumBin()), color[1], MarkerMult[1], "FT0M Multiplicity percentile", sdNdY, "", 0, 0, 0, 1.35, 1.25, 1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  StyleHisto(hIntegratedYieldRun2Syst, 0, 1.2 * hIntegratedYieldRun2Syst->GetBinContent(hIntegratedYieldRun2Syst->GetMaximumBin()), color[1], MarkerMult[1], "FT0M Multiplicity percentile", sdNdY, "", 0, 0, 0, 1.35, 1.25, 1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  hIntegratedYield->Draw("same EP");
  hIntegratedYieldSyst->SetFillStyle(0);
  hIntegratedYieldSyst->Draw("same e2");

  hIntegratedYieldRun2->Draw("same EP");
  hIntegratedYieldRun2Syst->SetFillStyle(0);
  hIntegratedYieldRun2Syst->Draw("same e2");

  legendIntegrated->AddEntry(hIntegratedYield, nameFit[typefit] + " fit", "pl");
  legendIntegrated->AddEntry(hIntegratedYieldRun2, "Eur.Phys.J.C.80 (2020) 167, pp, #sqrt{#it{s}} = 13 TeV", "pl");
  legendIntegrated->Draw("same");
  legendTitleIntegrated->Draw("same");
  outputfile->cd();
  canvasIntegratedYields->Write();

  // Integrated yield vs dN/dEta //
  TLegend *legendTitleIntegrateddNdEta = new TLegend(0.150, 0.669, 0.352, 0.868);
  StyleLegend(legendTitleIntegrateddNdEta, 0.0, 0.0);
  legendTitleIntegrateddNdEta->SetTextAlign(11);
  legendTitleIntegrateddNdEta->SetTextSize(0.04);
  legendTitleIntegrateddNdEta->SetTextFont(42);
  legendTitleIntegrateddNdEta->SetLineColorAlpha(0.,0.);
  legendTitleIntegrateddNdEta->SetFillColorAlpha(0.,0.);
  legendTitleIntegrateddNdEta->AddEntry("", "#bf{ALICE Work In Progress}", "");
  legendTitleIntegrateddNdEta->AddEntry("", "pp, #sqrt{#it{s}} = 13.6 TeV", "");
  legendTitleIntegrateddNdEta->AddEntry("", particleSymnbols[nParticle] + ", |y| < 0.5", "");

  TLegend *legendIntegrateddNdEta = new TLegend(0.191, 0.441, 0.392, 0.640);
  legendIntegrateddNdEta->SetFillStyle(0);
  legendIntegrateddNdEta->SetTextSize(0.03);

  TLegend * legUncert = new TLegend(0.604, 0.225, 0.914, 0.316);
  legUncert->SetFillColor(0);
  legUncert->SetTextAlign(11);
  legUncert->SetTextSize(0.04);

  TCanvas *canvasIntegratedYieldsdNdEta = new TCanvas("canvasIntegratedYieldsdNdEta","canvasIntegratedYieldsdNdEta", 0,66,857,708);
  canvasIntegratedYieldsdNdEta->cd();
  canvasIntegratedYieldsdNdEta->Draw();
  StyleCanvas(canvasIntegratedYieldsdNdEta, 0.14, 0.05, 0.11, 0.15);
  TMultiGraph *mg = new TMultiGraph("integratedYields", "");
  TGraphAsymmErrors *gdNdEta = new TGraphAsymmErrors(numMult);
  TGraphAsymmErrors *gdNdEtaSyst = new TGraphAsymmErrors(numMult);
  TGraphAsymmErrors *gdNdEtaMB = new TGraphAsymmErrors(1);
  TGraphAsymmErrors *gdNdEtaMBSyst = new TGraphAsymmErrors(1);
  TGraphAsymmErrors *gdNdEtaMB_run2 = new TGraphAsymmErrors(1);
  TGraphAsymmErrors *gdNdEtaMBSyst_run2 = new TGraphAsymmErrors(1);
  for (int i = 1; i < numMult + 1; ++i) {
    gdNdEta->SetPoint(i-1, etaDensity[i], integratedYield[i]);
    gdNdEta->SetPointError(i-1, etaDensityErrLow[i], etaDensityErrHigh[i], integratedYieldStat[i], integratedYieldStat[i]);
    gdNdEtaSyst->SetPoint(i-1, etaDensity[i], integratedYield[i]);
    gdNdEtaSyst->SetPointError(i-1, etaDensityErrLow[i], etaDensityErrHigh[i], integratedYieldSystLow[i], integratedYieldSystHigh[i]);
  }
  gdNdEtaMB->SetPoint(0, etaDensity[0], integratedYield[0]);
  gdNdEtaMB->SetPointError(0, etaDensityErrLow[0], etaDensityErrHigh[0], integratedYieldStat[0], integratedYieldStat[0]);
  gdNdEtaMBSyst->SetPoint(0, etaDensity[0], integratedYield[0]);
  gdNdEtaMBSyst->SetPointError(0, etaDensityErrLow[0], etaDensityErrHigh[0], integratedYieldSystLow[0], integratedYieldSystHigh[0]);

  gdNdEtaMB_run2->SetPoint(0, etaDensity_mb_run2[0], xiYield_mb_run2[0]);
  gdNdEtaMB_run2->SetPointError(0, etaDensity_errLow_mb_run2[0], etaDensity_errHigh_mb_run2[0], xiYield_stat_mb_run2[0], xiYield_stat_mb_run2[0]);
  gdNdEtaMBSyst_run2->SetPoint(0, etaDensity_mb_run2[0], xiYield_mb_run2[0]);
  gdNdEtaMBSyst_run2->SetPointError(0, etaDensityErrLow[0], etaDensityErrHigh[0], xiYield_syst_mb_run2[0], xiYield_syst_mb_run2[0]);

  TGraphAsymmErrors *gdNdEtaRun2 = new TGraphAsymmErrors(hHEPYieldInClasses->GetNbinsX());
  TGraphAsymmErrors *gdNdEtaRun2Syst = new TGraphAsymmErrors(hHEPYieldInClasses->GetNbinsX());
  TGraphAsymmErrors *gdNdEtaRun2SystUncorr = new TGraphAsymmErrors(hHEPYieldInClasses->GetNbinsX());

  TGraphAsymmErrors *gdNdEtaRun2Tracklets = new TGraphAsymmErrors(hHEPYieldInTrackletClasses->GetNbinsX());
  TGraphAsymmErrors *gdNdEtaRun2TrackletsSyst = new TGraphAsymmErrors(hHEPYieldInTrackletClasses->GetNbinsX());
  TGraphAsymmErrors *gdNdEtaRun2TrackletsSystUncorr = new TGraphAsymmErrors(hHEPYieldInTrackletClasses->GetNbinsX());

  for (int i = 1; i <= hHEPYieldInClasses->GetNbinsX(); ++i) {
    gdNdEtaRun2->SetPoint(i-1, hHEPYieldInClasses->GetBinCenter(i), hHEPYieldInClasses->GetBinContent(i));
    //std::cout << "x: " << hHEPYieldInClasses->GetBinCenter(i) << " y: " << hHEPYieldInClasses->GetBinContent(i) << std::endl;
    gdNdEtaRun2Syst->SetPoint(i-1, hHEPYieldInClasses->GetBinCenter(i), hHEPYieldInClasses->GetBinContent(i));
    gdNdEtaRun2SystUncorr->SetPoint(i-1, hHEPYieldInClasses->GetBinCenter(i), hHEPYieldInClasses->GetBinContent(i));

    Double_t xError = (hHEPYieldInClasses->GetXaxis()->GetBinUpEdge(i) - hHEPYieldInClasses->GetXaxis()->GetBinLowEdge(i)) / 2; // CHECK X ERROR CAUSE THERE'RE EMPTY BINS!
    gdNdEtaRun2->SetPointError(i-1, xError, xError, hHEPYieldInClassese1->GetBinContent(i), hHEPYieldInClassese1->GetBinContent(i));
    gdNdEtaRun2Syst->SetPointError(i-1, xError, xError, hHEPYieldInClassese2->GetBinContent(i), hHEPYieldInClassese2->GetBinContent(i));
    gdNdEtaRun2SystUncorr->SetPointError(i-1, xError, xError, hHEPYieldInClassese3->GetBinContent(i), hHEPYieldInClassese3->GetBinContent(i));

    gdNdEtaRun2Tracklets->SetPoint(i-1, hHEPYieldInTrackletClasses->GetBinCenter(i), hHEPYieldInTrackletClasses->GetBinContent(i));
    //std::cout << "x: " << hHEPYieldInClasses->GetBinCenter(i) << " y: " << hHEPYieldInClasses->GetBinContent(i) << std::endl;
    gdNdEtaRun2TrackletsSyst->SetPoint(i-1, hHEPYieldInTrackletClasses->GetBinCenter(i), hHEPYieldInTrackletClasses->GetBinContent(i));
    gdNdEtaRun2TrackletsSystUncorr->SetPoint(i-1, hHEPYieldInTrackletClasses->GetBinCenter(i), hHEPYieldInTrackletClasses->GetBinContent(i));

    xError = (hHEPYieldInTrackletClasses->GetXaxis()->GetBinUpEdge(i) - hHEPYieldInTrackletClasses->GetXaxis()->GetBinLowEdge(i)) / 2;
    gdNdEtaRun2Tracklets->SetPointError(i-1, xError, xError, hHEPYieldInTrackletClassese1->GetBinContent(i), hHEPYieldInTrackletClassese1->GetBinContent(i));
    gdNdEtaRun2TrackletsSyst->SetPointError(i-1, xError, xError, hHEPYieldInTrackletClassese2->GetBinContent(i), hHEPYieldInTrackletClassese2->GetBinContent(i));
    gdNdEtaRun2TrackletsSystUncorr->SetPointError(i-1, xError, xError, hHEPYieldInTrackletClassese3->GetBinContent(i), hHEPYieldInTrackletClassese3->GetBinContent(i));
  }

  // Run 3
  // In FT0M classes
  gdNdEta->SetMarkerColor(color[0]);
  gdNdEta->SetMarkerStyle(MarkerMult[0]);
  gdNdEta->SetMarkerSize(1.5);
  gdNdEtaSyst->SetMarkerColor(color[0]);
  gdNdEtaSyst->SetMarkerStyle(MarkerMult[0]);
  gdNdEtaSyst->SetFillStyle(0);
  gdNdEtaSyst->SetMarkerSize(1.5);

  // MB Run 3
  gdNdEtaMB->SetMarkerColor(color[1]);
  gdNdEtaMB->SetMarkerStyle(MarkerMult[1]);
  gdNdEtaMB->SetMarkerSize(2.0);
  gdNdEtaMBSyst->SetMarkerColor(color[1]);
  gdNdEtaMBSyst->SetMarkerStyle(MarkerMult[1]);
  gdNdEtaMBSyst->SetFillStyle(0);
  gdNdEtaMBSyst->SetMarkerSize(2.0);

  // MB Run 2
  gdNdEtaMB_run2->SetMarkerColor(color[3]);
  gdNdEtaMB_run2->SetMarkerStyle(MarkerMult[3]);
  gdNdEtaMB_run2->SetMarkerSize(2.0);
  gdNdEtaMBSyst_run2->SetMarkerColor(color[3]);
  gdNdEtaMBSyst_run2->SetMarkerStyle(MarkerMult[3]);
  gdNdEtaMBSyst_run2->SetFillStyle(0);
  gdNdEtaMBSyst_run2->SetMarkerSize(2.0); 

  // Run 2 V0M
  Int_t run2Style = 8;
  gdNdEtaRun2->SetMarkerColor(color[run2Style]);
  gdNdEtaRun2->SetMarkerStyle(MarkerMult[run2Style]);
  gdNdEtaRun2->SetMarkerSize(1.5);
  gdNdEtaRun2Syst->SetMarkerColor(color[run2Style]);
  gdNdEtaRun2Syst->SetMarkerStyle(MarkerMult[run2Style]);
  gdNdEtaRun2Syst->SetFillStyle(0);
  gdNdEtaRun2Syst->SetMarkerSize(1.5);
  gdNdEtaRun2SystUncorr->SetMarkerColor(color[run2Style]);
  gdNdEtaRun2SystUncorr->SetMarkerStyle(MarkerMult[run2Style]);
  gdNdEtaRun2SystUncorr->SetFillColor(color[run2Style]);
  gdNdEtaRun2SystUncorr->SetFillStyle(3001);
  gdNdEtaRun2SystUncorr->SetMarkerSize(1.5);

  // Run 2 Tracklets
  Int_t run2TrackletStyle = 6;
  gdNdEtaRun2Tracklets->SetMarkerColor(color[run2TrackletStyle]);
  gdNdEtaRun2Tracklets->SetMarkerStyle(MarkerMult[run2TrackletStyle]);
  gdNdEtaRun2Tracklets->SetMarkerSize(1.5);
  gdNdEtaRun2TrackletsSyst->SetMarkerColor(color[run2TrackletStyle]);
  gdNdEtaRun2TrackletsSyst->SetMarkerStyle(MarkerMult[run2TrackletStyle]);
  gdNdEtaRun2TrackletsSyst->SetFillStyle(0);
  gdNdEtaRun2TrackletsSyst->SetMarkerSize(1.5);
  gdNdEtaRun2TrackletsSystUncorr->SetMarkerColor(color[run2TrackletStyle]);
  gdNdEtaRun2TrackletsSystUncorr->SetMarkerStyle(MarkerMult[run2TrackletStyle]);
  gdNdEtaRun2TrackletsSystUncorr->SetFillColor(color[run2TrackletStyle]);
  gdNdEtaRun2TrackletsSystUncorr->SetFillStyle(3001);
  gdNdEtaRun2TrackletsSystUncorr->SetMarkerSize(1.5);

  legendIntegrateddNdEta->AddEntry(gdNdEta, "FT0M Classes", "epl");
  legendIntegrateddNdEta->AddEntry(gdNdEtaMB, "MB (Run 3)", "epl");
  legendIntegrateddNdEta->AddEntry(gdNdEtaRun2, "V0M Classes (Run 2)", "epl");
  legendIntegrateddNdEta->AddEntry(gdNdEtaMB_run2, "MB (Run 2)", "epl");
  legendIntegrateddNdEta->AddEntry(gdNdEtaRun2Tracklets, "#it{N}_{ tracklets}^{|#it{#eta}|<0.8} (Run 2)", "epl");
  mg->Add(gdNdEta, "P");
  mg->Add(gdNdEtaMB, "P");
  mg->Add(gdNdEtaMB_run2, "P");
  mg->Add(gdNdEtaSyst, "P 5");
  mg->Add(gdNdEtaMBSyst, "P 5");
  mg->Add(gdNdEtaMBSyst_run2, "P 5");
  mg->Add(gdNdEtaRun2, "P");
  mg->Add(gdNdEtaRun2Syst, "P 5");
  mg->Add(gdNdEtaRun2SystUncorr, "P 5");
  mg->Add(gdNdEtaRun2Tracklets, "P");
  mg->Add(gdNdEtaRun2TrackletsSyst, "P 5");
  mg->Add(gdNdEtaRun2TrackletsSystUncorr, "P 5");
  StyleMultGraph(mg, 1e-6, 0.165, sdNdEta, sdNdY, "", 1, 1, 35, 1.35, 1.25, 0.04, 0.04);
  mg->Draw("a");
  legUncert->AddEntry(gdNdEta, "stat.", "LE");
  legUncert->AddEntry(gdNdEtaMB, "syst.", "F");
  // Dummy histogram to display syst. ucorr. uncertainty gray box
  TGraphAsymmErrors *gDummy = new TGraphAsymmErrors(1);
  gDummy->SetMarkerColor(color[0]);
  gDummy->SetMarkerStyle(MarkerMult[0]);
  gDummy->SetFillStyle(3001);
  gDummy->SetFillColor(color[0]);
  gDummy->SetMarkerSize(1.5);
  mg->Add(gDummy, "P 5");
  legUncert->AddEntry(gDummy, "syst. uncorr.", "F");
  legUncert->Draw("same");
  legendTitleIntegrateddNdEta->Draw("same");
  legendIntegrateddNdEta->Draw("same");
  outputfile->cd();
  canvasIntegratedYieldsdNdEta->Write();
  //mg->Write();
  // End of Code
  std::cout << "\x1B[1;32m"; // Set text color to green
  std::cout << "\n************* Spectra in Multiplicity Classes is Fitted! *************\n";
  std::cout << "\x1B[0m"; // Reset text color
}

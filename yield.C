#include "help.h"

Bool_t reject = kTRUE;
Double_t fparab(Double_t *x, Double_t *par)
{
  const Int_t numPart = 2;
  // Signal region for the BG estimation
  Double_t liminf[numPart] = {1.308, 1.665};
  Double_t limsup[numPart] = {1.332, 1.677};
  Int_t part = par[3];
  if (reject && x[0] > liminf[part] && x[0] < limsup[part])
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
}

Double_t fpol1(Double_t *x, Double_t *par)
{
  const Int_t numPart = 2;
  // Signal region for the BG estimation
  Double_t liminf[numPart] = {1.308, 1.665};
  Double_t limsup[numPart] = {1.332, 1.677};
  Int_t part = par[2];
  if (reject && x[0] > liminf[part] && x[0] < limsup[part])
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0];
}

void yield(const Int_t nParticle = 2, // 0-2 : xi, 3-5 : omega
           const Int_t nMultInterval = 0, // 0 : MB, 1-9 - multiplicity intervals
           const Int_t inel = 0, // inel > N (0/1)
           const Bool_t isMC = 0, // 0 - data, 1 - MC
           const TString fileData = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisPostSQM/23jul-lhc22o-pass6-test/AnalysisResults.root", 
           const TString fileDataPP = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisPostSQM/run3_13tev/xi/lhc22o_pass6/2024-07-23/AnalysisResults.root",           
           const TString outputDir = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisPostSQM",
           const TString postFix = "")
{
  //gROOT->SetBatch(kTRUE);
  // Start of Code
  std::cout << "\x1B[1;33m"; // Set text color to yellow
  std::cout << "\n************ Starting Signal Extraction ************\n";
  std::cout << "\x1B[0m"; // Reset text color

  cout << "Signal extraction is performed for \e[35m" << particleNames[nParticle] << "\e[39m" <<  endl;
  if (nMultInterval > 0)
  {
    cout << "in centrality interval\e[35m " << multiplicityPerc[nMultInterval - 1] << " - " << multiplicityPerc[nMultInterval] << "\e[39m" << endl;
  } else {
    cout << "in\e[35m MB \e[39m" << endl;
  }
  cout << "for\e[35m INEL > " << inel << "\e[39m events" << endl;

  // Files with events statistics
  TFile* fileDataIn = TFile::Open(fileData);
  if (!fileDataIn || fileDataIn->IsZombie()) {
      std::cerr << "Error opening input data file!" << std::endl;
      return;
  }

  TFile* fileDataPPIn = TFile::Open(fileDataPP);
  if (!fileDataPPIn || fileDataPPIn->IsZombie()) {
      std::cerr << "Error opening input data PP file!" << std::endl;
      return;
  }

  // Output File
  TString outputfilePath;
  if (nMultInterval > 0)
  {
    outputfilePath = outputDir + "/yieldsOut/" + "yield_" + particleNames[nParticle] + "_" + multiplicityPerc[nMultInterval - 1] + "-" + multiplicityPerc[nMultInterval] + "_inel" + inel + postFix + ".root";
  } else {
    outputfilePath = outputDir + "/yieldsOut/" + "yield_" + particleNames[nParticle] + "_MB_inel" + inel + postFix + ".root";
  }
  TFile *outputfile = new TFile(outputfilePath, "RECREATE");
  TDirectory* dirBaseHists = outputfile->mkdir("baseHists");
  TDirectory* dirInvMassHists = outputfile->mkdir("invMassHists");
  TDirectory* dirFitParams = outputfile->mkdir("fitParams");
  fileDataIn->cd();

  TDirectory* cascqaanalysisDir = fileDataIn->GetDirectory("lf-cascqaanalysis");
  if (!cascqaanalysisDir)
  {
    std::cerr << "`lf-cascqaanalysis` directory is not found!" << std::endl;
    return;
  }

  cascqaanalysisDir->cd();

  TH1F* hNEvents = (TH1F *)cascqaanalysisDir->Get("hNEvents");
  if (!hNEvents)
  {
    std::cerr << "Histogram `hNEvents` is not found!" << std::endl;
    return;
  }

  // hNEvents histogram (event selection statistics for reconstructed part)
  dirBaseHists->cd();
  TCanvas *canvasNEvents = new TCanvas("canvasNEvents","canvasNEvents", 800, 600);
  canvasNEvents->cd();
  hNEvents->GetXaxis()->SetRange(1,12);
  hNEvents->SetTitle("");
  hNEvents->Draw("HIST TEXT0");
  //StyleHistoLight(hNEvents);
  canvasNEvents->Write();

  long int nEvents = 0;
  // Counts in multiplicity percentile
  TH2F* hCentFT0M_rec = (TH2F *)cascqaanalysisDir->Get("hCentFT0M_rec");
  if (!hCentFT0M_rec)
  {
    std::cerr << "Histogram `hCentFT0M_rec` is not found!" << std::endl;
    return;
  }

  hCentFT0M_rec->GetYaxis()->SetRange(2 + inel, 3); // INEL>0
  hCentFT0M_rec->GetXaxis()->SetRangeUser(0., 100.);

  TH1F* hCentFT0M_rec_1D = (TH1F*)hCentFT0M_rec->ProjectionX();

  Int_t multLowBin = 0;
  Int_t multUpBin = 0;

  if (nMultInterval == 0) { // MB
    multLowBin = hCentFT0M_rec_1D->GetXaxis()->FindBin(1e-6);
    multUpBin = hCentFT0M_rec_1D->GetXaxis()->FindBin(100 - 1e-6);
  } else { // in class
    multLowBin = hCentFT0M_rec_1D->GetXaxis()->FindBin(multiplicityPerc[nMultInterval - 1] + 1e-6);
    multUpBin = hCentFT0M_rec_1D->GetXaxis()->FindBin(multiplicityPerc[nMultInterval] - 1e-6);
  }

  // Count number of Events in mult. class (or MB)
  for (Int_t b = multLowBin; b <= multUpBin; b++)
  {
      nEvents += hCentFT0M_rec_1D->GetBinContent(b);
  }

  hCentFT0M_rec_1D->GetXaxis()->SetRange(multLowBin, multUpBin);
  dirBaseHists->cd();
  TCanvas *canvasFT0M = new TCanvas("canvasFT0M","canvasFT0M", 800, 600);
  canvasFT0M->cd();
  StyleCanvas(canvasFT0M, 0.14, 0.05, 0.11, 0.15);
  NormalizeHistogram(hCentFT0M_rec_1D); // MAKE SURE IT'S NOT USED AFTER NORMALIZATION
  StyleHisto(hCentFT0M_rec_1D, 0.5, 1.2 * hCentFT0M_rec_1D->GetBinContent(hCentFT0M_rec_1D->GetMaximumBin()), color[0], MarkerMult[0], "FT0M Multiplicity percentile", "Counts/Bin Width", "EventsFT0M", 0, 0, 0, 1.0, 1.25, 0.5, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  hCentFT0M_rec_1D->Draw();
  canvasFT0M->Write();

  cout << "number of counts in chosen mult percentile: " << nEvents << endl;

  // Setup pt binning
  Int_t numPtBins = 0;
  Int_t partType = -1;
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

  TDirectory* cascpostprocessingDir = fileDataPPIn->GetDirectory("lf-cascpostprocessing");
  if (!cascpostprocessingDir)
  {
    std::cerr << "`lf-cascpostprocessing` directory is not found!" << std::endl;
    return;
  }

  TH3F* hCascMinusInvMassvsPt_FT0M = (TH3F *)cascpostprocessingDir->Get("hCascMinusInvMassvsPt_FT0M");
  TH3F* hCascPlusInvMassvsPt_FT0M = (TH3F *)cascpostprocessingDir->Get("hCascPlusInvMassvsPt_FT0M");

  if (!hCascMinusInvMassvsPt_FT0M || !hCascPlusInvMassvsPt_FT0M)
  {
    std::cerr << "`hCascMinusInvMassvsPt_FT0M` or `hCascPlusInvMassvsPt_FT0M` histogram is not found!" << std::endl;
    return;
  }

  TH3F* hInvMass3D_minus_FT0M = (TH3F*)hCascMinusInvMassvsPt_FT0M->Clone("hInvMass3D_minus_FT0M");
  TH3F* hInvMass3D_plus_FT0M = (TH3F*)hCascPlusInvMassvsPt_FT0M->Clone("hInvMass3D_plus_FT0M");
  TH3F* hInvMass3D;

  if (nParticle == 2 || nParticle == 5) {
    hInvMass3D = (TH3F*)hInvMass3D_minus_FT0M->Clone("hInvMass3D_combined_FT0M");
    hInvMass3D->Add(hInvMass3D_plus_FT0M);
    hInvMass3D->SetTitle("hInvMass3D_combined_FT0M");
  } else if (nParticle == 0 || nParticle == 3) {
    hInvMass3D = (TH3F*)hInvMass3D_minus_FT0M->Clone("hInvMass3D_minus_FT0M_clone");
    hInvMass3D->SetTitle("hInvMass3D_minus_FT0M");
  } else if (nParticle == 1 || nParticle == 4) {
    hInvMass3D = (TH3F*)hInvMass3D_plus_FT0M->Clone("hInvMass3D_plus_FT0M_clone");
    hInvMass3D->SetTitle("hInvMass3D_plus_FT0M");
  }

  Double_t minRangeSignal[2] = {1.318, 1.668}; // to constrain mean
  Double_t maxRangeSignal[2] = {1.328, 1.677};
  Double_t minRange[2] =       {1.29, 1.655}; // range of total fit
  Double_t maxRange[2] =       {1.35, 1.690};
  Double_t minRangeDisplay[2] =       {1.28, 1.6}; // range of displayed inv mass
  Double_t maxRangeDisplay[2] =       {1.36, 1.7};
  Double_t minRangeDisplayTrue[2] =       {pdgMass[0] - 0.015, pdgMass[1] - 0.015}; // range of displayed true counts and fits
  Double_t maxRangeDisplayTrue[2] =       {pdgMass[0] + 0.015, pdgMass[1] + 0.015};

  TH2F* hCascadePlusTrue;
  TH2F* hCascadeMinusTrue;
  TH2F* hCascadeSumTrue;
  TH2F* hCascadeSumTrueClone;
  TH1D* hProjPtCascadeRecSumClone;
  TH1D* hProjPtCascadeRecSumCloneInPtBins[numPtBins];

  TH2F* hCascadePlusBG;
  TH2F* hCascadeMinusBG;
  TH2F* hCascadeBGSum;
  TH2F* hCascadeBGSumClone;
  TH1D* hProjPtCascadeBGSumClone;
  TH1D* hProjPtCascadeBGSumCloneInPtBins[numPtBins];

  TDirectory* dirMC;
  TCanvas *canvasInvMassFitsSummaryBG[2];
  TCanvas *canvasInvMassFitsSummaryTrueRec[2];

  // MC QC (BG/TRUE) //
  if (isMC) {
    dirMC = outputfile->mkdir("MC");
    dirMC->cd();
    
    for (Int_t c = 0; c < 2; c++)
    {
      canvasInvMassFitsSummaryBG[c] = new TCanvas(Form("canvasInvMassFitsSummaryBG_%i", c), Form("canvasInvMassFitsSummaryBG_%i", c), 1800, 1400);
      canvasInvMassFitsSummaryBG[c]->Divide(3, 3);
      StyleCanvas(canvasInvMassFitsSummaryBG[c], 0.15, 0.05, 0.05, 0.15);
      for (Int_t i = 1; i <= 9; i++) {
        canvasInvMassFitsSummaryBG[c]->cd(i)->SetMargin(0.12, 0.05, 0.11, 0.1);
      }
    }

    for (Int_t c = 0; c < 2; c++)
    {
      canvasInvMassFitsSummaryTrueRec[c] = new TCanvas(Form("canvasInvMassFitsSummaryTrueRec_%i", c), Form("canvasInvMassFitsSummaryTrueRec_%i", c), 1800, 1400);
      canvasInvMassFitsSummaryTrueRec[c]->Divide(3, 3);
      StyleCanvas(canvasInvMassFitsSummaryTrueRec[c], 0.15, 0.05, 0.05, 0.15);
      for (Int_t i = 1; i <= 9; i++) {
        canvasInvMassFitsSummaryTrueRec[c]->cd(i)->SetMargin(0.12, 0.05, 0.11, 0.1);
      }
    }

    Double_t leftRap = -0.5;
    Double_t rightRap = 0.5;

    // Reconstructed cascades (only true ones)
    hCascadePlusTrue = (TH2F*)fileDataPPIn->Get("lf-cascpostprocessing/hCascPlusMassvsPtTrueRec");
    hCascadeMinusTrue = (TH2F*)fileDataPPIn->Get("lf-cascpostprocessing/hCascMinusMassvsPtTrueRec");
    hCascadeSumTrue = (TH2F*)hCascadePlusTrue->Clone("SumRecTrue");
    hCascadeSumTrue->Add(hCascadeMinusTrue);
    // Sum //
    hCascadeSumTrue->GetXaxis()->SetRangeUser(binpt[0], binpt[numPtBins]); // Set pt interval
    hCascadeSumTrue->GetYaxis()->SetRangeUser(minRangeDisplay[partType], maxRangeDisplay[partType]);
    hCascadeSumTrue->Write();
    // Plus //
    hCascadePlusTrue->GetXaxis()->SetRangeUser(binpt[0], binpt[numPtBins]); // Set pt interval
    hCascadePlusTrue->GetYaxis()->SetRangeUser(minRangeDisplay[partType], maxRangeDisplay[partType]);
    hCascadePlusTrue->Write();
    // Minus //
    hCascadeMinusTrue->GetXaxis()->SetRangeUser(binpt[0], binpt[numPtBins]); // Set pt interval
    hCascadeMinusTrue->GetYaxis()->SetRangeUser(minRangeDisplay[partType], maxRangeDisplay[partType]);
    hCascadeMinusTrue->Write();
    // 1D Sum //
    hCascadeSumTrueClone = (TH2F* )hCascadeSumTrue->Clone("SumRecTrueClone");
    hProjPtCascadeRecSumClone = (hCascadeSumTrueClone->ProjectionX());
    hProjPtCascadeRecSumClone->Sumw2();
    hProjPtCascadeRecSumClone->Write();

    // BG counts
    hCascadePlusBG = (TH2F*)fileDataPPIn->Get("lf-cascpostprocessing/hCascPlusMassvsPtBG");
    hCascadeMinusBG = (TH2F*)fileDataPPIn->Get("lf-cascpostprocessing/hCascMinusMassvsPtBG");
    hCascadeBGSum = (TH2F*)hCascadePlusBG->Clone("SumBG");
    hCascadeBGSum->Add(hCascadeMinusBG);
    // Sum //
    hCascadeBGSum->GetXaxis()->SetRangeUser(binpt[0], binpt[numPtBins]); // Set pt interval
    hCascadeBGSum->GetYaxis()->SetRangeUser(minRangeDisplay[partType], maxRangeDisplay[partType]);
    hCascadeBGSum->Write();
    // Plus //
    hCascadePlusBG->GetXaxis()->SetRangeUser(binpt[0], binpt[numPtBins]); // Set pt interval
    hCascadePlusBG->GetYaxis()->SetRangeUser(minRangeDisplay[partType], maxRangeDisplay[partType]);
    hCascadePlusBG->Write();
    // Minus //
    hCascadeMinusBG->GetXaxis()->SetRangeUser(binpt[0], binpt[numPtBins]); // Set pt interval
    hCascadeMinusBG->GetYaxis()->SetRangeUser(minRangeDisplay[partType], maxRangeDisplay[partType]);
    hCascadeMinusBG->Write();
    // 1D Sum //
    hCascadeBGSumClone = (TH2F* )hCascadeBGSum->Clone("SumBGClone");
    hProjPtCascadeBGSumClone = (hCascadeBGSumClone->ProjectionX());
    hProjPtCascadeBGSumClone->Sumw2();
    hProjPtCascadeBGSumClone->Write();
  }

  Int_t multLowBinInvMass = 0;
  Int_t multUpBinInvMass = 0;
  if (nMultInterval == 0) {
    multLowBinInvMass = hInvMass3D->GetXaxis()->FindBin(1e-6);
    multUpBinInvMass = hInvMass3D->GetXaxis()->FindBin(100 - 1e-6);
  } else {
    multLowBinInvMass = hInvMass3D->GetXaxis()->FindBin(multiplicityPerc[nMultInterval - 1] + 1e-6);
    multUpBinInvMass = hInvMass3D->GetXaxis()->FindBin(multiplicityPerc[nMultInterval] - 1e-6);
  }
  hInvMass3D->GetXaxis()->SetRange(multLowBinInvMass, multUpBinInvMass);
  
  hInvMass3D->GetYaxis()->SetRangeUser(binpt[0], binpt[numPtBins]);
  hInvMass3D->GetZaxis()->SetRangeUser(minRangeDisplay[partType], maxRangeDisplay[partType]);
  TH2F* hInvMass2D = static_cast<TH2F*>(hInvMass3D->Project3D("yz"));
  TH2F* hInvMass2Dpt[numPtBins];
  TH1F* hInvMass1Dpt[numPtBins];

  dirBaseHists->cd();
  hInvMass3D->Write();
  hInvMass3D_minus_FT0M->Write();
  hInvMass3D_plus_FT0M->Write();
  hInvMass2D->Write();

  TCanvas *canvasInvMass[numPtBins];
  TF1 *total[numPtBins];
  TF1 *bkgfunc[numPtBins]; // bg function without signal region
  TF1 *signalFirstDraw[numPtBins]; // first gauss
  TF1 *signalSecondDraw[numPtBins]; // first gauss
  TF1 *bkgfuncDraw[numPtBins]; // bg function in full range to be drawn
  // fFitResultTotal is used for cov matrix
  TFitResultPtr fFitResultTotal[numPtBins];
  // fFitResultBKG is used for bg integral error
  TFitResultPtr fFitResultBKG[numPtBins];
  
  // Create mean, sigma and yield hsitograms
  TH1F *hMeans = new TH1F("Mean_" + particleNames[nParticle], "Mean_" + particleNames[nParticle], numPtBins, binpt);
  TH1F *hSigmas = new TH1F("Sigma_" + particleNames[nParticle], "Sigma_" + particleNames[nParticle], numPtBins, binpt);
  TH1F *hYields = new TH1F("Yield_" + particleNames[nParticle], "Yield_" + particleNames[nParticle], numPtBins, binpt);
  TH1F *hPurity = new TH1F("Purity_" + particleNames[nParticle], "Purity_" + particleNames[nParticle], numPtBins, binpt);

  // Reconstruction QC histograms
  TH1F *hRecCounts = new TH1F("hRecCounts", "hRecCounts", numPtBins, binpt);
  TH1F *hGenCounts = new TH1F("hGenCounts", "hGenCounts", numPtBins, binpt);

  Double_t meanYLow[2] = {1.315, 1.668};
  Double_t meanYUp[2] = {1.328, 1.677};
  Double_t sigmaYLow[2] = {0.0005, 0.0005};
  Double_t sigmaYUp[2] = {0.004, 0.01};
  Double_t yieldYLow[2] = {1e-9, 1e-9};
  Double_t yieldYUp[2] = {0.1, 0.01};

  Double_t mean[numPtBins];
  Double_t errMean[numPtBins];
  Double_t sigma[numPtBins];
  Double_t errSigma[numPtBins];
  Double_t yieldWithBG[numPtBins];
  Double_t yieldBG[numPtBins];
  Double_t errYieldBG[numPtBins];

  TCanvas *canvasInvMassFitsSummary[2];
  for (Int_t c = 0; c < 2; c++)
  {
    canvasInvMassFitsSummary[c] = new TCanvas(Form("canvasInvMassFitsSummary_%i", c), Form("canvasInvMassFitsSummary_%i", c), 1800, 1400);
    canvasInvMassFitsSummary[c]->Divide(3, 3);
    StyleCanvas(canvasInvMassFitsSummary[c], 0.15, 0.05, 0.05, 0.15);
    for (Int_t i = 1; i <= 9; i++) {
      canvasInvMassFitsSummary[c]->cd(i)->SetMargin(0.12, 0.05, 0.11, 0.1);
    }
  }

  Int_t nSigma = 6;

  Int_t fitFunc = 2; // 1 - pol1, 2 - pol2

  for (Int_t iPt = 0; iPt < numPtBins; iPt++) {
    hInvMass2Dpt[iPt] = (TH2F*)hInvMass2D->Clone(Form("2DHistInPtBin_%d", iPt));
    Int_t lowPtBin = hInvMass2Dpt[iPt]->GetYaxis()->FindBin(binpt[iPt] + 1e-6);
    Int_t upPtBin = hInvMass2Dpt[iPt]->GetYaxis()->FindBin(binpt[iPt+1] - 1e-6);
    hInvMass2Dpt[iPt]->GetYaxis()->SetRange(lowPtBin, upPtBin);
    hInvMass2Dpt[iPt]->GetXaxis()->SetRangeUser(minRangeDisplay[partType], maxRangeDisplay[partType]);
    hInvMass1Dpt[iPt] = (TH1F *)hInvMass2Dpt[iPt]->ProjectionX();
    hInvMass1Dpt[iPt]->Sumw2();
    //hInvMass1Dpt[iPt]->Rebin(4);
    hInvMass1Dpt[iPt]->GetYaxis()->SetTitle("Counts");
    hInvMass1Dpt[iPt]->SetTitle(Form("%.2f < #it{p}_{T} < %.2f GeV/#it{c}", binpt[iPt], binpt[iPt+1]));
    hInvMass1Dpt[iPt]->GetXaxis()->SetTitle(titleInvMass[nParticle] + " " + sInvMass);
    hInvMass1Dpt[iPt]->GetYaxis()->SetMaxDigits(3);
    hInvMass1Dpt[iPt]->GetYaxis()->SetDecimals(kTRUE);
    StyleHistoLight(hInvMass1Dpt[iPt]);
    // Fit inv mass histogram
    canvasInvMass[iPt] = new TCanvas(Form("%s_%d", hInvMass3D->GetName(), iPt), Form("%s_%d", hInvMass3D->GetName(), iPt), 1000, 800);
    StyleCanvas(canvasInvMass[iPt], 0.14, 0.05, 0.11, 0.15);
    canvasInvMass[iPt]->cd();

    if (fitFunc == 1) {
      bkgfunc[iPt] = new TF1(Form("pol1_%d", iPt), fpol1, minRange[partType], maxRange[partType], 3);
      bkgfunc[iPt]->FixParameter(2, partType); // Fix Particle type to reject signal range
    } else if (fitFunc == 2) {
      bkgfunc[iPt] = new TF1(Form("parab_%d", iPt), fparab, minRange[partType], maxRange[partType], 4);
      bkgfunc[iPt]->FixParameter(3, partType); // Fix Particle type to reject signal range
    }

    bkgfunc[iPt]->SetLineColor(kGreen);

    fFitResultBKG[iPt] = hInvMass1Dpt[iPt]->Fit(bkgfunc[iPt], "SRQL0");
    Double_t parBG[4]; // maximum number of parameters is 4
    bkgfunc[iPt]->GetParameters(parBG);
    if (fitFunc == 1) {
      total[iPt] = new TF1(Form("total_%d", iPt), "gaus(0) + gaus(3) + pol1(6)", minRange[partType], maxRange[partType]);
    } else if (fitFunc == 2) {
      total[iPt] = new TF1(Form("total_%d", iPt), "gaus(0) + gaus(3) + pol2(6)", minRange[partType], maxRange[partType]);
    }

    Double_t peakValue = hInvMass1Dpt[iPt]->GetBinContent(hInvMass1Dpt[iPt]->GetMaximumBin());
    hInvMass1Dpt[iPt]->GetYaxis()->SetRangeUser(0, 1.1 * peakValue);
    total[iPt]->SetLineColor(kRed);
    total[iPt]->SetParName(0, "norm");
    total[iPt]->SetParName(1, "mean");
    total[iPt]->SetParName(2, "sigma");
    total[iPt]->SetParName(3, "norm2");
    total[iPt]->SetParName(4, "mean2");
    total[iPt]->SetParName(5, "sigma2");
    total[iPt]->SetParName(6, "p0");
    total[iPt]->SetParName(7, "p1");
    if (fitFunc == 2) {
      total[iPt]->SetParName(8, "p2");
    }
    total[iPt]->SetParameters(peakValue, pdgMass[partType], 0.001, peakValue, pdgMass[partType], 0.001);
    total[iPt]->SetParLimits(0, 0., 1.2 * peakValue);
    total[iPt]->SetParLimits(1, minRangeSignal[partType], maxRangeSignal[partType]);
    total[iPt]->SetParLimits(2, 0.001, 0.01);
    total[iPt]->SetParLimits(3, 0., 1.2 * peakValue);
    total[iPt]->SetParLimits(4, minRangeSignal[partType], maxRangeSignal[partType]);
    total[iPt]->SetParLimits(5, 0.001, 0.05);
    //total[iPt]->FixParameter(6, 1.05*parBG[0]); // better reconstruction for MC
    total[iPt]->FixParameter(6, parBG[0]);
    total[iPt]->FixParameter(7, parBG[1]);

    if (fitFunc == 2) {
      total[iPt]->FixParameter(8, parBG[2]);
    }

    total[iPt]->SetNpx(1e4);
    total[iPt]->SetLineWidth(3);
    fFitResultTotal[iPt] = hInvMass1Dpt[iPt]->Fit(total[iPt], "SRQBL");
    Double_t parFinal[fitFunc + 7]; // maximum number of parameters is 9
    total[iPt]->GetParameters(parFinal);

    if (fitFunc == 1) {
      bkgfuncDraw[iPt] = new TF1(Form("pol1ToDraw_%d", iPt), "pol1", minRange[partType], maxRange[partType]);
    } else if (fitFunc == 2) {
      bkgfuncDraw[iPt] = new TF1(Form("parabToDraw_%d", iPt), "pol2", minRange[partType], maxRange[partType]);
      bkgfuncDraw[iPt]->FixParameter(2, parFinal[8]);
    }
    
    bkgfuncDraw[iPt]->FixParameter(0, parFinal[6]);
    bkgfuncDraw[iPt]->FixParameter(1, parFinal[7]);
    
    bkgfuncDraw[iPt]->SetLineStyle(2);
    bkgfuncDraw[iPt]->SetNpx(1e4);
    bkgfuncDraw[iPt]->SetLineWidth(3);
    bkgfuncDraw[iPt]->SetLineColor(kCyan+2); //{kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};

    signalFirstDraw[iPt] = new TF1(Form("signalToDraw1_%d", iPt), "gaus", minRange[partType], maxRange[partType]);
    signalFirstDraw[iPt]->FixParameter(0, total[iPt]->GetParameter(0));
    signalFirstDraw[iPt]->FixParameter(1, total[iPt]->GetParameter(1));
    signalFirstDraw[iPt]->FixParameter(2, total[iPt]->GetParameter(2));
    signalFirstDraw[iPt]->SetLineStyle(2);
    signalFirstDraw[iPt]->SetNpx(1e4);
    signalFirstDraw[iPt]->SetLineWidth(3);
    signalFirstDraw[iPt]->SetLineColor(kGreen+3); //{kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
    signalFirstDraw[iPt]->Draw("same");

    signalSecondDraw[iPt] = new TF1(Form("signalToDraw2_%d", iPt), "gaus", minRange[partType], maxRange[partType]);
    signalSecondDraw[iPt]->FixParameter(0, total[iPt]->GetParameter(3));
    signalSecondDraw[iPt]->FixParameter(1, total[iPt]->GetParameter(4));
    signalSecondDraw[iPt]->FixParameter(2, total[iPt]->GetParameter(5));
    signalSecondDraw[iPt]->SetLineStyle(2);
    signalSecondDraw[iPt]->SetNpx(1e4);
    signalSecondDraw[iPt]->SetLineWidth(3);
    signalSecondDraw[iPt]->SetLineColor(kBlack); //{kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
    signalSecondDraw[iPt]->Draw("same");

    // Uncertainties of double gauss fit
    Double_t N1a = total[iPt]->GetParameter(0);      // N1
    Double_t N2a = total[iPt]->GetParameter(3);      // N2
    Double_t mu1a = total[iPt]->GetParameter(1);     // mu1
    Double_t mu2a = total[iPt]->GetParameter(4);     // mu2
    Double_t sigma1a = total[iPt]->GetParameter(2);  // sigma1
    Double_t sigma2a = total[iPt]->GetParameter(5);  // sigma2
    TMatrixD cova = fFitResultTotal[iPt]->GetCovarianceMatrix();
    Double_t mu_wa = (N1a*mu1a + N2a*mu2a)/(N1a+N2a);
    Double_t sigma_wa = (N1a*sigma1a + N2a*sigma2a)/(N1a+N2a);

    Double_t sa = N1a + N2a;
    Double_t wa_mu = N1a*mu1a + N2a*mu2a;
    Double_t wa_sigma = N1a*sigma1a + N2a*sigma2a; 

    Double_t mu_wa_step = pow((mu1a - mu2a),2)*(pow(N1a,2)*cova(3,3) + pow(N2a,2)*cova(0,0))
                + 2*cova(0,3)*(wa_mu-sa*mu1a)*(wa_mu-sa*mu2a)
                + pow(sa,2)*(pow(N1a,2)*cova(1,1) + pow(N2a,2)*cova(4,4) + 2*N1a*N2a*cova(1,4))
                - 2*sa * (N1a * (cova(0,1)*(wa_mu-sa*mu1a) + cova(3,1)*(wa_mu-sa*mu2a)) + 
                      N2a * (cova(0,4)*(wa_mu-sa*mu1a) + cova(3,4)*(wa_mu-sa*mu2a)));

    mean[iPt] = mu_wa;
    errMean[iPt] = sqrt(mu_wa_step / pow(sa,4));

    Double_t sigma_wa_step = pow((sigma1a - sigma2a),2)*(pow(N1a,2)*cova(3,3) + pow(N2a,2)*cova(0,0))
                + 2*cova(0,3)*(wa_sigma-sa*sigma1a)*(wa_sigma-sa*sigma2a)
                + pow(sa,2)*(pow(N1a,2)*cova(2,2) + pow(N2a,2)*cova(5,5) + 2*N1a*N2a*cova(2,5))
                - 2*sa * (N1a * (cova(0,2)*(wa_sigma-sa*sigma1a) + cova(3,2)*(wa_sigma-sa*sigma2a)) + 
                      N2a * (cova(0,5)*(wa_sigma-sa*sigma1a) + cova(3,5)*(wa_sigma-sa*sigma2a)));

    sigma[iPt] = sigma_wa;
    errSigma[iPt] = sqrt(sigma_wa_step / pow(sa,4));

    // cout << "Analysed pt interval: " << binpt[iPt] << "-" << binpt[iPt+1] << endl;
    // cout << "sigma1: " << total[iPt]->GetParameter(2) << " err1: " << total[iPt]->GetParError(2) << " ampl1: " <<  total[iPt]->GetParameter(0) << std::endl;
    // cout << "sigma2: " << total[iPt]->GetParameter(5) << " err2: " << total[iPt]->GetParError(5) << " ampl2: " <<  total[iPt]->GetParameter(3) <<  std::endl;
    // cout << "total sigma: " << sigma[iPt] << " sigma err: " << errSigma[iPt] <<  std::endl;
    // cout << "total mean: " << mean[iPt] << " mean err: " << errMean[iPt] <<  std::endl;

    Double_t leftSignal = mean[iPt] - nSigma * sigma[iPt];
    Double_t rightSignal = mean[iPt] + nSigma * sigma[iPt];

    Double_t leftSignalBin = hInvMass1Dpt[iPt]->GetXaxis()->FindBin(leftSignal + 1e-6);
    Double_t rightSignalBin = hInvMass1Dpt[iPt]->GetXaxis()->FindBin(rightSignal - 1e-6);

    Double_t leftEdgeSignal = hInvMass1Dpt[iPt]->GetXaxis()->GetBinLowEdge(leftSignalBin);
    Double_t rightEdgeSignal = hInvMass1Dpt[iPt]->GetXaxis()->GetBinUpEdge(rightSignalBin);

    //std::cout << "leftSignal - leftSignal in bin " << leftSignal - leftEdgeSignal << std::endl;

    Double_t yaxisMin = hInvMass1Dpt[iPt]->GetMinimum();
    DrawVertLine(leftEdgeSignal, yaxisMin, peakValue, kBlack);
    DrawVertLine(rightEdgeSignal, yaxisMin, peakValue, kBlack);
    DrawVertLine(minRange[partType], yaxisMin, peakValue, kMagenta+1);
    DrawVertLine(maxRange[partType], yaxisMin, peakValue, kMagenta+1);

    gPad->Update();
    canvasInvMass[iPt]->Update();
    outputfile->cd();
    dirInvMassHists->cd();
    bkgfuncDraw[iPt]->Draw("same");
    canvasInvMass[iPt]->Write();
    delete canvasInvMass[iPt];

    if (iPt < 9)
      canvasInvMassFitsSummary[0]->cd(iPt + 1);
    else if (iPt < 18)
      canvasInvMassFitsSummary[1]->cd(iPt + 1 - 9);
    // else if (iPt < 27)
    //   canvasInvMassFitsSummary[2]->cd(iPt + 1 - 18);
    // else
    //   canvasInvMassFitsSummary[3]->cd(iPt + 1 - 27);
    hInvMass1Dpt[iPt]->Draw();
    DrawVertLine(leftEdgeSignal, yaxisMin, peakValue, kBlack);
    DrawVertLine(rightEdgeSignal, yaxisMin, peakValue, kBlack);
    DrawVertLine(minRange[partType], yaxisMin, peakValue, kMagenta+1);
    DrawVertLine(maxRange[partType], yaxisMin, peakValue, kMagenta+1);
    bkgfuncDraw[iPt]->Draw("same");
    signalFirstDraw[iPt]->Draw("same");
    signalSecondDraw[iPt]->Draw("same");

    hMeans->SetBinContent(iPt + 1, mean[iPt]);
    hMeans->SetBinError(iPt + 1, errMean[iPt]);

    hSigmas->SetBinContent(iPt + 1, sigma[iPt]);
    hSigmas->SetBinError(iPt + 1, errSigma[iPt]);

    yieldBG[iPt] = bkgfuncDraw[iPt]->Integral(leftEdgeSignal, rightEdgeSignal)/hInvMass1Dpt[iPt]->GetBinWidth(2);
    errYieldBG[iPt] = bkgfuncDraw[iPt]->IntegralError(leftEdgeSignal, rightEdgeSignal, fFitResultBKG[iPt]->GetParams(), (fFitResultBKG[iPt]->GetCovarianceMatrix()).GetMatrixArray());

    yieldWithBG[iPt] = 0;
    for (Int_t bin = leftSignalBin; bin <= rightSignalBin; bin++) {
      yieldWithBG[iPt] += hInvMass1Dpt[iPt]->GetBinContent(bin);
    }

    hYields->SetBinContent(iPt + 1, (yieldWithBG[iPt] - yieldBG[iPt]) / hYields->GetBinWidth(iPt + 1) / nEvents);
    hYields->SetBinError(iPt + 1, sqrt(yieldWithBG[iPt] + pow(errYieldBG[iPt], 2)) / hYields->GetBinWidth(iPt + 1) / nEvents);

    hPurity->SetBinContent(iPt + 1, (yieldWithBG[iPt] - yieldBG[iPt]) / (yieldWithBG[iPt]));
    Double_t num = yieldWithBG[iPt] - yieldBG[iPt];
    Double_t errNum = sqrt(yieldWithBG[iPt] + pow(errYieldBG[iPt], 2));
    hPurity->SetBinError(iPt + 1, hPurity->GetBinContent(iPt + 1) * sqrt(pow(errNum/num,2) + pow(errYieldBG[iPt]/yieldWithBG[iPt], 2)));

    if (isMC) {
      dirMC->cd();
      TH2F* hCascadeSumTrueClone2 = (TH2F* )hCascadeSumTrue->Clone(Form("SumTrue_%i", iPt));
      hCascadeSumTrueClone2->GetXaxis()->SetRangeUser(binpt[iPt], binpt[iPt+1]);
      hProjPtCascadeRecSumCloneInPtBins[iPt] = (TH1D* )hCascadeSumTrueClone2->ProjectionY();
      if (iPt < 9)
        canvasInvMassFitsSummaryTrueRec[0]->cd(iPt + 1);
      else if (iPt < 18)
        canvasInvMassFitsSummaryTrueRec[1]->cd(iPt + 1 - 9);
      hProjPtCascadeRecSumCloneInPtBins[iPt]->GetYaxis()->SetTitle("Counts");
      hProjPtCascadeRecSumCloneInPtBins[iPt]->SetTitle(Form("%.2f < #it{p}_{T} < %.2f GeV/#it{c}", binpt[iPt], binpt[iPt+1]));
      hProjPtCascadeRecSumCloneInPtBins[iPt]->GetXaxis()->SetTitle(titleInvMass[nParticle] + " " + sInvMass);
      StyleHistoLight(hProjPtCascadeRecSumCloneInPtBins[iPt]);
      hProjPtCascadeRecSumCloneInPtBins[iPt]->GetXaxis()->SetRangeUser(minRangeDisplayTrue[partType], maxRangeDisplayTrue[partType]);
      hProjPtCascadeRecSumCloneInPtBins[iPt]->Draw();
      signalFirstDraw[iPt]->Draw("same");
      signalSecondDraw[iPt]->Draw("same");
      TF1* sumGaussians = new TF1(Form("sumGaussians_%d", iPt), "gaus(0) + gaus(3)", minRange[partType], maxRange[partType]);
      sumGaussians->FixParameter(0, total[iPt]->GetParameter(0));
      sumGaussians->FixParameter(1, total[iPt]->GetParameter(1));
      sumGaussians->FixParameter(2, total[iPt]->GetParameter(2));
      sumGaussians->FixParameter(3, total[iPt]->GetParameter(3));
      sumGaussians->FixParameter(4, total[iPt]->GetParameter(4));
      sumGaussians->FixParameter(5, total[iPt]->GetParameter(5));
      sumGaussians->SetLineStyle(1);
      sumGaussians->SetNpx(1e4);
      sumGaussians->SetLineWidth(3);
      sumGaussians->SetLineColor(kRed+1); //{kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
      sumGaussians->Draw("same");

      Double_t genCountsInPtBin = hProjPtCascadeRecSumCloneInPtBins[iPt]->Integral();
      Double_t recYieldFit = sumGaussians->Integral(leftEdgeSignal, rightEdgeSignal)/hInvMass1Dpt[iPt]->GetBinWidth(2);
      //std::cout << "true counts: " << genCountsInPtBin << " rec. counts:" << yieldWithBG[iPt] - yieldBG[iPt] << " difference (true - rec): " << genCountsInPtBin - (yieldWithBG[iPt] - yieldBG[iPt]) << std::endl;
      // hRecCounts->SetBinContent(iPt + 1, recYieldFit);
      // hRecCounts->SetBinError(iPt + 1, 1e-6);
      hRecCounts->SetBinContent(iPt + 1, yieldWithBG[iPt] - yieldBG[iPt]);
      hRecCounts->SetBinError(iPt + 1, sqrt(yieldWithBG[iPt] - yieldBG[iPt]));
      hGenCounts->SetBinContent(iPt + 1, genCountsInPtBin);
      hGenCounts->SetBinError(iPt + 1, sqrt(genCountsInPtBin));

      TH2F* hCascadeBGSumClone2 = (TH2F* )hCascadeBGSum->Clone(Form("SumBG_%i", iPt));
      hCascadeBGSumClone2->GetXaxis()->SetRangeUser(binpt[iPt], binpt[iPt+1]);
      hProjPtCascadeBGSumCloneInPtBins[iPt] = (TH1D* )hCascadeBGSumClone2->ProjectionY();
      if (iPt < 9)
        canvasInvMassFitsSummaryBG[0]->cd(iPt + 1);
      else if (iPt < 18)
        canvasInvMassFitsSummaryBG[1]->cd(iPt + 1 - 9);
      // else if (iPt < 27)
      //   canvasInvMassFitsSummaryBG[2]->cd(iPt + 1 - 18);
      // else
      //   canvasInvMassFitsSummaryBG[3]->cd(iPt + 1 - 27);
      hProjPtCascadeBGSumCloneInPtBins[iPt]->GetYaxis()->SetTitle("Counts");
      hProjPtCascadeBGSumCloneInPtBins[iPt]->SetTitle(Form("%.2f < #it{p}_{T} < %.2f GeV/#it{c}", binpt[iPt], binpt[iPt+1]));
      hProjPtCascadeBGSumCloneInPtBins[iPt]->GetXaxis()->SetTitle(titleInvMass[nParticle] + " " + sInvMass);
      StyleHistoLight(hProjPtCascadeBGSumCloneInPtBins[iPt]);
      hProjPtCascadeBGSumCloneInPtBins[iPt]->Draw();
      bkgfuncDraw[iPt]->Draw("same");
    }
  }

  outputfile->cd();
  dirInvMassHists->cd();
  for (Int_t i = 0; i < 2; i++)
  {
    canvasInvMassFitsSummary[i]->Write();
  }

  if (isMC) {
    outputfile->cd();
    dirMC->cd();
    for (Int_t i = 0; i < 2; i++)
    {
      canvasInvMassFitsSummaryBG[i]->Write();
    }

    for (Int_t i = 0; i < 2; i++)
    {
      canvasInvMassFitsSummaryTrueRec[i]->Write();
    }
  }


  dirFitParams->cd();

  StyleHisto(hMeans, meanYLow[partType], meanYUp[partType], kBlack, 20, sPt, "Mean (GeV/#it{c}^{2})", "", 0, 0, 0, 1.0, 1.25, 1.1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  //hMeans->Draw("P");
  hMeans->Write();

  StyleHisto(hSigmas, sigmaYLow[partType], sigmaYUp[partType], kBlack, 20, sPt, "Sigma (GeV/#it{c}^{2})", "", 0, 0, 0, 1.0, 1.25, 1.1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  //hSigmas->Draw("P");
  hSigmas->Write();

  StyleHisto(hPurity, 0, 1, kBlack, 20, sPt, "S/(B+S)", "", 0, 0, 0, 1.0, 1.25, 1.1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  //hPurity->Draw("P");
  hPurity->Write();

  StyleHisto(hYields, yieldYLow[partType], yieldYUp[partType], kBlack, 20, sPt, sdNdPtdY, "", 0, 0, 0, 1.0, 1.25, 1.1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  //hYields->Draw("P");
  hYields->Write();

  if (isMC) {
    dirMC->cd();
    // QC of the reconstruction procedure //
    TCanvas* canvasRecQC = new TCanvas("RecQC", "RecQC", 0, 70, 620, 850);
    StyleCanvas(canvasRecQC, 0.15, 0.05, 0.05, 0.15);
    TPad *padRecQCUp = new TPad("padRecQCUp", "padRecQCUp", 0, 0.36, 1, 1);
    TPad *padRecQCLow = new TPad("padRecQCLow", "padRecQCLow", 0, 0, 1, 0.36);
    StylePad(padRecQCUp, 0.15, 0.05, 0.05, 0.);
    StylePad(padRecQCLow, 0.15, 0.05, 0.02, 0.2);
    canvasRecQC->cd();
    padRecQCUp->Draw();
    padRecQCLow->Draw();

    // Dummy histograms
    TH1F *hDummyRecQC = new TH1F("hDummyRecQC", "hDummyRecQC", 10000, 0., binpt[numPtBins] + 0.5);
    // Up
    for (Int_t i = 1; i <= hDummyRecQC->GetNbinsX(); i++)
      hDummyRecQC->SetBinContent(i, 1e-12);
    canvasRecQC->cd();
    padRecQCUp->cd();
    hDummyRecQC->GetYaxis()->SetMaxDigits(3);
    hDummyRecQC->GetYaxis()->SetDecimals(kTRUE);
    Double_t maxCounts = (hRecCounts->GetMaximum() > hGenCounts->GetMaximum()) ? hRecCounts->GetMaximum() : hGenCounts->GetMaximum();
    StyleHisto(hDummyRecQC, 0., maxCounts*1.2, 1, 1, sPt, "Raw Counts", "", 0, 0, 0, 1.5, 1.0, 0, 0.0, 0.05, 0.0, 0.035, 0.005, 0.005);
    SetTickLength(hDummyRecQC, 0.025, 0.03);
    TAxis *axisRecQCDummy = hDummyRecQC->GetYaxis();
    axisRecQCDummy->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
    hDummyRecQC->GetXaxis()->SetRangeUser(0., binpt[numPtBins] + 0.5);
    //padRecQCUp->SetLogy();
    hDummyRecQC->Draw("same");
    // Low
    TH1F *hDummyRecQCLow = new TH1F("hDummyRecQCLow", "hDummyRecQCLow", 10000, 0., binpt[numPtBins] + 0.5);
    for (Int_t i = 1; i <= hDummyRecQCLow->GetNbinsX(); i++)
      hDummyRecQCLow->SetBinContent(i, 1);
    padRecQCLow->cd();
    StyleHisto(hDummyRecQCLow, 0.5, 1.5, 1, 1, sPt, "rec./gen.", "", 0, 0, 0, 1.0, 0.7, 0, 0.08, 0.08, 0.08, 0.07, 0.01, 0.01);
    SetTickLength(hDummyRecQCLow, 0.025, 0.03);
    TAxis *axisDummyRecQCLow = hDummyRecQCLow->GetYaxis();
    axisDummyRecQCLow->ChangeLabel(-1, -1, -1, -1, -1, -1, " "); // last 
    axisDummyRecQCLow->ChangeLabel(1, -1, -1, -1, -1, -1, " "); // first
    hDummyRecQCLow->GetXaxis()->SetRangeUser(0., binpt[numPtBins] + 0.5);
    hDummyRecQCLow->GetYaxis()->CenterTitle();
    hDummyRecQCLow->Draw("same");

    // Plotting
    padRecQCUp->cd();
    StyleHisto(hRecCounts, 1e-6, hRecCounts->GetMaximum()*1e2, color[0], MarkerMult[0], "", hRecCounts->GetYaxis()->GetTitle(), "", 0, 0, 0, 1.5, 1.0, SizeMult[0], 0.0, 0.05, 0.0, 0.035, 0.005, 0.005);
    StyleHisto(hGenCounts, 1e-6, hGenCounts->GetMaximum()*1e2, color[3], MarkerMult[3], "", hRecCounts->GetYaxis()->GetTitle(), "", 0, 0, 0, 1.5, 1.0, SizeMult[3], 0.0, 0.05, 0.0, 0.035, 0.005, 0.005);
    hRecCounts->Draw("same");
    hGenCounts->Draw("same");

    TLegend *legRecQC = new TLegend(0.45, 0.56, 0.97, 0.82);
    legRecQC->SetHeader("Signal Extraction QC");
    legRecQC->SetFillStyle(0);
    legRecQC->AddEntry(hRecCounts, "Rec.", "p");
    legRecQC->AddEntry(hGenCounts, "Gen.", "p");
    StyleLegend(legRecQC, 0.0, 0.0);
    legRecQC->SetTextSize(0.04);
    legRecQC->Draw("same");

    padRecQCLow->cd();
    TH1F *hRatioCounts = (TH1F* )hRecCounts->Clone("hRatioCounts");
    hRatioCounts->Divide(hGenCounts);
    StyleHisto(hRatioCounts, 0.1 , 11, color[0], MarkerMult[0], "", "", "", 0, 0, 0, 1.0, 0.7, 0.7, 0.08, 0.08, 0.08, 0.07, 0.005, 0.01);
    hRatioCounts->SetStats(0);
    hRatioCounts->Draw("same");
    canvasRecQC->Write();

    dirFitParams->cd();
    TH1F* hGenCountsClone = (TH1F*)hGenCounts->Clone("hGenCountsClone");
    StyleHisto(hGenCountsClone, 1e-6, hGenCountsClone->GetMaximum()*1e2, kBlack, 20, sPt, sdNdPtdY, "", 0, 0, 0, 1.0, 1.1, 1, 0.04, 0.04, 0.08, 0.07, 0.005, 0.005);
    hGenCountsClone->Write();
  }

  delete[] binpt;

  // End of Code
  std::cout << "\x1B[1;32m"; // Set text color to green
  std::cout << "\n************* Signal Extraction Ended Successfully! *************\n";
  std::cout << "\x1B[0m"; // Reset text color

  //gROOT->SetBatch(kFALSE);
}
#include "../../help.h"
#include "../../effHelp.h"

void oobrej(const Int_t nParticle = 2) // 0-2 : xi, 3-5 : omega
{
  // Start of Code
  std::cout << "\x1B[1;33m"; // Set text color to yellow
  std::cout << "\n************ Starting OOB syst ************\n";
  std::cout << "\x1B[0m"; // Reset text color

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

  // output file
  TString outputfilePath;
  outputfilePath = "../sourcesOfSyst/OOBSyst-Xi.root";
  TFile *outputfile = new TFile(outputfilePath, "RECREATE");

  TFile *fdef = TFile::Open("data/yieldsOutEffCorr/yield_XiPm_MB_inel0_DEF.root");
  TFile *foob = TFile::Open("data/yieldsOutEffCorr/yield_XiPm_MB_inel0_OOB.root");

  TH1F *hOOB = (TH1F *)foob->Get("effCorrYield/Yield_XiPm");
  hOOB->SetName("hOOB");
  TH1F *hDEF = (TH1F *)fdef->Get("effCorrYield/Yield_XiPm");
  hDEF->SetName("hDEF");

  TCanvas *canvasSyst = new TCanvas("canvasSyst","canvasSyst", 800, 600);
  canvasSyst->cd();
  StyleCanvas(canvasSyst, 0.14, 0.05, 0.11, 0.15);
  gStyle->SetOptStat(0);

  // Efficiency multiplicity dependence
  TH1D *hSystOOB = (TH1D *)hPt->Clone("hSystOOB");
  for (int i = 1; i <= hSystOOB->GetNbinsX(); i++)
  {
    std::cout << "dev: " << hOOB->GetBinContent(i) / hDEF->GetBinContent(i) << std::endl;
    hSystOOB->SetBinContent(i, abs(1 - hOOB->GetBinContent(i) / hDEF->GetBinContent(i)));
  }

  outputfile->cd();
  StyleHisto(hSystOOB, -0.01, 0.01, color[0], MarkerMult[0], sPt, "1 - #it{Y}_{ITS||TOF}/#it{Y}_{def}", "", 0, 0, 0, 1.0, 1.25, SizeMult[0], 0.04, 0.04);
  hSystOOB->GetYaxis()->SetMaxDigits(2);
  hSystOOB->GetYaxis()->SetDecimals(kTRUE);
  hSystOOB->Draw();
  hSystOOB->Write();

  // End of Code
  std::cout << "\x1B[1;32m"; // Set text color to green
  std::cout << "\n************* OOB syst Finished Successfully! *************\n";
  std::cout << "\x1B[0m"; // Reset text color 
}
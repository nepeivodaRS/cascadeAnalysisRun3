#include "../../help.h"
#include "../../effHelp.h"

void bkgfunc(const Int_t nParticle = 2) // 0-2 : xi, 3-5 : omega
{
  // Start of Code
  std::cout << "\x1B[1;33m"; // Set text color to yellow
  std::cout << "\n************ Starting bkg variation syst ************\n";
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
  outputfilePath = "../sourcesOfSyst/SystBKG-Xi.root";
  TFile *outputfile = new TFile(outputfilePath, "RECREATE");

  TFile *fdef = TFile::Open("data/yieldsOut/yield_XiPm_MB_inel0_POL2.root");
  TFile *fpol1 = TFile::Open("data/yieldsOut/yield_XiPm_MB_inel0_POL1.root");

  TH1F *hPOL1 = (TH1F *)fpol1->Get("fitParams/Yield_XiPm");
  hPOL1->SetName("hPOL1");
  TH1F *hDEF = (TH1F *)fdef->Get("fitParams/Yield_XiPm");
  hDEF->SetName("hDEF");

  TCanvas *canvasSyst = new TCanvas("canvasSyst","canvasSyst", 800, 600);
  canvasSyst->cd();
  StyleCanvas(canvasSyst, 0.14, 0.05, 0.11, 0.15);
  gStyle->SetOptStat(0);

  // Efficiency multiplicity dependence
  TH1D *hSystBKG = (TH1D *)hPt->Clone("hSystBKG");
  for (int i = 1; i <= hSystBKG->GetNbinsX(); i++)
  {
    std::cout << "dev: " << hPOL1->GetBinContent(i) / hDEF->GetBinContent(i) << std::endl;
    hSystBKG->SetBinContent(i, abs(1 - hPOL1->GetBinContent(i) / hDEF->GetBinContent(i)));
  }

  outputfile->cd();
  StyleHisto(hSystBKG, 0., 0.1, color[0], MarkerMult[0], sPt, "1 - #it{Y}_{pol1}/#it{Y}_{pol2}", "", 0, 0, 0, 1.0, 1.25, SizeMult[0], 0.04, 0.04);
  hSystBKG->GetYaxis()->SetMaxDigits(2);
  hSystBKG->GetYaxis()->SetDecimals(kTRUE);
  hSystBKG->Draw("SAME");
  hSystBKG->Write();

  // End of Code
  std::cout << "\x1B[1;32m"; // Set text color to green
  std::cout << "\n************* Bkg variation Finished Successfully! *************\n";
  std::cout << "\x1B[0m"; // Reset text color 
}
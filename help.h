#include "TStyle.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"

void StyleCanvas(TCanvas *canvas, Double_t LMargin, Double_t RMargin, Double_t TMargin, Double_t BMargin)
{
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetLeftMargin(LMargin);
  canvas->SetRightMargin(RMargin);
  canvas->SetTopMargin(TMargin);
  canvas->SetBottomMargin(BMargin);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
}

void DrawHorLine(Double_t x, Double_t y){
  TLine *line = new TLine(0., y, x, y);
  line->SetLineStyle(2); // Set the line style to dashed (2)
  line->SetLineColor(kBlack); // Set the line color (kRed is a ROOT predefined color)
  line->Draw("same"); // Draw the line on the same canvas
}

void DrawHorLineColor(Double_t x, Double_t y, Color_t color){
  TLine *line = new TLine(0., y, x, y);
  line->SetLineStyle(2); // Set the line style to dashed (2)
  line->SetLineColor(color); // Set the line color (kRed is a ROOT predefined color)
  line->Draw("same"); // Draw the line on the same canvas
}

void DrawVertLine(Double_t x, Double_t yMin, Double_t yMax, Color_t color){
  TLine *line = new TLine(x, yMin, x, yMax);
  line->SetLineStyle(2); // Set the line style to dashed (2)
  line->SetLineWidth(2); // Set the line width
  line->SetLineColor(color); // Set the line color (kRed is a ROOT predefined color)
  line->Draw("same"); // Draw the line on the same canvas
}

void StylePad(TPad *pad, Double_t LMargin, Double_t RMargin, Double_t TMargin, Double_t BMargin)
{
  pad->SetFillColor(0);
  pad->SetTickx(1);
  pad->SetTicky(1);
  pad->SetLeftMargin(LMargin);
  pad->SetRightMargin(RMargin);
  pad->SetTopMargin(TMargin);
  pad->SetBottomMargin(BMargin);
}

void StyleLegend(TLegend *legend, Double_t bordersize, Double_t fillStyle)
{
  legend->SetBorderSize(bordersize);
  legend->SetFillStyle(fillStyle);
  legend->SetTextSize(0.03);
}

void SetTickLength(TH1F *histo, Double_t TickLengthX, Double_t TickLengthY)
{
  histo->GetXaxis()->SetTickLength(TickLengthX);
  histo->GetYaxis()->SetTickLength(TickLengthY);
}

template <typename THisto>
void NormalizeHistogram(THisto hist){
  const Int_t nBins = hist->GetXaxis()->GetNbins();
  for(Int_t bin = 1; bin <= nBins; bin++) {
    const Double_t binwidth = hist->GetXaxis()->GetBinWidth(bin);
    Double_t scale = 1.0/binwidth;
    hist->SetBinContent(bin, scale*hist->GetBinContent(bin));
    hist->SetBinError(bin, scale*hist->GetBinError(bin));
  }
}

template <typename THisto>
void StyleHisto(THisto *histo, Double_t Low, Double_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, 
                Bool_t XRange, Double_t XLow, Double_t XUp, Double_t xOffset, Double_t yOffset, Double_t mSize, 
                Double_t xTitleSize, Double_t yTitleSize, Double_t xLabelSize,  Double_t yLabelSize, Double_t xLabelOffset, Double_t yLabelOffset)
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
  histo->GetXaxis()->SetLabelOffset(yLabelOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(yTitleSize);
  histo->GetYaxis()->SetLabelSize(yLabelSize);
  histo->GetYaxis()->SetTitleOffset(yOffset);
  histo->GetYaxis()->SetLabelOffset(yLabelOffset);
  histo->SetTitle(title);
}

template <typename THisto>
void Style2DHisto(THisto *histo, TString titleX, TString titleY, TString title,
  Double_t xTitleSize, Double_t yTitleSize, 
  Double_t xLabelSize, Double_t yLabelSize, Double_t zLabelSize)
{
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetYaxis()->SetTitle(titleY);
  histo->SetTitle(title);
  histo->GetXaxis()->SetTitleSize(xTitleSize);
  histo->GetYaxis()->SetTitleSize(yTitleSize);
  histo->GetXaxis()->SetLabelSize(xLabelSize);
  histo->GetYaxis()->SetLabelSize(yLabelSize);
  histo->GetZaxis()->SetLabelSize(zLabelSize);
}

template <typename TMG>
void StyleMultGraph(TMG *mg, Double_t Low, Double_t Up, TString titleX, TString titleY, TString title, 
                Bool_t XRange, Double_t XLow, Double_t XUp, Double_t xOffset, Double_t yOffset, 
                Double_t xTitleSize, Double_t yTitleSize)
{
  mg->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)
    mg->GetXaxis()->SetLimits(XLow, XUp);
  mg->GetXaxis()->SetTitle(titleX);
  //histo->GetXaxis()->SetLabelSize(0.05);
  mg->GetXaxis()->SetTitleSize(xTitleSize);
  mg->GetXaxis()->SetTitleOffset(xOffset);
  mg->GetYaxis()->SetTitle(titleY);
  mg->GetYaxis()->SetTitleSize(yTitleSize);
  //histo->GetYaxis()->SetLabelSize(0.05);
  mg->GetYaxis()->SetTitleOffset(yOffset);
  mg->SetTitle(title);
}

template <typename THisto>
void StyleHistoLight(THisto *histo)
{
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetLabelOffset(0.02);
  histo->GetXaxis()->SetTitleOffset(1.1);
}

void SetErrorDivide(TH1D* histFinal, const TH1D* hist1, const TH1D* hist2)
{
  const Int_t nBins = hist1->GetXaxis()->GetNbins();
  for(Int_t bin = 1; bin <= nBins; bin++) {
    Int_t k = hist1->GetBinContent(bin);
    Int_t n = hist2->GetBinContent(bin);
    Double_t errorInBin = sqrt(((Double_t)k+1)*((Double_t)k+2)/(n+2)/(n+3) - pow((Double_t)(k+1),2)/pow(n+2,2));
    histFinal->SetBinError(bin, errorInBin);
  }
}

void ErrRatioCorr(TH1F* hNum, TH1F* hDenom, TH1F* hRatio, Int_t FullCorr){
  //FullCorr == 1 means ro = 1; //full correlation                                                           
  //FullCorr == 0 means ro = sigmaDenom/sigmaNum ; //data sample at numerator is a subset of denominator     
  //FullCorr == 2 :another possibility (followed by Fiorella's according to Roger Barlow): I think this is an approximation of Barlow's prescription, and my procedure is better. I found out that errors calculated in this way vary within +-10% from those calculated by me.
  Double_t Err1=0;
  Double_t Err2=0;
  Double_t ErrC=0;
  Double_t Err=0;

  for (Int_t b=1; b<=hNum->GetNbinsX();b++){
    if (hNum->GetBinContent(b)==0 ||hDenom->GetBinContent(b)==0){
      hRatio->SetBinError(b,0);
      continue;
    }
    Err1=pow(hNum->GetBinError(b)/hNum->GetBinContent(b),2);
    Err2=pow(hDenom->GetBinError(b)/hDenom->GetBinContent(b),2);
    if (FullCorr==0){
      if (hDenom->GetBinError(b)<hNum->GetBinError(b))      ErrC=pow(hDenom->GetBinError(b),2)/(hNum->GetBinContent(b)*hDenom->GetBinContent(b)); // num is a subsample of denom
      else       ErrC=pow(hNum->GetBinError(b),2)/(hNum->GetBinContent(b)*hDenom->GetBinContent(b)); // denom is a subsample of num
    }
    else if (FullCorr==1){
      ErrC=hDenom->GetBinError(b) * hNum->GetBinError(b)/(hNum->GetBinContent(b)*hDenom->GetBinContent(b));
    }
    Err=sqrt(Err1+Err2-2*ErrC); // are we sure there's the 2?
    if (Err1+Err2-ErrC<0) {
      cout << "Error not defined! (NaN)" << endl;
      Err = sqrt(Err1+Err2);
    }
    hRatio->SetBinError(b,Err*hRatio->GetBinContent(b));
    if (FullCorr==2){
      hRatio->SetBinError(b, sqrt(TMath::Abs(pow(hNum->GetBinError(b),2) - pow(hDenom->GetBinError(b),2)) ) /hDenom->GetBinContent(b));
    }
    //  cout << "bin: " << b << " err 1: " << Err*hRatio->GetBinContent(b) << " err 2: " <<  hRatio->GetBinError(b) << " 2/1: " << hRatio->GetBinError(b) / Err*hRatio->GetBinContent(b) << endl;
  }
}

void ListHistogramsInDirectory(TDirectory *dir) {
    if (!dir) return;

    TIter next(dir->GetListOfKeys());
    TKey *key;

    while ((key = (TKey*)next())) {
        TObject *obj = key->ReadObj();
        if (obj->InheritsFrom("TH1")) {
            TH1 *histogram = (TH1*)obj;
            std::cout << "Histogram Name: " << histogram->GetName() << std::endl;
        }
    }
}

Double_t pdgMass[2] = {1.32171, 1.67245};

//Int_t ColorPt[] = {634, 628, 807, kOrange - 4, 797, 815, 418, 429, 867, 856, 601, kViolet, kPink + 9, kPink + 1, 1};
Int_t MarkerMult[] = {20, 33, 34, 29, 20, 21, 33, 34, 29, 20, 21, 33, 34, 29, 21, 21, 21, 21, 21};
//Double_t SizeMult[] = {2, 2, 2.8, 2.5, 2.8, 2, 2, 2.8, 2.5, 2.8, 2, 2, 2.8, 2.5, 2.8};
Double_t SizeMult[] = {1.5, 1.8, 1.5, 1.8, 1.5, 1.5, 1.8, 1.5, 1.8, 1.5, 1.5, 1.8, 1.5, 1.8, 1.5};

const Int_t numFitFuncs = 6;
TString nameFit[numFitFuncs] = {"m_{T}-exponential", "Boltzmann", "Fermi-Dirac", "L#acute{e}vy", "Bose-Einstein", "Blast wave"};
TString nameFitFile[numFitFuncs] = {"mTexpo", "Boltzmann", "FermiDirac", "Levi", "BoseEinstein", "BW"};

TString particleNames[] = {"XiMinus", "XiPlus", "XiPm", "OmegaMinus", "OmegaPlus", "OmegaPm"};
TString particleSymnbols[] = {"#Xi^{-}", "#bar{#Xi}^{+}", "#Xi^{-}+#bar{#Xi}^{+}", "#Omega^{-}", "#bar{#Omega}^{+}", "#Omega^{-}+#bar{#Omega}^{+}"};

TString titleInvMass[] = {"(#Lambda, #pi^{-})", "(#bar{Lambda}, #pi^{+})", "(#Lambda, #pi)", "(#Lambda, K^{-})", "(#bar{#Lambda}, K^{+})", "(#Lambda, K)"};
TString sInvMass = "invariant mass (GeV/#it{c}^{2})";
TString sPt = "#it{p}_{T} (GeV/#it{c})";
TString sdNdPt = "1/#it{N}_{ev} d#it{N}/d#it{p}_{T} [(GeV/#it{c})^{-1}]";
TString sdNdPtdY = "1/#it{N}_{ev} d^{2}#it{N}/(d#it{p}_{T}d#it{y}) [(GeV/#it{c})^{-1}]";
TString sdNdY = "#LTd#it{N}/d#it{y}#GT";
TString sdNdEta = "#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}| < 0.5}";

// const Int_t numMult = 13; // number of multiplicity sub-intervals
// Double_t multiplicityPerc[numMult + 1] = {0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 50, 70, 85, 100};
const Int_t numMult = 10; // number of multiplicity sub-intervals
Double_t multiplicityPerc[numMult + 1] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
// const Int_t numMult = 1; // number of multiplicity sub-intervals
// Double_t multiplicityPerc[numMult + 1] = {0, 100};
// Double_t multiplicityPerc[numMult + 1] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
// const Int_t numMult = 5; // number of multiplicity sub-intervals for efficiency correction
// Double_t multiplicityPerc[numMult + 1] = {0, 1, 30, 50, 70, 100};
//const Int_t numMult = 11; // number of multiplicity sub-intervals
// Double_t multiplicityPerc[numMult + 1] = {0, 0.01, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
// const Int_t numMult = 3; // number of multiplicity sub-intervals
// Double_t multiplicityPerc[numMult + 1] = {0, 30, 70, 100};

// const Int_t numMultEff = 10; // number of multiplicity sub-intervals for efficiency correction
// Double_t multiplicityPercEff[numMultEff + 1] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
// const Int_t numMultEff = 6; // number of multiplicity sub-intervals for efficiency correction
// Double_t multiplicityPercEff[numMultEff + 1] = {0, 1, 5, 15, 30, 50, 100};
const Int_t numMultEff = 5; // number of multiplicity sub-intervals for efficiency correction
Double_t multiplicityPercEff[numMultEff + 1] = {0, 1, 30, 50, 70, 100};
// const Int_t numMultEff = 1; // number of multiplicity sub-intervals for efficiency correction
// Double_t multiplicityPercEff[numMultEff + 1] = {0, 100};
// 0, 0.001, 0.01, 1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0

Int_t color[] = {kBlack, kMagenta+1, kGreen+3,kCyan+2, 634, kBlue+1, kYellow+2, kViolet, kOrange-1, kRed+1};

// Color palette
const int NRGBs = 5;
double stops[NRGBs] = { 0.00, 0.24, 0.61, 0.84, 1.00 };
// double stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
// From blue to red 
double red[NRGBs]   = { 0.00, 0.00, 0.78, 1.00, 0.51 };
double green[NRGBs] = { 0.00, 0.81, 0.90, 0.20, 0.00 };
double blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
int FI = TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,numMult);
int FIEff = TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,numMultEff);

////// BINNING //////
// Xi binning
// const Int_t numPtXi = 18;
// Double_t binptXi[numPtXi + 1] = {0.6, 0.8, 1.0,
//                                 1.2, 1.4, 1.6, 1.8, 2.0,
//                                 2.2, 2.4, 2.6, 2.8, 3.0,
//                                 3.5, 4.0, 4.5, 5.0, 6.0, 8.0};

//Xi binning to compare with Run 2
const Int_t numPtXi = 13; // to compare with Run 2
Double_t binptXi[numPtXi + 1] = {0.6, 1.0,
                                1.2, 1.4, 1.6, 1.8, 2.0,
                                2.2, 2.5, 2.9, 3.4,
                                4.0, 5.0, 6.5};

// Topo study //
// const Int_t numPtXi = 5;
// Double_t binptXi[numPtXi + 1] = {0.6, 1.6, 2.6, 3.5, 4.5, 8.0};
// const Int_t numPtXi = 1;
// Double_t binptXi[numPtXi + 1] = {0.6, 8.0};                         



int FIptxi = TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,numPtXi);       

// Omega binning
const Int_t numPtOmega = 17;
Double_t binptOmega[numPtOmega + 1] = {0.8, 1.0,
                                      1.2, 1.4, 1.6, 1.8, 2.0,
                                      2.2, 2.4, 2.6, 2.8, 3.0,
                                      3.5, 4.0, 4.5, 5.0, 6.0, 8.0};

////// SPECTRA SCALING //////
// Omega
//Double_t ScaleFactor[] = {16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1};
//Double_t ScaleFactorMB = pow(2, 19);
// Xi
Double_t ScaleFactor[] = {4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1};
Float_t etaDensity[] = {7.15, 20.62,  17.62, 15.11, 13.37, 12.02, 10.43, 8.74,  7.34, 5.64, 3.42};
Float_t etaDensityErrLow[] = {0.12, 0.43, 0.28, 0.24, 0.21, 0.19, 0.19, 0.16, 0.13, 0.10, 0.07};
Float_t etaDensityErrHigh[] = {0.09, 0.37, 0.21, 0.18, 0.16, 0.15,0.13, 0.11, 0.09, 0.07, 0.05};
Float_t fractionsINEL0run2[] = {0.0, 0.90, 4.5, 8.9, 13.5, 18.0, 27.0, 36.1, 45.3, 64.5, 100.0};
// Float_t fractionsINEL0run2[] = {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
Float_t etaDensity_mb_run2[] = {6.89};
Float_t etaDensity_errLow_mb_run2[] = {0.11};
Float_t etaDensity_errHigh_mb_run2[] = {0.11};
Float_t xiYield_mb_run2[] = {0.0273555395};
Float_t xiYield_stat_mb_run2[] = {0.0001852769};
Float_t xiYield_syst_mb_run2[] = {0.0019195267};
Double_t ScaleFactorMB = pow(2, 16);
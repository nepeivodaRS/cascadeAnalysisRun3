#include "help.h"

void qc(const TString fileMC = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisPostSQM/data/4aug-lhc24b1b-tight/AnalysisResults.root",
        const TString fileData = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisPostSQM/data/6aug-lhc22o-pass6-tight/AnalysisResults.root",
        const TString outputDir = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisPostSQM",
        const TString postFix = "8aug",
        const Bool_t isMultQC = 0)
{
  //gROOT->SetBatch(kTRUE);
  // Start of Code
  std::cout << "\x1B[1;33m"; // Set text color to yellow
  std::cout << "\n************ Starting QC ************\n";
  std::cout << "\x1B[0m"; // Reset text color

  gStyle->SetOptStat(0);

  // Get Data File //
  TFile* fileDataIn = TFile::Open(fileData);
  if (!fileDataIn || fileDataIn->IsZombie()) {
      std::cerr << "Error opening input data file!" << std::endl;
      return;
  }

  // Get MC File //
  TFile* fileMCin = TFile::Open(fileMC);
  if (!fileMCin || fileMCin->IsZombie()) {
      std::cerr << "Error opening `fileMC` data file!" << std::endl;
      return;
  }

  // Output File //
  TString outputfilePath = outputDir + "/qc/" + "qc_" + postFix + ".root";
  TFile *outputfile = new TFile(outputfilePath, "RECREATE");
  outputfile->cd();

  // Centrality calibration (DATA)//
  TH2F* hCentFT0M_rec_data = (TH2F *)fileDataIn->Get("lf-cascqaanalysis/hCentFT0M_rec");
  if (!hCentFT0M_rec_data)
  {
    std::cerr << "Histogram `hCentFT0M_rec` is not found!" << std::endl;
    return;
  }

  TH2F* hCentFT0M_rec_data_clone = (TH2F*)hCentFT0M_rec_data->Clone("hCentFT0M_rec_data_clone");
  const Int_t inel = 0; // INEL>0
  hCentFT0M_rec_data_clone->GetYaxis()->SetRange(2 + inel, 3);

  TCanvas *canvasCentrality = new TCanvas("canvasCentrality","canvasCentrality", 800, 600);
  canvasCentrality->cd();
  StyleCanvas(canvasCentrality, 0.14, 0.05, 0.11, 0.15);

  TH1F* hCentFT0M_rec_data_1D = (TH1F*)hCentFT0M_rec_data_clone->ProjectionX();
  hCentFT0M_rec_data_1D = (TH1F*)hCentFT0M_rec_data_1D->Rebin(numMult, "hCentFT0M_rec_data_1D_rebinned", multiplicityPerc);
  NormalizeHistogram(hCentFT0M_rec_data_1D);
  StyleHisto(hCentFT0M_rec_data_1D, 0.5, 1.4 * hCentFT0M_rec_data_1D->GetBinContent(hCentFT0M_rec_data_1D->GetMaximumBin()), color[0], MarkerMult[0], "FT0M Multiplicity percentile", "Counts / Bin width", "", 0, 0, 0, 1.0, 1.25, 1.1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  hCentFT0M_rec_data_1D->Draw("SAME");

  TLatex labelCentr_data;
  labelCentr_data.SetTextSize(0.05);
  labelCentr_data.SetTextAlign(12);
  labelCentr_data.DrawLatexNDC(0.56,0.78,"LHC22o_pass6");

  canvasCentrality->Write();

  // Centrality calibration (MC)//
  TH2F* hCentFT0M_rec_mc = (TH2F *)fileMCin->Get("lf-cascqaanalysis/hCentFT0M_rec");
  if (!hCentFT0M_rec_mc)
  {
    std::cerr << "Histogram `hCentFT0M_rec` is not found!" << std::endl;
    return;
  }

  TH2F* hCentFT0M_rec_mc_norm = (TH2F*)hCentFT0M_rec_mc->Clone("hCentFT0M_rec_mc_norm");
  hCentFT0M_rec_mc_norm->GetYaxis()->SetRange(2 + inel, 3);

  TCanvas *canvasCentrality_mc_norm = new TCanvas("canvasCentrality_mc_norm","canvasCentrality_mc_norm", 800, 600);
  canvasCentrality_mc_norm->cd();
  StyleCanvas(canvasCentrality_mc_norm, 0.14, 0.05, 0.11, 0.15);

  TH1F* hCentFT0M_rec_mc_norm_1D = (TH1F*)hCentFT0M_rec_mc_norm->ProjectionX();
  hCentFT0M_rec_mc_norm_1D = (TH1F*)hCentFT0M_rec_mc_norm_1D->Rebin(numMult, "hCentFT0M_rec_mc_norm_1D_rebinned", multiplicityPerc);
  NormalizeHistogram(hCentFT0M_rec_mc_norm_1D);
  StyleHisto(hCentFT0M_rec_mc_norm_1D, 0.5, 1.4 * hCentFT0M_rec_mc_norm_1D->GetBinContent(hCentFT0M_rec_mc_norm_1D->GetMaximumBin()), color[0], MarkerMult[0], "FT0M Multiplicity percentile", "Counts", "", 0, 0, 0, 1.0, 1.25, 1.1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  hCentFT0M_rec_mc_norm_1D->Draw("SAME");

  TLatex labelCentr_mc_norm;
  labelCentr_mc_norm.SetTextSize(0.05);
  labelCentr_mc_norm.SetTextAlign(12);
  labelCentr_mc_norm.DrawLatexNDC(0.56,0.78,"LHC24b1b");

  canvasCentrality_mc_norm->Write();

  // Centrality calibration (MC) not normalized//
  TH2F* hCentFT0M_rec_mc_def = (TH2F*)hCentFT0M_rec_mc->Clone("hCentFT0M_rec_mc_def");

  hCentFT0M_rec_mc_def->GetYaxis()->SetRange(2 + inel, 3);

  TCanvas *canvasCentrality_mc_def = new TCanvas("canvasCentrality_mc_def","canvasCentrality_mc_def", 800, 600);
  canvasCentrality_mc_def->cd();
  StyleCanvas(canvasCentrality_mc_def, 0.14, 0.05, 0.11, 0.15);

  TH1F* hCentFT0M_rec_mc_def_1D = (TH1F*)hCentFT0M_rec_mc_def->ProjectionX();
  hCentFT0M_rec_mc_def_1D = (TH1F*)hCentFT0M_rec_mc_def_1D->Rebin(numMult, "hCentFT0M_rec_mc_def_1D_rebinned", multiplicityPerc);
  StyleHisto(hCentFT0M_rec_mc_def_1D, 0.5, 1.4 * hCentFT0M_rec_mc_def_1D->GetBinContent(hCentFT0M_rec_mc_def_1D->GetMaximumBin()), color[0], MarkerMult[0], "FT0M Multiplicity percentile", "Counts", "", 0, 0, 0, 1.0, 1.25, 1.1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  hCentFT0M_rec_mc_def_1D->Draw("SAME");

  TLatex labelCentr_mc_def;
  labelCentr_mc_def.SetTextSize(0.05);
  labelCentr_mc_def.SetTextAlign(12);
  labelCentr_mc_def.DrawLatexNDC(0.56,0.78,"LHC24b1b");

  canvasCentrality_mc_def->Write();

  // Centrality calibration (MC) - all events //
  TH1F* hCentFT0M_mc_all = (TH1F *)fileMCin->Get("mc-centrality/FT0M/percentile");
  if (!hCentFT0M_mc_all)
  {
    std::cerr << "Histogram `hCentFT0M_mc_all` is not found!" << std::endl;
    return;
  }

  TCanvas *canvasCentrality_mc_all = new TCanvas("canvasCentrality_mc_all","canvasCentrality_mc_all", 800, 600);
  canvasCentrality_mc_all->cd();
  StyleCanvas(canvasCentrality_mc_all, 0.14, 0.05, 0.11, 0.15);

  hCentFT0M_mc_all = (TH1F*)hCentFT0M_mc_all->Rebin(numMult, "hCentFT0M_mc_all_rebinned", multiplicityPerc);
  StyleHisto(hCentFT0M_mc_all, 0.5, 1.4 * hCentFT0M_mc_all->GetBinContent(hCentFT0M_mc_all->GetMaximumBin()), color[0], MarkerMult[0], "FT0M Multiplicity percentile", "Counts", "", 0, 0, 0, 1.0, 1.25, 1.1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  hCentFT0M_mc_all->Draw("SAME");

  TLatex labelCentr_mc_all;
  labelCentr_mc_all.SetTextSize(0.05);
  labelCentr_mc_all.SetTextAlign(12);
  labelCentr_mc_all.DrawLatexNDC(0.56,0.78,"LHC24b1b");

  canvasCentrality_mc_all->Write();

  // Centrality calibration (MC) - gen //
  TH1F* hCentFT0M_mc_gen = (TH1F *)fileMCin->Get("lf-cascqaanalysis/hCentFT0M_genMC");
  if (!hCentFT0M_mc_gen)
  {
    std::cerr << "Histogram `hCentFT0M_mc_gen` is not found!" << std::endl;
    return;
  }

  TH2F* hCentFT0M_rec_mc_gen = (TH2F*)hCentFT0M_mc_gen->Clone("hCentFT0M_rec_mc_gen");
  hCentFT0M_rec_mc_gen->GetYaxis()->SetRange(2 + inel, 3);
  TH1F* hCentFT0M_rec_mc_gen_1D = (TH1F*)hCentFT0M_rec_mc_gen->ProjectionX();

  TCanvas *canvasCentrality_mc_gen = new TCanvas("canvasCentrality_mc_gen","canvasCentrality_mc_gen", 800, 600);
  canvasCentrality_mc_gen->cd();
  StyleCanvas(canvasCentrality_mc_gen, 0.14, 0.05, 0.11, 0.15);

  hCentFT0M_rec_mc_gen_1D = (TH1F*)hCentFT0M_rec_mc_gen_1D->Rebin(numMult, "hCentFT0M_rec_mc_gen_1D_rebinned", multiplicityPerc);
  StyleHisto(hCentFT0M_rec_mc_gen_1D, 0.5, 1.4 * hCentFT0M_rec_mc_gen_1D->GetBinContent(hCentFT0M_rec_mc_gen_1D->GetMaximumBin()), color[0], MarkerMult[0], "FT0M Multiplicity percentile", "Counts", "", 0, 0, 0, 1.0, 1.25, 1.1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  hCentFT0M_rec_mc_gen_1D->Draw("SAME");

  TLatex labelCentr_mc_gen;
  labelCentr_mc_gen.SetTextSize(0.05);
  labelCentr_mc_gen.SetTextAlign(12);
  labelCentr_mc_gen.DrawLatexNDC(0.56,0.78,"LHC24b1b, gen.");

  canvasCentrality_mc_gen->Write();

  // Event Selection Stats (DATA)//
  TCanvas *canvasEventSel_data = new TCanvas("canvasEventSel_data","canvasEventSel_data", 800, 600);
  canvasEventSel_data->cd();
  StyleCanvas(canvasEventSel_data, 0.14, 0.05, 0.11, 0.15);

  TH1F* hNEvents_data = (TH1F *)fileDataIn->Get("lf-cascqaanalysis/hNEvents");
  if (!hNEvents_data)
  {
    std::cerr << "Histogram `hNEvents_data` is not found!" << std::endl;
    return;
  }

  StyleHisto(hNEvents_data, 0.5, 1.4 * hNEvents_data->GetBinContent(hNEvents_data->GetMaximumBin()), color[0], MarkerMult[0], "", "Events", "", 0, 0, 0, 1.0, 1.25, 1.1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  hNEvents_data->Draw("SAME");

  TLatex labelEvSel_data;
  labelEvSel_data.SetTextSize(0.05);
  labelEvSel_data.SetTextAlign(12);
  labelEvSel_data.DrawLatexNDC(0.56,0.78,"LHC22o_pass6");

  canvasEventSel_data->Write();

  // Event Selection Stats (MC)//
  TCanvas *canvasEventSel_mc = new TCanvas("canvasEventSel_mc","canvasEventSel_mc", 800, 600);
  canvasEventSel_mc->cd();
  StyleCanvas(canvasEventSel_mc, 0.14, 0.05, 0.11, 0.15);

  TH1F* hNEvents_mc = (TH1F *)fileMCin->Get("lf-cascqaanalysis/hNEvents");
  if (!hNEvents_mc)
  {
    std::cerr << "Histogram `hNEvents_mc` is not found!" << std::endl;
    return;
  }

  StyleHisto(hNEvents_mc, 0.5, 1.4 * hNEvents_mc->GetBinContent(hNEvents_mc->GetMaximumBin()), color[0], MarkerMult[0], "", "Events", "", 0, 0, 0, 1.0, 1.25, 1.1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  hNEvents_mc->Draw("SAME");

  TLatex labelEvSel_mc;
  labelEvSel_mc.SetTextSize(0.05);
  labelEvSel_mc.SetTextAlign(12);
  labelEvSel_mc.DrawLatexNDC(0.71,0.8,"LHC24b1b");

  canvasEventSel_mc->Write();

  // hNchFT0MPVContr MC histogram //
  Int_t nNtracksLimit = 220;
  Int_t nNchLimit = 220;

  TH3F* hNchFT0MPVContr;
  TH2F* hProjPVContrNchFT0M;
  TH2F* hProjPVContrFT0Mperc;
  if (isMultQC == kTRUE) {
    hNchFT0MPVContr = (TH3F*)fileMCin->Get("lf-cascqaanalysis/hNchFT0MPVContr");
    TH3F* hNchFT0MPVContrClone = (TH3F*)hNchFT0MPVContr->Clone("hNchFT0MPVContrClone");
    hNchFT0MPVContrClone->GetZaxis()->SetRange(2 + inel, 3);

    hProjPVContrNchFT0M = static_cast<TH2F*>(hNchFT0MPVContrClone->Project3D("xy"));
    TH1D* hProjPVContrMC = (hProjPVContrNchFT0M->ProjectionX());
    Double_t meanPVcontrMC = hProjPVContrMC->GetMean();
    std::cout << "\n********* mean number of PV contributors in MC: " << meanPVcontrMC << " *********\n";
    hProjPVContrNchFT0M->SetTitle("");
    hProjPVContrNchFT0M->GetXaxis()->SetRangeUser(0., nNtracksLimit); // Set Ntracks axis range
    hProjPVContrNchFT0M->GetYaxis()->SetRangeUser(0., nNchLimit); // Set Nch in FT0M region axis range
    hProjPVContrNchFT0M->Write();

    // hFT0MpvContr DATA histogram
    TH3F* hFT0MPVContr = (TH3F*)fileDataIn->Get("lf-cascqaanalysis/hFT0MpvContr");
    TH3F* hFT0MPVContrClone = (TH3F*)hFT0MPVContr->Clone("hFT0MPVContrClone");
    hFT0MPVContrClone->GetZaxis()->SetRange(2 + inel, 3);
    hProjPVContrFT0Mperc = static_cast<TH2F*>(hFT0MPVContrClone->Project3D("xy"));
    TH1D* hProjPVContrData = (hProjPVContrFT0Mperc->ProjectionX());
    Double_t meanPVcontrDATA = hProjPVContrData->GetMean();
    std::cout << "********* mean number of PV contributors in DATA: " << meanPVcontrDATA << " *********\n";
    hProjPVContrFT0Mperc->GetXaxis()->SetRangeUser(0., nNtracksLimit); // Set Ntracks axis range
    hProjPVContrFT0Mperc->GetYaxis()->SetRangeUser(0., 105.5); // Set FT0M axis range
    hProjPVContrFT0Mperc->GetYaxis()->SetTitle("FT0M Multiplicity percentile");
    hProjPVContrFT0Mperc->SetTitle("");
    hProjPVContrFT0Mperc->Write();
  }

  // hNcandidates (check number of candidates) //
  TCanvas *canvasNcandidates = new TCanvas("canvasNcandidates","canvasNcandidates", 800, 600);
  canvasNcandidates->cd();
  StyleCanvas(canvasNcandidates, 0.14, 0.05, 0.11, 0.15);
  TH3F* hNcandidates = (TH3F*)fileDataIn->Get("lf-cascqaanalysis/hNcandidates");
  // range for qc
  Double_t dLeftCentr = 0.;
  Double_t dRightCentr = 100.;
  hNcandidates->GetYaxis()->SetRangeUser(dLeftCentr, dRightCentr);

  TH2F* hCentFT0M_rec_data_ncand = (TH2F*)hCentFT0M_rec_data->Clone("hCentFT0M_rec_data_ncand");
  hCentFT0M_rec_data_ncand->GetYaxis()->SetRange(1, 3); // take all events (including inel)
  TH1F* hCentFT0M_rec_data_1D_inel = (TH1F*)hCentFT0M_rec_data_ncand->ProjectionX();

  Double_t nEvents_data = 0;
  for (Int_t i = hCentFT0M_rec_data_1D_inel->FindBin(dLeftCentr + 1e-6); i <= hCentFT0M_rec_data_1D_inel->FindBin(dRightCentr - 1e-6); i++) {
    nEvents_data += hCentFT0M_rec_data_1D_inel->GetBinContent(i);
  }
  std::cout << "Number of events in data is: " << nEvents_data << std::endl;

  // Not pre-selected
  TH3F* hNcandidatesClone = (TH3F*)hNcandidates->Clone("hNcandidatesClone");
  hNcandidatesClone->GetZaxis()->SetRange(1, 1);
  TH2F* hNcandidatesClone_2D = static_cast<TH2F*>(hNcandidatesClone->Project3D("yx"));
  Style2DHisto(hNcandidatesClone_2D, "N_{cand.}", "FT0M %", "", 0.05, 0.05, 0.05, 0.05, 0.05);
  TH1D* hNcandidatesClone_1D = (hNcandidatesClone_2D->ProjectionX());
  canvasNcandidates->cd()->SetLogy();
  //hNcandidatesClone_1D->Scale(1. / nEvents_data);
  //StyleHisto(hNcandidatesClone_1D, 1e-1, 1e1 * hNcandidatesClone_1D->GetBinContent(hNcandidatesClone_1D->GetMaximumBin()), color[0], MarkerMult[0], "N_{cand.}", "Frequency", "", 1, 0, 10, 1.0, 1.25, 1.1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  StyleHisto(hNcandidatesClone_1D, 1e-1, 1e1 * hNcandidatesClone_1D->GetBinContent(hNcandidatesClone_1D->GetMaximumBin()), color[0], MarkerMult[0], "N_{cand.}", "Events", "", 1, 0, 10, 1.0, 1.25, 1.1, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  hNcandidatesClone_1D->Draw("P same");

  // Fit with Poisson
  TF1 *poissonFit = new TF1("poissonFit", "[0]*TMath::Poisson(x, [1])", 0.5, 3.5);
  poissonFit->SetParameters(1e8, 0.01); // norm, mean
  hNcandidatesClone_1D->Fit("poissonFit", "RL");
  poissonFit->SetNpx(1e4);
  poissonFit->SetLineWidth(2);
  poissonFit->SetLineColor(kRed+1); //{kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
  poissonFit->Draw("same");
  // Fit with Exponential
  // TF1 *expFit = new TF1("expFit", "[0]*exp([1]*x) + [2]*log(x)", 2., 6.);
  // expFit->SetParameters(1e8, -1, 0.005); // norm, power
  // hNcandidatesClone_1D->Fit("expFit", "RL");
  // expFit->SetNpx(1e4);
  // expFit->SetLineWidth(2);
  // expFit->SetLineColor(kRed+1); //{kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
  // expFit->Draw("same");

  // Pre-selected 
  TH3F* hNcandidatesClone_presel = (TH3F*)hNcandidates->Clone("hNcandidatesClone_presel");
  hNcandidatesClone_presel->GetZaxis()->SetRange(2, 2);
  TH2F* hNcandidatesClone_presel_2D = static_cast<TH2F*>(hNcandidatesClone_presel->Project3D("yx"));
  Style2DHisto(hNcandidatesClone_presel_2D, "N_{cand.}", "FT0M %", "", 0.05, 0.05, 0.05, 0.05, 0.05);
  TH1D* hNcandidatesClone_presel_1D = (hNcandidatesClone_presel_2D->ProjectionX());
  //hNcandidatesClone_presel_1D->Scale(1. / nEvents_data);
  StyleHisto(hNcandidatesClone_presel_1D, 1e-10, 1e1 * hNcandidatesClone_presel_1D->GetBinContent(hNcandidatesClone_presel_1D->GetMaximumBin()), color[1], MarkerMult[1], "N_{cand.}", "Frequency", "", 1, 0, 10, 1.0, 1.25, 1.4, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  hNcandidatesClone_presel_1D->Draw("P same");

  // MC not pre-selected
  TH3F* hNcandidates_mc = (TH3F*)fileMCin->Get("lf-cascqaanalysis/hNcandidates");
  hNcandidates_mc->GetYaxis()->SetRangeUser(dLeftCentr, dRightCentr);

  TH2F* hCentFT0M_rec_mc_ncand = (TH2F*)hCentFT0M_rec_mc->Clone("hCentFT0M_rec_mc_clone");
  hCentFT0M_rec_mc_ncand->GetYaxis()->SetRange(1, 3); // take all events (including inel)
  TH1F* hCentFT0M_rec_mc_1D_inel = (TH1F*)hCentFT0M_rec_mc_ncand->ProjectionX();

  Double_t nEvents_mc = 0;
  for (Int_t i = hCentFT0M_rec_mc_1D_inel->FindBin(dLeftCentr + 1e-6); i <= hCentFT0M_rec_mc_1D_inel->FindBin(dRightCentr - 1e-6); i++) {
    nEvents_mc += hCentFT0M_rec_mc_1D_inel->GetBinContent(i);
  }

  std::cout << "Number of events in mc is: " << nEvents_mc << std::endl;

  TH3F* hNcandidatesClone_mc = (TH3F*)hNcandidates_mc->Clone("hNcandidatesClone_mc");
  hNcandidatesClone_mc->GetZaxis()->SetRange(1, 1);
  TH2F* hNcandidatesClone_mc_2D = static_cast<TH2F*>(hNcandidatesClone_mc->Project3D("yx"));
  Style2DHisto(hNcandidatesClone_mc_2D, "N_{cand.}", "FT0M %", "", 0.05, 0.05, 0.05, 0.05, 0.05);
  TH1D* hNcandidatesClone_mc_1D = (hNcandidatesClone_mc_2D->ProjectionX());
  hNcandidatesClone_mc_1D->Scale(1. / nEvents_mc);
  StyleHisto(hNcandidatesClone_mc_1D, 1e-10, 1e1 * hNcandidatesClone_mc_1D->GetBinContent(hNcandidatesClone_mc_1D->GetMaximumBin()), color[2], MarkerMult[2], "N_{cand.}", "Frequency", "", 1, 0, 10, 1.0, 1.25, 1.4, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  //hNcandidatesClone_mc_1D->Draw("P same");

  // MC pre-selected
  TH3F* hNcandidatesClone_mc_presel = (TH3F*)hNcandidates_mc->Clone("hNcandidatesClone_mc_presel");
  hNcandidatesClone_mc_presel->GetZaxis()->SetRange(2, 2);
  TH2F* hNcandidatesClone_mc_presel_2D = static_cast<TH2F*>(hNcandidatesClone_mc_presel->Project3D("yx"));
  Style2DHisto(hNcandidatesClone_mc_presel_2D, "N_{cand.}", "FT0M %", "", 0.05, 0.05, 0.05, 0.05, 0.05);
  TH1D* hNcandidatesClone_mc_presel_1D = (hNcandidatesClone_mc_presel_2D->ProjectionX());
  hNcandidatesClone_mc_presel_1D->Scale(1. / nEvents_mc);
  StyleHisto(hNcandidatesClone_mc_presel_1D, 1e-10, 1e1 * hNcandidatesClone_mc_presel_1D->GetBinContent(hNcandidatesClone_mc_presel_1D->GetMaximumBin()), color[3], MarkerMult[3], "N_{cand.}", "Frequency", "", 1, 0, 10, 1.0, 1.25, 1.4, 0.05, 0.05, 0.05, 0.05, 0.007, 0.007);
  //hNcandidatesClone_mc_presel_1D->Draw("P same");

  TLegend *legend_Ncand = new TLegend(0.354, 0.63, 0.875, 0.875);
  legend_Ncand->SetHeader("");
  legend_Ncand->SetFillStyle(0);
  legend_Ncand->AddEntry(hNcandidatesClone_1D, "LHC22o_pass6", "p");
  legend_Ncand->AddEntry(hNcandidatesClone_presel_1D, "LHC22o_pass6, pre-selected", "p");
  // legend_Ncand->AddEntry(hNcandidatesClone_mc_1D, "LHC24b1b", "p");
  // legend_Ncand->AddEntry(hNcandidatesClone_mc_presel_1D, "LHC24b1b, pre-selected", "p");
  StyleLegend(legend_Ncand, 0.0, 0.0);
  legend_Ncand->SetTextSize(0.04);
  legend_Ncand->Draw("same");

  TLatex labelEvSel;
  labelEvSel.SetTextSize(0.04);
  labelEvSel.SetTextAlign(12);
  labelEvSel.DrawLatexNDC(0.18, 0.23,"Tight event sel. (inel)");

  canvasNcandidates->Write();

  if (isMultQC == kTRUE) {
    // FV0A vs FT0M (DATA) //
    TCanvas *canvasFV0AFT0M_data = new TCanvas("canvasFV0AFT0M_data","canvasFV0AFT0M_data", 800, 600);
    canvasFV0AFT0M_data->cd();
    StyleCanvas(canvasFV0AFT0M_data, 0.14, 0.14, 0.11, 0.15);
    TH3F* hFV0AFT0M = (TH3F*)fileDataIn->Get("lf-cascqaanalysis/hFV0AFT0M");
    hFV0AFT0M->GetZaxis()->SetRange(2 + inel, 3);
    TH2F* hFV0AFT0M_2D = static_cast<TH2F*>(hFV0AFT0M->Project3D("xy"));
    hFV0AFT0M_2D->GetXaxis()->SetRangeUser(0., 100.);
    hFV0AFT0M_2D->GetYaxis()->SetRangeUser(0., 100.);
    Style2DHisto(hFV0AFT0M_2D, "FV0A %", "FT0M %", "", 0.05, 0.05, 0.05, 0.05, 0.05);
    canvasFV0AFT0M_data->cd()->SetLogz();
    hFV0AFT0M_2D->Draw("COLZ");

    // Compare hFT0Msignal vs PVcontr in DATA and MC
    TH3F* hFT0MsignalPVContrData = (TH3F*)fileDataIn->Get("lf-cascqaanalysis/hFT0MsignalPVContr");
    TH3F* hFT0MsignalPVContrDataClone = (TH3F*)hFT0MsignalPVContrData->Clone("hFT0MsignalPVContrDataClone");
    hFT0MsignalPVContrDataClone->GetZaxis()->SetRange(2 + inel, 3);
    TH2F* hProjPVContrFT0MsignalData = static_cast<TH2F*>(hFT0MsignalPVContrDataClone->Project3D("xy"));
    hProjPVContrFT0MsignalData->GetXaxis()->SetRangeUser(0., nNtracksLimit); // Set Ntracks axis range

    TH3F* hFT0MsignalPVContrMC = (TH3F*)fileMCin->Get("lf-cascqaanalysis/hFT0MsignalPVContr");
    TH3F* hFT0MsignalPVContrMCClone = (TH3F*)hFT0MsignalPVContrMC->Clone("hFT0MsignalPVContrMCClone");
    hFT0MsignalPVContrMCClone->GetZaxis()->SetRange(2 + inel, 3);
    TH2F* hProjPVContrFT0MsignalMC = static_cast<TH2F*>(hFT0MsignalPVContrMCClone->Project3D("xy"));
    hProjPVContrFT0MsignalMC->GetXaxis()->SetRangeUser(0., nNtracksLimit); // Set Ntracks axis range

    TCanvas* canvas2D_PVcontrFT0Msignal = new TCanvas("canvas2D_PVcontrFT0Msignal", "canvas2D_PVcontrFT0Msignal", 1000,800);
    canvas2D_PVcontrFT0Msignal->Divide(2, 2);
    StyleCanvas(canvas2D_PVcontrFT0Msignal, 0.15, 0.05, 0.05, 0.15);
    for (Int_t i = 1; i <= 4; i++) {
      canvas2D_PVcontrFT0Msignal->cd(i)->SetMargin(0.16, 0.15, 0.11, 0.1); 
    }
    canvas2D_PVcontrFT0Msignal->cd(3)->SetLogz();
    canvas2D_PVcontrFT0Msignal->cd(4)->SetLogz();

    TLegend *LegendFT0MsignalData = new TLegend(0.40, 0.64, 0.83, 0.84);
    LegendFT0MsignalData->SetFillStyle(0);
    LegendFT0MsignalData->SetTextAlign(33);
    LegendFT0MsignalData->SetTextSize(0.04);
    LegendFT0MsignalData->SetTextFont(42);
    LegendFT0MsignalData->SetLineColorAlpha(0.,0.);
    LegendFT0MsignalData->SetFillColorAlpha(0.,0.);
    LegendFT0MsignalData->SetBorderSize(0.);
    LegendFT0MsignalData->AddEntry("", "#bf{ALICE Work In Progress}", "");
    LegendFT0MsignalData->AddEntry("", "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    LegendFT0MsignalData->AddEntry("", "LHC22o_pass6, INEL > 0", "");

    canvas2D_PVcontrFT0Msignal->cd(1);
    Style2DHisto(hProjPVContrFT0MsignalData, "#it{N}_{tracks}", "FT0M signal", "DATA", 0.05, 0.05, 0.05, 0.05, 0.05);
    hProjPVContrFT0MsignalData->GetYaxis()->SetTitleOffset(1.7);
    hProjPVContrFT0MsignalData->Scale(1.0 / hProjPVContrFT0MsignalData->GetEntries());
    hProjPVContrFT0MsignalData->SetMinimum(0.0);  // Set the minimum Z value (e.g., 0)
    hProjPVContrFT0MsignalData->SetMaximum(0.65e-3);  // Set the maximum Z value (e.g., 1 for a normalized histogram)
    hProjPVContrFT0MsignalData->Draw("COLZ");
    LegendFT0MsignalData->Draw("same");

    TLegend *LegendFT0MsignalMC = new TLegend(0.40, 0.64, 0.83, 0.84);
    LegendFT0MsignalMC->SetFillStyle(0);
    LegendFT0MsignalMC->SetTextAlign(33);
    LegendFT0MsignalMC->SetTextSize(0.04);
    LegendFT0MsignalMC->SetTextFont(42);
    LegendFT0MsignalMC->SetLineColorAlpha(0.,0.);
    LegendFT0MsignalMC->SetFillColorAlpha(0.,0.);
    LegendFT0MsignalMC->SetBorderSize(0.);
    LegendFT0MsignalMC->AddEntry("", "LHC24b1b, rec. true INEL > 0", "");

    canvas2D_PVcontrFT0Msignal->cd(2);
    Style2DHisto(hProjPVContrFT0MsignalMC, "#it{N}_{tracks}", "FT0M signal", "MC", 0.05, 0.05, 0.05, 0.05, 0.05);
    hProjPVContrFT0MsignalMC->GetYaxis()->SetTitleOffset(1.7);
    hProjPVContrFT0MsignalMC->Scale(1.0 / hProjPVContrFT0MsignalMC->GetEntries());
    hProjPVContrFT0MsignalMC->SetMinimum(0.0);  // Set the minimum Z value (e.g., 0)
    hProjPVContrFT0MsignalMC->SetMaximum(0.65e-3);  // Set the maximum Z value (e.g., 1 for a normalized histogram)
    hProjPVContrFT0MsignalMC->Draw("COLZ");
    LegendFT0MsignalMC->Draw();

    canvas2D_PVcontrFT0Msignal->cd(3);
    TH2F* hProjPVContrFT0MsignalDataZoomed = (TH2F*)hProjPVContrFT0MsignalData->Clone("hProjPVContrFT0MsignalDataZoomed");
    Int_t nFT0MsignalZoomedLimit = 4500;
    Int_t nNtracksZoomedLimit = 70;
    hProjPVContrFT0MsignalDataZoomed->GetYaxis()->SetRangeUser(0., nFT0MsignalZoomedLimit); // Set FT0M raw signal axis range
    hProjPVContrFT0MsignalDataZoomed->GetXaxis()->SetRangeUser(0., nNtracksZoomedLimit); // Set Ntracks axis range
    hProjPVContrFT0MsignalDataZoomed->SetTitle("Zoomed DATA");
    hProjPVContrFT0MsignalDataZoomed->Draw("COLZ");

    canvas2D_PVcontrFT0Msignal->cd(4);
    TH2F* hProjPVContrFT0MsignalMCZoomed = (TH2F*)hProjPVContrFT0MsignalMC->Clone("hProjPVContrFT0MsignalMCZoomed");
    hProjPVContrFT0MsignalMCZoomed->GetYaxis()->SetRangeUser(0., nFT0MsignalZoomedLimit); // Set FT0M raw signal axis range
    hProjPVContrFT0MsignalMCZoomed->GetXaxis()->SetRangeUser(0., nNtracksZoomedLimit); // Set Ntracks axis range
    hProjPVContrFT0MsignalMCZoomed->SetTitle("Zoomed MC");
    hProjPVContrFT0MsignalMCZoomed->Draw("COLZ");

    canvas2D_PVcontrFT0Msignal->Write();

    // Calc. the fraction of BG events in DATA
    Double_t nWeirdEventsData = 0;
    Double_t nAllEventsData = 0;
    for (Int_t binx = 1; binx <= hProjPVContrFT0MsignalData->GetNbinsX(); binx++) {
      for (Int_t biny = 1; biny <= hProjPVContrFT0MsignalData->GetNbinsY(); biny++) {
        Double_t binContent = hProjPVContrFT0MsignalData->GetBinContent(binx, biny);
        nAllEventsData += binContent;
        Double_t binValueY = hProjPVContrFT0MsignalData->GetYaxis()->GetBinCenter(biny);
        if(binValueY > 5000){
          nWeirdEventsData += binContent;
        }
      }
    }
    std::cout << "********* nWeirdEventsData (FT0M signal > 5k): " << nWeirdEventsData << " nAllEventsData: " <<  nAllEventsData << " Fraction of weird events in DATA: " << nWeirdEventsData/nAllEventsData << " *********\n";

    // Just for QC, check how mean of PV contributors depends on Ngen in FT0M //
    TH2F* hProjPVContrNchFT0MClone = (TH2F*)hProjPVContrNchFT0M->Clone("hProjPVContrNchFT0MClone");
    TProfile *hProfileNFT0M = hProjPVContrNchFT0MClone->ProfileY("ProfileNgen", 1, hProjPVContrNchFT0MClone->GetXaxis()->GetNbins());
    hProfileNFT0M->SetName("TProfile_plot_for_NgenFT0M");
    hProfileNFT0M->GetYaxis()->SetTitle("<#it{N}_{tracks}>");
    hProfileNFT0M->Write();
    // Just for QC, check number of counts for each Ngen in FT0M bin //
    TH2F* hProjPVContrNchFT0MClone2 = (TH2F*)hProjPVContrNchFT0M->Clone("hProjPVContrNchFT0MClone2");
    TH1D* hProjNFT0M = (hProjPVContrNchFT0MClone2->ProjectionY());
    hProjNFT0M->GetYaxis()->SetTitle("Counts");
    hProjNFT0M->SetTitle("");
    hProjNFT0M->Write();

    std::cout << "\n****** Try to estimate the best mult. % division for MC to have a flat distribution (just for QC) *********\n";
    // Just for QC, try to evaluate best mult. % division for MC to have a flat distribution //
    std::cout << "****** Total number of events in GP MC: " << hProjNFT0M->Integral() << " *********\n";
    Double_t evIn1Perc = hProjNFT0M->Integral()/100.0;
    Double_t classWidth = 0;
    long int numEventsFLAT = 0;
    Int_t classCount = 0;
    do {
      classCount = 0;
      classWidth += 1;
      std::cout << "****** Iteration number " << classWidth << " *********\n";
      std::cout << "* Should be in " << classWidth << " % wide mult. class (GP MC): " << classWidth * evIn1Perc << " *\n";
      for (Int_t bin = 1; bin <= hProjNFT0M->GetNbinsX(); ++bin) {
        numEventsFLAT += hProjNFT0M->GetBinContent(bin);
        if(numEventsFLAT >= classWidth*evIn1Perc - evIn1Perc/3.){
          std::cout << "NumEvents: " << numEventsFLAT << " Ngen: " << hProjNFT0M->GetXaxis()->GetBinUpEdge(bin) << " Class: " << 100. - (classCount+1)*classWidth << " -- " << 100. - classCount*classWidth << std::endl;
          classCount++;
          numEventsFLAT = 0;
        }
      }

    } while (classCount != 100/int(classWidth));
    std::cout << "\n****** Class width to make calibration flat using GP MC is: " << classWidth << " %" <<  " *********\n";
  }
  if (isMultQC == kTRUE) {
    // Compare f(Ntracks) distributions in data and GP MC //
    TCanvas* canvasPVContrComp = new TCanvas("PVContrComp", "PVContrComp", 0, 70, 620, 850);
    StyleCanvas(canvasPVContrComp, 0.15, 0.05, 0.05, 0.15);
    TPad *padPVContrCompUp = new TPad("padPVContrCompUp", "padPVContrCompUp", 0, 0.36, 1, 1);
    TPad *padPVContrCompLow = new TPad("padPVContrCompLow", "padPVContrCompLow", 0, 0, 1, 0.36);
    StylePad(padPVContrCompUp, 0.15, 0.05, 0.05, 0.);
    StylePad(padPVContrCompLow, 0.15, 0.05, 0.02, 0.2);
    canvasPVContrComp->cd();
    padPVContrCompUp->Draw();
    padPVContrCompLow->Draw();

    // Dummy histograms
    TH1F *hDummyPVContrComp = new TH1F("hDummyPVContrComp", "hDummyPVContrComp", 10000, 0, nNtracksLimit);
    Int_t nNtracksLimitOffset = 100;
    // Up
    for (Int_t i = 1; i <= hDummyPVContrComp->GetNbinsX(); i++)
      hDummyPVContrComp->SetBinContent(i, 1e-12);
    canvasPVContrComp->cd();
    padPVContrCompUp->cd();
    StyleHisto(hDummyPVContrComp, 1e-6, 2, 1, 1, "#it{N}_{tracks}", "f(#it{N}_{tracks})", "", 0, 0, 0, 1.5, 1.0, 0, 0.0, 0.05, 0.0, 0.035, 0.005, 0.005);
    SetTickLength(hDummyPVContrComp, 0.025, 0.03);
    TAxis *axisPVContrCompDummy = hDummyPVContrComp->GetYaxis();
    axisPVContrCompDummy->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
    hDummyPVContrComp->GetXaxis()->SetRangeUser(0, nNtracksLimit - nNtracksLimitOffset);
    padPVContrCompUp->SetLogy();
    hDummyPVContrComp->Draw("same");
    // Low
    TH1F *hDummyPVContrLow = new TH1F("hDummyPVContrLow", "hDummyPVContrLow", 10000, 0, nNtracksLimit);
    for (Int_t i = 1; i <= hDummyPVContrLow->GetNbinsX(); i++)
      hDummyPVContrLow->SetBinContent(i, 1);
    padPVContrCompLow->cd();
    StyleHisto(hDummyPVContrLow, 0.09, 11, 1, 1, "#it{N}_{tracks}", "Data/MC", "", 0, 0, 0, 1.0, 0.7, 0, 0.08, 0.08, 0.08, 0.07, 0.01, 0.01);
    SetTickLength(hDummyPVContrLow, 0.025, 0.03);
    TAxis *axisDummyPVContrLow = hDummyPVContrLow->GetYaxis();
    axisDummyPVContrLow->ChangeLabel(-1, -1, -1, -1, -1, -1, " "); // last 
    //axisDummyPVContrLow->ChangeLabel(1, -1, -1, -1, -1, -1, " "); // first
    hDummyPVContrLow->GetXaxis()->SetRangeUser(0, nNtracksLimit - nNtracksLimitOffset);
    hDummyPVContrLow->GetYaxis()->CenterTitle();
    hDummyPVContrLow->Draw("same");

    // Plotting
    padPVContrCompUp->cd();
    Int_t rebinFactor = 4;
    TH2F* hProjPVContrNchFT0MClone3 = (TH2F*)hProjPVContrNchFT0M->Clone("NtracksProjMC");
    TH1D* hProjPVContrMC2 = hProjPVContrNchFT0MClone3->ProjectionX();
    hProjPVContrMC2->Rebin(rebinFactor);
    Int_t markTypeMC = 3;
    hProjPVContrMC2->SetTitle("");
    hProjPVContrMC2->SetMarkerColor(color[markTypeMC]);
    hProjPVContrMC2->SetLineColor(color[markTypeMC]);
    hProjPVContrMC2->SetMarkerStyle(MarkerMult[markTypeMC]);
    hProjPVContrMC2->SetMarkerSize(SizeMult[markTypeMC]);
    hProjPVContrMC2->Scale(1.0/hProjPVContrMC2->Integral());
    hProjPVContrMC2->GetXaxis()->SetRangeUser(0, nNtracksLimit); // Ntracks axis range
    hProjPVContrMC2->GetYaxis()->SetRangeUser(1e-8, hProjPVContrMC2->GetMaximum() * 10); // f(Ntracks)
    hProjPVContrMC2->GetYaxis()->SetTitle("f(N_{tracks})");
    hProjPVContrMC2->Draw("same");

    TH2F* hProjPVContrFT0MpercClone = (TH2F*)hProjPVContrFT0Mperc->Clone("NtracksProjData");
    TH1D* hProjPVContrData2 = hProjPVContrFT0MpercClone->ProjectionX();
    hProjPVContrData2->Rebin(rebinFactor);
    Int_t markTypeData = 0;
    hProjPVContrData2->SetMarkerColor(color[markTypeData]);
    hProjPVContrData2->SetLineColor(color[markTypeData]);
    hProjPVContrData2->SetMarkerStyle(MarkerMult[markTypeData]);
    hProjPVContrData2->SetMarkerSize(SizeMult[4]);
    hProjPVContrData2->Scale(1.0/hProjPVContrData2->Integral());
    hProjPVContrData2->GetXaxis()->SetRangeUser(0, nNtracksLimit); // Set Ntracks axis range
    hProjPVContrData2->Draw("same");

    TLegend *legPVContrComp = new TLegend(0.42, 0.63, 0.94, 0.88);
    legPVContrComp->SetHeader("Comparison of N_{tracks} distributions");
    legPVContrComp->SetFillStyle(0);
    legPVContrComp->AddEntry(hProjPVContrData2, "Data", "p");
    legPVContrComp->AddEntry(hProjPVContrMC2, "MC", "p");
    StyleLegend(legPVContrComp, 0.0, 0.0);
    legPVContrComp->SetTextSize(0.04);
    legPVContrComp->Draw("same");

    padPVContrCompLow->cd();
    TH1D* hPVContrRatio = (TH1D*)hProjPVContrData2->Clone("PVContrData");
    hPVContrRatio->Divide(hProjPVContrMC2);
    StyleHisto(hPVContrRatio, 0.1 , 11, color[markTypeData], MarkerMult[markTypeData], "#it{N}_{tracks}", "Data/MC", "", 0, 0, 0, 1.0, 0.7, 0.7, 0.08, 0.08, 0.08, 0.07, 0.01, 0.01);
    hPVContrRatio->SetStats(0);
    hPVContrRatio->Draw("same");

    canvasPVContrComp->Write();
  }
  // End of Code
  std::cout << "\x1B[1;32m"; // Set text color to green
  std::cout << "\n************* QC is Done Successfully! *************\n";
  std::cout << "\x1B[0m"; // Reset text color
}
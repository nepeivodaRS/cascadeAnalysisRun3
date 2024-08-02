#include "../../help.h"
#include "../../effHelp.h"

void topoStudy1D(const Int_t nParticle = 2, // 0-2 : xi, 3-5 : omega
                         TString fileListDATA = "casccospaDATA.txt",
                         TString fileListMC = "casccospaMC.txt",
                         const TString xLabelIn = "cos(PA)",
                         Double_t cutMin = 0.970,
                         Double_t cutMax = 0.999,
                         Double_t yMin = 0.1,
                         Double_t yMax = 1.5,
                         Double_t yRatioMin = 0.96,
                         Double_t yRatioMax = 1.04,
                         const bool isGeq = true, 
                         const std::string workingDir = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/topoStudy/casccospa",
                         const std::string postFix = "_casccospa")
{
  // Start of Code
  std::cout << "\x1B[1;33m"; // Set text color to yellow
  std::cout << "\n************ Starting Topo/Kinem Study ************\n";
  std::cout << "\x1B[0m"; // Reset text color

  // output file
  TString outputfilePath;
  outputfilePath = workingDir + "/results/" + "topoStudy" + postFix + ".root";
  TFile *outputfile = new TFile(outputfilePath, "RECREATE");

  // Files (DATA) //
  std::vector<std::string> nameDATA;
  std::vector<std::string> nameLegendDATA;
  std::ifstream fileDATA(Form("%s", fileListDATA.Data()));

  std::string removeDATA = workingDir + "/data/yieldsOut/yield_XiPm_MB_inel0";
  std::string remove2 = ".root";
  std::string remove3 = "_";
  cout << "Files:" << endl;

  if (fileDATA.is_open())
  {
    std::string line;
    while (std::getline(fileDATA, line))
    {
      nameDATA.push_back(line);
      cout << line << endl;
      size_t pos = line.find(removeDATA);
      if (pos != std::string::npos)
      {
        line.erase(pos, removeDATA.length());
      }

      size_t pos2 = line.find(remove2);
      if (pos2 != std::string::npos)
      {
        line.erase(pos2, remove2.length());
      }
      //nameLegendDATA.push_back(line);
      size_t pos3 = line.find(remove3);
      if (pos3 != std::string::npos)
      {
        line.erase(pos3, remove3.length());
      }
      nameLegendDATA.push_back(line);
    }
    fileDATA.close();
  }
  else
  {
    std::cerr << "Unable to open fileListDATA!" << std::endl;
  }

  cout << "List:" << endl;
  cout << fileListDATA.Data() << endl;

  const Int_t numFilesDATA = nameDATA.size();
  cout << "Number of files: " << numFilesDATA << endl;

  // Files (MC) //
  std::vector<std::string> nameMC;
  std::vector<std::string> nameLegendMC;
  std::ifstream fileMC(Form("%s", fileListMC.Data()));

  std::string removeMC = workingDir + "/mc/yieldsOut/yield_XiPm_MB_inel0";
  cout << "Files:" << endl;

  if (fileMC.is_open())
  {
    std::string line;
    while (std::getline(fileMC, line))
    {
      nameMC.push_back(line);
      cout << line << endl;
      size_t pos = line.find(removeMC);
      if (pos != std::string::npos)
      {
        line.erase(pos, removeMC.length());
      }

      size_t pos2 = line.find(remove2);
      if (pos2 != std::string::npos)
      {
        line.erase(pos2, remove2.length());
      }
      //nameLegendMC.push_back(line);
      size_t pos3 = line.find(remove3);
      if (pos3 != std::string::npos)
      {
        line.erase(pos3, remove3.length());
      }
      nameLegendMC.push_back(line);
    }
    fileMC.close();
  }
  else
  {
    std::cerr << "Unable to open fileListMC!" << std::endl;
  }

  cout << "List:" << endl;
  cout << fileListMC.Data() << endl;

  const Int_t numFilesMC = nameMC.size();
  cout << "Number of files: " << numFilesMC << endl;

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

  TFile *fileInMC[numFilesMC];
  TFile *fileInDATA[numFilesDATA];

  // DATA //
  std::vector<std::tuple<Double_t, std::string, Int_t>> yieldDATA; // value, config #, bin number
  std::vector<std::tuple<Double_t, std::string, Int_t>> varsDATA; // value, config #, bin number
  for (Int_t iFile = 0; iFile < numFilesDATA; iFile++)
  {
    // Open input file
    fileInDATA[iFile] = TFile::Open(Form("%s", nameDATA[iFile].c_str()));
    if (!fileInDATA[iFile] || fileInDATA[iFile]->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
        return;
    }

    fileInDATA[iFile]->cd();

    TDirectory* dirFitParams = fileInDATA[iFile]->GetDirectory("fitParams");
    if (!dirFitParams)
    {
      std::cerr << "`fitParams` directory is not found!" << std::endl;
      return;
    }

    dirFitParams->cd();

    TH1F* hYield = (TH1F *)dirFitParams->Get("Yield_" + particleNames[nParticle]);
    if (!hYield)
    {
      std::cerr << "Histogram `hYield` is not found!" << std::endl;
      return;
    }

    for (Int_t b = 1; b <= numPtBins; b++) {
      Double_t yieldValue = hYield->GetBinContent(b);
      yieldDATA.push_back(std::make_tuple(yieldValue, nameLegendDATA[iFile], b));
    }
  }

  // Find the loosest config in DATA
  std::string sLoosestNumber = isGeq ? "1" : std::to_string(numFilesDATA - 1);
  std::vector<std::pair<Double_t, Int_t>> yieldLoosestDATA; // yield, bin number
  for (const auto& tuple : yieldDATA) {
      if (std::get<1>(tuple) == sLoosestNumber) {
        yieldLoosestDATA.push_back(std::make_pair(std::get<0>(tuple), std::get<2>(tuple)));
      }
  }

  // Make signal loss fraction
  for (const auto& tuple : yieldDATA) {
    for (const auto& pair : yieldLoosestDATA) {
      if (pair.second == std::get<2>(tuple)) { // match bin number
        Double_t recYield = std::get<0>(tuple);
        Double_t yieldVariation = recYield / pair.first;
        varsDATA.push_back(std::make_tuple(yieldVariation, std::get<1>(tuple), std::get<2>(tuple)));
      }
    }
  }

  // MC //
  std::vector<std::tuple<Double_t, std::string, Int_t>> yieldMC; // value, config #, bin number
  std::vector<std::tuple<Double_t, std::string, Int_t>> varsMC; // value, config #, bin number
  for (Int_t iFile = 0; iFile < numFilesMC; iFile++)
  {
    // Open input file
    fileInMC[iFile] = TFile::Open(Form("%s", nameMC[iFile].c_str()));
    if (!fileInMC[iFile] || fileInMC[iFile]->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
        return;
    }

    fileInMC[iFile]->cd();

    TDirectory* dirFitParams = fileInMC[iFile]->GetDirectory("fitParams");
    if (!dirFitParams)
    {
      std::cerr << "`fitParams` directory is not found!" << std::endl;
      return;
    }

    dirFitParams->cd();

    TH1F* hYield = (TH1F *)dirFitParams->Get("hGenCountsClone"); // take true gen. number of rec. cascades
    //TH1F* hYield = (TH1F *)dirFitParams->Get("Yield_" + particleNames[nParticle]); // take rec. number of cascades from signal extraction procedure
    if (!hYield)
    {
      std::cerr << "Histogram `hYield` is not found!" << std::endl;
      return;
    }

    for (Int_t b = 1; b <= numPtBins; b++) {
      Double_t yieldValue = hYield->GetBinContent(b);
      yieldMC.push_back(std::make_tuple(yieldValue, nameLegendMC[iFile], b));
    }
  }

  // Find the loosest config in MC
  std::vector<std::pair<Double_t, Int_t>> yieldLoosestMC; // yield, bin number
  for (const auto& tuple : yieldMC) {
      if (std::get<1>(tuple) == sLoosestNumber) {
        yieldLoosestMC.push_back(std::make_pair(std::get<0>(tuple), std::get<2>(tuple)));
      }
  }

  // Make signal loss fraction
  for (const auto& tuple : yieldMC) {
    for (const auto& pair : yieldLoosestMC) {
      if (pair.second == std::get<2>(tuple)) { // match bin number
        Double_t recYield = std::get<0>(tuple);
        Double_t yieldVariation = recYield / pair.first;
        varsMC.push_back(std::make_tuple(yieldVariation, std::get<1>(tuple), std::get<2>(tuple)));
      }
    }
  }

  TH1F* signaLossDATA[numPtBins];
  TH1F* signaLossMC[numPtBins];
  TH1F* signaLossRatio[numPtBins];
  for (Int_t b = 1; b <= numPtBins; b++) {
    signaLossDATA[b - 1] = new TH1F(Form("signalLossFractionDATA_%i", b), Form("signalLossFractionDATA_%i", b), numFilesDATA - 1, cutMin - (cutMax - cutMin) / (numFilesDATA - 2) / 2, cutMax + (cutMax - cutMin) / (numFilesDATA - 2) / 2);
    signaLossMC[b - 1] = new TH1F(Form("signalLossFractionMC_%i", b), Form("signalLossFractionMC_%i", b), numFilesMC - 1, cutMin - (cutMax - cutMin) / (numFilesDATA - 2) / 2, cutMax + (cutMax - cutMin) / (numFilesDATA - 2) / 2);
  }

  // DATA
  for (Int_t b = 1; b <= numPtBins; b++) {
    for (auto& tuple : varsDATA) {
      if (b == std::get<2>(tuple)) { // match bin number
        if (std::get<1>(tuple) != "DEF") {
          signaLossDATA[b - 1]->SetBinContent(std::stoi(std::get<1>(tuple)), std::get<0>(tuple));
        }
      }
    }
  }
  // MC
  for (Int_t b = 1; b <= numPtBins; b++) {
    for (auto& tuple : varsMC) {
      if (b == std::get<2>(tuple)) { // match bin number
        if (std::get<1>(tuple) != "DEF") {
          signaLossMC[b - 1]->SetBinContent(std::stoi(std::get<1>(tuple)), std::get<0>(tuple));

        }
      }
    }
  }

  {
    Double_t xAxisLimits[2] = {cutMin - (cutMax - cutMin)*0.1, cutMax + (cutMax - cutMin)*0.1};
    Double_t yAxisLimitsUp[2] = {yMin, yMax};
    Double_t yAxisLimitsLow[2] = {yRatioMin, yRatioMax};
    TString xLabel = xLabelIn;

    Int_t colorPt[numPtBins];
    // Define colors
    for(int iClass = 0; iClass < numPtBins; iClass++) { 
      colorPt[iClass] = numPtBins - 1 + FIptxi - iClass;
    }

    TString legBits[numPtBins];
    for (Int_t m = 1; m <= numPtBins; m++) {
      legBits[m - 1] += Form("%.1f#minus%.1f%s", binpt[m-1], binpt[m], " (GeV/#it{c})");
    }

    // Plot The Signal Loss Fraction and DATA/MC ratios //
    TCanvas* canvasSignalLoss; 
    Double_t xLegend;
    Double_t padRightOffset;
    if (xLabel == "#Delta#it{M}_{#Lambda}, GeV/#it{c}^{2}" || xLabel == "#Delta#it{M}_{#Omega}, GeV/#it{c}^{2}") {
      canvasSignalLoss = new TCanvas("SignalLoss", "SignalLoss", 0, 70, 650, 850);
      padRightOffset = 0.1;
      xLegend = 0.503;
    } else {
      padRightOffset = 0.05;
      xLegend = 0.609;
      canvasSignalLoss = new TCanvas("SignalLoss", "SignalLoss", 0, 70, 620, 850);
    }
    TPad *padSignalLossUp = new TPad("padSignalLossUp", "padSignalLossUp", 0, 0.36, 1, 1);
    TPad *padSignalLossLow = new TPad("padSignalLossLow", "padSignalLossLow", 0, 0, 1, 0.36);
    StylePad(padSignalLossUp, 0.15, padRightOffset, 0.05, 0.);
    StylePad(padSignalLossLow, 0.15, padRightOffset, 0.02, 0.2);
    StyleCanvas(canvasSignalLoss, 0.15, 0.05, 0.05, 0.15);
    canvasSignalLoss->cd();
    padSignalLossUp->Draw();
    padSignalLossLow->Draw();

    // Dummy histograms
    TH1F *hDummySignalLoss = new TH1F("hDummySignalLoss", "hDummySignalLoss", 10000, xAxisLimits[0], xAxisLimits[1]);
    // Up
    for (Int_t i = 1; i <= hDummySignalLoss->GetNbinsX(); i++)
      hDummySignalLoss->SetBinContent(i, 10000);
    canvasSignalLoss->cd();
    padSignalLossUp->cd();
    hDummySignalLoss->GetYaxis()->SetMaxDigits(3);
    hDummySignalLoss->GetYaxis()->SetDecimals(kTRUE);
    StyleHistoMultiPlot(hDummySignalLoss, yAxisLimitsUp[0], yAxisLimitsUp[1], 1, 1, xLabel, "Signal(x)/Signal(loosest)", "", 0, 0, 0, 1.5, 1.2, 0, 0.0, 0.05, 0.0, 0.035, 0.005);
    SetTickLength(hDummySignalLoss, 0.025, 0.03);
    TAxis *axisSignalLossDummy = hDummySignalLoss->GetYaxis();
    axisSignalLossDummy->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
    hDummySignalLoss->GetXaxis()->SetRangeUser(xAxisLimits[0], xAxisLimits[1]);
    //padSignalLossUp->SetLogy();
    hDummySignalLoss->Draw("same");
    // Low
    TH1F *hDummySignalLossLow = new TH1F("hDummySignalLossLow", "hDummySignalLossLow", 10000, xAxisLimits[0], xAxisLimits[1]);
    for (Int_t i = 1; i <= hDummySignalLossLow->GetNbinsX(); i++)
      hDummySignalLossLow->SetBinContent(i, 1);
    padSignalLossLow->cd();
    StyleHistoMultiPlot(hDummySignalLossLow, yAxisLimitsLow[0], yAxisLimitsLow[1], 1, 1, xLabel, "DATA/MC", "", 0, 0, 0, 1.0, 0.9, 0, 0.08, 0.08, 0.08, 0.07, 0.01);
    SetTickLength(hDummySignalLossLow, 0.025, 0.03);
    hDummySignalLossLow->GetXaxis()->SetMaxDigits(3);
    hDummySignalLossLow->GetXaxis()->SetDecimals(kTRUE);
    TAxis *axisDummySignalLossLow = hDummySignalLossLow->GetYaxis();
    axisDummySignalLossLow->ChangeLabel(-1, -1, -1, -1, -1, -1, " "); // last 
    axisDummySignalLossLow->ChangeLabel(1, -1, -1, -1, -1, -1, " "); // first
    hDummySignalLossLow->GetXaxis()->SetRangeUser(xAxisLimits[0], xAxisLimits[1]);
    hDummySignalLossLow->GetYaxis()->CenterTitle();
    hDummySignalLossLow->Draw("same");

    // Yield Legend
    TLegend *legSignalLoss = new TLegend(0.186, 0.0625, 0.893, 0.3125);
    legSignalLoss->SetHeader("Open Circles: MC; Full ones: DATA ");
    legSignalLoss->SetNColumns(3);
    TLegendEntry *legYieldHeader = (TLegendEntry *)legSignalLoss->GetListOfPrimitives()->First();
    legYieldHeader->SetTextSize(0.04);
    //legSignalLoss->SetTextAlign(13); // adjust text wrt markers
    StyleLegend(legSignalLoss, 0.0, 0.0);

    TLegend *legendTitle = new TLegend(xLegend, 0.693913, 0.908, 0.893);
    StyleLegend(legendTitle, 0.0, 0.0);
    legendTitle->SetTextAlign(33);
    legendTitle->SetTextSize(0.05);
    legendTitle->SetTextFont(42);
    legendTitle->SetLineColorAlpha(0.,0.);
    legendTitle->SetFillColorAlpha(0.,0.);
    legendTitle->AddEntry("", "#bf{ALICE Work In Progress}", "");
    legendTitle->AddEntry("", "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    legendTitle->AddEntry("", particleSymnbols[nParticle] + ", |y| < 0.5", "");

    // Plotting
    for (Int_t b = 1; b <= numPtBins; b++) {
      padSignalLossUp->cd();
      // Data
      StyleHistoMultiPlot(signaLossDATA[b - 1], 0., 1., colorPt[b - 1], MarkerMult[0], "", "", "", 0, 0, 0, 1.5, 1.0, SizeMult[0], 0.0, 0.05, 0.0, 0.035, 0.005);
      signaLossDATA[b - 1]->Draw("same p");
      legSignalLoss->AddEntry(signaLossDATA[b - 1], legBits[b - 1], "pef");
      // MC
      StyleHistoMultiPlot(signaLossMC[b - 1], 0., 1., colorPt[b - 1], 4, "", "", "", 0, 0, 0, 1.5, 1.0, SizeMult[0], 0.0, 0.05, 0.0, 0.035, 0.005);
      signaLossMC[b - 1]->Draw("same p");
      //legSignalLoss->AddEntry(signaLossMC[b - 1], legBits[b - 1], "pef");

      padSignalLossLow->cd();
      signaLossRatio[b - 1] = (TH1F*)signaLossDATA[b - 1]->Clone(Form("Ratio-in-%i", b));
      signaLossRatio[b - 1]->Divide(signaLossMC[b - 1]);
      StyleHistoMultiPlot(signaLossRatio[b - 1], 0., 10, colorPt[b - 1], MarkerMult[0], "", "", "", 0, 0, 0, 1.5, 1.0, SizeMult[0], 0.0, 0.05, 0.0, 0.035, 0.005);
      signaLossRatio[b - 1]->Draw("same p");
    }

    DrawHorLineSyst(xAxisLimits[0], xAxisLimits[1], 1.02 , "2%");
    DrawHorLineSyst(xAxisLimits[0], xAxisLimits[1], 0.98 , "2%");

    padSignalLossUp->cd();
    legSignalLoss->Draw("same");
    legendTitle->Draw("same");
    outputfile->cd();
    canvasSignalLoss->Write();
  }

  // End of Code
  std::cout << "\x1B[1;32m"; // Set text color to green
  std::cout << "\n************* Topo/Kinem Study Finished Successfully! *************\n";
  std::cout << "\x1B[0m"; // Reset text color
}

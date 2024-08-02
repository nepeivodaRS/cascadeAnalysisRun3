#include "../../help.h"
#include "../../effHelp.h"

void topoStudy(const Int_t nParticle = 2, // 0-2 : xi, 3-5 : omega
               TString fileListDATA = "DATA.txt",
               TString fileListMC = "MC.txt",
               const std::string storageDir = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/topoStudy",
               const std::string postFix = "") // list of yields files
{
  // Start of Code
  std::cout << "\x1B[1;33m"; // Set text color to yellow
  std::cout << "\n************ Starting Topo/Kinem Study ************\n";
  std::cout << "\x1B[0m"; // Reset text color

  // Files (DATA) //
  std::vector<std::string> nameDATA;
  std::vector<std::string> nameLegendDATA;
  std::ifstream fileDATA(Form("%s", fileListDATA.Data()));

  std::string removeDATA = storageDir + "/dataPurity/yieldsOut/yield_XiPm_MB_inel0";
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

  std::string removeMC = storageDir + "/mcPurity/yieldsOut/yield_XiPm_MB_inel0";
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

  TFile *fileInMC[numFilesMC];
  TFile *fileInDATA[numFilesDATA];

  // output file
  TString outputfilePath;
  outputfilePath = storageDir + "/results/" + "topoStudy" + postFix + ".root";
  TFile *outputfile = new TFile(outputfilePath, "RECREATE");

  // Create vectors for DATA
  std::vector<std::tuple<Double_t, std::string, Int_t>> yieldDATA; // value, config #, bin number
  std::vector<std::tuple<Double_t, std::string, Int_t>> varsDATA; // value, config #, bin number
  std::vector<std::pair<Double_t, std::string>> purityDATA; // value, config #, bin number

  // Create vectors for MC
  std::vector<std::tuple<Double_t, std::string, Int_t>> yieldMC; // value, config #, bin number
  std::vector<std::tuple<Double_t, std::string, Int_t>> varsMC; // value, config #, bin number

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

  // DATA //
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

    TH1F* hPurity = (TH1F *)dirFitParams->Get("Purity_" + particleNames[nParticle]);
    if (!hPurity)
    {
      std::cerr << "Histogram `hPurity` is not found!" << std::endl;
      return;
    }

    TH1F* hYield = (TH1F *)dirFitParams->Get("Yield_" + particleNames[nParticle]);
    if (!hYield)
    {
      std::cerr << "Histogram `hYield` is not found!" << std::endl;
      return;
    }

    Double_t purityValue = 0;
    for (Int_t b = 1; b <= numPtBins; b++) {
      Double_t yieldValue = hYield->GetBinContent(b);
      purityValue += hPurity->GetBinContent(b);
      yieldDATA.push_back(std::make_tuple(yieldValue, nameLegendDATA[iFile], b));
    }
    purityDATA.push_back(std::make_pair(purityValue/numPtBins, nameLegendMC[iFile]));
  }

  // Find the loosest config in DATA
  std::vector<std::pair<Double_t, Int_t>> yieldLoosestDATA; // yield, bin number
  for (const auto& tuple : yieldDATA) {
      if (std::get<1>(tuple) == "LOOSEST") {
        yieldLoosestDATA.push_back(std::make_pair(std::get<0>(tuple), std::get<2>(tuple)));
      }
  }

  // Make signal loss fraction
  for (const auto& tuple : yieldDATA) {
    for (const auto& pair : yieldLoosestDATA) {
      if (pair.second == std::get<2>(tuple)) { // match bin number
        Double_t recYield = std::get<0>(tuple);
        Double_t yeildVariation = recYield / pair.first;
        varsDATA.push_back(std::make_tuple(yeildVariation, std::get<1>(tuple), std::get<2>(tuple)));
      }
    }
  }

  // MC //
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

    TH1F* hPurity = (TH1F *)dirFitParams->Get("Purity_" + particleNames[nParticle]);
    if (!hPurity)
    {
      std::cerr << "Histogram `hPurity` is not found!" << std::endl;
      return;
    }

    //TH1F* hYield = (TH1F *)dirFitParams->Get("Yield_" + particleNames[nParticle]);
    TH1F* hYield = (TH1F *)dirFitParams->Get("hGenCountsClone"); // take true gen. number of rec. cascades
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
      if (std::get<1>(tuple) == "LOOSEST") {
        yieldLoosestMC.push_back(std::make_pair(std::get<0>(tuple), std::get<2>(tuple)));
      }
  }

  // Make signal loss fraction
  for (const auto& tuple : yieldMC) {
    for (const auto& pair : yieldLoosestMC) {
      if (pair.second == std::get<2>(tuple)) { // match bin number
        Double_t recYield = std::get<0>(tuple);
        Double_t yeildVariation = recYield / pair.first;
        //std::cout << "yield fraction: " << yeildVariation << " config: " << std::get<1>(tuple) << std::endl;
        varsMC.push_back(std::make_tuple(yeildVariation, std::get<1>(tuple), std::get<2>(tuple)));
      }
    }
  }

  TH1F* signaLossDATA[numPtBins];
  TH1F* signaLossMC[numPtBins];
  TH1F* signaLossRatio[numPtBins];
  for (Int_t b = 1; b <= numPtBins; b++) {
    signaLossDATA[b - 1] = new TH1F(Form("signalLossFractionDATA_%i", b), Form("signalLossFractionDATA_%i", b), numFilesDATA, 1. / (numFilesDATA - 1) / 2, 1 + 1. / (numFilesDATA - 1) / 2);
    signaLossMC[b - 1] = new TH1F(Form("signalLossFractionMC_%i", b), Form("signalLossFractionMC_%i", b), numFilesMC, 1. / (numFilesMC - 1) / 2, 1 + 1. / (numFilesMC - 1) / 2);
  }

  Double_t signalLossBestPurity[numPtBins];
  Double_t signalLossBestPurityMC[numPtBins];
  TString sBestPurity = "428";

  std::vector<std::tuple<Double_t, std::string, Int_t>> varsDATAsorted; // value, config #, bin number
  for (Int_t b = 1; b <= numPtBins; b++) {
    // DATA //
    std::vector<std::pair<Double_t, std::string>> varsInBinDATA; // yield, config #
    for (const auto& tuple : varsDATA) {
      if (b == std::get<2>(tuple)) { // match bin number
        varsInBinDATA.push_back(std::make_pair(std::get<0>(tuple), std::get<1>(tuple)));
        if (sBestPurity == std::get<1>(tuple)) {
          signalLossBestPurity[b - 1] = std::get<0>(tuple);
        }
      }
    }
    // Sort the vector by keys using a lambda comparator
    std::sort(varsInBinDATA.begin(), varsInBinDATA.end(), [](const auto& a, const auto& b) {
        return a.first > b.first;
    });


    for (const auto& pair : varsInBinDATA) {
      varsDATAsorted.push_back(std::make_tuple(pair.first, pair.second, b));
    }

    Int_t bestPurityNumber = 0;
    Int_t bestPurityCounter = 0;
    for (const auto& pair : varsInBinDATA) {
      bestPurityCounter += 1;
      if (sBestPurity == std::get<1>(pair)) {
        bestPurityNumber = bestPurityCounter;
      }
    }

    std::cout << "Conseq. number of the best purity config in signalLoss-sorted bin number " << b << " is " << bestPurityNumber << std::endl;

    for (const auto& tuple : varsMC) {
      if (sBestPurity == std::get<1>(tuple)) {
        signalLossBestPurityMC[std::get<2>(tuple) - 1] = std::get<0>(tuple); // to track signal loss in MC in the config corresponding to the best purity in DATA
      }
    }

    for (Int_t bin = 1; bin <= signaLossDATA[0]->GetNbinsX(); bin++) {
      signaLossDATA[b - 1]->SetBinContent(bin, varsInBinDATA[bin-1].first);
      signaLossDATA[b - 1]->SetBinError(bin, 0.);
      // MC //
      for (const auto& tuple : varsMC) {
        if (varsInBinDATA[bin-1].second == std::get<1>(tuple) && b == std::get<2>(tuple)) { // match the tag and pt bin
          signaLossMC[b - 1]->SetBinContent(bin, std::get<0>(tuple));
          signaLossMC[b - 1]->SetBinError(bin, 0.);
        }
      }
    }
  }

  // Sort the vector by keys using a lambda comparator
  std::sort(purityDATA.begin(), purityDATA.end(), [](const auto& a, const auto& b) {
      return a.first > b.first;
  });

  for (const auto& pair : purityDATA) {
      std::cout << pair.first << ": " << pair.second << std::endl;
  }

  for (const auto& value : signalLossBestPurity) {
      std::cout << "signal loss in purity-best config in DATA: " << value << std::endl;
  }

  for (const auto& value : signalLossBestPurityMC) {
      std::cout << "signal loss in purity-best config in MC: " << value << std::endl;
  }

  {
    Double_t cutMin = 0;
    Double_t cutMax = 1;
    Double_t xAxisLimits[2] = {cutMin - (cutMax - cutMin)*0.1, cutMax + (cutMax - cutMin)*0.1};
    Double_t yAxisLimitsUp[2] = {-0.1, 1.5};
    Double_t yAxisLimitsLow[2] = {0.8, 1.2};
    TString xLabel = "Configs (bin-by-bin data-sorted by signal loss)";

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
    TCanvas* canvasSignalLoss = new TCanvas("SignalLoss", "SignalLoss", 0, 70, 620, 850);
    StyleCanvas(canvasSignalLoss, 0.15, 0.05, 0.05, 0.15);
    TPad *padSignalLossUp = new TPad("padSignalLossUp", "padSignalLossUp", 0, 0.36, 1, 1);
    TPad *padSignalLossLow = new TPad("padSignalLossLow", "padSignalLossLow", 0, 0, 1, 0.36);
    StylePad(padSignalLossUp, 0.15, 0.05, 0.05, 0.);
    StylePad(padSignalLossLow, 0.15, 0.05, 0.02, 0.2);
    canvasSignalLoss->cd();
    padSignalLossUp->Draw();
    padSignalLossLow->Draw();

    // Dummy histograms
    TH1F *hDummySignalLoss = new TH1F("hDummySignalLoss", "hDummySignalLoss", 10000, xAxisLimits[0], xAxisLimits[1]);
    // Up
    for (Int_t i = 1; i <= hDummySignalLoss->GetNbinsX(); i++)
      hDummySignalLoss->SetBinContent(i, 1000);
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

    Int_t scaling = 20;

    TLegend *legendTitleLow = new TLegend(0.609, 0.693913, 0.908, 0.893);
    StyleLegend(legendTitleLow, 0.0, 0.0);
    legendTitleLow->SetTextAlign(33);
    legendTitleLow->SetTextSize(0.09);
    legendTitleLow->SetTextFont(42);
    legendTitleLow->SetLineColorAlpha(0.,0.);
    legendTitleLow->SetFillColorAlpha(0.,0.);
    legendTitleLow->AddEntry("", "#bf{1000} Configs (15 varied cuts)", "");
    legendTitleLow->AddEntry("", Form("Rebin factor: %d", scaling), "");

    // Plotting
    for (Int_t b = 1; b <= numPtBins; b++) {
      padSignalLossUp->cd();

      // Rebin for a better representation
      signaLossDATA[b - 1]->Rebin(scaling);
      signaLossMC[b - 1]->Rebin(scaling);
      signaLossDATA[b - 1]->Scale(1./scaling);
      signaLossMC[b - 1]->Scale(1./scaling);

      // Data
      StyleHistoMultiPlot(signaLossDATA[b - 1], 0., 1., colorPt[b - 1], MarkerMult[0], "", "", "", 0, 0, 0, 1.5, 1.0, SizeMult[0], 0.0, 0.05, 0.0, 0.035, 0.005);
      signaLossDATA[b - 1]->Draw("same p");
      legSignalLoss->AddEntry(signaLossDATA[b - 1], legBits[b - 1], "pef");
      // MC
      StyleHistoMultiPlot(signaLossMC[b - 1], 0., 1., colorPt[b - 1], 4, "", "", "", 0, 0, 0, 1.5, 1.0, SizeMult[0], 0.0, 0.05, 0.0, 0.035, 0.005);
      signaLossMC[b - 1]->Draw("same p");
      //legSignalLoss->AddEntry(signaLossMC[b - 1], legBits[b - 1], "pef");

      //DrawHorLineColor(1.0, signalLossBestPurity[b - 1], colorPt[b - 1]);

      padSignalLossLow->cd();
      signaLossRatio[b - 1] = (TH1F*)signaLossDATA[b - 1]->Clone(Form("Ratio-in-%i", b));
      signaLossRatio[b - 1]->Divide(signaLossMC[b - 1]);
      if (b == 1) {
        std::cout << std::endl << "STUDY OF RELIABLE CONFIGS" << std::endl << std::endl;
        for (Int_t i = 1; i <= signaLossRatio[b - 1]->GetNbinsX(); i++) {
          if(signaLossRatio[b - 1]->GetBinContent(i) > 0.995) {
            std::cout << "************************" << std::endl;
            std::cout << "conseq. number of reliable config: " << i << " ratio is: " << signaLossRatio[b - 1]->GetBinContent(i) << std::endl;
            Double_t iter = 0;
            for (const auto& tuple : varsDATAsorted) {
              if (std::get<2>(tuple) == b) {
                iter++;
                if (iter == i) {
                  std::cout << "config number of conseq. config number " << i << " is " << std::get<1>(tuple) << std::endl;
                  for (const auto& pair : purityDATA) {
                    if (pair.second == std::get<1>(tuple)) {
                      std::cout << "purity for the config " << pair.second << " is " << pair.first << std::endl;
                      std::cout << "************************" << std::endl << std::endl;
                    }
                  }
                }
              }
            }
          }
        }
      }

      StyleHistoMultiPlot(signaLossRatio[b - 1], 0., 10, colorPt[b - 1], MarkerMult[0], "", "", "", 0, 0, 0, 1.5, 1.0, SizeMult[0], 0.0, 0.05, 0.0, 0.035, 0.005);
      signaLossRatio[b - 1]->Draw("same p");
    }

    DrawHorLineSyst(xAxisLimits[0], xAxisLimits[1], 1.02 , "2%");
    DrawHorLineSyst(xAxisLimits[0], xAxisLimits[1], 0.98 , "2%");

    padSignalLossUp->cd();
    //DrawVertLine(signaLossDATA[0]->GetBinCenter(double(bestPurityNumber)/scaling), 0.5, 1.0, kBlack);
    legSignalLoss->Draw("same");
    legendTitle->Draw("same");
    padSignalLossLow->cd();
    legendTitleLow->Draw("same");
    outputfile->cd();
    canvasSignalLoss->Write();
  }

  std::cout << "                Purity-Best config " << std::endl;
  for (Int_t b = 1; b <= numPtBins; b++) {
    std::cout << "DATA/MC signal loss: " << signalLossBestPurity[b-1]/signalLossBestPurityMC[b-1] << " in " << binpt[b-1] << " < p_T < " << binpt[b] << std::endl;
  }



  // End of Code
  std::cout << "\x1B[1;32m"; // Set text color to green
  std::cout << "\n************* Topo/Kinem Study Finished Successfully! *************\n";
  std::cout << "\x1B[0m"; // Reset text color
}

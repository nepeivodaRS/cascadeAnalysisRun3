#include "../../help.h"
#include "../../effHelp.h"

void topoStudy2D(const Int_t nParticle = 2, // 0-2 : xi, 3-5 : omega
                         TString fileListDATA = "DATA.txt",
                         TString fileListMC = "MC.txt",
                         const TString xLabelIn = "cos(PA)", const TString yLabelIn = "cos(PA)",
                         Double_t cutMin1 = 0.970, Double_t cutMax1 = 0.999,
                         Double_t cutMin2 = 0.970, Double_t cutMax2 = 0.999,
                         const Int_t N1 = 2, const Int_t N2 = 3,
                         const Int_t compOrder = 1, // see code for comments
                         const std::string workingDir = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/topoStudy/casccospa",
                         const std::string postFix = "_casccospa")
{
  // Start of Code
  std::cout << "\x1B[1;33m"; // Set text color to yellow
  std::cout << "\n************ Starting Topo/Kinem 2D Study ************\n";
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
  std::vector<std::pair<Double_t, std::string>> yieldIntDATA; // value, config #
  std::vector<std::pair<Double_t, std::string>> varsIntDATA; // value, config #
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

    Double_t totYield = 0;
    for (Int_t b = 1; b <= numPtBins; b++) {
      Double_t yieldValue = hYield->GetBinContent(b);
      totYield += hYield->GetBinContent(b);
      yieldDATA.push_back(std::make_tuple(yieldValue, nameLegendDATA[iFile], b));
    }
    yieldIntDATA.push_back(std::make_pair(totYield, nameLegendDATA[iFile]));
  }

  // Find the loosest config in DATA
  std::string sLoosestNumber;
  switch (compOrder) {
    case 1:
      sLoosestNumber = "1"; // > cut1, > cut2
      break;
    case 2:
      sLoosestNumber = std::to_string(N2); // > cut1, < cut2
      break;
    case 3:
      sLoosestNumber = std::to_string(N1*N2 - (N2 - 1)); // < cut1, > cut2
      break;
    case 4:
      sLoosestNumber = std::to_string(N1*N2); // < cut1, < cut2
      break;
    default:
      std::cerr << "Wrong order type" << std::endl;
      break;
  }

  // in pt bins
  std::vector<std::pair<Double_t, Int_t>> yieldLoosestDATA; // yield, bin number
  for (const auto& tuple : yieldDATA) {
      if (std::get<1>(tuple) == sLoosestNumber) {
        yieldLoosestDATA.push_back(std::make_pair(std::get<0>(tuple), std::get<2>(tuple)));
      }
  }

  // integrated
  Double_t yieldLoosestIntDATA;
  for (const auto& pair : yieldIntDATA) {
      if (std::get<1>(pair) == sLoosestNumber) {
        yieldLoosestIntDATA = std::get<0>(pair);
      }
  }

  // Make signal loss fraction
  // in pt bins
  for (const auto& tuple : yieldDATA) {
    for (const auto& pair : yieldLoosestDATA) {
      if (pair.second == std::get<2>(tuple)) { // match bin number
        Double_t recYield = std::get<0>(tuple);
        Double_t yieldVariation = recYield / pair.first;
        varsDATA.push_back(std::make_tuple(yieldVariation, std::get<1>(tuple), std::get<2>(tuple)));
      }
    }
  }

  // integrated
  for (const auto& pair : yieldIntDATA) {
    Double_t recYield = std::get<0>(pair);
    Double_t yieldVariation = recYield / yieldLoosestIntDATA;
    varsIntDATA.push_back(std::make_pair(yieldVariation, std::get<1>(pair)));
  }

  // MC //
  std::vector<std::tuple<Double_t, std::string, Int_t>> yieldMC; // value, config #, bin number
  std::vector<std::tuple<Double_t, std::string, Int_t>> varsMC; // value, config #, bin number
  std::vector<std::pair<Double_t, std::string>> yieldIntMC; // value, config #
  std::vector<std::pair<Double_t, std::string>> varsIntMC; // value, config #
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
    Double_t totYield = 0;
    for (Int_t b = 1; b <= numPtBins; b++) {
      Double_t yieldValue = hYield->GetBinContent(b);
      totYield += hYield->GetBinContent(b);
      yieldMC.push_back(std::make_tuple(yieldValue, nameLegendMC[iFile], b));
    }
    yieldIntMC.push_back(std::make_pair(totYield, nameLegendMC[iFile]));
  }

  // Find the loosest config in MC
  // in pt bins
  std::vector<std::pair<Double_t, Int_t>> yieldLoosestMC; // yield, bin number
  for (const auto& tuple : yieldMC) {
      if (std::get<1>(tuple) == sLoosestNumber) {
        yieldLoosestMC.push_back(std::make_pair(std::get<0>(tuple), std::get<2>(tuple)));
      }
  }

  // integrated
  Double_t yieldLoosestIntMC;
  for (const auto& pair : yieldIntMC) {
      if (std::get<1>(pair) == sLoosestNumber) {
        yieldLoosestIntMC = std::get<0>(pair);
      }
  }  

  // Make signal loss fraction
  // in pt bins
  for (const auto& tuple : yieldMC) {
    for (const auto& pair : yieldLoosestMC) {
      if (pair.second == std::get<2>(tuple)) { // match bin number
        Double_t recYield = std::get<0>(tuple);
        Double_t yieldVariation = recYield / pair.first;
        varsMC.push_back(std::make_tuple(yieldVariation, std::get<1>(tuple), std::get<2>(tuple)));
      }
    }
  }

  // integrated
  for (const auto& pair : yieldIntMC) {
    Double_t recYield = std::get<0>(pair);
    Double_t yieldVariation = recYield / yieldLoosestIntMC;
    varsIntMC.push_back(std::make_pair(yieldVariation, std::get<1>(pair)));
  }

  TH2F* signaLossDATA[numPtBins];
  TH2F* signaLossDATAInt = new TH2F("signalLossFractionDATA_int", "signalLossFractionDATA_int", N1, cutMin1 - (cutMax1 - cutMin1) / (numFilesDATA - 2) / 2, cutMax1 + (cutMax1 - cutMin1) / (numFilesDATA - 2) / 2, 
                                                                                                N2, cutMin2 - (cutMax2 - cutMin2) / (numFilesDATA - 2) / 2, cutMax2 + (cutMax2 - cutMin2) / (numFilesDATA - 2) / 2);
  TH2F* signaLossMC[numPtBins];
  TH2F* signaLossMCInt = new TH2F("signalLossFractionMC_int", "signalLossFractionMC_int", N1, cutMin1 - (cutMax1 - cutMin1) / (numFilesMC - 2) / 2, cutMax1 + (cutMax1 - cutMin1) / (numFilesMC - 2) / 2, 
                                                                                          N2, cutMin2 - (cutMax2 - cutMin2) / (numFilesMC - 2) / 2, cutMax2 + (cutMax2 - cutMin2) / (numFilesMC - 2) / 2);  
  TH2F* signaLossRatio[numPtBins];
  for (Int_t b = 1; b <= numPtBins; b++) {
    signaLossDATA[b - 1] = new TH2F(Form("signalLossFractionDATA_%i", b), Form("signalLossFractionDATA_%i", b), N1, cutMin1 - (cutMax1 - cutMin1) / (numFilesDATA - 2) / 2, cutMax1 + (cutMax1 - cutMin1) / (numFilesDATA - 2) / 2, 
                                                                                                                N2, cutMin2 - (cutMax2 - cutMin2) / (numFilesDATA - 2) / 2, cutMax2 + (cutMax2 - cutMin2) / (numFilesDATA - 2) / 2);
    signaLossMC[b - 1] = new TH2F(Form("signalLossFractionMC_%i", b), Form("signalLossFractionMC_%i", b), N1, cutMin1 - (cutMax1 - cutMin1) / (numFilesMC - 2) / 2, cutMax1 + (cutMax1 - cutMin1) / (numFilesMC - 2) / 2, 
                                                                                                          N2, cutMin2 - (cutMax2 - cutMin2) / (numFilesMC - 2) / 2, cutMax2 + (cutMax2 - cutMin2) / (numFilesMC - 2) / 2);
  }

  // DATA
  // in pt bins
  for (Int_t b = 1; b <= numPtBins; b++) {
    for (auto& tuple : varsDATA) {
      if (b == std::get<2>(tuple)) { // match bin number
        if (std::get<1>(tuple) != "DEF") {
          Int_t varValue = std::stoi(std::get<1>(tuple));
          signaLossDATA[b - 1]->SetBinContent(varValue / N2 + 1, varValue % N2 + 1, std::get<0>(tuple));
        }
      }
    }
  }

  // integrated
  for (auto& pair : varsIntDATA) {
    if (std::get<1>(pair) != "DEF") {
      Int_t varValue = std::stoi(std::get<1>(pair));
      //std::cout << "config: " << varValue << " x: " << (varValue - 1) / N2 + 1 << " y: " << (varValue - 1) % N2 + 1 << " z: " << std::get<0>(pair) << std::endl;
      signaLossDATAInt->SetBinContent( (varValue - 1) / N2 + 1, (varValue - 1) % N2 + 1, std::get<0>(pair));
    }
  }

  // MC
  // in pt bins
  for (Int_t b = 1; b <= numPtBins; b++) {
    for (auto& tuple : varsMC) {
      if (b == std::get<2>(tuple)) { // match bin number
        if (std::get<1>(tuple) != "DEF") {
          Int_t varValue = std::stoi(std::get<1>(tuple));
          signaLossMC[b - 1]->SetBinContent(varValue / N2 + 1, varValue % N2 + 1, std::get<0>(tuple));
        }
      }
    }
  }

  // integrated
  for (auto& pair : varsIntMC) {
    if (std::get<1>(pair) != "DEF") {
      Int_t varValue = std::stoi(std::get<1>(pair));
      //std::cout << "config: " << varValue << " x: " << (varValue - 1) / N2 + 1 << " y: " << (varValue - 1) % N2 + 1 << " z: " << std::get<0>(pair) << std::endl;
      signaLossMCInt->SetBinContent( (varValue - 1) / N2 + 1, (varValue - 1) % N2 + 1, std::get<0>(pair));
    }
  }

  // for (Int_t b = 1; b <= numPtBins; b++) {
  //   signaLossDATAInt->Add(signaLossDATA[b - 1]);
  // }

  TCanvas* canvas2D = new TCanvas("canvas2D", "canvas2D", 1000,800);
  canvas2D->Divide(2, 2);
  StyleCanvas(canvas2D, 0.15, 0.05, 0.05, 0.15);
  for (Int_t i = 1; i <= 4; i++) {
    canvas2D->cd(i)->SetMargin(0.18, 0.15, 0.11, 0.1);
  }

  canvas2D->cd(1);
  gStyle->SetOptStat(0);
  Style2DHisto(signaLossDATAInt, xLabelIn, yLabelIn, "Signal Fraction in DATA wrt loosest cut", 0.05, 0.05, 0.05, 0.05, 0.05);
  signaLossDATAInt->GetYaxis()->SetTitleOffset(1.4);
  signaLossDATAInt->Draw("COLZ");

  canvas2D->cd(2);
  Style2DHisto(signaLossMCInt, xLabelIn, yLabelIn, "Signal Fraction in MC wrt loosest cut", 0.05, 0.05, 0.05, 0.05, 0.05);
  signaLossMCInt->GetYaxis()->SetTitleOffset(1.4);
  signaLossMCInt->Draw("COLZ");

  canvas2D->cd(3);
  TH2F* signaLossRatioInt = (TH2F*)signaLossDATAInt->Clone("RatioInt");
  signaLossRatioInt->Divide(signaLossMCInt);
  Style2DHisto(signaLossRatioInt, xLabelIn, yLabelIn, "DATA/MC", 0.05, 0.05, 0.05, 0.05, 0.05);
  signaLossRatioInt->Draw("COLZ");

  TH2F* signaLossRatioInt2Perc = (TH2F*)signaLossRatioInt->Clone("Ratio2perc");
  TH1F* signaLossRatioInt2PercXaxis = new TH1F("signaLossRatioInt2PercXaxis", "signaLossRatioInt2PercXaxis", N1, cutMin1 - (cutMax1 - cutMin1) / (numFilesMC - 2) / 2, cutMax1 + (cutMax1 - cutMin1) / (numFilesMC - 2) / 2);
  signaLossRatioInt2Perc->Reset();
  for (Int_t i = 1; i <= signaLossRatioInt->GetNbinsX(); i++) {
    Double_t averageY = 0;
    Int_t counter = 0;
    for (Int_t j = 1; j <= signaLossRatioInt->GetNbinsY(); j++) {
      //std::cout << "bincontent: " << signaLossRatioInt->GetBinContent(i, j) << std::endl;
      if (signaLossRatioInt->GetBinContent(i, j) < 0.978 || signaLossRatioInt->GetBinContent(i, j) > 0.982) {
      } else {
        Double_t content = signaLossRatioInt->GetBinContent(i, j);
        averageY += signaLossRatioInt->GetYaxis()->GetBinCenter(j);
        counter++;
        signaLossRatioInt2Perc->SetBinContent(i, j, content);
      }
    }
    if(counter != 0) {
      signaLossRatioInt2PercXaxis->SetBinContent(i, averageY/counter);
    }
  }

  canvas2D->cd(4);
  Style2DHisto(signaLossRatioInt2Perc, xLabelIn, yLabelIn, "Region of DATA/MC = 0.980 #pm 0.002 (2% border)", 0.05, 0.05, 0.05, 0.05, 0.05);
  signaLossRatioInt2Perc->SetMinimum(signaLossRatioInt->GetMinimum());
  signaLossRatioInt2Perc->SetMaximum(signaLossRatioInt->GetMaximum());
  signaLossRatioInt2Perc->Draw("COLZ");

  if (postFix == "_casccospa-ntpccrrows") {
    gROOT->SetBatch(kTRUE);
    TCanvas* testcanvas = new TCanvas("test", "test", 800, 600);
    testcanvas->cd();
    TF1 *fitExp = new TF1("fitExp", "[0]*exp([1]*x-[2]) + [3]", 0.97, 0.999);
    fitExp->SetParameters(-1.5, 0.9, 0.5, 63);
    signaLossRatioInt2PercXaxis->Fit(fitExp, "QR");
    Double_t fitPars[4];
    fitExp->GetParameters(fitPars);
    gROOT->SetBatch(kTRUE);
    fitExp->SetNpx(1e5);
    fitExp->SetLineWidth(3);
    fitExp->SetLineColor(kCyan+2);
    gROOT->SetBatch(kFALSE);
    canvas2D->cd(4);
    fitExp->Draw("same");
    TLatex *label = new TLatex(cutMin1 + (cutMax1 - cutMin1)/20., cutMin2 + (cutMax2 - cutMin2)/1.2, Form("y = %.1e*exp(%.1f*x-%.1f) + %.1f", fitPars[0], fitPars[1], fitPars[2], fitPars[3]));
    label->SetTextSize(0.05); // Set the text size
    label->SetTextAngle(0);
    label->SetTextColor(kBlack); // Set the text color to match the line
    label->Draw("same");    
  }


  // End of Code
  std::cout << "\x1B[1;32m"; // Set text color to green
  std::cout << "\n************* Topo/Kinem 2D Study Finished Successfully! *************\n";
  std::cout << "\x1B[0m"; // Reset text color
}

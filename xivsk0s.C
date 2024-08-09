#include "help.h"

void xivsk0s(const TString workingDir = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisPostSQM",
             const TString postFix = "")
{
  // Start of Code
  std::cout << "\x1B[1;33m"; // Set text color to yellow
  std::cout << "\n************ Starting Plotting Integrated k0s vs xi Yields ************\n";
  std::cout << "\x1B[0m"; // Reset text color

  // Define colors
  for(int iClass = 0; iClass < numMult + 1; iClass++) { 
    if(iClass==0) {
      color[iClass] = kBlack;
    } else {
      color[iClass] = numMult + FI - iClass;
    }
  }

  // Open k0s File //
  //TFile* file_k0s = TFile::Open(workingDir + "/data/k0s/YieldsIntegrated_NoTOF.root");
  TFile* file_k0s = TFile::Open(workingDir + "/data/k0s/YieldsIntegrated_newsel.root");
  if (!file_k0s || file_k0s->IsZombie()) {
      std::cerr << "Error opening k0s file" << std::endl;
      return;
  }

  // Retrieve the TGraphErrors object from the file
  TGraphErrors* graph_k0s = dynamic_cast<TGraphErrors*>(file_k0s->Get("YieldsNchStat"));
  if (!graph_k0s) {
      std::cerr << "Error: Could not find graph_k0s in k0s file " << std::endl;
      file_k0s->Close();
      return;
  }

  // Open xi File //
  TFile* file_xi = TFile::Open(workingDir + "/yieldInMultFitted/yieldInMultFitted_XiPm_inel0-medium-tight.root");
  if (!file_xi || file_xi->IsZombie()) {
      std::cerr << "Error opening k0s file" << std::endl;
      return;
  }

  // Retrieve the canvas_xi by name
  TCanvas* canvas_xi = dynamic_cast<TCanvas*>(file_xi->Get("canvasIntegratedYieldsdNdEta"));
  if (!canvas_xi) {
      std::cerr << "Error retrieving canvas_xi" << std::endl;
      file_xi->Close();
      return;
  }

  // Get the list of primitives_xi in the canvas_xi
  TList* primitives_xi = canvas_xi->GetListOfPrimitives();
  if (!primitives_xi) {
      std::cerr << "No primitives_xi found in the canvas_xi" << std::endl;
      file_xi->Close();
      return;
  }

  // Iterate over the primitives and print their names and types; find multigraph and take it
  TIter next(primitives_xi);
  TObject* obj = nullptr;
  TMultiGraph* multiGraph = nullptr;
  std::cout << "Primitives in Xi canvas: " << std::endl;
  while ((obj = next())) {
      std::cout << "  Name: " << obj->GetName() << ", Type: " << obj->ClassName() << std::endl;
      if(strcmp(obj->ClassName(), "TMultiGraph") == 0){
        multiGraph = dynamic_cast<TMultiGraph*>(obj);
      }
  }

  if (!multiGraph) {
      std::cerr << "Error: Could not find TMultiGraph in canvas" << std::endl;
      file_xi->Close();
      return;
  }

  std::vector<Double_t> yield_k0s;
  std::vector<Double_t> yield_xi;

  std::vector<Double_t> error_k0s;
  std::vector<Double_t> error_xi;

  // Get Xi yields //
  // Loop over all graphs in the TMultiGraph and print their points and errors
  TList* graphList = multiGraph->GetListOfGraphs();
  TIter next_multiGraph(graphList);
  TObject* graphObj;
  int graphIndex = 0;

  while ((graphObj = next_multiGraph())) {
      TGraphAsymmErrors* graph = dynamic_cast<TGraphAsymmErrors*>(graphObj);
      if (graph) {
          //std::cout << "Graph " << graphIndex << " Graph Name: " << graph->GetName() << "\n";
          if(graphIndex == 0){
            std::cout << "Xi yields: " << std::endl;
            int nPoints = graph->GetN();
            for (int i = 0; i < nPoints; ++i) {
              Double_t x, y;
              graph->GetPoint(i, x, y);
              Double_t ex = graph->GetErrorX(i);
              Double_t ey = graph->GetErrorY(i);
              yield_xi.push_back(y);
              error_xi.push_back(ey);
              std::cout << "  Point " << i << ": x = " << x << ", y = " << y << ", ex = " << ex << ", ey = " << ey << std::endl;
            }
          }
          ++graphIndex;
      }
  }

  // Get k0s yields //
  // Access the number of points in the graph
  int nPoints_k0s = graph_k0s->GetN();
  std::cout << "k0s yields: " << std::endl;
  for (int i = 0; i < nPoints_k0s; ++i) {
    Double_t x, y;
    graph_k0s->GetPoint(i, x, y);
    Double_t ex = graph_k0s->GetErrorX(i);
    Double_t ey = graph_k0s->GetErrorY(i);
    yield_k0s.push_back(y);
    error_k0s.push_back(ey);
    std::cout << "  Point " << i << ": x = " << x << ", y = " << y << ", ex = " << ex << ", ey = " << ey << std::endl;
  }

  // Create TGraph k0svsxi Run3 //
  Double_t* arr_yield_k0s = new Double_t[yield_k0s.size()];
  std::copy(yield_k0s.begin(), yield_k0s.end(), arr_yield_k0s);
  Double_t* arr_yield_xi = new Double_t[yield_xi.size()];
  std::copy(yield_xi.begin(), yield_xi.end(), arr_yield_xi);

  Double_t* arr_err_k0s = new Double_t[error_k0s.size()];
  std::copy(error_k0s.begin(), error_k0s.end(), arr_err_k0s);
  Double_t* arr_err_xi = new Double_t[error_xi.size()];
  std::copy(error_xi.begin(), error_xi.end(), arr_err_xi);

  TGraphErrors *g_k0svsxi = new TGraphErrors(nPoints_k0s, arr_yield_k0s, arr_yield_xi, arr_err_k0s, arr_err_xi);
  g_k0svsxi->SetMarkerColor(color[0]);
  g_k0svsxi->SetMarkerStyle(MarkerMult[0]);
  g_k0svsxi->SetMarkerSize(1.5);

  // Get Published Data (V0M) //
  TString inputHEPDataPath = workingDir + "/data/published/HEPData-ins1748157-v2-Table_8c.root";
  TFile* fileHEPDataIn = TFile::Open(inputHEPDataPath);
  if (!fileHEPDataIn || fileHEPDataIn->IsZombie()) {
      std::cerr << "Error opening HEP data file!" << std::endl;
      return;
  }

  TH1F* hHEP_yield_k0s = (TH1F *)fileHEPDataIn->Get("Table 8c/Hist1D_y1");
  TH1F* hHEP_yield_xi = (TH1F *)fileHEPDataIn->Get("Table 8c/Hist1D_y3");

  TH1F* hHEP_e1_k0s = (TH1F *)fileHEPDataIn->Get("Table 8c/Hist1D_y1_e1");
  TH1F* hHEP_e1_xi = (TH1F *)fileHEPDataIn->Get("Table 8c/Hist1D_y3_e1");

  TH1F* hHEP_e2_k0s = (TH1F *)fileHEPDataIn->Get("Table 8c/Hist1D_y1_e2");
  TH1F* hHEP_e2_xi = (TH1F *)fileHEPDataIn->Get("Table 8c/Hist1D_y3_e2");

  TH1F* hHEP_e3_k0s = (TH1F *)fileHEPDataIn->Get("Table 8c/Hist1D_y1_e3");
  TH1F* hHEP_e3_xi = (TH1F *)fileHEPDataIn->Get("Table 8c/Hist1D_y3_e3");

  if (!hHEP_yield_k0s || !hHEP_e1_k0s || !hHEP_yield_xi || !hHEP_e1_xi || !hHEP_e2_k0s || !hHEP_e2_xi || !hHEP_e3_k0s || !hHEP_e3_xi)
  {
    std::cerr << "Histogram HEPdata is not found!" << std::endl;
    return;
  }

  Double_t* arr_yield_k0s_published_v0m = new Double_t[hHEP_yield_k0s->GetEntries()];
  Double_t* arr_yield_xi_published_v0m = new Double_t[hHEP_yield_xi->GetEntries()];

  Double_t* arr_err1_k0s_published_v0m = new Double_t[hHEP_e1_k0s->GetEntries()];
  Double_t* arr_err1_xi_published_v0m = new Double_t[hHEP_e1_xi->GetEntries()];

  Double_t* arr_err2_k0s_published_v0m = new Double_t[hHEP_e2_k0s->GetEntries()];
  Double_t* arr_err2_xi_published_v0m = new Double_t[hHEP_e2_xi->GetEntries()];

  Double_t* arr_err3_k0s_published_v0m = new Double_t[hHEP_e3_k0s->GetEntries()];
  Double_t* arr_err3_xi_published_v0m = new Double_t[hHEP_e3_xi->GetEntries()];

  std::cout << "Published Run 2 (V0M) Yields:" << std::endl;
  std::cout << "Number of data points: " << hHEP_yield_k0s->GetEntries() << std::endl;
  for(Int_t i = 0; i < hHEP_yield_k0s->GetEntries(); i++) {
    arr_yield_k0s_published_v0m[i] = hHEP_yield_k0s->GetBinContent(2*i+1);
    arr_yield_xi_published_v0m[i] = hHEP_yield_xi->GetBinContent(2*i+1);
    std::cout << "  k0s: " << hHEP_yield_k0s->GetBinContent(2*i+1) << " xi: " << hHEP_yield_xi->GetBinContent(2*i+1) << std::endl;

    arr_err1_k0s_published_v0m[i] = hHEP_e1_k0s->GetBinContent(2*i+1);
    arr_err1_xi_published_v0m[i] = hHEP_e1_xi->GetBinContent(2*i+1);

    arr_err2_k0s_published_v0m[i] = hHEP_e2_k0s->GetBinContent(2*i+1);
    arr_err2_xi_published_v0m[i] = hHEP_e2_xi->GetBinContent(2*i+1);

    arr_err3_k0s_published_v0m[i] = hHEP_e3_k0s->GetBinContent(2*i+1);
    arr_err3_xi_published_v0m[i] = hHEP_e3_xi->GetBinContent(2*i+1);
  }

  TGraphErrors *g_k0svsxi_run2_v0m = new TGraphErrors(hHEP_yield_k0s->GetEntries(), arr_yield_k0s_published_v0m, arr_yield_xi_published_v0m, arr_err1_k0s_published_v0m, arr_err1_xi_published_v0m);
  Int_t style_run2_v0m = 8;
  g_k0svsxi_run2_v0m->SetMarkerColor(color[style_run2_v0m]);
  g_k0svsxi_run2_v0m->SetMarkerStyle(MarkerMult[style_run2_v0m]);
  g_k0svsxi_run2_v0m->SetMarkerSize(1.5);

  TGraphErrors *g_k0svsxi_run2_v0m_e2 = new TGraphErrors(hHEP_yield_k0s->GetEntries(), arr_yield_k0s_published_v0m, arr_yield_xi_published_v0m, arr_err2_k0s_published_v0m, arr_err2_xi_published_v0m);
  g_k0svsxi_run2_v0m_e2->SetMarkerColor(color[style_run2_v0m]);
  g_k0svsxi_run2_v0m_e2->SetMarkerStyle(MarkerMult[style_run2_v0m]);
  g_k0svsxi_run2_v0m_e2->SetFillStyle(0);
  g_k0svsxi_run2_v0m_e2->SetMarkerSize(1.5);

  TGraphErrors *g_k0svsxi_run2_v0m_e3 = new TGraphErrors(hHEP_yield_k0s->GetEntries(), arr_yield_k0s_published_v0m, arr_yield_xi_published_v0m, arr_err3_k0s_published_v0m, arr_err3_xi_published_v0m);
  g_k0svsxi_run2_v0m_e3->SetMarkerColor(color[style_run2_v0m]);
  g_k0svsxi_run2_v0m_e3->SetMarkerStyle(MarkerMult[style_run2_v0m]);
  g_k0svsxi_run2_v0m_e3->SetFillColor(color[style_run2_v0m]);
  g_k0svsxi_run2_v0m_e3->SetFillStyle(3001);
  g_k0svsxi_run2_v0m_e3->SetMarkerSize(1.5);

  // Get Published Data (SPDtracklets0815) //
  TString inputHEPDataPath_trck0815 = workingDir + "/data/published/HEPData-ins1748157-v2-Table_8b.root";
  TFile* fileHEPDataIn_trck0815 = TFile::Open(inputHEPDataPath_trck0815);
  if (!fileHEPDataIn_trck0815 || fileHEPDataIn_trck0815->IsZombie()) {
      std::cerr << "Error opening HEP data file!" << std::endl;
      return;
  }

  TH1F* hHEP_yield_k0s_trck0815 = (TH1F *)fileHEPDataIn_trck0815->Get("Table 8b/Hist1D_y1");
  TH1F* hHEP_yield_xi_trck0815 = (TH1F *)fileHEPDataIn_trck0815->Get("Table 8b/Hist1D_y3");

  TH1F* hHEP_e1_k0s_trck0815 = (TH1F *)fileHEPDataIn_trck0815->Get("Table 8b/Hist1D_y1_e1");
  TH1F* hHEP_e1_xi_trck0815 = (TH1F *)fileHEPDataIn_trck0815->Get("Table 8b/Hist1D_y3_e1");

  TH1F* hHEP_e2_k0s_trck0815 = (TH1F *)fileHEPDataIn_trck0815->Get("Table 8b/Hist1D_y1_e2");
  TH1F* hHEP_e2_xi_trck0815 = (TH1F *)fileHEPDataIn_trck0815->Get("Table 8b/Hist1D_y3_e2");

  TH1F* hHEP_e3_k0s_trck0815 = (TH1F *)fileHEPDataIn_trck0815->Get("Table 8b/Hist1D_y1_e3");
  TH1F* hHEP_e3_xi_trck0815 = (TH1F *)fileHEPDataIn_trck0815->Get("Table 8b/Hist1D_y3_e3");

  if (!hHEP_yield_k0s_trck0815 || !hHEP_e1_k0s_trck0815 || !hHEP_yield_xi_trck0815 || !hHEP_e1_xi_trck0815 || !hHEP_e2_k0s_trck0815 || !hHEP_e2_xi_trck0815 || !hHEP_e3_k0s_trck0815 || !hHEP_e3_xi_trck0815)
  {
    std::cerr << "Histogram HEPdata is not found!" << std::endl;
    return;
  }

  Double_t* arr_yield_k0s_published_trck0815 = new Double_t[hHEP_yield_k0s_trck0815->GetEntries()];
  Double_t* arr_yield_xi_published_trck0815 = new Double_t[hHEP_yield_xi_trck0815->GetEntries()];

  Double_t* arr_err1_k0s_published_trck0815 = new Double_t[hHEP_e1_k0s_trck0815->GetEntries()];
  Double_t* arr_err1_xi_published_trck0815 = new Double_t[hHEP_e1_xi_trck0815->GetEntries()];

  Double_t* arr_err2_k0s_published_trck0815 = new Double_t[hHEP_e2_k0s_trck0815->GetEntries()];
  Double_t* arr_err2_xi_published_trck0815 = new Double_t[hHEP_e2_xi_trck0815->GetEntries()];

  Double_t* arr_err3_k0s_published_trck0815 = new Double_t[hHEP_e3_k0s_trck0815->GetEntries()];
  Double_t* arr_err3_xi_published_trck0815 = new Double_t[hHEP_e3_xi_trck0815->GetEntries()];

  std::cout << "Published Run 2 (SPDtracklets0815) Yields:" << std::endl;
  std::cout << "Number of data points: " << hHEP_yield_k0s_trck0815->GetEntries() << std::endl;
  for(Int_t i = 0; i < hHEP_yield_k0s_trck0815->GetEntries(); i++) {
    arr_yield_k0s_published_trck0815[i] = hHEP_yield_k0s_trck0815->GetBinContent(2*i+1);
    arr_yield_xi_published_trck0815[i] = hHEP_yield_xi_trck0815->GetBinContent(2*i+1);
    std::cout << "  k0s: " << hHEP_yield_k0s_trck0815->GetBinContent(2*i+1) << " xi: " << hHEP_yield_xi_trck0815->GetBinContent(2*i+1) << std::endl;

    arr_err1_k0s_published_trck0815[i] = hHEP_e1_k0s_trck0815->GetBinContent(2*i+1);
    arr_err1_xi_published_trck0815[i] = hHEP_e1_xi_trck0815->GetBinContent(2*i+1);

    arr_err2_k0s_published_trck0815[i] = hHEP_e2_k0s_trck0815->GetBinContent(2*i+1);
    arr_err2_xi_published_trck0815[i] = hHEP_e2_xi_trck0815->GetBinContent(2*i+1);

    arr_err3_k0s_published_trck0815[i] = hHEP_e3_k0s_trck0815->GetBinContent(2*i+1);
    arr_err3_xi_published_trck0815[i] = hHEP_e3_xi_trck0815->GetBinContent(2*i+1);
  }

  TGraphErrors *g_k0svsxi_run2_trck0815 = new TGraphErrors(hHEP_yield_k0s_trck0815->GetEntries(), arr_yield_k0s_published_trck0815, arr_yield_xi_published_trck0815, arr_err1_k0s_published_trck0815, arr_err1_xi_published_trck0815);
  Int_t style_run2_trck0815 = 6;
  g_k0svsxi_run2_trck0815->SetMarkerColor(color[style_run2_trck0815]);
  g_k0svsxi_run2_trck0815->SetMarkerStyle(MarkerMult[style_run2_trck0815]);
  g_k0svsxi_run2_trck0815->SetMarkerSize(1.5);

  TGraphErrors *g_k0svsxi_run2_trck0815_e2 = new TGraphErrors(hHEP_yield_k0s_trck0815->GetEntries(), arr_yield_k0s_published_trck0815, arr_yield_xi_published_trck0815, arr_err2_k0s_published_trck0815, arr_err2_xi_published_trck0815);
  g_k0svsxi_run2_trck0815_e2->SetMarkerColor(color[style_run2_trck0815]);
  g_k0svsxi_run2_trck0815_e2->SetMarkerStyle(MarkerMult[style_run2_trck0815]);
  g_k0svsxi_run2_trck0815_e2->SetFillStyle(0);
  g_k0svsxi_run2_trck0815_e2->SetMarkerSize(1.5);

  TGraphErrors *g_k0svsxi_run2_trck0815_e3 = new TGraphErrors(hHEP_yield_k0s_trck0815->GetEntries(), arr_yield_k0s_published_trck0815, arr_yield_xi_published_trck0815, arr_err3_k0s_published_trck0815, arr_err3_xi_published_trck0815);
  g_k0svsxi_run2_trck0815_e3->SetMarkerColor(color[style_run2_trck0815]);
  g_k0svsxi_run2_trck0815_e3->SetMarkerStyle(MarkerMult[style_run2_trck0815]);
  g_k0svsxi_run2_trck0815_e3->SetFillColor(color[style_run2_trck0815]);
  g_k0svsxi_run2_trck0815_e3->SetFillStyle(3001);
  g_k0svsxi_run2_trck0815_e3->SetMarkerSize(1.5);

  // Get Published Data (SPDtracklets08) //
  TString inputHEPDataPath_trck08 = workingDir + "/data/published/HEPData-ins1748157-v2-Table_8a.root";
  TFile* fileHEPDataIn_trck08 = TFile::Open(inputHEPDataPath_trck08);
  if (!fileHEPDataIn_trck08 || fileHEPDataIn_trck08->IsZombie()) {
      std::cerr << "Error opening HEP data file!" << std::endl;
      return;
  }

  TH1F* hHEP_yield_k0s_trck08 = (TH1F *)fileHEPDataIn_trck08->Get("Table 8a/Hist1D_y1");
  TH1F* hHEP_yield_xi_trck08 = (TH1F *)fileHEPDataIn_trck08->Get("Table 8a/Hist1D_y3");

  TH1F* hHEP_e1_k0s_trck08 = (TH1F *)fileHEPDataIn_trck08->Get("Table 8a/Hist1D_y1_e1");
  TH1F* hHEP_e1_xi_trck08 = (TH1F *)fileHEPDataIn_trck08->Get("Table 8a/Hist1D_y3_e1");

  TH1F* hHEP_e2_k0s_trck08 = (TH1F *)fileHEPDataIn_trck08->Get("Table 8a/Hist1D_y1_e2");
  TH1F* hHEP_e2_xi_trck08 = (TH1F *)fileHEPDataIn_trck08->Get("Table 8a/Hist1D_y3_e2");

  TH1F* hHEP_e3_k0s_trck08 = (TH1F *)fileHEPDataIn_trck08->Get("Table 8a/Hist1D_y1_e3");
  TH1F* hHEP_e3_xi_trck08 = (TH1F *)fileHEPDataIn_trck08->Get("Table 8a/Hist1D_y3_e3");

  if (!hHEP_yield_k0s_trck08 || !hHEP_e1_k0s_trck08 || !hHEP_yield_xi_trck08 || !hHEP_e1_xi_trck08 || !hHEP_e2_k0s_trck08 || !hHEP_e2_xi_trck08 || !hHEP_e3_k0s_trck08 || !hHEP_e3_xi_trck08)
  {
    std::cerr << "Histogram HEPdata is not found!" << std::endl;
    return;
  }

  Double_t* arr_yield_k0s_published_trck08 = new Double_t[hHEP_yield_k0s_trck08->GetEntries()];
  Double_t* arr_yield_xi_published_trck08 = new Double_t[hHEP_yield_xi_trck08->GetEntries()];

  Double_t* arr_err1_k0s_published_trck08 = new Double_t[hHEP_e1_k0s_trck08->GetEntries()];
  Double_t* arr_err1_xi_published_trck08 = new Double_t[hHEP_e1_xi_trck08->GetEntries()];

  Double_t* arr_err2_k0s_published_trck08 = new Double_t[hHEP_e2_k0s_trck08->GetEntries()];
  Double_t* arr_err2_xi_published_trck08 = new Double_t[hHEP_e2_xi_trck08->GetEntries()];

  Double_t* arr_err3_k0s_published_trck08 = new Double_t[hHEP_e3_k0s_trck08->GetEntries()];
  Double_t* arr_err3_xi_published_trck08 = new Double_t[hHEP_e3_xi_trck08->GetEntries()];

  std::cout << "Published Run 2 (SPDtracklets08) Yields:" << std::endl;
  std::cout << "Number of data points: " << hHEP_yield_k0s_trck08->GetEntries() << std::endl;
  for(Int_t i = 0; i < hHEP_yield_k0s_trck08->GetEntries(); i++) {
    arr_yield_k0s_published_trck08[i] = hHEP_yield_k0s_trck08->GetBinContent(2*i+1);
    arr_yield_xi_published_trck08[i] = hHEP_yield_xi_trck08->GetBinContent(2*i+1);
    std::cout << "  k0s: " << hHEP_yield_k0s_trck08->GetBinContent(2*i+1) << " xi: " << hHEP_yield_xi_trck08->GetBinContent(2*i+1) << std::endl;

    arr_err1_k0s_published_trck08[i] = hHEP_e1_k0s_trck08->GetBinContent(2*i+1);
    arr_err1_xi_published_trck08[i] = hHEP_e1_xi_trck08->GetBinContent(2*i+1);

    arr_err2_k0s_published_trck08[i] = hHEP_e2_k0s_trck08->GetBinContent(2*i+1);
    arr_err2_xi_published_trck08[i] = hHEP_e2_xi_trck08->GetBinContent(2*i+1);

    arr_err3_k0s_published_trck08[i] = hHEP_e3_k0s_trck08->GetBinContent(2*i+1);
    arr_err3_xi_published_trck08[i] = hHEP_e3_xi_trck08->GetBinContent(2*i+1);
  }

  TGraphErrors *g_k0svsxi_run2_trck08 = new TGraphErrors(hHEP_yield_k0s_trck08->GetEntries(), arr_yield_k0s_published_trck08, arr_yield_xi_published_trck08, arr_err1_k0s_published_trck08, arr_err1_xi_published_trck08);
  Int_t style_run2_trck08 = 3;
  g_k0svsxi_run2_trck08->SetMarkerColor(color[style_run2_trck08]);
  g_k0svsxi_run2_trck08->SetMarkerStyle(MarkerMult[style_run2_trck08]);
  g_k0svsxi_run2_trck08->SetMarkerSize(1.5);

  TGraphErrors *g_k0svsxi_run2_trck08_e2 = new TGraphErrors(hHEP_yield_k0s_trck08->GetEntries(), arr_yield_k0s_published_trck08, arr_yield_xi_published_trck08, arr_err2_k0s_published_trck08, arr_err2_xi_published_trck08);
  g_k0svsxi_run2_trck08_e2->SetMarkerColor(color[style_run2_trck08]);
  g_k0svsxi_run2_trck08_e2->SetMarkerStyle(MarkerMult[style_run2_trck08]);
  g_k0svsxi_run2_trck08_e2->SetFillStyle(0);
  g_k0svsxi_run2_trck08_e2->SetMarkerSize(1.5);

  TGraphErrors *g_k0svsxi_run2_trck08_e3 = new TGraphErrors(hHEP_yield_k0s_trck08->GetEntries(), arr_yield_k0s_published_trck08, arr_yield_xi_published_trck08, arr_err3_k0s_published_trck08, arr_err3_xi_published_trck08);
  g_k0svsxi_run2_trck08_e3->SetMarkerColor(color[style_run2_trck08]);
  g_k0svsxi_run2_trck08_e3->SetMarkerStyle(MarkerMult[style_run2_trck08]);
  g_k0svsxi_run2_trck08_e3->SetFillColor(color[style_run2_trck08]);
  g_k0svsxi_run2_trck08_e3->SetFillStyle(3001);
  g_k0svsxi_run2_trck08_e3->SetMarkerSize(1.5);

  // Plot k0svsxi //
  // Integrated yield vs dN/dEta //
  TLegend *legendTitle_k0svsxi = new TLegend(0.143, 0.670, 0.341, 0.849);
  StyleLegend(legendTitle_k0svsxi, 0.0, 0.0);
  legendTitle_k0svsxi->SetTextAlign(11);
  legendTitle_k0svsxi->SetTextSize(0.04);
  legendTitle_k0svsxi->SetTextFont(42);
  legendTitle_k0svsxi->SetLineColorAlpha(0.,0.);
  legendTitle_k0svsxi->SetFillColorAlpha(0.,0.);
  legendTitle_k0svsxi->AddEntry("", "#bf{ALICE Work In Progress}", "");
  legendTitle_k0svsxi->AddEntry("", "pp, #sqrt{#it{s}} = 13.6 TeV", "");
  legendTitle_k0svsxi->AddEntry("", "|y| < 0.5", "");

  TLegend *legend_k0svsxi = new TLegend(0.191, 0.461, 0.392, 0.639);
  legend_k0svsxi->SetFillStyle(0);
  legend_k0svsxi->SetTextSize(0.03);

  TLegend * legUncert = new TLegend(0.65, 0.225, 0.960, 0.316);
  StyleLegend(legUncert, 0.0, 0.0);
  legUncert->SetFillColor(0);
  legUncert->SetTextSize(0.04);

  TCanvas *canvas_k0svsxi = new TCanvas("canvas_k0svsxi","canvas_k0svsxi", 0,66,857,708);
  canvas_k0svsxi->cd();
  canvas_k0svsxi->Draw();
  StyleCanvas(canvas_k0svsxi, 0.14, 0.05, 0.14, 0.15);

  TMultiGraph *mg_k0svsxi = new TMultiGraph("k0svsxi", "");

  mg_k0svsxi->Add(g_k0svsxi, "P");

  mg_k0svsxi->Add(g_k0svsxi_run2_v0m, "P");
  mg_k0svsxi->Add(g_k0svsxi_run2_v0m_e2, "P 5");
  mg_k0svsxi->Add(g_k0svsxi_run2_v0m_e3, "P 5");

  mg_k0svsxi->Add(g_k0svsxi_run2_trck0815, "P");
  mg_k0svsxi->Add(g_k0svsxi_run2_trck0815_e2, "P 5");
  mg_k0svsxi->Add(g_k0svsxi_run2_trck0815_e3, "P 5");

  mg_k0svsxi->Add(g_k0svsxi_run2_trck08, "P");
  mg_k0svsxi->Add(g_k0svsxi_run2_trck08_e2, "P 5");
  mg_k0svsxi->Add(g_k0svsxi_run2_trck08_e3, "P 5");

  mg_k0svsxi->Draw("A");
  // Update the canvas to reflect changes
  canvas_k0svsxi->Update();
  gPad->Modified();
  StyleMultGraph(mg_k0svsxi, 0, 0.16, "K_{S}^{0}", particleSymnbols[2], "", 1, 0, 2., 1.15, 1.25, 0.05, 0.05);

  legendTitle_k0svsxi->Draw("SAME");

  legend_k0svsxi->AddEntry(g_k0svsxi, "FT0M Classes (stat. only)", "epl");
  legend_k0svsxi->AddEntry(g_k0svsxi_run2_v0m, "V0M Classes (Run 2)", "epl");
  legend_k0svsxi->AddEntry(g_k0svsxi_run2_trck0815, "#it{N}_{ tracklets}^{0.8<|#it{#eta}|<1.5} (Run 2)", "epl");
  legend_k0svsxi->AddEntry(g_k0svsxi_run2_trck08, "#it{N}_{ tracklets}^{|#it{#eta}|<0.8} (Run 2)", "epl");
  StyleLegend(legend_k0svsxi, 0.0, 0.0);
  legend_k0svsxi->Draw("SAME");

  legUncert->AddEntry(g_k0svsxi, "stat.", "LE");
  legUncert->AddEntry(g_k0svsxi_run2_v0m_e2, "syst.", "F");
  // Dummy histogram to display syst. ucorr. uncertainty gray box
  TGraphAsymmErrors *gDummy = new TGraphAsymmErrors(1);
  gDummy->SetPoint(0, -1, -1);
  gDummy->SetMarkerColor(0);
  gDummy->SetMarkerStyle(MarkerMult[0]);
  gDummy->SetFillStyle(3001);
  gDummy->SetFillColor(color[0]);
  gDummy->SetMarkerSize(1.5);
  mg_k0svsxi->Add(gDummy, "P 5");
  legUncert->AddEntry(gDummy, "syst. uncorr.", "F");
  legUncert->Draw("same");

  // End of Code
  std::cout << "\x1B[1;32m"; // Set text color to green
  std::cout << "\n************* Integrated k0s vs xi Yields are Plotted Successfully! *************\n";
  std::cout << "\x1B[0m"; // Reset text color
}
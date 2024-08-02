void rootlogon() {
  const Int_t numMult = 10; // number of multiplicity sub-intervals
// Color palette
  const int NRGBs = 5;
  double stops[NRGBs] = { 0.00, 0.24, 0.61, 0.84, 1.00 };
  // double stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  // From blue to red 
  double red[NRGBs]   = { 0.00, 0.00, 0.78, 1.00, 0.51 };
  double green[NRGBs] = { 0.00, 0.81, 0.90, 0.20, 0.00 };
  double blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  int FI = TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,numMult);
}
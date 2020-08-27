#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../utils/Estimator.hh"
#include "../utils/HistGetter.hh"
#include "../SUSYNano19/SRParameters.hh"
#include "../utils/json.hpp"

#include <TRandom.h>
#include <TH2Poly.h>
#include "TPaletteAxis.h"
#include <math.h> 
#include <map>

using namespace std;
using namespace EstTools;
using json = nlohmann::json;

TString outputDir = "";
TString baseDir = "moduleTolerances";
vector<Color_t> colors;


void print2DPlots(TH2Poly *hc, TString geometry, TString BinLatex = "", TString name = "testHoneycomb", double width = 0.2);
void moduleTolerances();
#endif

void ModuleStudies(){
  moduleTolerances();
}

double GetWidthToA(double width){
  //width/2 == a/2*tan(pi/6)
  return width*TMath::Tan(TMath::Pi()/6);
}

double GetAToWidth(double a){
  return a/TMath::Tan(TMath::Pi()/6);
}

double GetSideFive(double width){
  //width/2 == hyp/TMath::Sqrt(3)/2
  //
  //double hyp = (width/2)*TMath::Sqrt(3)/2;
  return GetWidthToA(width);
}

vector<double> GetPolar(double x1, double y1, double x0, double y0){
  vector<double> polar;
  double r = TMath::Sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
  polar.push_back(r);
  double phi = TMath::ATan((y1 - y0)/(x1 - x0));
  //cout << "(x, y) = " << "(" << (x1 - x0) << ", " << (y1 - y0) << ")" << endl;
  //cout << "(r, phi) = " << "(" << r << ", " << phi << endl;
  if((x1 - x0) < 0) phi += TMath::Pi();
  if((x1 - x0) >= 0 && (y1 - y0) < 0) phi += 2*TMath::Pi();
  //cout << "(r, phi) = " << "(" << r << ", " << phi << endl;

  polar.push_back(phi);
  return polar;
}

vector<double> GetCartesian(double r, double phi, double x0, double y0){
  vector<double> cart;
  double x1 = r*TMath::Cos(phi) + x0;
  double y1 = r*TMath::Sin(phi) + y0;
  //cout << "(r, phi) = " << "(" << r << ", " << phi << endl;
  //cout << "(x, y) = " << "(" << x1 << ", " << y1 << ")" << endl;
  
  cart.push_back(x1);
  cart.push_back(y1);
  return cart;
}

double integrateFrom(TH1D* h, double value = 0.){
  int xbin = h->FindFixBin(value);
  double prob = h->Integral(0,xbin+1);
  return prob;
}

vector<double> findHighXbin(TH1D* h){
  double mod = -999., belowZero = 0., width = 0.;;
  vector<double> output;
  for(auto bin = h->GetNbinsX() + 1; bin > 0; --bin){
    mod = h->GetBinContent(bin) > 0 ? h->GetBinLowEdge(bin) : mod;
    if(mod != -999.){
      width = h->GetBinWidth(bin);
    }
    if(h->GetBinLowEdge(bin) < 0){
      belowZero += h->GetBinContent(bin);
    }
  }
  output.push_back((mod + width/2)*1000);
  output.push_back(width*1000);
  output.push_back(belowZero);
  return output;
}

vector<TH1*> IntegrateHex(TH2Poly* h, TString geo, TString type, pair < double, double > center, double width, double axis, double binSize, int nPeaks){
  int nbins = (axis - 0.1)/binSize*4;
  vector<TH1*> hist;
  for(auto iP = 0; iP != nPeaks + 1; iP++){ 
    TH1* hist_1D = nullptr;
    if(TString(h->GetName()).Contains("Gold Bond") || TString(h->GetName()).Contains("Guard Bond"))
      hist_1D = new TH1D(TString(h->GetName()) + "_Peak" + to_string(iP), ";" + type + " Diff Width [mm]; Modules", 750, 0.0, 1.5);
    else
      hist_1D = new TH1D(TString(h->GetName()) + "_Peak" + to_string(iP), ";" + type + " Diff Width [mm]; Modules", 425, -0.15, 0.7);
    hist.push_back(hist_1D);
  }

  vector<double> angle = {0., TMath::Pi()/3, TMath::Pi()*2/3, TMath::Pi(), TMath::Pi()*4/3, TMath::Pi()*5/3};
  if(geo == "Five" || geo == "Semi") angle = {0., TMath::Pi()/2, TMath::Pi(), TMath::Pi()*4/3, TMath::Pi()*5/3};
  for(auto iA = 0; iA != angle.size(); iA++){
    int index = 1;
         if(nPeaks > 1 && iA == 1) index = 2;
    else if(nPeaks > 1 && (iA == 3 || iA == 4)) index = 3;
    for(auto iB = 0; iB != nbins; iB++){
      double x_ = (binSize*iB + (width/2 - 300*binSize))*TMath::Cos(angle[iA]);
      double y_ = (binSize*iB + (width/2 - 300*binSize))*TMath::Sin(angle[iA]);
      int bin = h->FindBin(x_, y_);
      double binValue = h->GetBinContent(bin);

      //get distance to the edge
      double r = TMath::Sqrt(x_*x_ + y_*y_) - width/2 + 0.0009;
      double dr_center = center.first*TMath::Cos(angle[iA]) + center.second*TMath::Sin(angle[iA]);
      if(geo == "Five" && iA == 1) r += width/2 - GetWidthToA(width)/2;
      if(geo == "Semi" && iA == 1) r += width/2;
      r += dr_center;

      double bin1d = hist[index]->GetXaxis()->FindBin(r);
      double tot = hist[index]->GetBinContent(bin1d);
      hist[index]->SetBinContent(bin1d, tot + binValue);
      //fill plot with total
      bin1d = hist[0]->GetXaxis()->FindBin(r);
      tot = hist[0]->GetBinContent(bin1d);
      hist[0]->SetBinContent(bin1d, tot + binValue);
    }
  }

  return hist;
}

vector < pair <double, pair <double, double> > > newShapeInput(){
  double baseplate_w = 166.94;
  vector<double> baseplateInput = {6.572557, 6.573377, 6.573227, 6.568969, 6.56959, 6.569731, 
				  6.568373, 6.569033, 6.56842, 6.572988, 6.572592, 6.57396, 
				  6.573593, 6.573748, 6.573184, 6.572817, 6.572886, 6.573836, 
				  6.573114, 6.573073, 6.573021, 6.568787, 6.56881, 6.569536, 
				  6.569422, 6.569377, 6.569091, 6.573068, 6.573528, 6.573849, 
				  6.573138, 6.573199, 6.573384, 6.572994, 6.573926, 6.573588, 
				  6.568979, 6.569386, 6.569742};

  vector<double> kaptonInput = {6.573537, 6.571026, 6.569803, 6.569896, 6.572914, 6.570429,
         			6.571285, 6.569886, 6.572366, 6.570863, 6.569884, 6.57224,
                  		6.573037, 6.571025, 6.569289, 6.571091, 6.569678, 6.572536,
 				6.573057, 6.570624, 6.569357, 6.573133, 6.570727, 6.569481, 
   				6.57096,  6.569481, 6.572649, 6.571264, 6.569811, 6.572641
                               };

  vector<double> convert;
  TH1* baseplate_meas_hist = new TH1D("Baseplate", ";Measured Width [mm]; Modules", 50, 166.76, 167.09);
  TH1* kapton_meas_hist = new TH1D("Kapton", ";Measured Width [mm]; Modules", 50, 166.76, 167.09);
  for(auto c: baseplateInput){
    baseplate_meas_hist->Fill(c*25.4);
  }
  for(auto c: kaptonInput){
    kapton_meas_hist->Fill(c*25.4);
  }

  baseplate_meas_hist->Fit("gausn");
  TF1 *baseplate_meas_avg_param = (TF1*)baseplate_meas_hist->GetFunction("gausn");
  kapton_meas_hist->Fit("gausn");
  TF1 *kapton_meas_avg_param = (TF1*)kapton_meas_hist->GetFunction("gausn");

  prepHists({baseplate_meas_hist, kapton_meas_hist}, false, false, false, {kRed, kBlue});
  auto leg = prepLegends({baseplate_meas_hist, kapton_meas_hist}, {baseplate_meas_hist->GetName(), kapton_meas_hist->GetName()}, "L");
  leg->SetTextSize(0.03);
  leg->SetY1NDC(leg->GetY2NDC() - 0.3);
  TCanvas* c = drawCompMatt({baseplate_meas_hist, kapton_meas_hist}, leg);
  gStyle->SetOptStat(0);
  TString name = "measured_dist";
  c->SetTitle(name);
  c->Print(outputDir+"/"+name+".pdf");

  vector < pair <double, pair <double, double> > > avg_param;
  avg_param.push_back(make_pair(baseplate_meas_avg_param->GetParameter(0), make_pair(baseplate_meas_hist->GetMean(), baseplate_meas_hist->GetStdDev())));
  avg_param.push_back(make_pair(kapton_meas_avg_param->GetParameter(0), make_pair(kapton_meas_hist->GetMean(), kapton_meas_hist->GetStdDev())));

  return avg_param;
}

json makeJSONModule(vector< pair< pair< string, string>, pair< double, double > > >& fit, vector< pair< pair< string, string>, pair< double, double > > >& worst){
  json j;
  int binnum = 0;
  for ( const auto &it : fit ) {
    j["Fit"][it.first.first][it.first.second] = {it.second.first, it.second.second};
    if((TString(it.first.first).Contains("_PCB") || TString(it.first.first).Contains("Nominal")) && TString(it.first.second).Contains("Gold")){
      j["PCB_Dist"][it.first.first] = it.first.first;
    }
    if(TString(it.first.first).Contains("Gaussian") && j["BinNum"].find(it.first.second) == j["BinNum"].end()){
      j["BinNum"][it.first.second] = binnum;
      binnum++;
    }
    if((TString(it.first.first).Contains("Gaussian") || TString(it.first.first).Contains("Nominal")) && TString(it.first.second).Contains("Gold")){
      j["Sensor"][it.first.first] = it.first.first;
    }
    if (!TString(it.first.first).Contains("_PCB") && !TString(it.first.first).Contains("Custom") && TString(it.first.second).Contains("Gold")){
      j["Distribution"][it.first.first] = it.first.first;
    }
    if ((TString(it.first.first).Contains("Custom") || TString(it.first.first).Contains("Nominal")) && TString(it.first.second).Contains("Gold")){
      j["Custom"][it.first.first] = it.first.first;
    }
  }
  for ( const auto & it : worst ) {
    j["Worst"][it.first.first][it.first.second] = {it.second.first, it.second.second};
  }

  return j;
}

double getRandomValue(double width, double error, TString type){
  double rand = 0.;
  if(type.Contains("Gaus"))
    rand = gRandom->Gaus(width, error/3);
  else if(type.Contains("Land"))
    rand = gRandom->Landau(width, error/3);
  else if(type.Contains("Flat"))
    rand = gRandom->Uniform(width - error, width + error);

  return rand;
}

pair < double, double >  getPolyCenter(vector<pair<double, double> > p){
  double xf = 0., yf = 0., div = 0.;
  for(auto i = 0; i < p.size(); i++){
    int n0 = i, n1 = (i + 1) % p.size();
    xf += (p[n0].first + p[n1].first)*(p[n0].first*p[n1].second - p[n1].first*p[n0].second);
    yf += (p[n0].second + p[n1].second)*(p[n0].first*p[n1].second - p[n1].first*p[n0].second);
    div += (p[n0].first*p[n1].second - p[n1].first*p[n0].second);
  }

  return {xf/(3*div), yf/(3*div)};
}

vector<pair<double, double> > GetPoints(TString geo, double a_in, double width){
  double x[6], y[6];
  if(geo == "Full"){
    x[0] = -1*width/2;
    y[0] = -1*a_in/2;
    x[1] = x[0];
    y[1] = y[0] + a_in;
    x[2] = x[1] + a_in*TMath::Sqrt(3)/2.0;
    y[2] = y[1] + a_in/2.0;
    x[3] = x[2] + a_in*TMath::Sqrt(3)/2.0;
    y[3] = y[1];
    x[4] = x[3];
    y[4] = y[0];
    x[5] = x[2];
    y[5] = y[4] - a_in/2.0;
  } else if(geo == "Five"){
    x[0] = -1*width/2;
    y[0] = -1*a_in/2;
    x[1] = x[0];
    y[1] = y[0] + a_in;
    x[2] = x[1] + a_in*TMath::Sqrt(3);
    y[2] = y[1];
    x[3] = x[2];
    y[3] = y[0];
    x[4] = x[1] + a_in*TMath::Sqrt(3)/2.0;
    y[4] = y[3] - a_in/2.0;
  } else if(geo == "Semi"){
    x[0] = -1*width/2;
    y[0] = -1*a_in/2;
    x[1] = x[0];
    y[1] = y[0] + a_in/2;
    x[2] = x[1] + a_in*TMath::Sqrt(3);
    y[2] = y[1];
    x[3] = x[2];
    y[3] = y[0];
    x[4] = x[1] + a_in*TMath::Sqrt(3)/2.0;
    y[4] = y[3] - a_in/2.0;
  } else if(geo == "Half"){
    x[0] = -1*width/2;
    y[0] = -1*a_in/2;
    x[1] = x[0] + a_in*TMath::Sqrt(3);
    y[1] = y[0] + a_in;
    x[2] = x[1];
    y[2] = y[0];
    x[3] = x[0] + a_in*TMath::Sqrt(3)/2.0;
    y[3] = y[2] - a_in/2.0;
  }

  vector<pair<double, double> > points;
  int nPoints = 6;
  if(geo == "Five" or geo == "Semi") nPoints = 5;
  else if(geo == "Half") nPoints = 4;
  for(auto i = 0; i < nPoints; i++){
    points.push_back(make_pair(x[i], y[i]));
  }

  return points;
}

pair < double, double > HoneycombCustom(TString geo, TH2Poly* hc, Double_t a, double width, double a_in = 96.67, int steps = 10){
  // Add the bins
  vector<pair<double, double> >  points = GetPoints(geo, a_in, width);
  pair < double, double > center = getPolyCenter(points);
  Double_t xloop[4], yloop[4];
  //Radius of points
  vector<double> radius, angle, ratio;
  float rMin = 9999999.;
  for(auto i = 0; i < points.size(); i++){
    points.at(i).first = points.at(i).first - center.first;
    points.at(i).second = points.at(i).second - center.second;
    vector<double> polar = GetPolar(points.at(i).first, points.at(i).second, 0., 0.);
    radius.push_back(polar[0]);
    angle.push_back(polar[1]);
    if(polar[0] < rMin) rMin = polar[0];
  }
  for(auto i = 0; i < points.size(); i++) ratio.push_back(radius[i]/rMin);
  for (auto theta = 0; theta != angle.size(); theta++){
    int first = theta % angle.size();
    int second = (theta + 1) % angle.size();
    for(int neg = 0; neg < 2; neg++){
      xloop[0] = points.at(first).first;
      yloop[0] = points.at(first).second;
      xloop[1] = points.at(second).first;
      yloop[1] = points.at(second).second;
      double sign = neg ? -1. : 1.;
      int count = 0;
      for(auto ibin = 1; ibin < steps; ibin++){
        vector<double> cart = GetCartesian(radius[first] + sign*a*ibin*ratio[first], angle[first], 0., 0.);
        xloop[3] = cart[0];
        yloop[3] = cart[1];
        cart = GetCartesian(radius[second] + sign*a*ibin*ratio[second], angle[second], 0., 0.);
        xloop[2] = cart[0];
        yloop[2] = cart[1];
        if(count == 300 && neg == 1) break;
        if(geo == "Semi" && count == 200 && neg == 1) break;
        hc->AddBin(4, xloop, yloop);
        //if(theta == 0 && neg == 0){
        //  cout << theta << endl;
        //  cout << "(" << xloop[0] << ", " << yloop[0] << ")" << endl;
        //  cout << "(" << xloop[1] << ", " << yloop[1] << ")" << endl;
        //  cout << "(" << xloop[2] << ", " << yloop[2] << ")" << endl;
        //  cout << "(" << xloop[3] << ", " << yloop[3] << ")" << endl;
        //}
        xloop[0] = xloop[3];
        yloop[0] = yloop[3];
        xloop[1] = xloop[2];
        yloop[1] = yloop[2];
        count++;
      }
    }
  }
  return center;
}

void FillAllSides(TString geo, double max, pair < double, double > center, TH2Poly* hc, vector<double> rand, double nom, double shift_x, double shift_y, double width_new, vector<double> comp = {}, double comp_shift_x = 0., double comp_shift_y = 0.){
  //Widths
  double r_out[3];
  double secondary = nom;
  double dx = shift_x - comp_shift_x, dy = shift_y - comp_shift_y;
  double dx_cen = 0., dy_cen = 0.;
  int iR = 0;
  vector<double> forward = {0., TMath::Pi()/3, TMath::Pi()*2/3, TMath::Pi(), TMath::Pi()*4/3, TMath::Pi()*5/3};
  if(geo == "Five" || geo == "Semi") forward = {0., TMath::Pi()/2, TMath::Pi(), TMath::Pi()*4/3, TMath::Pi()*5/3};
  for(auto iF = 0; iF != forward.size(); iF++){
    if(forward[iF] == TMath::Pi()) iR = 0;
    if(comp.size() > 0) secondary = comp[iR];

    r_out[iR] = (rand[iR] - secondary + width_new)/2;
    double dr = dx*TMath::Cos(forward[iF]) + dy*TMath::Sin(forward[iF]);
    double dr_cen = TMath::Abs(center.first)*TMath::Cos(forward[iF]) + TMath::Abs(center.second)*TMath::Sin(forward[iF]);

    if((geo == "Five" || geo == "Semi") && (iF == 1)){ 
      r_out[iR] = (GetWidthToA(rand[iR])/2 - GetWidthToA(secondary)/2 + GetWidthToA(width_new)/2);
      //cout << TString(hc->GetName()) 
      //     << " a_new: " << GetWidthToA(rand[iR])/2 
      //     << ", sub: " << GetWidthToA(secondary)/2 
      //     << ", nom: " << GetWidthToA(width_new)/2 
      //     << ", out: " << (r_out[iR] + dr - dr_cen)*TMath::Sin(forward[iF]) 
      //     << ", dr: " << dr 
      //     << ", dr_Center: " << dr_cen << endl;
    }
    hc->Fill((r_out[iR] + dr + dr_cen)*TMath::Cos(forward[iF]), (r_out[iR] + dr + dr_cen)*TMath::Sin(forward[iF]));
    iR++;
  }
}

void AddLineNomHex(double width, TString geo = "Full"){
   double a_in = GetWidthToA(width);
   // Add the bins
  vector<pair<double, double> >  points = GetPoints(geo, a_in, width);
  pair < double, double > center = getPolyCenter(points);
  for(auto i = 0; i <  points.size(); i++){
    points.at(i).first = points.at(i).first - center.first;
    points.at(i).second = points.at(i).second - center.second;
  }

  for(auto iS = 0; iS < points.size(); iS++){
    int first = iS % points.size();
    int second = (iS + 1) % points.size();
    TLine *line = new TLine(points.at(first).first, points.at(first).second, points.at(second).first, points.at(second).second);
    line->SetLineColor(kBlack);
    line->Draw();
  }
}

void moduleTolerances(){
  lumistr = "";
  // Input Parameters units are Millimeters
  double baseplate_center = 3.05, baseplate_err = 0.05;
  double pcb_w = 166.64, sensor_w = 166.57, kapton_w = 166.94, baseplate_w = 166.94;
  double pcb_w_err = 0.1, sensor_w_err = 0.02, kapton_w_err = 0.05, baseplate_w_err = 0.05;
  double pcb_shift_r = 0.04, pcb_shift_theta = TMath::Pi();
  double sensor_shift_x = 0.500, sensor_shift_y = 0.050, sensor_shift_theta = 0.050; //theta is radians
  double baseValue = sensor_w;
  double baseValue_total = baseValue;
  double guard_ring_min = 0.25, guard_ring_size = 1.500;
  double gold_min = 0.25, gold_size = 1.500;
  int max = 100000;
  max = 30000;

  int nbins = 550;
  //nbins = 5;
  double step = 0.002;
  double width_new = 1.4;
  double axis = width_new/2 + nbins*step + 0.1;
  double a_new = GetWidthToA(width_new);
  vector < pair <double, pair<double, double> > > custom = newShapeInput();
  double meas_sigma = custom[0].second.second;
  double meas_mean = custom[0].second.first;
  double kapton_meas_sigma = custom[1].second.second;
  double kapton_meas_mean = custom[1].second.first;

  double mean = meas_mean/baseplate_w;
  double sigma = meas_sigma/baseplate_w_err;
  double kapton_mean = kapton_meas_mean/kapton_w;
  double kapton_sigma = kapton_meas_sigma/kapton_w_err;

  cout << "baseplate: " << mean << " +/- " << sigma << endl;
  cout << "kapton: " << kapton_mean << " +/- " << kapton_sigma << endl;

  //Colors
  colors.push_back(kAzure+6);
  colors.push_back(kMagenta-10);
  colors.push_back(kRed-9);
  colors.push_back(kGreen+3);
  colors.push_back(kRed-7);
  colors.push_back(kViolet-4);
  colors.push_back(kSpring-9);

  map<TString, pair<double, double> > center;
  //vector<string> Dist = {"Gaussian_PCBplus25", "Gaussian_newSensor", "Gaussian_PCBplus25_newSensor"};
  vector<string> Dist = {"Gaussian"};
  //vector<string> Geometry = {"Full", "Five", "Semi", "Half"};
  vector<string> Geometry = {"Five"};
  //vector<string> Dist = {"Gaussian", "Landau", "Flat", "CustomGaus", "CustomLandau", "CustomFlat", "Gaussian_PCBplus25", "Gaussian_PCBplus50", "Gaussian_PCBplus75", "Gaussian_PCBminus25", "Gaussian_PCBminus50", "Gaussian_PCBminus75"};
  for(auto &geo_str : Geometry){
    vector< pair< pair< string, string>, pair< double, double > > > worst_values;
    vector< pair< pair< string, string>, pair< double, double > > > fit_values;
    TString geo = TString(geo_str);
    outputDir = baseDir + "/Test" + geo;
    gSystem->mkdir(outputDir, true);

    vector < pair < double, double > > points = GetPoints(geo, GetWidthToA(baseplate_w - baseValue + width_new), baseplate_w - baseValue + width_new);
    center.emplace("baseplate", getPolyCenter(points));
    points = GetPoints(geo, GetWidthToA(kapton_w - baseValue + width_new), kapton_w - baseValue + width_new);
    center.emplace("kapton", getPolyCenter(points));
    points = GetPoints(geo, GetWidthToA(sensor_w - baseValue + width_new), sensor_w - baseValue + width_new);
    center.emplace("sensor", getPolyCenter(points));
    points = GetPoints(geo, GetWidthToA(pcb_w - baseValue + width_new), pcb_w - baseValue + width_new);
    center.emplace("pcb", getPolyCenter(points));
    points = GetPoints(geo, GetWidthToA(pcb_w - sensor_w + width_new), pcb_w - sensor_w + width_new);
    center.emplace("sen_pcb_stack_hist", getPolyCenter(points));
    points = GetPoints(geo, GetWidthToA(baseplate_w - pcb_w + width_new), baseplate_w - pcb_w + width_new);
    center.emplace("pcb_base_stack_hist", getPolyCenter(points));
    points = GetPoints(geo, GetWidthToA(baseplate_w - sensor_w + width_new), baseplate_w - sensor_w + width_new);
    center.emplace("sen_base_stack_hist", getPolyCenter(points));
    points = GetPoints(geo, GetWidthToA(kapton_w - pcb_w + width_new), kapton_w - pcb_w + width_new);
    center.emplace("pcb_kap_stack_hist", getPolyCenter(points));
    points = GetPoints(geo, GetWidthToA(kapton_w - sensor_w + width_new), kapton_w - sensor_w + width_new);
    center.emplace("sen_kap_stack_hist", getPolyCenter(points));
    points = GetPoints(geo, GetWidthToA((kapton_w  - (gold_size - 0.650)*2) - (pcb_w - gold_size*2) + width_new), (kapton_w  - (gold_size - 0.650)*2) - (pcb_w - gold_size*2) + width_new);
    center.emplace("kap_pcb_hist", getPolyCenter(points));
    points = GetPoints(geo, GetWidthToA((sensor_w - (guard_ring_size - 0.580 - 0.150)*2) - (pcb_w - gold_size*2) + width_new), (sensor_w - (guard_ring_size - 0.580 - 0.150)*2) - (pcb_w - gold_size*2) + width_new);
    center.emplace("sen_pcb_hist", getPolyCenter(points));

    int nPeaks = 1;
         if(geo == "Five" || geo == "Semi") nPeaks = 3;
    else if(geo == "Half") nPeaks = 2;
    for(auto &type_str: Dist){
      TString type = TString(type_str);
      if(type.Contains("PCB")){
        pcb_w = 166.64;
        TString buffer = type;
        TString pcb_new_shift_str = buffer.ReplaceAll("Gaussian_PCBplus", "");
        pcb_new_shift_str.ReplaceAll("Gaussian_PCBminus", "");
        float pcb_new_shift = pcb_new_shift_str.Atof()/1000;
       
        if(type.Contains("plus")) pcb_w  = pcb_w + pcb_new_shift;
        if(type.Contains("minus")) pcb_w = pcb_w -  pcb_new_shift;
        cout << "pcb_w: " << pcb_w << "pcb_new_shift: " << pcb_new_shift << endl;
      } else pcb_w = 166.64;
      if(type.Contains("newSensor")){ 
        sensor_shift_x = 0.050; 
        sensor_shift_y = 0.050;
      } else {
	sensor_shift_x = 0.500;
	sensor_shift_y = 0.050;
      }
      vector<double> pcb, kapton, baseplate, sensor, pcb_bond, kapton_bond, sensor_bond;
      double x_shift = 0., s_x_shift = 0., y_shift = 0., s_y_shift = 0.;
      map<TString, TH2Poly*> components {
      	{"pcb", new TH2Poly("pcb_2d", ";" + geo + " " + type + " Width [mm];Width [mm]", -1*axis, axis, -1*axis, axis)},
      	{"kapton", new TH2Poly("kapton_2d", ";" + geo + " " + type + " Width [mm];Width [mm]", -1*axis, axis, -1*axis, axis)},
      	{"sensor", new TH2Poly("sensor_2d", ";" + geo + " " + type + " Width [mm];Width [mm]", -1*axis, axis, -1*axis, axis)},
      	{"baseplate", new TH2Poly("baseplate_2d", ";" + geo + " " + type + " Width [mm];Width [mm]", -1*axis, axis, -1*axis, axis)},
      };
      map<TString, TH2Poly*> overlap {
      	{"kap_pcb_hist",        new TH2Poly("Gold Bond",         ";" + geo + " " + type + " Diff Width [mm];Diff Width [mm]", -1*axis, 1*axis, -1*axis, 1*axis)}, 			
      	{"sen_pcb_hist",        new TH2Poly("Guard Bond",        ";" + geo + " " + type + " Diff Width [mm];Diff Width [mm]", -1*axis, 1*axis, -1*axis, 1*axis)}, 			
      	{"pcb_base_stack_hist", new TH2Poly("Overlap pcb - base",";" + geo + " " + type + " Diff Width [mm];Diff Width [mm]", -1*axis, 1*axis, -1*axis, 1*axis)},    
      	{"sen_base_stack_hist", new TH2Poly("Overlap sen - base",";" + geo + " " + type + " Diff Width [mm];Diff Width [mm]", -1*axis, 1*axis, -1*axis, 1*axis)},
      	{"pcb_kap_stack_hist",  new TH2Poly("Overlap pcb - kap", ";" + geo + " " + type + " Diff Width [mm];Diff Width [mm]", -1*axis, 1*axis, -1*axis, 1*axis)},
      	{"sen_kap_stack_hist",  new TH2Poly("Overlap sen - kap", ";" + geo + " " + type + " Diff Width [mm];Diff Width [mm]", -1*axis, 1*axis, -1*axis, 1*axis)},
      	{"sen_pcb_stack_hist",  new TH2Poly("Overlap sen - pcb", ";" + geo + " " + type + " Diff Width [mm];Diff Width [mm]", -1*axis, 1*axis, -1*axis, 1*axis)},
      };
     
      for(std::map<TString,TH2Poly*>::iterator iter = components.begin(); iter != components.end(); ++iter)
        HoneycombCustom(geo, iter->second, step, width_new, a_new, nbins);
      for(std::map<TString,TH2Poly*>::iterator iter = overlap.begin(); iter != overlap.end(); ++iter)
        HoneycombCustom(geo, iter->second, step, width_new, a_new, nbins);

      cout << type << endl;
      for(auto i = 0; i != max; i++){
        double mean_ = 1., sigma_ = 1., kapton_mean_ = 1., kapton_sigma_ = 1.;
        pcb.clear();
        kapton.clear();
        baseplate.clear();
        sensor.clear();
        if(type.Contains("CustomGaus")){
          baseValue_total = baseValue*mean;
          mean_ = mean;
          sigma_ = sigma;
          kapton_mean_ = kapton_mean;
          kapton_sigma_ = kapton_sigma;
        }
        if(i == 0){
          cout << "baseplate: " << baseplate_w*mean_ << " +/- " << baseplate_w_err*sigma_ << endl;
          cout << "pcb: " << pcb_w*kapton_mean_ << " +/- " << pcb_w_err*kapton_sigma_ << endl;
          cout << "kapton: " << kapton_w*kapton_mean_ << " +/- " << kapton_w_err*kapton_sigma_ << endl;
          cout << "sensor: " << sensor_w*kapton_mean_ << " +/- " << sensor_w_err*kapton_sigma_ << endl;
          cout << "center: " << baseplate_err << endl;
          cout << "PCB shift: " << pcb_shift_r<< endl;
          cout << "Sensor shift: " << sensor_shift_x << endl;
	}
	//Get widths for each component and each side
        for(auto iS = 0; iS != 3; iS++){
          pcb.push_back(getRandomValue(pcb_w*kapton_mean_, pcb_w_err*kapton_sigma_, type));
          kapton.push_back(getRandomValue(kapton_w*kapton_mean_, kapton_w_err*kapton_sigma_, type));
          baseplate.push_back(getRandomValue(baseplate_w*mean_, baseplate_w_err*sigma_, type));
          sensor.push_back(getRandomValue(sensor_w*kapton_mean_, sensor_w_err*kapton_sigma_, type));
          kapton_bond.push_back(kapton[iS] - (gold_size - 0.650)*2);
          sensor_bond.push_back(sensor[iS] - (guard_ring_size - 0.580 - 0.150)*2);
          pcb_bond.push_back(pcb[iS] - gold_size*2);
        }
        //Get shift in x and y, soon add theta also
        double 
        x_shift = getRandomValue(0, baseplate_err, type);
        y_shift = getRandomValue(0, baseplate_err, type);
        vector<double> pcb_shift = GetCartesian(getRandomValue(0, pcb_shift_r, type), getRandomValue(0, pcb_shift_theta, type), 0., 0.);
        s_x_shift = getRandomValue(0, sensor_shift_x, type);
        s_y_shift = getRandomValue(0, sensor_shift_y, type);

        FillAllSides(geo, double(max), center["baseplate"], components["baseplate"], baseplate, baseValue_total, x_shift, y_shift, width_new);
        FillAllSides(geo, double(max), center["kapton"],    components["kapton"],    kapton, baseValue_total, x_shift, y_shift, width_new);
        FillAllSides(geo, double(max), center["sensor"],    components["sensor"],    sensor, baseValue_total, s_x_shift, s_y_shift, width_new);
        FillAllSides(geo, double(max), center["pcb"],       components["pcb"], 	pcb, baseValue_total, pcb_shift[0], pcb_shift[1], width_new);
        FillAllSides(geo, double(max), center["sen_pcb_stack_hist"], overlap["sen_pcb_stack_hist"], 	pcb, baseValue_total, pcb_shift[0] + s_x_shift, pcb_shift[1] + s_y_shift, width_new, 
										     	    		sensor, s_x_shift, s_y_shift);
        FillAllSides(geo, double(max), center["pcb_base_stack_hist"], overlap["pcb_base_stack_hist"], baseplate, baseValue_total, x_shift + s_x_shift, y_shift + s_y_shift, width_new, 
										            		   pcb, pcb_shift[0] + s_x_shift, pcb_shift[1] + s_y_shift);
        FillAllSides(geo, double(max), center["sen_base_stack_hist"], overlap["sen_base_stack_hist"], baseplate, baseValue_total, x_shift + s_x_shift, y_shift + s_y_shift, width_new, 
											    		      sensor, s_x_shift, s_y_shift);
        FillAllSides(geo, double(max), center["pcb_kap_stack_hist"], overlap["pcb_kap_stack_hist"],     kapton, baseValue_total, x_shift + s_x_shift, y_shift + s_y_shift, width_new, 
										            		   pcb, pcb_shift[0] + s_x_shift, pcb_shift[1] + s_y_shift);
        FillAllSides(geo, double(max), center["sen_kap_stack_hist"], overlap["sen_kap_stack_hist"],  kapton, baseValue_total, 0.0, 0.0, width_new, 
										            		   sensor, s_x_shift, s_y_shift);
        FillAllSides(geo, double(max), center["kap_pcb_hist"], overlap["kap_pcb_hist"], 	   kapton_bond, baseValue_total, x_shift + s_x_shift, y_shift + s_y_shift, width_new, 
										  	    		   pcb_bond, pcb_shift[0] + s_x_shift, pcb_shift[1] + s_y_shift);
        FillAllSides(geo, double(max), center["sen_pcb_hist"], overlap["sen_pcb_hist"], 	   sensor_bond, baseValue_total, s_x_shift, s_y_shift, width_new, 
										  	    		   pcb_bond, pcb_shift[0] + s_x_shift, pcb_shift[1] + s_y_shift);
      }

      for(std::map<TString,TH2Poly*>::iterator iter = components.begin(); iter != components.end(); ++iter){
        cout << iter->first << ": (x,y) = (" << center[iter->first].first << "," << center[iter->first].second << ")" << endl;
        print2DPlots(iter->second,      geo,	iter->first + " width",       	geo + "_" + type + "_Width_"+iter->first, width_new);
      }
      print2DPlots(overlap["kap_pcb_hist"], 		geo,	"Gold Bond", 		geo + "_" + type + "_Diff_Width_guard_bond", width_new);
      print2DPlots(overlap["sen_pcb_hist"], 		geo,	"Guard Bond", 		geo + "_" + type + "_Diff_Width_gold_bond", width_new);
      print2DPlots(overlap["sen_pcb_stack_hist"], 	geo,	"(PCB - Sensor)", 	geo + "_" + type + "_Diff_Width_sen_pcb_stack", width_new);
      print2DPlots(overlap["pcb_base_stack_hist"], 	geo,	"(Baseplate - PCB)", 	geo + "_" + type + "_Diff_Width_pcb_base_stack", width_new);
      print2DPlots(overlap["sen_base_stack_hist"], 	geo,	"(Baseplate - Sensor)", geo + "_" + type + "_Diff_Width_sen_base_stack", width_new);
      print2DPlots(overlap["pcb_kap_stack_hist"], 	geo,	"(Kapton - PCB)", 	geo + "_" + type + "_Diff_Width_pcb_kap_stack", width_new);
      print2DPlots(overlap["sen_kap_stack_hist"], 	geo,	"(Kapton - Sensor)", 	geo + "_" + type + "_Diff_Width_sen_kap_stack", width_new);

      map<TString, TH1D*> overlap_1D;
      int iC = 0;
      for(std::map<TString,TH2Poly*>::iterator iter = overlap.begin(); iter != overlap.end(); ++iter){
        cout << iter->first << ": (x,y) = (" << center[iter->first].first << "," << center[iter->first].second << ")" << endl;
        vector<TH1*> hist = IntegrateHex(iter->second, geo, type, center[iter->first], width_new, axis, step, nPeaks);
        normalize(hist);
        prepHists(hist, false, false, false, {colors[iC], colors[iC], colors[iC], colors[iC]});
        for(auto iP = 0; iP != nPeaks + 1; iP++)
          overlap_1D.insert(pair<TString, TH1D*>(iter->first + "_Peak" + to_string(iP), (TH1D*)hist[iP]));
        iC++;
      }

      TString name = "";
      for(auto iP = 0; iP != nPeaks + 1; iP++){
        TString pName = "_Peak" + to_string(iP);
        auto leg_diff = prepLegends({overlap_1D["pcb_base_stack_hist" + pName], overlap_1D["pcb_kap_stack_hist" + pName]}, 
                                    {overlap_1D["pcb_base_stack_hist" + pName]->GetName(), overlap_1D["pcb_kap_stack_hist" + pName]->GetName()}, "P");
        leg_diff->SetTextSize(0.03);
        leg_diff->SetY1NDC(leg_diff->GetY2NDC() - 0.3);
        TCanvas* diff = drawCompMatt({overlap_1D["pcb_base_stack_hist" + pName], overlap_1D["pcb_kap_stack_hist" + pName]}, leg_diff, 0.001);
        name = geo + "_" + type + "_diff" + pName;
        diff->Update();
        diff->SetTitle(name);
        diff->Print(outputDir+"/"+name+".pdf");

        auto leg_sen_diff = prepLegends({overlap_1D["sen_base_stack_hist" + pName], overlap_1D["sen_pcb_stack_hist" + pName], overlap_1D["sen_kap_stack_hist" + pName]}, 
                                    {overlap_1D["sen_base_stack_hist" + pName]->GetName(), overlap_1D["sen_pcb_stack_hist" + pName]->GetName(), overlap_1D["sen_kap_stack_hist" + pName]->GetName()}, "P");
        leg_sen_diff->SetTextSize(0.03);
        leg_sen_diff->SetY1NDC(leg_sen_diff->GetY2NDC() - 0.3);
        TCanvas* sen_diff = drawCompMatt({overlap_1D["sen_base_stack_hist" + pName], overlap_1D["sen_pcb_stack_hist" + pName], overlap_1D["sen_kap_stack_hist" + pName]}, leg_sen_diff, 0.001);
        name = geo + "_" + type + "_sen_diff" + pName;
        sen_diff->Update();
        sen_diff->SetTitle(name);
        sen_diff->Print(outputDir+"/"+name+".pdf");

        auto leg_bond_diff = prepLegends({overlap_1D["sen_pcb_hist" + pName], overlap_1D["kap_pcb_hist" + pName]},
                                    {overlap_1D["sen_pcb_hist" + pName]->GetName(), overlap_1D["kap_pcb_hist" + pName]->GetName()}, "P");
        leg_bond_diff->SetTextSize(0.03);
        leg_bond_diff->SetY1NDC(leg_bond_diff->GetY2NDC() - 0.3);
        TCanvas* bond_diff = drawCompMatt({overlap_1D["sen_pcb_hist" + pName], overlap_1D["kap_pcb_hist" + pName]}, leg_bond_diff, 0.001);
        name = geo + "_" + type + "_bond_diff" + pName;
        bond_diff->Update();
        bond_diff->SetTitle(name);
        bond_diff->Print(outputDir+"/"+name+".pdf");
      }

      if(fit_values.size() == 0){
        double base_nom = (baseplate_w - baseValue_total)/2;
        double kap_nom = (kapton_w - baseValue_total)/2;
        double sen_nom = (sensor_w - baseValue_total)/2;
        double pcb_nom = (pcb_w - baseValue_total)/2;
        double guard_ring_nom = (sensor_w/2 - (guard_ring_size - 0.580 - 0.150)) - (pcb_w/2 - guard_ring_size);
        double gold_nom = (kapton_w/2 - (gold_size - 0.650)) - (pcb_w/2 - gold_size);

        fit_values.push_back(make_pair(make_pair("Nominal", "Guard Bond"), 	make_pair((guard_ring_nom - guard_ring_min), 0.0))); 
        fit_values.push_back(make_pair(make_pair("Nominal", "Gold Bond"), 	make_pair((gold_nom - gold_min), 0.0)));
        fit_values.push_back(make_pair(make_pair("Nominal", "PCB to Base"), 	make_pair((base_nom - pcb_nom), 0.0)));
        fit_values.push_back(make_pair(make_pair("Nominal", "Sen to Base"),	make_pair((base_nom - sen_nom), 0.0)));
        fit_values.push_back(make_pair(make_pair("Nominal", "PCB to Kap"), 	make_pair((kap_nom - pcb_nom), 0.0)));
        fit_values.push_back(make_pair(make_pair("Nominal", "Sen to Kap"), 	make_pair((kap_nom - sen_nom), 0.0)));
        fit_values.push_back(make_pair(make_pair("Nominal", "Sen to PCB"), 	make_pair((pcb_nom - sen_nom), 0.0)));
      }

      for(std::map<TString,TH1D*>::iterator iter = overlap_1D.begin(); iter != overlap_1D.end(); ++iter){
        cout << iter->first << endl;
        if(TString(iter->first).Contains("_Peak0")) continue;

        vector<double> max_tol = findHighXbin(iter->second);
        string name = "";
             if(iter->first.Contains("pcb_base_stack")) name = "PCB to Base";
        else if(iter->first.Contains("sen_pcb_stack"))  name = "Sen to PCB";
        else if(iter->first.Contains("sen_base_stack")) name = "Sen to Base";
        else if(iter->first.Contains("pcb_kap_stack"))  name = "PCB to Kap";
        else if(iter->first.Contains("sen_kap_stack"))  name = "Sen to Kap";
        else if(iter->first.Contains("kap_pcb_hist"))   name = "Gold Bond";
        else if(iter->first.Contains("sen_pcb_hist"))   name = "Guard Bond";
        cout << name << " lowest xbin: " << max_tol[0] << " +/- " << max_tol[1] << " P(<0) =  " << max_tol[2] << endl;
        worst_values.push_back(make_pair(make_pair(type_str, name), make_pair(max_tol[0]/1000, max_tol[1]/1000)));

        //Get Fit
        TString dist_type = "gausn";
        if(type.Contains("Land")) dist_type = "landau";

        iter->second->Fit(dist_type);
        TF1* avg_param = (TF1*)iter->second->GetFunction(dist_type);
        fit_values.push_back(make_pair(make_pair(type_str, name), make_pair(avg_param->GetParameter(1), avg_param->GetParameter(2))));
        cout << "Fit to " << name << " = " << avg_param->GetParameter(1)*1000 << " +/- " << avg_param->GetParameter(2)*1000 << endl;
      }
      
      TFile *outFile = new TFile(outputDir+"/"+name+".root", "RECREATE");
      for(std::map<TString,TH2Poly*>::iterator iter = components.begin(); iter != components.end(); ++iter){
        iter->second->Write();
        delete iter->second;
      }
      for(std::map<TString,TH1D*>::iterator iter = overlap_1D.begin(); iter != overlap_1D.end(); ++iter){
        iter->second->Write();
        delete iter->second;
      }
      outFile->Close();
    }

    std::ofstream jout;
    jout.open(outputDir+"/" + geo + "_moduleTolerances.json");
    json jtot = makeJSONModule(fit_values, worst_values);
    jout << jtot.dump(3);
    jout.close();

    for (json::iterator d = jtot.begin(); d != jtot.end(); ++d){
      if(d.key() == "BinNum") continue;
      if(d.key() == "Fit") continue;
      if(d.key() == "Worst") continue;
      vector<TH1*> hFit, hWorst;
      vector<TString> binLabels(7);
      vector<Quantity> hNom_val(7);
      for (json::iterator comp = jtot[d.key()].begin(); comp != jtot[d.key()].end(); ++comp) {
        vector<Quantity> hFit_val(jtot["Fit"][comp.key()].size());
        vector<Quantity> hWorst_val(jtot["Fit"][comp.key()].size());

        for (json::iterator dist = jtot["Fit"][comp.key()].begin(); dist != jtot["Fit"][comp.key()].end(); ++dist) {
          int binnum = jtot["BinNum"][dist.key()];
          binLabels.at(binnum) = dist.key();
          if(TString(comp.key()) == "Nominal"){
            hNom_val.at(binnum).value = jtot["Fit"][comp.key()][dist.key()][0];
            hNom_val.at(binnum).error = 0.0;
          } else {
            hFit_val.at(binnum).value = jtot["Fit"][comp.key()][dist.key()][0];
            hFit_val.at(binnum).error = jtot["Fit"][comp.key()][dist.key()][1];
          }
        }
        for (json::iterator dist = jtot["Worst"][comp.key()].begin(); dist != jtot["Worst"][comp.key()].end(); ++dist) {
          int binnum = jtot["BinNum"][dist.key()];
          hWorst_val.at(binnum).value = jtot["Worst"][comp.key()][dist.key()][0];
          hWorst_val.at(binnum).error = jtot["Worst"][comp.key()][dist.key()][1];
        }

        if(TString(comp.key()) != "Nominal"){ 
          hFit.push_back(convertToHist(hFit_val, "Fit" + TString(d.key()) + TString(comp.key()), ";; Fitted Value", nullptr));
          hWorst.push_back(convertToHist(hWorst_val, "Worst" + TString(d.key()) + TString(comp.key()), ";; Worst Value", nullptr));
        }
      }
      TH1* hNominal = convertToHist(hNom_val, "Nominal", ";; Nominal Value", nullptr);

      prepHists({hNominal}, false, false, false, {kBlack});
      prepHists(hFit, false, false, false, colors);
      prepHists(hWorst, false, false, false, colors);

      setBinLabels(hNominal, binLabels);
      for(auto *f : hFit) setBinLabels(f, binLabels);
      for(auto *f : hWorst) setBinLabels(f, binLabels);

      vector<TH1*> hDiv;
      for(auto *f : hFit){
        TH1* hratio = (TH1*)f->Clone();
        hratio->Divide(hNominal);
        hDiv.push_back(hratio);
      }
      vector<TH1*> hDiv_w;
      for(auto *f : hWorst){
        TH1* hratio = (TH1*)f->Clone();
        hratio->Divide(hNominal);
        hDiv_w.push_back(hratio);
      }

      auto leg = prepLegends({}, {""}, "L");
      appendLegends(leg, {hNominal}, {"Nominal"}, "L");
      for(auto h = 0; h != hFit.size(); h++){
        TString legName = hFit[h]->GetName();
        legName.ReplaceAll("FitDistribution", "");
        legName.ReplaceAll("FitPCB_Dist", "");
        legName.ReplaceAll("Gaussian_PCBminus", "PCB - ");
        legName.ReplaceAll("Gaussian_PCBplus", "PCB + ");
        legName.ReplaceAll("Gaussian_newSensor", "sensor [0.05 mm]");
        legName.ReplaceAll("_newSensor", " [0.05 mm]");
        legName.ReplaceAll("FitCustom", "");
        legName.ReplaceAll("FitSensor", "");
        appendLegends(leg, {hFit[h]}, {legName}, "L");
      }
      hFit.push_back(hNominal);
      leg->SetTextSize(0.04);
      leg->SetY1NDC(leg->GetY2NDC() - 0.2);
      if(TString(d.key()) == "PCB_Dist") setLegend(leg, 2, 0.45, 0.6, 0.92, 0.87);
      TCanvas* c = drawCompAndRatio(hFit, hDiv, leg, "Fit/Nominal", 0.499, 1.501, true, 0.01, 2., false, nullptr, false);
      TString outputName = outputDir + "/Fit_" + geo + "_" + TString(d.key());
      c->SetTitle(outputName);
      c->Print(outputName+".pdf");

      leg = prepLegends({}, {""}, "L");
      appendLegends(leg, {hNominal}, {"Nominal"}, "L");
      for(auto h = 0; h != hWorst.size(); h++){
        TString legName = hWorst[h]->GetName();
        legName.ReplaceAll("WorstDistribution", "");
        legName.ReplaceAll("WorstPCB_Dist", "");
        legName.ReplaceAll("Gaussian_PCBminus", "PCB - ");
        legName.ReplaceAll("Gaussian_PCBplus", "PCB + ");
        legName.ReplaceAll("Gaussian_newSensor", "sensor [0.05 mm]");
        legName.ReplaceAll("_newSensor", " [0.05 mm]");
        legName.ReplaceAll("WorstCustom", "");
        legName.ReplaceAll("WorstSensor", "");
        appendLegends(leg, {hWorst[h]}, {legName}, "L");
      }
      hWorst.push_back(hNominal);
      leg->SetTextSize(0.04);
      leg->SetY1NDC(leg->GetY2NDC() - 0.2);
      if(TString(d.key()) == "PCB_Dist") setLegend(leg, 2, 0.45, 0.6, 0.92, 0.87);
      c = drawCompAndRatio(hWorst, hDiv_w, leg, "Worst/Nominal", 0.001, 2.001, true, -1., 0.7, false, nullptr, false, -0.3);
      outputName = outputDir + "/Worst_" + geo + "_" + TString(d.key());
      c->SetTitle(outputName);
      c->Print(outputName+".pdf");

      TFile *outFile = new TFile(outputName+".root", "RECREATE");
      for(auto *f : hFit) f->Write();
      for(auto *w : hWorst) w->Write();
      outFile->Close();

    }

    map<TString, pair<double, double> > averages;
    for ( const pair< pair< string, string>, pair< double, double > > &avg : fit_values ){
      if(avg.first.first == "Nominal") continue;
      double average = 0., error = 0.;
      map<TString, pair<double, double> >::iterator it = averages.find(avg.first.second);
      if(it != averages.end()){
        average = it->second.first + avg.second.first;
        error = it->second.second + avg.second.second;
        averages.erase(it);
      } else {
        average = avg.second.first;
        error = avg.second.second;
      }
      averages.insert(pair<TString, pair<double, double> >(avg.first.second, pair<double, double>(average, error)));
    }
    for(std::map<TString,pair<double, double> >::iterator iter = averages.begin(); iter != averages.end(); ++iter){
      cout << "Average of " << iter->first << " = " << iter->second.first/Dist.size() << " +/- " << iter->second.second/Dist.size() << endl;
    }
  }
}

float getMinimum(TH1F* p_Hist){
  float min = p_Hist->GetMinimum(0);
  float temp_min = 99999.;
  return min;
}

float histYield(TH1F* hist){
  return hist->Integral(0, hist->GetNbinsX() + 1);
}

void setTitleOffset(TCanvas *c, double xOff = .950, double yOff = 1.400){
  TList * list = c->GetListOfPrimitives();
  for(auto iP = 0; iP < list->GetSize(); ++iP){
    TH1 * h = dynamic_cast<TH1*>(list->At(iP));
    if(h == 0) continue;
   h->GetXaxis()->SetTitleOffset(xOff);
    h->GetYaxis()->SetTitleOffset(yOff);
  }
}

TH2D *create2Dhisto( pair <TTree*, TString>& tree,TString intLumi,TString cuts,
                     pair <TString, TString>& branchX,int binsX,float xmin,float xmax,
                     pair <TString, TString>& branchY,int binsY,float ymin,float ymax) {

  //TH1::SetDefaultSumw2(kTRUE);
  delete gROOT->FindObject("hTemp");

  TString cut ="("+cuts+")";
  TH2D *hTemp = new TH2D("hTemp","hTemp",binsX,xmin,xmax,binsY,ymin,ymax);
  tree.first->Project("hTemp",branchY.first+":"+branchX.first,cut);

  hTemp->GetXaxis()->SetTitle(branchX.second);
  hTemp->GetYaxis()->SetTitle(branchY.second);

  return hTemp;
} //~ end of create2Dhisto 

void print2DPlots(TH2Poly *hc, TString geometry, TString BinLatex, TString name, double width){
  TString hName = TString(hc->GetName()) + "_c";
  TCanvas* canvas;
  delete gROOT->FindObject(hName);

  canvas = new TCanvas(hName, "2Dhist", 500, 500);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetRightMargin(0.19);
  canvas->SetTopMargin(0.08);
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.14);
  int NCont = 255;
  double stops[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
  double red[5]= {0.50, 0.50, 1.00, 1.00, 1.00};
  double green[5] = {0.50, 1.00, 1.00, 0.60, 0.50};
  double blue[5] = {1.00, 1.00, 0.50, 0.40, 0.50};
  TColor::CreateGradientColorTable(5, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  canvas->cd();
  hc->Draw("colz");

  gPad->Update();
  TPaletteAxis* palette = (TPaletteAxis*)hc->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(1.-0.18);
  palette->SetY1NDC(0.14);
  palette->SetX2NDC(1.-0.13);
  palette->SetY2NDC(1.-0.08);
  palette->SetLabelFont(42);
  palette->SetLabelSize(0.035);
  canvas->Modified();
  if(BinLatex != ""){
    TLatex *tl = new TLatex();
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.17, 0.87, BinLatex);
  }

  AddLineNomHex(width, geometry);

  setTitleOffset(canvas, .950, .950);
  canvas->Update();
  canvas->SaveAs(outputDir + "/" + name + ".pdf");
}

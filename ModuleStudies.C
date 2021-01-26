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
#include <iomanip>
#include<bits/stdc++.h> 

using namespace std;
using namespace EstTools;
using json = nlohmann::json;

TString outputDir = "";
TString baseDir = "moduleTolerances_complete_012521_backside";
vector<Color_t> colors;
bool debug = false;
TString whichOverlap = "Backside Bond";

vector<TString> Order = {
			 "Gaussian_Kaptonminus170_newSensor", "Gaussian_Kaptonminus170_midSensor", "Gaussian_Kaptonminus170_oldSensor",  
			 "Gaussian_Kaptonminus150_newSensor", "Gaussian_Kaptonminus150_midSensor", "Gaussian_Kaptonminus150_oldSensor",  
			 "Gaussian_Kaptonminus75_newSensor", "Gaussian_Kaptonminus75_midSensor", "Gaussian_Kaptonminus75_oldSensor",  
			 "Gaussian_Kaptonminus50_newSensor", "Gaussian_Kaptonminus50_midSensor", "Gaussian_Kaptonminus50_oldSensor",  
			 "Gaussian_Kaptonminus25_newSensor", "Gaussian_Kaptonminus25_midSensor", "Gaussian_Kaptonminus25_oldSensor",  
			 "Gaussian_Kaptonminus0_newSensor", "Gaussian_Kaptonminus0_midSensor", "Gaussian_Kaptonminus0_oldSensor", 
			 "Gaussian_Kaptonplus25_newSensor", "Gaussian_Kaptonplus25_midSensor", "Gaussian_Kaptonplus25_oldSensor",  
			 "Gaussian_Kaptonplus50_newSensor", "Gaussian_Kaptonplus50_midSensor", "Gaussian_Kaptonplus50_oldSensor",  
			 "Gaussian_Kaptonplus75_newSensor", "Gaussian_Kaptonplus75_midSensor", "Gaussian_Kaptonplus75_oldSensor"};

std::map<TString, TString> latexMap{
  {"Custom", ""},
  {"Distribution", ""},
  {"PCBDist", ""},
  {"KaptonMultiDist", ""},
  {"TotalBadModules", ""},
  {"Sensor", ""},
  {"otherCenter", R"(center $= + 29 \mu m$)"},
  {"PCBminus25", R"(PCB - 25 $\mu m$)"},
  {"PCBplus25", R"(PCB + 25 $\mu m$)"},
  {"PCBminus50", R"(PCB - 50 $\mu m$)"},
  {"PCBplus50", R"(PCB + 50 $\mu m$)"},
  {"PCBminus75", R"(PCB - 75 $\mu m$)"},
  {"PCBplus75", R"(PCB + 75 $\mu m$)"},
  {"Kaptonminus0", R"(Kapton - 0 $\mu m$)"},
  {"Kaptonminus25", R"(Kapton - 25 $\mu m$)"},
  {"Kaptonplus25", R"(Kapton + 25 $\mu m$)"},
  {"Kaptonminus50", R"(Kapton - 50 $\mu m$)"},
  {"Kaptonplus50", R"(Kapton + 50 $\mu m$)"},
  {"Kaptonminus150", R"(Kapton - 150 $\mu m$)"},
  {"Kaptonminus170", R"(Kapton - 170 $\mu m$)"},
  {"Kaptonminus75", R"(Kapton - 75 $\mu m$)"},
  {"Kaptonplus75", R"(Kapton + 75 $\mu m$)"},
  {"Nominal", "Nominal"},
  {"Gaussian", "Gaussian"},
  {"Landau", "Landau"},
  {"Flat", "Flat"},
  {"CustomGaus", "Custom Gaussian"},
  {"CustomLandau", "Custom Landau"},
  {"CustomFlat", "Custom Flat"},
  {"newSensor", R"([50 $\mu m$])"},
  {"midSensor", R"([100 $\mu m$])"},
  {"oldSensor", R"([200 $\mu m$])"},
  {"newSplitSensor", R"( & [50 $\mu m$])"},
  {"midSplitSensor", R"( & [100 $\mu m$])"},
  {"oldSplitSensor", R"( & [200 $\mu m$])"},
  {"Fit", ""},
  {"Worst", ""},
  {"Peak1", "Peak 1"},
  {"Peak2", "Peak 2"},
  {"Peak3", "Peak 3"},
};

std::map<TString, TString> plotMap{
  {"Custom", ""},
  {"Distribution", ""},
  {"PCBDist", ""},
  {"KaptonMultiDist", ""},
  {"TotalBadModules", ""},
  {"Sensor", ""},
  {"otherCenter", R"(center $= + 29 #mum)"},
  {"PCBminus25", R"(PCB - 25 #mum)"},
  {"PCBplus25", R"(PCB + 25 #mum)"},
  {"PCBminus50", R"(PCB - 50 #mum)"},
  {"PCBplus50", R"(PCB + 50 #mum)"},
  {"PCBminus75", R"(PCB - 75 #mum)"},
  {"PCBplus75", R"(PCB + 75 #mum)"},
  {"Kaptonminus0", R"(Kapton - 0 #mum)"},
  {"Kaptonminus25", R"(Kapton - 25 #mum)"},
  {"Kaptonplus25", R"(Kapton + 25 #mum)"},
  {"Kaptonminus50", R"(Kapton - 50 #mum)"},
  {"Kaptonplus50", R"(Kapton + 50 #mum)"},
  {"Kaptonminus150", R"(Kapton - 150 #mum)"},
  {"Kaptonminus170", R"(Kapton - 170 #mum)"},
  {"Kaptonminus75", R"(Kapton - 75 #mum)"},
  {"Kaptonplus75", R"(Kapton + 75 #mum)"},
  {"Nominal", "Nominal"},
  {"Gaussian", "Gaussian"},
  {"Landau", "Landau"},
  {"Flat", "Flat"},
  {"CustomGaus", "Custom Gaussian"},
  {"CustomLandau", "Custom Landau"},
  {"CustomFlat", "Custom Flat"},
  {"newSensor", R"([50 #mum])"},
  {"midSensor", R"([100 #mum])"},
  {"oldSensor", R"([200 #mum])"},
  {"Fit", "Fit"},
  {"Worst", "Worst"},
  {"Peak1", "Peak 1"},
  {"Peak2", "Peak 2"},
  {"Peak3", "Peak 3"},
};

void print2DPlots(TH2Poly *hc, TString geometry, TString BinLatex = "", TString name = "testHoneycomb", double width = 0.2);
void print2DPlots(TH2D *hc, TString geometry, TString BinLatex = "", TString name = "testHoneycomb");
void moduleTolerances();
#endif

void ModuleStudies(){
  moduleTolerances();
}

float Round(float var, float decimal = 1000.){ 
    float value = (int)(var * decimal + .5); 
    return (float)value / decimal; 
} 

void printToleranceTableLatex(json jtot, TString outputfile="/tmp/yields.tex"){
  ofstream outfile(outputfile);
  Quantity::printStyle = Quantity::LATEX;

  outfile << R"(\documentclass[12pt,notitlepage]{revtex4-1})" << endl;
  outfile << R"(\usepackage{amsmath,mathrsfs})" << endl;
  outfile << R"(\usepackage{graphicx})" << endl;
  outfile << R"(\usepackage{amsmath})" << endl;
  outfile << R"(\usepackage{siunitx})" << endl;
  outfile << R"(\usepackage{graphicx})" << endl;
  outfile << R"(\usepackage{indentfirst})" << endl;
  outfile << R"(\usepackage{amssymb})" << endl;
  outfile << R"(\usepackage{natbib})" << endl;
  outfile << R"(\usepackage{longtable})" << endl;
  outfile << R"(\usepackage[hidelinks]{hyperref})" << endl;
  outfile << R"(\usepackage{listings})" << endl;
  outfile << R"(\usepackage{rotating})" << endl;
  outfile << R"(\usepackage[utf8]{inputenc})" << endl;
  outfile << R"(\usepackage[english]{babel})" << endl;
  outfile << R"(\usepackage[usenames, dvipsnames]{color})" << endl;
  outfile << R"(\pagenumbering{gobble})" << endl;
  outfile << R"(\def\arraystretch{0.5})" << endl;
  outfile << R"(\usepackage[margin=1cm]{geometry})" << endl;
  outfile << R"(\begin{document})" << endl;

  for (json::iterator type = jtot.begin(); type != jtot.end(); ++type){
    if(type.key() == "BinNum") continue;
    if(type.key() == "Fit") continue;
    if(type.key() == "Worst") continue;
    outfile << R"(\newpage)" << endl;
    outfile << R"(\begin{center})" << endl;
         if(type.key() == "KaptonMultiDist") outfile << R"(\begin{longtable}{| c | c | c | c | c |})" << endl;
    else if(type.key() == "TotalBadModules") outfile << R"(\begin{longtable}{| c | c |})" << endl;
    else                                outfile << R"(\begin{longtable}{| c | c | c | c |})" << endl;
    outfile << R"(\hline )" << endl;

    int ncols = (type.key() == "KaptonMultiDist") ? 5 : 4;
    if(type.key() == "TotalBadModules") ncols = 2;
    //key for value type
         if(type.key() == "KaptonMultiDist") outfile << "Component Overlaps & \t Sensor Placement & \t Nominal " << R"([$\mu m]$)" << " & \t Fitted " << R"($[\mu m]$)" << " & \t Worst " << R"($[\mu m]$)" << R"( \\)" << endl;
    else if(type.key() == "TotalBadModules") outfile << "Component Overlaps & \t Total Number of Bad Modules" << endl;
    else                                outfile << "Component Overlaps & \t Nominal " << R"([$\mu m]$)" << " & \t Fitted " << R"($[\mu m]$)" << " & \t Worst " << R"($[\mu m]$)" << R"( \\)" << endl;
    if(type.key() == "TotalBadModules"){
      for (json::iterator comp = jtot[type.key()].begin(); comp != jtot[type.key()].end(); ++comp) {
        if(TString(comp.key()).Contains("Nominal")) continue;
        TString buff = TString(comp.key());
        buff.ReplaceAll("Gaussian_", "");
        if(outputfile.Contains("Full")) buff.ReplaceAll("_Peak1", "");

        outfile << translateString(buff, latexMap, "_", ", ") << " & \t " << jtot["Worst"]["Total Modules"][comp.key()][2] << R"( \\)" << endl;
      }//for comparison
    } else{
      for (json::iterator dist = jtot["BinNum"].begin(); dist != jtot["BinNum"].end(); ++dist) {
        //key for overlap comparison
        outfile << R"(\hline)" << endl << R"(\multicolumn{)"+to_string(ncols)+R"(}{c}{$\vcenter{)" + dist.key() + R"(}$} \\)" << endl << R"(\hline)" << endl;
        if(type.key() == "KaptonMultiDist"){
          for (int i = 0; i != Order.size(); i++){
            for (json::iterator comp = jtot[type.key()].begin(); comp != jtot[type.key()].end(); ++comp) {
              if(TString(comp.key()).Contains("Nominal")) continue;
              if(!TString(comp.key()).Contains(Order[i])) continue;
              TString buff = TString(comp.key());
              buff.ReplaceAll("Sensor", "SplitSensor");
              buff.ReplaceAll("Gaussian_", "");
              if(outputfile.Contains("Full")) buff.ReplaceAll("_Peak1", "");

              outfile << translateString(buff, latexMap, "_", ", "); 
              outfile << " & \t " << Round(jtot["Fit"][dist.key()][comp.key() + "_Nominal"][0]) 
                      << " & \t " << Round(jtot["Fit"][dist.key()][comp.key()][0]) << R"( $\pm$ )"   << Round(jtot["Fit"][dist.key()][comp.key()][1]);
              outfile << " & \t " << Round(jtot["Worst"][dist.key()][comp.key()][0], 1.0) << R"( $\pm$ )" << Round(jtot["Worst"][dist.key()][comp.key()][1], 1.0)
                	         << R"( $N_{mod} = $)" << Round(jtot["Worst"][dist.key()][comp.key()][2], 100000.)*30000  << R"( \\)" << endl;

            }//for comparison
          if(i % 3 == 2) outfile << R"(\hline)" << endl;
          }//for Order
        } else {
          for (json::iterator comp = jtot[type.key()].begin(); comp != jtot[type.key()].end(); ++comp) {
            if(TString(comp.key()).Contains("Nominal")) continue;

            outfile << translateString(comp.key(), latexMap, "_", ", "); 
            outfile << " & \t " << Round(jtot["Fit"][dist.key()][comp.key() + "_Nominal"][0]) 
                    << " & \t " << Round(jtot["Fit"][dist.key()][comp.key()][0]) << R"( $\pm$ )"   << Round(jtot["Fit"][dist.key()][comp.key()][1]);
            outfile << " & \t " << Round(jtot["Worst"][dist.key()][comp.key()][0], 1.0) << R"( $\pm$ )" << Round(jtot["Worst"][dist.key()][comp.key()][1], 1.0)
              	         << R"( $N_{mod} = $)" << Round(jtot["Worst"][dist.key()][comp.key()][2], 100000.)*30000  << R"( \\)" << endl;

          }//for comparison
        }
      }//for overlap
    }//for total

    outfile << R"(\hline)" << endl;
    outfile << R"(\end{longtable})" << endl;
    outfile << R"(\end{center})" << endl;
  }
  outfile << R"(\end{document})" << endl;

  outfile.close();
}//for country!

void printWorstTableLatex(json jtot, TString outputfile="/tmp/yields.tex"){
  ofstream outfile(outputfile);
  Quantity::printStyle = Quantity::LATEX;

  outfile << R"(\documentclass[12pt,notitlepage]{revtex4-1})" << endl;
  outfile << R"(\usepackage{amsmath,mathrsfs})" << endl;
  outfile << R"(\usepackage{graphicx})" << endl;
  outfile << R"(\usepackage{amsmath})" << endl;
  outfile << R"(\usepackage{siunitx})" << endl;
  outfile << R"(\usepackage{graphicx})" << endl;
  outfile << R"(\usepackage{indentfirst})" << endl;
  outfile << R"(\usepackage{amssymb})" << endl;
  outfile << R"(\usepackage{natbib})" << endl;
  outfile << R"(\usepackage{longtable})" << endl;
  outfile << R"(\usepackage[hidelinks]{hyperref})" << endl;
  outfile << R"(\usepackage{listings})" << endl;
  outfile << R"(\usepackage{rotating})" << endl;
  outfile << R"(\usepackage[utf8]{inputenc})" << endl;
  outfile << R"(\usepackage[english]{babel})" << endl;
  outfile << R"(\usepackage[usenames, dvipsnames]{color})" << endl;
  outfile << R"(\pagenumbering{gobble})" << endl;
  outfile << R"(\def\arraystretch{0.5})" << endl;
  outfile << R"(\usepackage[margin=1cm]{geometry})" << endl;
  outfile << R"(\begin{document})" << endl;

  for (json::iterator type = jtot.begin(); type != jtot.end(); ++type){
    if(type.key() == "BinNum") continue;
    if(type.key() == "Fit") continue;
    if(type.key() == "Worst") continue;
    outfile << R"(\newpage)" << endl;
    outfile << R"(\begin{center})" << endl;
    outfile << R"(\begin{longtable}{| c | c | c | c |})" << endl;
    outfile << R"(\hline )" << endl;

    int ncols = 4;
    //key for value type
    outfile << "Component Overlaps" << " & \t" << "Nominal " << R"([$\mu m]$)" << " & \t" << "Fitted " << R"($[\mu m]$)" << " & \t" << "Worst " << R"($[\mu m]$)" << R"( \\)" << endl;;
    for (json::iterator dist = jtot["BinNum"].begin(); dist != jtot["BinNum"].end(); ++dist) {
      bool isFirst = true;
      //key for overlap comparison
      for (json::iterator comp = jtot[type.key()].begin(); comp != jtot[type.key()].end(); ++comp) {
        if(TString(comp.key()).Contains("Nominal")) continue;
        if((double)jtot["Worst"][dist.key()][comp.key()][2] < 0.001) continue;

        if(isFirst) outfile << R"(\hline)" << endl << R"(\multicolumn{)"+to_string(ncols)+R"(}{c}{$\vcenter{)" + dist.key() + R"(}$} \\)" << endl << R"(\hline)" << endl;
        isFirst = false;

        outfile << translateString(comp.key(), latexMap, "_", ", "); 
        outfile << " & \t " << Round(jtot["Fit"][dist.key()][comp.key() + "_Nominal"][0]) 
                << " & \t " << Round(jtot["Fit"][dist.key()][comp.key()][0]) << R"( $\pm$ )"   << Round(jtot["Fit"][dist.key()][comp.key()][1]);
        outfile << " & \t " << Round(jtot["Worst"][dist.key()][comp.key()][0], 1.0) << R"( $\pm$ )" << Round(jtot["Worst"][dist.key()][comp.key()][1], 1.0)
          	         << R"( $N_{mod}= $)" << Round(jtot["Worst"][dist.key()][comp.key()][2], 100000.)*30000  << R"( \\)" << endl;

      }//for comparison
    }//for overlap

    outfile << R"(\hline)" << endl;
    outfile << R"(\end{longtable})" << endl;
    outfile << R"(\end{center})" << endl;
  }
  outfile << R"(\end{document})" << endl;

  outfile.close();
}//for country!

double GetWidthToA(double width){
  //width/2 == a/2*tan(pi/6)
  return width*TMath::Tan(TMath::Pi()/6);
}

double GetAToWidth(double a){
  return a/TMath::Tan(TMath::Pi()/6);
}

double getCutoffWidth(double width, double cutoff, double delta = 0.){
  //B = sqrt(3)/3*(6*width - cutoff) + 2*delta
  //return (TMath::Sqrt(3)/3)*(6*width - cutoff) + 2*delta;
  double cut = cutoff - 2*TMath::Sqrt(3)*delta;
  return width - (cut/2);
}

double getCutoffInBetween(double width, double cutoff){
  return width - (cutoff/2 - 0.250)/TMath::Sqrt(3);
}

pair<double, double> GetPolar(TString name, double x1, double y1, double x0, double y0){
  double r = TMath::Sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
  double phi = TMath::ATan((y1 - y0)/(x1 - x0));
  if((x1 - x0) < 0) phi += TMath::Pi();
  if((x1 - x0) >= 0 && (y1 - y0) < 0) phi += 2*TMath::Pi();

  return make_pair(r, phi);
}

pair<double, double> GetCartesian(TString name, double r, double phi, double x0, double y0){
  double x1 = r*TMath::Cos(phi) + x0;
  double y1 = r*TMath::Sin(phi) + y0;
  
  return make_pair(x1, y1);
}

vector<double> findHighXbin(TH1D* h){
  double mod = -999., belowZero = 0., width = 0.;;
  vector<double> output;
  int nbins = h->GetNbinsX();
  bool isFirst = true;
  for(int bin = 1; bin <= nbins; bin++) {
    if(h->GetBinLowEdge(bin) > 0 && !isFirst) break;
    if(h->GetBinContent(bin) > 0 && isFirst){
      width = h->GetBinWidth(bin);
      mod = h->GetBinLowEdge(bin);
      isFirst = false;
    }
    if(h->GetBinLowEdge(bin) < 0){
      belowZero = h->GetBinContent(bin);
    }
  }
  output.push_back((mod + width/2)*1000);
  output.push_back(width*1000);
  output.push_back(mod == -999. ? 0. : belowZero);
  return output;
}

vector<TH1*> IntegrateHex(TH2Poly* h, TString geo, TString type, pair < double, double > center, double width, double axis, double binSize, int nPeaks){
  float min = -1. * binSize * 300;
  if((geo == "Semi" || geo == "Half")) min = -1. * binSize * 200;
  if(geo == "Three") min = -1. * binSize * 100;
  float max = binSize * 650;
  if(geo == "Semi" || geo == "Three") max += binSize * 50;
  if(geo == "Half") max += binSize * 150;
  int nbins = (axis - 0.1)/binSize*4;
  int nbins_1D = int((max - min)/binSize);
  vector<TH1*> hist;
  for(auto iP = 0; iP != nPeaks + 1; iP++) 
    hist.push_back(new TH1D(TString(h->GetName()) + "_Peak" + to_string(iP), ";" + type + " Diff Width [mm]; Modules", nbins_1D, min, max));

  vector<double> angle = {0., TMath::Pi()/3, TMath::Pi()*2/3, TMath::Pi(), TMath::Pi()*4/3, TMath::Pi()*5/3};
  if(geo == "Five" || geo == "Semi") angle = {0., TMath::Pi()/2, TMath::Pi(), TMath::Pi()*4/3, TMath::Pi()*5/3};
  else if(geo == "Half") angle = {0., TMath::Pi()*2/3, TMath::Pi()*4/3, TMath::Pi()*5/3};
  else if(geo == "Three") angle = {TMath::Pi()/2, TMath::Pi()*4/3, TMath::Pi()*5/3};
  vector<int> index = {1, 1, 1, 1, 1, 1};
  if(geo == "Five" || geo == "Semi") index = {1, 2, 1, 3, 3};
  else if(geo == "Half") index = {1, 2, 1, 2};
  else if(geo == "Three") index = {1, 2, 2};
  for(int iA = 0; iA != (signed)angle.size(); iA++){
    for(auto iB = 0; iB != nbins; iB++){
      double shift = 0.;
      pair<double, double> cart = GetCartesian("", (binSize*iB + (width/2 - 300*binSize) + shift), angle[iA], center.first, center.second);
      pair<double, double> polar = GetPolar("", cart.first, cart.second, center.first, center.second);
      int bin = h->FindBin(cart.first, cart.second);
      if(bin == -5) continue;
      double binValue = h->GetBinContent(bin);

      //get distance to the edge
      double r = TMath::Sqrt(cart.first*cart.first + cart.second*cart.second) - width/2 + shift;
      double dr_center = center.first*TMath::Cos(angle[iA]) + center.second*TMath::Sin(angle[iA]);
           if(geo == "Five" && iA == 1) r += width/2 - GetWidthToA(width)/2;
      else if(geo == "Semi" && iA == 1) r += width/2;
      else if(geo == "half" && iA == 1) r += width/2;
      else if(geo == "Three" && iA == 0) r += width/2 + GetWidthToA(width)/2;
      r += dr_center;

      double bin1d = hist[index[iA]]->GetXaxis()->FindBin(r);
      double tot = hist[index[iA]]->GetBinContent(bin1d);
      hist[index[iA]]->SetBinContent(bin1d, tot + binValue);
      //if(bin1d == 427){
      //  cout << "First " << h->GetName() << endl;
      //  cout << "Poly Bin: " << bin << endl;
      //  cout << "(" << cart.first << "," << cart.second << ")" << " = " << binValue << endl;
      //  cout << "bin: " << bin1d << ", r = " << r << ", value = " << tot << " + " << binValue << " = " << tot + binValue << endl;
      //}
      //fill plot with total
      bin1d = hist[0]->GetXaxis()->FindBin(r);
      tot = hist[0]->GetBinContent(bin1d);
      hist[0]->SetBinContent(bin1d, tot + binValue);
      //if(TString(h->GetName()).Contains(whichOverlap) && r > 0.255 && r < 0.257){
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
  vector<double> PCBBaseplateInput = {6.556401, 6.555577, 6.555663, 6.555603, 6.556096, 6.555644,
				6.558309, 6.558489, 6.558413, 6.557552, 6.557295, 6.557938,
				6.554244, 6.555199, 6.5551,   6.555635, 6.555705, 6.555156,
				6.556504, 6.556046, 6.556088, 6.554441, 6.5542,   6.55515,
				6.559327, 6.555097, 6.559529, 6.558453, 6.555274, 6.557107,
				6.56,     6.555561, 6.559438, 6.559253, 6.555536, 6.559403,
				6.559719, 6.555871, 6.558826
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
  TCanvas* c = drawCompMatt({baseplate_meas_hist, kapton_meas_hist}, leg, -1., nullptr, "hist", true);
  gStyle->SetOptStat(0);
  TString name = "measured_dist";
  c->SetTitle(name);
  c->Print(outputDir+"/"+name+".pdf");

  vector < pair <double, pair <double, double> > > avg_param;
  avg_param.push_back(make_pair(baseplate_meas_avg_param->GetParameter(0), make_pair(baseplate_meas_hist->GetMean(), baseplate_meas_hist->GetStdDev())));
  avg_param.push_back(make_pair(kapton_meas_avg_param->GetParameter(0), make_pair(kapton_meas_hist->GetMean(), kapton_meas_hist->GetStdDev())));

  return avg_param;
}

json makeJSONModuleLatex(vector< pair< pair< string, string>, pair< double, double > > >& fit, vector< pair< pair< string, string>, pair< pair< double, double >, double > > >& worst){
  json j;
  int binnum = 0;
  for ( const auto &it : fit ) {
    double fit = it.second.first;
    double error = it.second.second;
    //if(fit == NULL) fit = 0.0;
    //if(error == NULL) error = 0.0;
    
    j["Fit"][it.first.second][it.first.first] = {fit*1000, error*1000};
    if((TString(it.first.first).Contains("_PCB")) && TString(it.first.second).Contains("Guard") && !TString(it.first.first).Contains("Sensor")){
      j["PCBDist"][it.first.first] = it.first.first;
    }
    if((TString(it.first.first).Contains("_mid") || TString(it.first.first).Contains("_old") || TString(it.first.first).Contains("_new")) && TString(it.first.second).Contains("Guard") && TString(it.first.first).Contains("_Kapton")){
      j["KaptonMultiDist"][it.first.first] = it.first.first;
    }
    if(TString(it.first.first).Contains("Gaussian") && j["BinNum"].find(it.first.second) == j["BinNum"].end()){
      j["BinNum"][it.first.second] = binnum;
      binnum++;
    }
    if((TString(it.first.first).Contains("Sensor")) && TString(it.first.second).Contains("Guard")){
      j["Sensor"][it.first.first] = it.first.first;
    }
    if (!TString(it.first.first).Contains("_Kapton") && !TString(it.first.first).Contains("_PCB") && !TString(it.first.first).Contains("Custom") && !TString(it.first.first).Contains("Sensor") && TString(it.first.second).Contains("Guard")){
      j["Distribution"][it.first.first] = it.first.first;
    }
    if ((TString(it.first.first).Contains("Custom")) && TString(it.first.second).Contains("Guard")){
      j["Custom"][it.first.first] = it.first.first;
    }
  }
  for ( const auto & it : worst ) {
    double worst = it.second.first.first;
    double error = it.second.first.second;
    double prob  = it.second.second;

    //if(worst == NULL) worst = 0.0;
    //if(error == NULL) error = 0.0;
    //if(prob == NULL) prob = 0.0;
    j["Worst"][it.first.second][it.first.first] = {worst*1000, error*1000, prob};
    if(TString(it.first.second).Contains("Total")){
      j["TotalBadModules"][it.first.first] = it.first.first;
    }
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
  for(int i = 0; i < (signed)p.size(); i++){
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
  } else if(geo == "Three"){
    x[0] = -1*width/2;
    y[0] = -1*a_in/2;
    x[1] = x[0] + a_in*TMath::Sqrt(3);
    y[1] = y[0];
    x[2] = x[0] + a_in*TMath::Sqrt(3)/2.0;
    y[2] = y[0] - a_in/2.0;
  }

  vector<pair<double, double> > points;
  int nPoints = 6;
  if(geo == "Five" or geo == "Semi") nPoints = 5;
  else if(geo == "Half") nPoints = 4;
  else if(geo == "Three") nPoints = 3;
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
  for(int i = 0; i < (signed)points.size(); i++){
    points.at(i).first = points.at(i).first - center.first;
    points.at(i).second = points.at(i).second - center.second;
    pair<double, double> polar = GetPolar(TString(hc->GetName()), points.at(i).first, points.at(i).second, 0., 0.);
    radius.push_back(polar.first);
    angle.push_back(polar.second);
    if(polar.first < rMin) rMin = polar.first;
  }
  for(int i = 0; i < (signed)points.size(); i++) ratio.push_back(radius[i]/rMin);
  for(int theta = 0; theta != (signed)angle.size(); theta++){
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
        pair<double, double> cart = GetCartesian(TString(hc->GetName()), radius[first] + sign*a*ibin*ratio[first], angle[first], 0., 0.);
        xloop[3] = cart.first;
        yloop[3] = cart.second;
        cart = GetCartesian(TString(hc->GetName()), radius[second] + sign*a*ibin*ratio[second], angle[second], 0., 0.);
        xloop[2] = cart.first;
        yloop[2] = cart.second;
        if(count == 300 && neg == 1) break;
        if((geo == "Semi" || geo == "Half" || geo == "Three") && count == 200 && neg == 1) break;
        if(geo == "Three" && count == 100 && neg == 1) break;
        hc->AddBin(4, xloop, yloop);
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

double ScaleSide(TString geo, int side, double width){
  double newSide = width;

       if(geo == "Five"  && side == 1) newSide = width*TMath::Sqrt(3)/2;
  else if(geo == "Five"  && side == 2) newSide = width*3/4;
  else if(geo == "Semi"  && side == 1) newSide = GetWidthToA(width);
  else if(geo == "Semi"  && side == 2) newSide = GetWidthToA(width)*TMath::Sqrt(3)/2;
  else if(geo == "Half"  && side == 0) newSide = GetWidthToA(width)*3/2;
  else if(geo == "Half"  && side == 1) newSide = width*3/4;
  else if(geo == "Three" && side == 0) newSide = GetWidthToA(width)/2;
  else if(geo == "Three" && side > 0)  newSide = (GetWidthToA(width)*TMath::Sqrt(3))/4; 

  return newSide;
}

bool FillAllSides(TString geo, double max, TH2Poly* hc, double nom, 
		vector<double> rand, double shift_x, double shift_y, 
		double width_new, pair < double, double > cen_base,
		vector<double> comp = {}, double min = 0.,
		vector<double> comp2 = {}, 
                vector<double> comp_between = {}, pair<double, double> backside_x_err = make_pair(0.0, 0.0), pair<double, double> backside_y_err = make_pair(0.0, 0.0), TString coordAxis = ""){
  //Widths
  vector<pair<pair<double, double>, double>> location;
  double dx = (shift_x), dy = (shift_y);
  bool isNegative = false;
  pair<double, double> rOut;
  vector<double> forward = {0., TMath::Pi()/3, TMath::Pi()*2/3, TMath::Pi(), TMath::Pi()*4/3, TMath::Pi()*5/3};
  if(geo == "Five" || geo == "Semi") forward = {0., TMath::Pi()/2, TMath::Pi(), TMath::Pi()*4/3, TMath::Pi()*5/3};
  else if(geo == "Half") forward = {0., TMath::Pi()*2/3, TMath::Pi()*4/3, TMath::Pi()*5/3};
  else if(geo == "Three") forward = {TMath::Pi()/2, TMath::Pi()*4/3, TMath::Pi()*5/3};
  vector<int> index = {0, 1, 2, 0, 1, 2};
  if(geo == "Five" || geo == "Semi") index = {0, 1, 0, 2, 2};
  else if(geo == "Half") index = {0, 1, 0, 1};
  else if(geo == "Three") index = {0, 1, 2};
  double r_out = 0.;
  for(int iF = 0; iF != (signed)forward.size(); iF++){
    double secondary = ScaleSide(geo, index[iF], nom);
    double third = 0.;
    if(comp.size() > 0) secondary = comp[index[iF]];
    if(comp2.size() > 0) third = comp2[index[iF]];
 
    double r_out = 0., r_out1 = 0., r_out2 = 0., dr = 0.;
    if(TString(hc->GetName()).Contains("Backside Bond")){
      // Get test values
      r_out1 = (rand[index[iF]] - secondary)/2;
      if(debug && TString(hc->GetName()).Contains(whichOverlap)) cout << hc->GetName() << ": " << r_out1 << " = (" << rand[index[iF]] << " - " << secondary << ")/2" << endl;
      r_out2 = (third - rand[index[iF]])/2;
      if(debug && TString(hc->GetName()).Contains(whichOverlap)) cout << hc->GetName() << ": " << r_out2 << " = (" << third << " - " << rand[index[iF]] << ")/2" << endl;
      double r_out3 = (secondary - third)/2;
      if(debug && TString(hc->GetName()).Contains(whichOverlap)) cout << hc->GetName() << ": " << r_out3 << " = (" << secondary << " - " << third << ")/2" << endl;


      r_out = (rand[index[iF]] - secondary)/2;
      if(debug && TString(hc->GetName()).Contains(whichOverlap)) cout << hc->GetName() << ": " << r_out << " = (" << rand[index[iF]] << " - " << secondary << ")/2" << endl;
      dr = dx*TMath::Cos(forward[iF]) + dy*TMath::Sin(forward[iF]);
      double r_btw = (2*r_out + comp_between[index[iF]]);
      r_out = (r_btw - third)/2;
      if(debug && TString(hc->GetName()).Contains(whichOverlap)) cout << r_out << " = (" << r_out << " + " << comp_between[index[iF]] << " - " << third << ")/2" << endl;
    } else if(TString(hc->GetName()).Contains("base - kap")){
      // Get test values
      r_out1 = (rand[index[iF]] - secondary)/2;
      r_out2 = (nom - third)/2;
      //if(debug && TString(hc->GetName()).Contains(whichOverlap)){
      //  cout << hc->GetName() << ": " << r_out1 << " = (" << rand[index[iF]] << " - " << secondary << ")/2" << endl;
      //  cout << hc->GetName() << ": " << r_out2 << " = (" << rand[index[iF]] << " - " << third << ")/2" << endl;
      //}
      r_out = r_out1 + r_out2;
      dr = dx*TMath::Cos(forward[iF]) + dy*TMath::Sin(forward[iF]);
    } else {
      r_out = (rand[index[iF]] - secondary)/2;
      if(debug && TString(hc->GetName()).Contains(whichOverlap)) cout << hc->GetName() << ": " << r_out << " = (" << rand[index[iF]] << " - " << secondary << ")/2" << endl;
      dr = dx*TMath::Cos(forward[iF]) + dy*TMath::Sin(forward[iF]);
    }
    double dr_cen = cen_base.first*TMath::Cos(forward[iF]) + cen_base.second*TMath::Sin(forward[iF]);

    if(geo == "Five" && index[iF] == 1) {
      r_out += GetWidthToA(width_new)/2;
    } else if(geo == "Three" && index[iF] == 0) {
      r_out -= GetWidthToA(width_new)/2;
    } else if(!(geo == "Semi" && index[iF] == 1) && !(geo == "Half" && forward[iF] == TMath::Pi()*2/3)){
      r_out += width_new/2;
    }
    if(debug && TString(hc->GetName()).Contains(whichOverlap) && (r_out + dr - dr_cen - min - 0.7) > 0.254 && (r_out + dr - dr_cen - min - 0.7) < 0.258){ 
      cout << hc->GetName() << ": " << r_out + dr - dr_cen - min << " = " << r_out << " + " << dr << " + " << dr_cen << endl;
      cout << "shifts = (" << dx << ", " << dy << ")" << endl;
    }   

    rOut = GetCartesian(TString(hc->GetName()), r_out + dr - dr_cen - min, forward[iF], 0., 0.);
    hc->Fill(rOut.first, rOut.second);

    //Check Negative
    double r = TMath::Sqrt(rOut.first*rOut.first + rOut.second*rOut.second) - width_new/2;
    double dr_cen_base = cen_base.first*TMath::Cos(forward[iF]) + cen_base.second*TMath::Sin(forward[iF]);
         if(geo == "Five" && iF == 1) r += width_new/2 - GetWidthToA(width_new)/2;
    else if(geo == "Semi" && iF == 1) r += width_new/2;
    else if(geo == "half" && iF == 1) r += width_new/2;
    else if(geo == "Three" && iF == 0) r += width_new/2 + GetWidthToA(width_new)/2;
    r += dr_cen_base;

    if(r < 0) isNegative = true;
  }
  return isNegative;
}

void AddLineNomHex(double width, TString geo = "Full"){
   double a_in = GetWidthToA(width);
   // Add the bins
  vector<pair<double, double> >  points = GetPoints(geo, a_in, width);
  pair < double, double > center = getPolyCenter(points);
  for(int i = 0; i < (signed)points.size(); i++){
    points.at(i).first = points.at(i).first - center.first;
    points.at(i).second = points.at(i).second - center.second;
  }

  for(int iS = 0; iS < (signed)points.size(); iS++){
    int first = iS % points.size();
    int second = (iS + 1) % points.size();
    TLine *line = new TLine(points.at(first).first, points.at(first).second, points.at(second).first, points.at(second).second);
    line->SetLineColor(kBlack);
    line->Draw();
  }
}

void moduleTolerances(){
  lumistr = "";
  SetStyle();
  TDR_EXTRA_LABEL_ = "";
  TDR_EXTRA_LABEL_2 = "";

  //Colors
  colors.push_back(kBlue);
  colors.push_back(kGreen +3);
  colors.push_back(kRed);
  colors.push_back(kCyan);
  colors.push_back(kViolet-4);
  colors.push_back(kSpring-9);
  colors.push_back(kMagenta);
  colors.push_back(kAzure);

  vector<string> Geometry = {"Full", "Five", "Semi", "Half", "Three"};
  vector<string> Dist = {"Gaussian_Kaptonminus0_oldSensor", "Landau", "Flat", "Gaussian_Kaptonminus0_newSensor", "Gaussian_Kaptonminus0_midSensor", 
			 "CustomGaus", "CustomLandau", "CustomFlat", 
			 "Gaussian_PCBplus25", "Gaussian_PCBplus50", "Gaussian_PCBplus75", "Gaussian_PCBminus25", "Gaussian_PCBminus50", "Gaussian_PCBminus75", 
			 "Gaussian_Kaptonplus25_oldSensor", "Gaussian_Kaptonplus50_oldSensor", "Gaussian_Kaptonplus75_oldSensor", 
			 "Gaussian_Kaptonminus25_oldSensor", "Gaussian_Kaptonminus50_oldSensor", "Gaussian_Kaptonminus75_oldSensor", "Gaussian_Kaptonminus150_oldSensor",  "Gaussian_Kaptonminus170_oldSensor",
			 "Gaussian_Kaptonplus25_newSensor", "Gaussian_Kaptonplus50_newSensor", "Gaussian_Kaptonplus75_newSensor",                                         
			 "Gaussian_Kaptonminus25_newSensor", "Gaussian_Kaptonminus50_newSensor", "Gaussian_Kaptonminus75_newSensor", "Gaussian_Kaptonminus150_newSensor",  "Gaussian_Kaptonminus170_newSensor",
			 "Gaussian_Kaptonplus25_midSensor", "Gaussian_Kaptonplus50_midSensor", "Gaussian_Kaptonplus75_midSensor",                                         
			 "Gaussian_Kaptonminus25_midSensor", "Gaussian_Kaptonminus50_midSensor", "Gaussian_Kaptonminus75_midSensor", "Gaussian_Kaptonminus150_midSensor",  "Gaussian_Kaptonminus170_midSensor",
			 "Gaussian_PCBplus25_oldSensor", "Gaussian_newSensor", "Gaussian_PCBplus25_newSensor", };
//  Dist = {"Gaussian_Kaptonminus0_oldSensor", "Gaussian_Kaptonminus0_newSensor", "Gaussian_Kaptonminus0_midSensor", 
//			 "Gaussian_Kaptonminus25_oldSensor", "Gaussian_Kaptonminus25_newSensor", "Gaussian_Kaptonminus25_midSensor",  
//			 "Gaussian_Kaptonminus50_oldSensor", "Gaussian_Kaptonminus50_newSensor", "Gaussian_Kaptonminus50_midSensor",  
//			 "Gaussian_Kaptonminus75_oldSensor", "Gaussian_Kaptonminus75_newSensor", "Gaussian_Kaptonminus75_midSensor",  
//			 "Gaussian_Kaptonminus150_oldSensor", "Gaussian_Kaptonminus150_newSensor", "Gaussian_Kaptonminus150_midSensor",  
//			 "Gaussian_Kaptonminus170_oldSensor", "Gaussian_Kaptonminus170_newSensor", "Gaussian_Kaptonminus170_midSensor",  
//			 "Gaussian_Kaptonplus25_oldSensor", "Gaussian_Kaptonplus25_newSensor", "Gaussian_Kaptonplus25_midSensor",  
//			 "Gaussian_Kaptonplus50_oldSensor", "Gaussian_Kaptonplus50_newSensor", "Gaussian_Kaptonplus50_midSensor",  
//			 "Gaussian_Kaptonplus75_oldSensor", "Gaussian_Kaptonplus75_newSensor", "Gaussian_Kaptonplus75_midSensor"};

//  vector<string> Dist = {"Gaussian", "Gaussian_newSensor", "Gaussian_midSensor", 
//			 "Gaussian_Kaptonplus25", "Gaussian_Kaptonplus50", "Gaussian_Kaptonplus75", 
//			 "Gaussian_Kaptonminus25", "Gaussian_Kaptonminus50", "Gaussian_Kaptonminus75", 
//			 "Gaussian_Kaptonplus25_newSensor", "Gaussian_Kaptonplus50_newSensor", "Gaussian_Kaptonplus75_newSensor", 
//			 "Gaussian_Kaptonminus25_newSensor", "Gaussian_Kaptonminus50_newSensor", "Gaussian_Kaptonminus75_newSensor",
//			 "Gaussian_Kaptonplus25_midSensor", "Gaussian_Kaptonplus50_midSensor", "Gaussian_Kaptonplus75_midSensor", 
//			 "Gaussian_Kaptonminus25_midSensor", "Gaussian_Kaptonminus50_midSensor", "Gaussian_Kaptonminus75_midSensor", };
  Dist = {"Gaussian"};
  Geometry = {"Full"};
  for(auto &geo_str : Geometry){
    vector< pair< pair< string, string>, pair< pair< double, double >, double > > > worst_values;
    vector< pair< pair< string, string>, pair< double, double > > > fit_values;
    TString geo = TString(geo_str);
    outputDir = baseDir + "/Test" + geo;
    gSystem->mkdir(outputDir, true);

    int nPeaks = 1;
         if(geo == "Five" || geo == "Semi") nPeaks = 3;
    else if(geo == "Half" || geo == "Three") nPeaks = 2;

    map<TString, TH1D*> overlap_1D_integrate_minus;
    for(auto &type_str: Dist){
      TString type = TString(type_str);
      gSystem->mkdir(outputDir+"/"+type, true);
      //Define Constant Widths
      double pcb_w_const = 166.64, sensor_w_const = 166.57, kapton_w_const = 166.94, baseplate_w_const = 166.94;

      // Input Parameters units are Millimeters
      double baseplate_center = 3.05, baseplate_err = 0.05;
      double baseplate_pinhole_err = 0.100, baseplate_pinhole_theta = 2*TMath::Pi();
      double pcb_w = pcb_w_const, sensor_w = sensor_w_const, kapton_w = kapton_w_const, baseplate_w = baseplate_w_const;
      double pcb_w_err = 0.1, sensor_w_err = 0.02, kapton_w_err = 0.05, baseplate_w_err = 0.05;
      double PcbToSen_shift_r = 0.04, PcbToSen_shift_theta = 2*TMath::Pi();
      double sensor_shift_x= 0.500, sensor_shift_y = 0.050, sensor_shift_theta = 0.050; //theta is radians
      double baseValue_total = sensor_w;
      double base_kap_min = 0.100;
      double gold_min = 0.25, gold_size = 1.500, backside_min = 0.100;
      pair<double, double> backside_x_err = make_pair(0.2133, 0.1769);
      pair<double, double> backside_y_err = make_pair(0.500, 0.000);

      double nBadModules = 0.;

      int max = 300000;
      //max = 30000;

      int nbins = 650;
      if(geo == "Half" || geo == "Semi") nbins = 650;
      
      //nbins = 5;
      double step = 0.002;
      double width_new = 1.4;
      double axis = width_new/2 + nbins*step + 0.1;
      double a_new = GetWidthToA(width_new);

      double sigma = 1., mean = 1., kapton_sigma = 1., kapton_mean = 1.;
      if(type.Contains("Custom")){
        vector < pair <double, pair<double, double> > > custom = newShapeInput();
        sigma = custom[0].second.second/baseplate_w_err;
        mean = custom[0].second.first/baseplate_w;
        kapton_sigma = custom[1].second.second/kapton_w_err;
        kapton_mean = custom[1].second.first/kapton_w;
        baseValue_total *= mean;
      }

      cout << "baseplate: " << mean << " +/- " << sigma << endl;
      cout << "kapton: " << kapton_mean << " +/- " << kapton_sigma << endl;

      map<TString, pair<double, double> > center;
      vector < pair < double, double > > points = GetPoints(geo, GetWidthToA(width_new), width_new);
      center.emplace("new", getPolyCenter(points));

      if(geo == "Semi"){
        axis += 0.5;
        nbins += 50;
      } else if(geo == "Half"){
        axis += 1.5;
        nbins += 150;
      } else if(geo == "Three"){
        axis += 1.4;
        nbins += 50;
      }

       
      if(type.Contains("otherCenter")){
        baseplate_err = 0.029;
      }
      if(type.Contains("PCB")){
        pcb_w = pcb_w_const;
        TString buffer = type;
        TString pcb_new_shift_str = buffer.ReplaceAll("Gaussian_PCBplus", "");
        pcb_new_shift_str.ReplaceAll("Gaussian_PCBminus", "");
        float pcb_new_shift = pcb_new_shift_str.Atof()/1000;
       
        if(type.Contains("plus")) pcb_w  = pcb_w + pcb_new_shift;
        if(type.Contains("minus")) pcb_w = pcb_w -  pcb_new_shift;
        cout << "pcb_w: " << pcb_w << "pcb_new_shift: " << pcb_new_shift << endl;
      } else pcb_w = pcb_w_const;
      if(type.Contains("Kapton")){
        kapton_w = kapton_w_const;
        TString buffer = type;
        TString kapton_new_shift_str = buffer.ReplaceAll("Gaussian_Kaptonplus", "");
        kapton_new_shift_str.ReplaceAll("Gaussian_Kaptonminus", "");
        float kapton_new_shift = kapton_new_shift_str.Atof()/1000;
       
        if(type.Contains("plus")) kapton_w  = kapton_w + kapton_new_shift;
        if(type.Contains("minus")) kapton_w = kapton_w -  kapton_new_shift;
        cout << "kapton_w: " << kapton_w << "kapton_new_shift: " << kapton_new_shift << endl;
      } else kapton_w = kapton_w_const;
      if(type.Contains("newSensor")){ 
        sensor_shift_x = 0.050; 
        sensor_shift_y = 0.050;
      } else if(type.Contains("midSensor")) {
	sensor_shift_x = 0.100;
	sensor_shift_y = 0.050;
      } else if(type.Contains("oldSensor")) {
	sensor_shift_x = 0.200;
	sensor_shift_y = 0.050;
      }
      vector<double> pcb, kapton, baseplate, sensor, pcb_bond, kapton_bond, sensor_bond, pcb_backside, kapton_backside, sensor_backside, sensor_backside_between;
      double KapToBas_x_shift = 0., SenToKap_x_shift = 0., KapToBas_y_shift = 0., SenToKap_y_shift = 0.;
      map<TString, TH2Poly*> components {
      	{"pcb", new TH2Poly("pcb_2d", ";" + geo + " " + type + " Width [mm];Width [mm];Components", -1*axis, axis, -1*axis, axis)},
      	{"kap", new TH2Poly("kap_2d", ";" + geo + " " + type + " Width [mm];Width [mm];Components", -1*axis, axis, -1*axis, axis)},
      	{"sen", new TH2Poly("sen_2d", ";" + geo + " " + type + " Width [mm];Width [mm];Components", -1*axis, axis, -1*axis, axis)},
      	{"bas", new TH2Poly("bas_2d", ";" + geo + " " + type + " Width [mm];Width [mm];Components", -1*axis, axis, -1*axis, axis)},
      };
      map<TString, TH2Poly*> overlap {
      	{"kap_pcb_hist",        new TH2Poly("Shield Bond",       ";" + geo + " " + type + " Diff Width [mm];Diff Width [mm];Overlaps", -1*axis, 1*axis, -1*axis, 1*axis)}, 			
      	{"sen_pcb_hist",        new TH2Poly("Guard Bond",        ";" + geo + " " + type + " Diff Width [mm];Diff Width [mm];Overlaps", -1*axis, 1*axis, -1*axis, 1*axis)}, 			
      	{"sen_pcb_kap_x_hist",  new TH2Poly("Backside Bond X Pos",";" + geo + " " + type + " Diff Width [mm];Diff Width [mm];Overlaps", -1*axis, 1*axis, -1*axis, 1*axis)}, 			
      	{"sen_pcb_kap_y_hist",  new TH2Poly("Backside Bond Y Pos",";" + geo + " " + type + " Diff Width [mm];Diff Width [mm];Overlaps", -1*axis, 1*axis, -1*axis, 1*axis)}, 			
      	{"pcb_bas_stack_hist",  new TH2Poly("Overlap pcb - base",";" + geo + " " + type + " Diff Width [mm];Diff Width [mm];Overlaps", -1*axis, 1*axis, -1*axis, 1*axis)},    
      	{"sen_bas_stack_hist",  new TH2Poly("Overlap sen - base",";" + geo + " " + type + " Diff Width [mm];Diff Width [mm];Overlaps", -1*axis, 1*axis, -1*axis, 1*axis)},
      	{"pcb_kap_stack_hist",  new TH2Poly("Overlap pcb - kap", ";" + geo + " " + type + " Diff Width [mm];Diff Width [mm];Overlaps", -1*axis, 1*axis, -1*axis, 1*axis)},
      	{"bas_kap_stack_hist",  new TH2Poly("Overlap base - kap",";" + geo + " " + type + " Diff Width [mm];Diff Width [mm];Overlaps", -1*axis, 1*axis, -1*axis, 1*axis)},
      	{"sen_kap_stack_hist",  new TH2Poly("Overlap sen - kap", ";" + geo + " " + type + " Diff Width [mm];Diff Width [mm];Overlaps", -1*axis, 1*axis, -1*axis, 1*axis)},
      	{"sen_pcb_stack_hist",  new TH2Poly("Overlap sen - pcb", ";" + geo + " " + type + " Diff Width [mm];Diff Width [mm];Overlaps", -1*axis, 1*axis, -1*axis, 1*axis)},
      };
     
      TH2D* corr = new TH2D("correlation", "correlation" + geo, overlap.size(), 0, overlap.size(), overlap.size(), 0, overlap.size());

      map<TString, string> nameMap {
        {"pcb",                 "PCB Width"},
        {"kap",                 "Kapton Width"},
        {"sen",                 "Sensor Width"},
        {"bas",                 "Baseplate Width"},
        {"kap_pcb_hist",        "Shield Bond"},
        {"sen_pcb_hist",        "Guard Bond"},
        {"sen_pcb_kap_x_hist",  "Backside Bond X Position"},
        {"sen_pcb_kap_y_hist",  "Backside Bond Y Position"},
        {"pcb_bas_stack_hist",  "PCB to Base"},
        {"sen_bas_stack_hist",  "Sen to Base"},
        {"pcb_kap_stack_hist",  "PCB to Kap"},
        {"bas_kap_stack_hist",  "Base to Kap"},
        {"sen_kap_stack_hist",  "Sen to Kap"},
        {"sen_pcb_stack_hist",  "Sen to PCB"},
      };

      for(std::map<TString,TH2Poly*>::iterator iter = components.begin(); iter != components.end(); ++iter)
        HoneycombCustom(geo, iter->second, step, width_new, a_new, nbins);
      for(std::map<TString,TH2Poly*>::iterator iter = overlap.begin(); iter != overlap.end(); ++iter)
        HoneycombCustom(geo, iter->second, step, width_new, a_new, nbins);

      cout << type << endl;
      for(auto i = 0; i != max; i++){
        pcb.clear();
        kapton.clear();
        baseplate.clear();
        sensor.clear();
        if(i == 0){
          cout << "baseplate: " << baseplate_w*mean << " +/- " << baseplate_w_err*sigma << endl;
          cout << "pcb: " << pcb_w*kapton_mean << " +/- " << pcb_w_err*kapton_sigma << endl;
          cout << "kapton: " << kapton_w*kapton_mean << " +/- " << kapton_w_err*kapton_sigma << endl;
          cout << "sensor: " << sensor_w*kapton_mean << " +/- " << sensor_w_err*kapton_sigma << endl;
          cout << "center: " << baseplate_err << endl;
          cout << "PCB shift: " << PcbToSen_shift_r<< endl;
          cout << "Sensor shift: " << sensor_shift_x << endl;
	}
	//Get widths for each component and each side
        int sides = 3;
        if(geo == "Half") sides = 2;
        for(auto iS = 0; iS != sides; iS++){
          double pcb_w_new = ScaleSide(geo, iS, pcb_w*kapton_mean);
          double kapton_w_new = ScaleSide(geo, iS, kapton_w*kapton_mean);
          double baseplate_w_new = ScaleSide(geo, iS, baseplate_w*mean);
          double sensor_w_new = ScaleSide(geo, iS, sensor_w*kapton_mean);
          pcb.push_back(getRandomValue(pcb_w_new, pcb_w_err*kapton_sigma, type));
          kapton.push_back(getRandomValue(kapton_w_new, kapton_w_err*kapton_sigma, type));
          baseplate.push_back(getRandomValue(baseplate_w_new, baseplate_w_err*sigma, type));
          sensor.push_back(getRandomValue(sensor_w_new, sensor_w_err*kapton_sigma, type));
          kapton_bond.push_back(kapton.at(iS) - (gold_size - 0.650)*2);
          sensor_bond.push_back(sensor.at(iS) - (gold_size - 0.730)*2);
          pcb_bond.push_back(pcb.at(iS) - gold_size*2);
          sensor_backside.push_back(getCutoffWidth(sensor.at(iS), 15.590));
          pcb_backside.push_back(getCutoffWidth(pcb.at(iS), 15.590, 0.540));
          kapton_backside.push_back(kapton.at(iS) - 0.101*2);
          sensor_backside_between.push_back(getCutoffInBetween(sensor.at(iS), 15.590));
        }//get three random sides
        //Get shift in x and y, soon add theta also
        pair<double, double> KapToBas_shift = GetCartesian("", getRandomValue(0, baseplate_err, type), getRandomValue(0, PcbToSen_shift_theta, type), 0., 0.);
        pair<double, double> PcbToSen_shift = GetCartesian("", getRandomValue(0, PcbToSen_shift_r, type), getRandomValue(0, PcbToSen_shift_theta, type), 0., 0.);
        SenToKap_x_shift = getRandomValue(0, sensor_shift_x, type);
        SenToKap_y_shift = getRandomValue(0, sensor_shift_y, type);

        //Component values
        int isNegative = 0x0;
        FillAllSides(geo, double(max), components["bas"], baseValue_total, baseplate, KapToBas_shift.first,     KapToBas_shift.second,      width_new, center["new"]);
        FillAllSides(geo, double(max), components["kap"], baseValue_total,    kapton, KapToBas_shift.first,     KapToBas_shift.second,      width_new, center["new"]);
        FillAllSides(geo, double(max), components["sen"], baseValue_total,    sensor, SenToKap_x_shift,     SenToKap_y_shift,      width_new, center["new"]);
        FillAllSides(geo, double(max), components["pcb"], baseValue_total, 	 pcb, PcbToSen_shift.first, PcbToSen_shift.second, width_new, center["new"]);
        //Overlap Values
        isNegative |= FillAllSides(geo, double(max), overlap["sen_pcb_stack_hist"], baseValue_total,       pcb, PcbToSen_shift.first, PcbToSen_shift.second, width_new, center["new"], sensor) << 0;
        isNegative |= FillAllSides(geo, double(max), overlap["pcb_bas_stack_hist"], baseValue_total, baseplate, PcbToSen_shift.first - (KapToBas_shift.first + SenToKap_x_shift), PcbToSen_shift.second - (KapToBas_shift.second + SenToKap_y_shift), width_new, center["new"], pcb) << 1;
        isNegative |= FillAllSides(geo, double(max), overlap["sen_bas_stack_hist"], baseValue_total, baseplate, KapToBas_shift.first - SenToKap_x_shift, KapToBas_shift.second - SenToKap_y_shift, width_new, center["new"], sensor) << 2;
        isNegative |= FillAllSides(geo, double(max), overlap["pcb_kap_stack_hist"], baseValue_total,    kapton, SenToKap_x_shift - PcbToSen_shift.first, SenToKap_y_shift - PcbToSen_shift.second, width_new, center["new"], pcb) << 3;
        isNegative |= FillAllSides(geo, double(max), overlap["bas_kap_stack_hist"], kapton_w_const,     kapton, KapToBas_shift.first + SenToKap_x_shift, KapToBas_shift.second + SenToKap_y_shift, width_new, center["new"], baseplate, base_kap_min, sensor) << 4;
        isNegative |= FillAllSides(geo, double(max), overlap["sen_kap_stack_hist"], baseValue_total,    kapton, SenToKap_x_shift, SenToKap_y_shift, width_new, center["new"], sensor) << 5;
        isNegative |= FillAllSides(geo, double(max), overlap["kap_pcb_hist"],       baseValue_total, kapton_bond, SenToKap_x_shift - PcbToSen_shift.first, SenToKap_y_shift - PcbToSen_shift.second, width_new, center["new"], pcb_bond, gold_min) << 6;
        isNegative |= FillAllSides(geo, double(max), overlap["sen_pcb_hist"],       baseValue_total, sensor_bond, PcbToSen_shift.first, PcbToSen_shift.second, width_new, center["new"], pcb_bond, gold_min) << 7;
        isNegative |= FillAllSides(geo, double(max), overlap["sen_pcb_kap_x_hist"],   baseValue_total,    pcb_backside, SenToKap_x_shift - PcbToSen_shift.first, SenToKap_y_shift - PcbToSen_shift.second, width_new, center["new"], sensor_backside, backside_min, kapton_backside, sensor_backside_between, backside_x_err, backside_y_err, "X") << 8;
        isNegative |= FillAllSides(geo, double(max), overlap["sen_pcb_kap_y_hist"],   baseValue_total,    pcb_backside, SenToKap_x_shift - PcbToSen_shift.first, SenToKap_y_shift - PcbToSen_shift.second, width_new, center["new"], sensor_backside, backside_min, kapton_backside, sensor_backside_between, backside_x_err, backside_y_err, "Y") << 9;

        if(isNegative > 0){
          nBadModules += 30000/max;
          for(int iX = 0; iX != overlap.size(); iX++){
            for(int iY = 0; iY != overlap.size() && iY <= iX; iY++){
              int value = ((isNegative >> iX) & 1) & ((isNegative >> iY) & 1);
              if(value > 0){
                corr->Fill(iX, iY, double(value*30000/max));
                if(iY != iX) corr->Fill(iY, iX, double(value*30000/max));
              }
            }
          }
        }
      }//for all modules

      vector<TString> binLabels = {
        TString(nameMap["sen_pcb_stack_hist"]), 
        TString(nameMap["pcb_bas_stack_hist"]), 
        TString(nameMap["sen_bas_stack_hist"]), 
        TString(nameMap["pcb_kap_stack_hist"]), 
        TString(nameMap["bas_kap_stack_hist"]), 
        TString(nameMap["sen_kap_stack_hist"]), 
        TString(nameMap["kap_pcb_hist"]),       
        TString(nameMap["sen_pcb_hist"]),       
        TString(nameMap["sen_pcb_kap_x_hist"]), 
        TString(nameMap["sen_pcb_kap_y_hist"])
        };
      for(std::map<TString,TH2Poly*>::iterator iter = components.begin(); iter != components.end(); ++iter)
        print2DPlots(iter->second, geo, TString(nameMap[iter->first]), type+"/"+geo + "_" + type + "_Width_" + iter->first, width_new);
      for(std::map<TString,TH2Poly*>::iterator iter = overlap.begin(); iter != overlap.end(); ++iter)
        print2DPlots(iter->second, geo,	TString(nameMap[iter->first]), type+"/"+geo + "_" + type + "_Diff_Width_" + iter->first, width_new);

      setBinLabelsYX(corr, binLabels);
      print2DPlots(corr, geo, "Bad Module Overlaps", type+"/"+geo + "_" + type + "_corr");

      map<TString, TH1D*> overlap_1D, overlap_1D_integrate;
      int iC = 0;
      for(std::map<TString,TH2Poly*>::iterator iter = overlap.begin(); iter != overlap.end(); ++iter){
        vector<TH1*> hist = IntegrateHex(iter->second, geo, type, center["new"], width_new, axis, step, nPeaks);
        prepHists(hist, true, false, false, {colors[iC], colors[iC], colors[iC], colors[iC]});
        vector<TH1*> intHist = getIntegratedHist(hist, false, false, false);
        for(auto iP = 0; iP != nPeaks + 1; iP++)
          overlap_1D.insert(pair<TString, TH1D*>(iter->first + "_Peak" + to_string(iP), (TH1D*)hist[iP]));
        for(auto iP = 0; iP != nPeaks + 1; iP++){
          overlap_1D_integrate.insert(pair<TString, TH1D*>(iter->first + "_integrate_Peak" + to_string(iP), (TH1D*)intHist[iP]));
          if(TString(iter->first).Contains("bas_kap_stack_hist") && (type.Contains("minus") || type.Contains("plus"))) 
            overlap_1D_integrate_minus.insert(pair<TString, TH1D*>(iter->first + type + "_integrate_Peak" + to_string(iP), (TH1D*)intHist[iP]));
        }
        iC++;
      }

      TString name = "";
      for(auto iP = 0; iP != nPeaks + 1; iP++){
        TString pName = "_Peak" + to_string(iP);
        auto leg_diff = prepLegends({overlap_1D["pcb_bas_stack_hist" + pName], overlap_1D["pcb_kap_stack_hist" + pName], overlap_1D["bas_kap_stack_hist" + pName]}, 
                                    {overlap_1D["pcb_bas_stack_hist" + pName]->GetName(), overlap_1D["pcb_kap_stack_hist" + pName]->GetName(), overlap_1D["bas_kap_stack_hist" + pName]->GetName()}, "P");
        leg_diff->SetTextSize(0.03);
        leg_diff->SetY1NDC(leg_diff->GetY2NDC() - 0.3);
        TCanvas* diff = drawCompMatt({overlap_1D["pcb_bas_stack_hist" + pName], overlap_1D["pcb_kap_stack_hist" + pName], overlap_1D["bas_kap_stack_hist" + pName]}, leg_diff, 0.00001, nullptr, "hist", true);
        name = geo + "_" + type + "_diff" + pName;
        diff->Update();
        diff->SetTitle(name);
        diff->Print(outputDir+"/"+type+"/"+name+".pdf");
        //integral
        diff = drawCompMatt({overlap_1D_integrate["pcb_bas_stack_hist_integrate" + pName], overlap_1D_integrate["pcb_kap_stack_hist_integrate" + pName], overlap_1D_integrate["bas_kap_stack_hist_integrate" + pName]}, leg_diff, -1., nullptr, "hist", true);
        name = geo + "_" + type + "_diff_integral" + pName;
        diff->Update();
        diff->SetTitle(name);
        diff->Print(outputDir+"/"+type+"/"+name+".pdf");

        auto leg_sen_diff = prepLegends({overlap_1D["sen_bas_stack_hist" + pName], overlap_1D["sen_pcb_stack_hist" + pName], overlap_1D["sen_kap_stack_hist" + pName]}, 
                                    {overlap_1D["sen_bas_stack_hist" + pName]->GetName(), overlap_1D["sen_pcb_stack_hist" + pName]->GetName(), overlap_1D["sen_kap_stack_hist" + pName]->GetName()}, "P");
        leg_sen_diff->SetTextSize(0.03);
        leg_sen_diff->SetY1NDC(leg_sen_diff->GetY2NDC() - 0.3);
        TCanvas* sen_diff = drawCompMatt({overlap_1D["sen_bas_stack_hist" + pName], overlap_1D["sen_pcb_stack_hist" + pName], overlap_1D["sen_kap_stack_hist" + pName]}, leg_sen_diff, 0.001, nullptr, "hist", true);
        name = geo + "_" + type + "_sen_diff" + pName;
        sen_diff->Update();
        sen_diff->SetTitle(name);
        sen_diff->Print(outputDir+"/"+type+"/"+name+".pdf");
        //integral
        sen_diff = drawCompMatt({overlap_1D_integrate["sen_bas_stack_hist_integrate" + pName], overlap_1D_integrate["sen_pcb_stack_hist_integrate" + pName], overlap_1D_integrate["sen_kap_stack_hist_integrate" + pName]}, leg_sen_diff, -1., nullptr, "hist", true);
        name = geo + "_" + type + "_sen_diff_integral" + pName;
        sen_diff->Update();
        sen_diff->SetTitle(name);
        sen_diff->Print(outputDir+"/"+type+"/"+name+".pdf");

        auto leg_bond_diff = prepLegends({overlap_1D["sen_pcb_hist" + pName], overlap_1D["kap_pcb_hist" + pName], overlap_1D["sen_pcb_kap_x_hist" + pName], overlap_1D["sen_pcb_kap_y_hist" + pName]},
                                    {overlap_1D["sen_pcb_hist" + pName]->GetName(), overlap_1D["kap_pcb_hist" + pName]->GetName(), overlap_1D["sen_pcb_kap_x_hist" + pName]->GetName(), overlap_1D["sen_pcb_kap_y_hist" + pName]->GetName()}, "P");
        leg_bond_diff->SetTextSize(0.03);
        leg_bond_diff->SetY1NDC(leg_bond_diff->GetY2NDC() - 0.3);
        TCanvas* bond_diff = drawCompMatt({overlap_1D["sen_pcb_hist" + pName], overlap_1D["kap_pcb_hist" + pName], overlap_1D["sen_pcb_kap_x_hist" + pName], overlap_1D["sen_pcb_kap_y_hist" + pName]}, leg_bond_diff, 0.001, nullptr, "hist", true);
        name = geo + "_" + type + "_bond_diff" + pName;
        bond_diff->Update();
        bond_diff->SetTitle(name);
        bond_diff->Print(outputDir+"/"+type+"/"+name+".pdf");
        //integral
        bond_diff = drawCompMatt({overlap_1D_integrate["sen_pcb_hist_integrate" + pName], overlap_1D_integrate["kap_pcb_hist_integrate" + pName], overlap_1D_integrate["sen_pcb_kap_x_hist_integrate" + pName], overlap_1D_integrate["sen_pcb_kap_y_hist_integrate" + pName]}, leg_bond_diff, 0.001, nullptr, "hist", true);
        name = geo + "_" + type + "_bond_diff_integral" + pName;
        bond_diff->Update();
        bond_diff->SetTitle(name);
        bond_diff->Print(outputDir+"/"+type+"/"+name+".pdf");
      }


      int sides = 1;
      vector<double> pcb_nom, kapton_nom, baseplate_nom, sensor_nom, pcb_bond_nom, kapton_bond_nom, sensor_bond_nom, pcb_backside_nom, sensor_backside_nom, kapton_backside_nom, sensor_backside_between_nom;
           if(geo == "Half" || geo == "Three") sides = 2;                                                               
      else if(geo == "Five" || geo == "Semi")  sides = 3;                                                               
      for(auto iS = 0; iS != sides; iS++){                                                                              
        string peak = "";
             if(iS == 0) peak = "_Peak1";
        else if(iS == 1) peak = "_Peak2";
        else if(iS == 2) peak = "_Peak3";

        double pcb_w_new = ScaleSide(geo, iS, pcb_w);
        double kapton_w_new = ScaleSide(geo, iS, kapton_w);
        double baseplate_w_new = ScaleSide(geo, iS, baseplate_w);
        double sensor_w_new = ScaleSide(geo, iS, sensor_w);
        pcb_nom.push_back(pcb_w_new);
        kapton_nom.push_back(kapton_w_new);
        baseplate_nom.push_back(baseplate_w_new);
        sensor_nom.push_back(sensor_w_new);
        kapton_bond_nom.push_back(kapton_w_new - (gold_size - 0.650)*2);
        sensor_bond_nom.push_back(sensor_w_new - (gold_size - 0.730)*2);
        pcb_bond_nom.push_back(pcb_w_new - gold_size*2);
        sensor_backside_nom.push_back(getCutoffWidth(sensor.at(iS), 15.590));
        pcb_backside_nom.push_back(getCutoffWidth(pcb.at(iS), 15.590, 0.540));
        kapton_backside_nom.push_back(kapton.at(iS) - 0.101*2);
        sensor_backside_between_nom.push_back(getCutoffInBetween(sensor.at(iS), 15.590));


        fit_values.push_back(make_pair(make_pair(type_str + peak + "_Nominal", "Guard Bond"), 	make_pair((sensor_bond_nom[iS] - pcb_bond_nom[iS])/2 - gold_min, 0.0))); 
        fit_values.push_back(make_pair(make_pair(type_str + peak + "_Nominal", "Shield Bond"), 	make_pair((kapton_bond_nom[iS] - pcb_bond_nom[iS])/2 - gold_min, 0.0)));
        //fit_values.push_back(make_pair(make_pair(type_str + peak + "_Nominal", "Backside Bond"),
        //    make_pair((pcb_bond_nom[iS] - sensor_bond_nom[iS] + sensor_backside_between_nom[iS] - kapton_backside_nom[iS])/2, 0.0)));
        fit_values.push_back(make_pair(make_pair(type_str + peak + "_Nominal", "PCB to Base"), 	make_pair((baseplate_nom[iS] - pcb_nom[iS])/2, 0.0)));
        fit_values.push_back(make_pair(make_pair(type_str + peak + "_Nominal", "Sen to Base"),	make_pair((baseplate_nom[iS] - sensor_nom[iS])/2, 0.0)));
        fit_values.push_back(make_pair(make_pair(type_str + peak + "_Nominal", "PCB to Kap"), 	make_pair((kapton_nom[iS] - pcb_nom[iS])/2, 0.0)));
        fit_values.push_back(make_pair(make_pair(type_str + peak + "_Nominal", "Base to Kap"), 	make_pair((kapton_nom[iS] - baseplate_nom[iS])/2 + 0.185 - base_kap_min, 0.0)));
        fit_values.push_back(make_pair(make_pair(type_str + peak + "_Nominal", "Sen to Kap"), 	make_pair((kapton_nom[iS] - sensor_nom[iS])/2, 0.0)));
        fit_values.push_back(make_pair(make_pair(type_str + peak + "_Nominal", "Sen to PCB"), 	make_pair((pcb_nom[iS] - sensor_nom[iS])/2, 0.0)));
      }

      for(std::map<TString,TH1D*>::iterator iter = overlap_1D.begin(); iter != overlap_1D.end(); ++iter){
        if(TString(iter->first).Contains("_Peak0")) continue;

        TString name_buff = iter->first;
        name_buff.ReplaceAll("_Peak1", "");
        name_buff.ReplaceAll("_Peak2", "");
        name_buff.ReplaceAll("_Peak3", "");
        string name = nameMap[name_buff];
        string peak = "";
             if(iter->first.Contains("Peak1")) peak = "_Peak1";
        else if(iter->first.Contains("Peak2")) peak = "_Peak2";
        else if(iter->first.Contains("Peak3")) peak = "_Peak3";

        //Get Fit
        TString dist_type = "gausn";
        if(type.Contains("Land")) dist_type = "landau";

        iter->second->Fit(dist_type);
        TF1* avg_param = (TF1*)iter->second->GetFunction(dist_type);
        if(!isnan(avg_param->GetParameter(1))) fit_values.push_back(make_pair(make_pair(type_str + peak, name), make_pair(avg_param->GetParameter(1), avg_param->GetParameter(2))));
        else fit_values.push_back(make_pair(make_pair(type_str + peak, name), make_pair(0., 0.)));
      }

      for(std::map<TString,TH1D*>::iterator iter = overlap_1D_integrate.begin(); iter != overlap_1D_integrate.end(); ++iter){
        if(TString(iter->first).Contains("_Peak0")) continue;

        vector<double> max_tol = findHighXbin(iter->second);
        TString name_buff = iter->first;
        name_buff.ReplaceAll("_integrate_Peak1", "");
        name_buff.ReplaceAll("_integrate_Peak2", "");
        name_buff.ReplaceAll("_integrate_Peak3", "");
        string name = nameMap[name_buff];
        string peak = "";
             if(iter->first.Contains("Peak1")) peak = "_Peak1";
        else if(iter->first.Contains("Peak2")) peak = "_Peak2";
        else if(iter->first.Contains("Peak3")) peak = "_Peak3";

        worst_values.push_back(make_pair(make_pair(type_str + peak, name), make_pair(make_pair(max_tol[0]/1000, max_tol[1]/1000), max_tol[2])));
      }
      worst_values.push_back(make_pair(make_pair(type_str, "Total Modules"), make_pair(make_pair(0., 0.), nBadModules)));

      TFile *outFile = new TFile(outputDir+"/"+type+"/"+name+".root", "RECREATE");
      for(std::map<TString,TH2Poly*>::iterator iter = components.begin(); iter != components.end(); ++iter){
        iter->second->Write();
        delete iter->second;
      }
      for(std::map<TString,TH1D*>::iterator iter = overlap_1D.begin(); iter != overlap_1D.end(); ++iter){
        iter->second->Write();
        delete iter->second;
      }
      outFile->Close();
    }//Different distribution type

    for(auto iP = 0; iP != nPeaks + 1; iP++){
      TString pName = "_Peak" + to_string(iP);
      vector<TH1*> plots;
      auto leg = prepLegends({}, {}, "P");
      for(std::map<TString,TH1D*>::iterator iter = overlap_1D_integrate_minus.begin(); iter != overlap_1D_integrate_minus.end(); ++iter){
        addLegendEntry(leg, iter->second, iter->first, "P");
        plots.push_back(iter->second);
      }
      leg->SetTextSize(0.03);
      leg->SetY1NDC(leg->GetY2NDC() - 0.3);
      //integral
      TCanvas* diff = drawCompMatt(plots, leg, -1., nullptr, "hist", true);
      TString name = geo + "_" + "_diff_kaptonAll" + pName;
      diff->Update();
      diff->SetTitle(name);
      diff->Print(outputDir+"/"+name+".pdf");
    } 
 

    std::ofstream jout;
    jout.open(outputDir+"/" + geo + "_moduleTolerances.json");
    json jtot = makeJSONModuleLatex(fit_values, worst_values);
    jout << jtot.dump(3);
    jout.close();
    printToleranceTableLatex(jtot, outputDir+"/" + geo + "_moduleTolerances.tex");
    printWorstTableLatex(jtot, outputDir+"/" + geo + "_worst_moduleTolerances.tex");

    for (int iP = 1; iP < nPeaks + 1; iP++){
      TString pName = "_Peak" + to_string(iP);
      for (json::iterator d = jtot.begin(); d != jtot.end(); ++d){
        if(d.key() == "BinNum") continue;
        if(d.key() == "Fit") continue;
        if(d.key() == "Worst") continue;
        vector<TH1*> hFit, hWorst, hNominal;
        vector<TString> binLabels(jtot["BinNum"].size());
        for (json::iterator comp = jtot[d.key()].begin(); comp != jtot[d.key()].end(); ++comp) {
          if(!TString(comp.key()).Contains(pName)) continue;
          vector<Quantity> hNom_val(jtot["BinNum"].size());
          vector<Quantity> hFit_val(jtot["BinNum"].size());
          vector<Quantity> hWorst_val(jtot["BinNum"].size());

          for (json::iterator dist = jtot["BinNum"].begin(); dist != jtot["BinNum"].end(); ++dist) {
            int binnum = jtot["BinNum"][dist.key()];
            binLabels.at(binnum) = dist.key();
            if(TString(comp.key()).Contains("Nominal")){
              hNom_val.at(binnum).value = jtot["Fit"][dist.key()][comp.key()][0];
              hNom_val.at(binnum).error = jtot["Fit"][dist.key()][comp.key()][1];
            } else {
              hFit_val.at(binnum).value = jtot["Fit"][dist.key()][comp.key()][0];
              hFit_val.at(binnum).error = jtot["Fit"][dist.key()][comp.key()][1];
              hWorst_val.at(binnum).value = jtot["Worst"][dist.key()][comp.key()][0];
              hWorst_val.at(binnum).error = jtot["Worst"][dist.key()][comp.key()][1];
            }
          }

          if(!TString(comp.key()).Contains("Nominal")){ 
            hFit.push_back(convertToHist(hFit_val, "Fit_" + TString(d.key()) + "_" + TString(comp.key()), ";; " + geo + " Fitted Value [#mum]", nullptr));
            hWorst.push_back(convertToHist(hWorst_val, "Worst_" + TString(d.key()) + "_" + TString(comp.key()), ";; " + geo + " Worst Value [#mum]", nullptr));
          }
          else hNominal.push_back(convertToHist(hNom_val, "Nominal_" + TString(d.key()) + "_" + TString(comp.key()), ";; " + geo + " Nominal Value [#mum]", nullptr));
        }

        if(hFit.size() == 0 || hWorst.size() == 0) continue;

        prepHists(hNominal, false, false, false, colors);
        prepHists(hFit, false, false, false, colors);
        prepHists(hWorst, false, false, false, colors);

        for(auto *f : hNominal) setBinLabels(f, binLabels);
        for(auto *f : hFit) setBinLabels(f, binLabels);
        for(auto *f : hWorst) setBinLabels(f, binLabels);

        vector<TH1*> hDiv;
        int i = 0;
        for(auto *f : hFit){
          TH1* hratio = (TH1*)f->Clone();
          hratio->Divide(hNominal[i]);
          hDiv.push_back(hratio);
          i++;
        }
        i = 0;
        vector<TH1*> hDiv_w;
        for(auto *f : hWorst){
          TH1* hratio = (TH1*)f->Clone();
          hratio->Divide(hNominal[i]);
          hDiv_w.push_back(hratio);
          i++;
        }

        auto leg = prepLegends({}, {""}, "L");
        for(int h = 0; h != (signed)hFit.size(); h++){
          cout << hFit[h]->GetName() << endl;
          TString legName = TString(hFit[h]->GetName());
          legName = translateString(legName, plotMap, "_", ", ");
          legName.ReplaceAll("Fit, ", "");
          if(legName.Contains("PCB ")) legName.ReplaceAll("Gaussian, ", "");
          if(legName.Contains(" mm")) legName.ReplaceAll("Gaussian, ", "");
          appendLegends(leg, {hFit[h]}, {legName}, "L");
        }
        leg->SetTextSize(0.04);
        leg->SetY1NDC(leg->GetY2NDC() - 0.2);
        if(TString(d.key()) == "PCBDist") setLegend(leg, 2, 0.3, 0.6, 0.92, 0.87);
        if(TString(d.key()) == "KaptonMultiDist") setLegend(leg, 2, 0.3, 0.6, 0.92, 0.87);
        if(TString(d.key()) == "Sensor") setLegend(leg, 1, 0.55, 0.6, 0.92, 0.87);
        TCanvas* c = drawCompAndRatio(hFit, hDiv, leg, "Fit/Nominal", 0.499, 1.501, true, 1., 2000., false, nullptr, false, 9999., true);
        TString outputName = outputDir + "/Fit_" + geo + "_" + TString(d.key()) + pName;
        c->SetTitle(outputName);
        c->Print(outputName+".pdf");

        leg = prepLegends({}, {""}, "L");
        for(int h = 0; h != (signed)hWorst.size(); h++){
          cout << hWorst[h]->GetName() << endl;
          TString legName = TString(hWorst[h]->GetName());
          legName = translateString(legName, plotMap, "_", ", ");
          legName.ReplaceAll("Worst, ", "");
          if(legName.Contains("PCB ")) legName.ReplaceAll("Gaussian, ", "");
          if(legName.Contains(" mm")) legName.ReplaceAll("Gaussian, ", "");
          appendLegends(leg, {hWorst[h]}, {legName}, "L");
        }
        leg->SetTextSize(0.04);
        leg->SetY1NDC(leg->GetY2NDC() - 0.2);
        if(TString(d.key()) == "PCBDist") setLegend(leg, 2, 0.3, 0.6, 0.92, 0.87);
        if(TString(d.key()) == "KaptonMultiDist") setLegend(leg, 2, 0.3, 0.6, 0.92, 0.87);
        if(TString(d.key()) == "Sensor") setLegend(leg, 1, 0.55, 0.6, 0.92, 0.87);
        c = drawCompAndRatio(hWorst, hDiv_w, leg, "Worst/Nominal", -1.001, 1.001, true, -1., 600, false, nullptr, false, -500, true);
        outputName = outputDir + "/Worst_" + geo + "_" + TString(d.key()) + pName;
        c->SetTitle(outputName);
        c->Print(outputName+".pdf");

        outputName.ReplaceAll("Worst_", "");
        TFile *outFile = new TFile(outputName+".root", "RECREATE");
        for(auto *f : hFit) f->Write();
        for(auto *f : hWorst) f->Write();
        for(auto *f : hNominal) f->Write();
        outFile->Close();

      }
    }// for different module shape

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
  TGaxis::SetMaxDigits(2);
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

void print2DPlots(TH2D *hc, TString geometry, TString BinLatex, TString name){
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
  TGaxis::SetMaxDigits(4);
  TPaletteAxis* palette = (TPaletteAxis*)hc->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(1.-0.18);
  palette->SetY1NDC(0.14);
  palette->SetX2NDC(1.-0.13);
  palette->SetY2NDC(1.-0.08);
  palette->SetLabelFont(42);
  palette->SetLabelSize(0.030);
  canvas->Modified();
  if(BinLatex != ""){
    TLatex *tl = new TLatex();
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.17, 0.87, BinLatex);
  }

  setTitleOffset(canvas, .950, .950);
  canvas->Update();
  canvas->SaveAs(outputDir + "/" + name + ".pdf");
}

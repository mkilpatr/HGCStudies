/*
 *
 *  Created on: Nov 11, 2019
 *      Author: mkilpatr
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <functional>
#include <numeric>
#include <glob.h>

#include "../utils/json.hpp"
using json = nlohmann::json;

#include "../utils/EstHelper.hh"

using namespace EstTools;

//files format
//time     bias voltage           meas current               meas voltage

map<TString, string> nLeg {
  {"Module805 Original", "Module805 Original"},
  {"16Apr21", "1 Cycle"},
  {"19Apr21", "2 Cycles"},
  {"19Apr21_n2", "3 Cycles"},
  {"21Apr21", "4 Cycles"},
  {"26Apr21_3", "5 Cycles"},
  {"28Apr21_10thermalcycles", "10 Cycles"},
  {"30Apr21_15thermalcycles", "15 Cycles"},
  {"JMay2821_20thermalcycles", "20 Cycles"},
  {"Jun0121_30thermalcycles", "30 Cycles"},
  {"Jun0221_40thermalcycles", "40 Cycles"},
  {"Jun0321_50thermalcycles", "50 Cycles"},
  {"1Cycle_07Jul21", "1 Cycle"},
  {"1Cycle_11Aug21", "1 Cycle"},
  {"2Cycle_08Jul21", "2 Cycles"},
  {"2Cycle_16Aug21", "2 Cycles"},
  {"2Cycle_16Aug21_attempt2", "2 Cycles"},
  {"3Cycle_09Jul21", "3 Cycles"},
  {"4Cycle_13Jul21", "4 Cycles"},
  {"5Cycle_14Jul21", "5 Cycles"},
  {"10Cycle_19Jul21", "10 Cycles, Bad IV"},
  {"11Cycle_10Aug21_Fix", "11 Cycles, Fix HV Connection"},
  {"1Cycle_24Aug21", "1 Cycle"},
  {"1_Cycle_24Aug21_attempt2", "1 Cycle v2"},
  {"2Cycles_25Aug21", "2 Cycle"},
  {"2Cycles_25Aug21_attempt2", "2 Cycle v2"},
  {"2Cycles_25Aug21_attempt3", "2 Cycle v3"},
};

map<TString, Color_t> nColor {
  {"Module805 Original", kBlack},
  {"16Apr21", kRed-9},
  {"19Apr21", kAzure+6},
  {"19Apr21_n2", kRed},
  {"21Apr21", kGreen},
  {"26Apr21_3", 606},
  {"28Apr21_10thermalcycles", kViolet+2},
  {"30Apr21_15thermalcycles", kSpring-8},
  {"JMay2821_20thermalcycles", kOrange-3},
  {"Jun0121_30thermalcycles",  kYellow-9},
  {"Jun0221_40thermalcycles", kGreen+3},
  {"Jun0321_50thermalcycles", kOrange},
  {"1Cycle_07Jul21", kRed-9},
  {"2Cycle_08Jul21", kAzure+6},
  {"2Cycle_16Aug21", kAzure+6},
  {"2Cycle_16Aug21_attempt2", kAzure+6},
  {"3Cycle_09Jul21", kRed},
  {"4Cycle_13Jul21", kGreen},
  {"5Cycle_14Jul21", 606},
  {"1Cycle_11Aug21", kRed},
  {"10Cycle_19Jul21", COLOR_MAP["color_comp1"]},
  {"11Cycle_10Aug21_Fix", COLOR_MAP["color_comp2"]},
  {"1Cycle_24Aug21", kRed -9},
  {"1_Cycle_24Aug21_attempt2", kAzure+6},
  {"2Cycles_25Aug21", kRed},
  {"2Cycles_25Aug21_attempt2", kGreen},
  {"2Cycles_25Aug21_attempt3", 606},
};

double strToDouble(string s){
  std::string::size_type sz;
  return std::stod(s,&sz);
}

json readFile(std::string FILENAME){
  std::ifstream file(FILENAME);
  string arr[4];
  int i = 0;
  json j;
  double biasV = 0., measV = 0., measC = 0.;
  double nAvg = 0.;
  if (file.is_open()) {
    std::string line;
    double prevVol = 0.;
    while (getline(file, line)) {
      string arr[4];
      int i = 0;
      stringstream ssin(line);
      std::string::size_type sz;
      while (ssin.good() && i < 4){
      	ssin >> arr[i];
      	++i;
      }

      int isEnd = ssin.rdbuf()->in_avail();

      if(isEnd != 0) biasV = strToDouble(arr[1]);
      if((biasV != prevVol && nAvg > 0) || isEnd == 0){ 
        if(isnan(measV/nAvg) || isnan(measC/nAvg)){
          cout << "Measured Current or Voltage was NULL in file " << FILENAME << "!!!" << endl;
          cout << "(MeasV, MeasC, nAvg) = (" << measV << ", " << measC << ", " << nAvg << ")" << endl;
        } else j[FILENAME][to_string(prevVol)] = {measV/nAvg, measC/nAvg};
        measC = 0.;
        measV = 0.;
        nAvg  = 0.;
      }
      if(biasV == prevVol && isEnd != 0){
        measC += strToDouble(arr[2]);
        measV += strToDouble(arr[3]);
        nAvg++;
      }
      prevVol = biasV;
    }
    file.close();
  }
  return j;
}

std::vector<std::string> readFileTotal(const string& pattern){
    glob_t glob_result;
    glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> files;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        files.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return files;
}

void IVCurvePlotter(std::string indir_ = "testDir", TString suffix = "module805", TString label = "Module X"){
  gROOT->SetBatch(1);

  std::string indir = indir_ + "/*";
  std::vector<std::string> files = readFileTotal(indir);
  for(auto s : files) cout << s << endl;
  std::vector<double> val_up_total, val_down_total;
  json j, jtot, j_orig;
  std::ofstream jout;
  jout.open(indir_ + "/" + suffix + ".json");
  for(unsigned i = 0; i != files.size(); i++){
    TString name = TString(files[i]);
    if(!name.Contains("txt")) continue;
    cout << files[i] << endl;
    j = readFile(files[i]);
    for (json::iterator it = j.begin(); it != j.end(); ++it) {
      jtot[it.key()] = it.value();
      
    }
  }
  jout << jtot.dump(3);
  jout.close();

  std::ifstream orig(suffix + "Master.json");
  orig >> j_orig;

  vector<Color_t> colors;
  vector<TH1*> hTotal;
  for (json::iterator unc = j_orig.begin(); unc != j_orig.end(); ++unc) {
    TString type = TString(unc.key());
    TH1D* hIV = new TH1D(type, ";Voltage (V);Current (#muA)", 32, 0, 800);
    cout << type << endl;
    TAxis *xaxis = hIV->GetXaxis();
    for (json::iterator bin = j_orig[unc.key()].begin(); bin != j_orig[unc.key()].end(); ++bin) {
      int binnum = xaxis->FindBin(strToDouble(bin.key()));
      hIV->SetBinContent(binnum, double(j_orig[unc.key()][bin.key()][1])*1000000);
      hIV->SetBinError(binnum, 0.);
    }
    hTotal.push_back(hIV);
    TString colName = hIV->GetName();
    colName.ReplaceAll(indir_ + "/data_IVCurve_Neg30_25Vper10sec_", "");
    colName.ReplaceAll("sec", "s ");
    colName.ReplaceAll(".txt", "");
    colName.ReplaceAll("_25Vper10s", "");
    cout << colName << endl;
    colors.push_back(nColor[colName]);
  }

  for (json::iterator unc = jtot.begin(); unc != jtot.end(); ++unc) {
    TString type = TString(unc.key());
    TH1D* hIV = new TH1D(type, ";Voltage (V);Current (#muA)", 32, 0, 800);
    cout << type << endl;
    TAxis *xaxis = hIV->GetXaxis();
    for (json::iterator bin = jtot[unc.key()].begin(); bin != jtot[unc.key()].end(); ++bin) {
      int binnum = xaxis->FindBin(strToDouble(bin.key()));
      hIV->SetBinContent(binnum, double(jtot[unc.key()][bin.key()][1])*1000000);
      hIV->SetBinError(binnum, 0.);
    }
    hTotal.push_back(hIV);
    TString colName = hIV->GetName();
    colName.ReplaceAll(indir_ + "/data_", "");
    colName.ReplaceAll("IVCurve_Neg30_25Vper10sec_", "");
    colName.ReplaceAll("IVCurve_30C_", "");
    colName.ReplaceAll("sec", "s ");
    colName.ReplaceAll(".txt", "");
    colName.ReplaceAll("_25Vper10s", "");
    cout << colName << endl;
    colors.push_back(nColor[colName]);

  }

  prepHists(hTotal, false, false, false, colors);

  for(auto* h : hTotal){
    h->SetMarkerStyle(8);
  }

  auto leg = prepLegends({}, {""}, "l");
  for(unsigned h = 0; h != hTotal.size(); h++){
    TString legName = hTotal[h]->GetName();
    cout << legName << endl;
    legName.ReplaceAll(indir_ + "/data_", "");
    legName.ReplaceAll("IVCurve_Neg30_25Vper10sec_", "");
    legName.ReplaceAll("IVCurve_30C_", "");
    legName.ReplaceAll("sec", "s ");
    legName.ReplaceAll(".txt", "");
    legName.ReplaceAll("_25Vper10s", "");
    legName.ReplaceAll(legName, nLeg[legName]);
    appendLegends(leg, {hTotal[h]}, {legName}, "P");
  }

  std::function<void(TCanvas*)> plotextra = [&](TCanvas *c){ c->cd(); drawTLatexNDC(label, 0.2, 0.94, 0.04); };

  setLegend(leg, 1, 0.2, 0.57, 0.94, 0.90);
  leg->SetTextSize(0.04);
  leg->SetY1NDC(leg->GetY2NDC() - 0.2);
  TCanvas* c = drawCompMatt(hTotal, leg, -1., &plotextra, "P", true);
  TString typeName = "IVCurves_"+suffix;
  c->SetTitle(typeName);
  c->Print(indir_+"/"+typeName+".pdf");

  c = drawCompMatt(hTotal, leg, 0.001, nullptr, "P", true);
  c->SetTitle(typeName + "_log");
  c->Print(indir_+"/"+typeName+"_log.pdf");

  TFile *outFile = new TFile(indir_+"/"+typeName+".root", "RECREATE");
  for(unsigned h = 0; h != hTotal.size(); h++){
    hTotal[h]->Write();
  }
  outFile->Close();

}

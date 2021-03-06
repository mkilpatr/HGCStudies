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

vector<TString> hNames = {"First Cycle", "Second Cycle", "Third Cycle"};

double strToDouble(string s){
  std::string::size_type sz;
  return std::stod(s,&sz);
}

json readFile(std::string FILENAME){
//std::pair< std::pair<vector<TString>, vector<double> >, std::pair<vector<TString>, vector<double> > > readFile(std::string FILENAME){
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
      //cout << line << endl;
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
      if(biasV != prevVol || isEnd == 0){ 
        if(isnan(measV/nAvg) || isnan(measC/nAvg)) cout << "Measured Current or Voltage was NULL!!!" << endl;
        else j[FILENAME][to_string(prevVol)] = {measV/nAvg, measC/nAvg};
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

void IVCurvePlotter(std::string indir_ = "testDir", TString suffix = "module805"){

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
    j = readFile(files[i]);
    for (json::iterator it = j.begin(); it != j.end(); ++it) {
      jtot[it.key()] = it.value();
      
    }
  }
  jout << jtot.dump(3);
  jout.close();

  std::ifstream orig("module805Master.json");
  orig >> j_orig;

  vector<TH1*> hTotal;
  for (json::iterator unc = j_orig.begin(); unc != j_orig.end(); ++unc) {
    TString type = TString(unc.key());
    TH1D* hIV = new TH1D(type, ";Voltage (V);Current (#muA)", 22, 0, 1100);
    cout << type << endl;
    TAxis *xaxis = hIV->GetXaxis();
    for (json::iterator bin = j_orig[unc.key()].begin(); bin != j_orig[unc.key()].end(); ++bin) {
      int binnum = xaxis->FindBin(strToDouble(bin.key()));
      hIV->SetBinContent(binnum, double(j_orig[unc.key()][bin.key()][1])*1000000);
      hIV->SetBinError(binnum, 0.);
    }
    hTotal.push_back(hIV);
  }

  cout << "Made it past master" << endl;
  for (json::iterator unc = jtot.begin(); unc != jtot.end(); ++unc) {
    TString type = TString(unc.key());
    TH1D* hIV = new TH1D(type, ";Voltage (V);Current (#muA)", 44, 0, 1100);
    cout << type << endl;
    TAxis *xaxis = hIV->GetXaxis();
    for (json::iterator bin = jtot[unc.key()].begin(); bin != jtot[unc.key()].end(); ++bin) {
      int binnum = xaxis->FindBin(strToDouble(bin.key()));
      hIV->SetBinContent(binnum, double(jtot[unc.key()][bin.key()][1])*1000000);
      hIV->SetBinError(binnum, 0.);
    }
    hTotal.push_back(hIV);
  }

  prepHists(hTotal, false, false, false, {kBlack, kRed-9, kAzure+6, kRed, kGreen, 606, kViolet+2});

  for(auto* h : hTotal){
    h->SetMarkerStyle(8);
  }

  auto leg = prepLegends({}, {""}, "l");
  for(unsigned h = 0; h != hTotal.size(); h++){
    TString legName = hTotal[h]->GetName();//hNames[h];
    legName.ReplaceAll(indir_ + "/data_IVCurve_Neg30_", "");
    legName.ReplaceAll("sec", "s ");
    legName.ReplaceAll(".txt", "");
    legName.ReplaceAll("_25Vper10s", "");
    legName.ReplaceAll("_10thermalcycles", "10 Thermal Cycles");
    appendLegends(leg, {hTotal[h]}, {legName}, "P");
  }
  setLegend(leg, 1, 0.2, 0.60, 0.94, 0.87);
  leg->SetTextSize(0.04);
  leg->SetY1NDC(leg->GetY2NDC() - 0.2);
  TCanvas* c = drawCompMatt(hTotal, leg, -1., nullptr, "P", true);
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

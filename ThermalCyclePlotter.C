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
#include "TVectorD.h" 
#include "TGraph.h"

#include "../utils/json.hpp"
using json = nlohmann::json;

#include "../utils/EstHelper.hh"

using namespace EstTools;

//files format
//time     bias voltage           meas current               meas voltage

map<TString, string> nLeg {
  {"Module805 Original", "Module805 Original"},
};

map<TString, Color_t> nColor {
  {"Corner", kBlack},
  {"Edges",  kRed},
  {"Middle", kMagenta},
  {"Point",  kAzure},
  {"Center", kGreen},
};

double strToDouble(string s){
  std::string::size_type sz;
  return std::stod(s,&sz);
}


vector<string> ModuleHeight = {"Corner", "Edges", "Middle", "Point", "Center"};
//vector<string> ModuleCycle_805 = {"0 Cycles", "1 Cycle", "2 Cycles", "3 Cycles", "4 Cycles", "5 Cycles", "10 Cycles", "After Fix", "11 Cycles", "15 Cycles"};
//vector<string> ModuleCycle = {"0 Cycles", "1 Cycle", "2 Cycles", "3 Cycles", "4 Cycles", "5 Cycles", "10 Cycles", "After Fix", "11 Cycles", "15 Cycles"};
vector<string> ModuleCycle = {"OGP", "0 Cycles", "1 Cycle", "2 Cycles", "3 Cycles", "Remount", "4 Cycles", "5 Cycles", "6 Cycles"};
vector<int> PadCenter = {88};
vector<int> PadCorner = {1, 18, 111, 197, 190, 81, 8, 95, 189, 191, 96, 9};
vector<int> PadEdge = {3, 39, 140, 195, 169, 52, 5, 65, 168, 193, 141, 29};
vector<int> PadMiddle = {22, 49, 124, 175, 158, 68, 24, 78, 152, 173, 129, 42};
vector<int> PadPoint = {45, 76, 122, 148, 131, 71};
vector<int> ChanCenter = {46};
vector<int> ChanCorner = {33, 61, 32, 57, 31, 70, 6, 33, 5, 29, 11, 35};
vector<int> ChanEdge   = {21, 67, 21, 68, 23, 67, 14, 29, 16, 66, 21, 29};
vector<int> ChanMiddle = {25, 57, 25, 50, 34, 57, 10, 71, 9, 59, 15, 59};
vector<int> ChanPoint  = {36, 52, 40, 48, 36, 52};

json readFile(std::string FILENAME){
  std::ifstream file(FILENAME);
  string arr[7];
  json j;
  ifstream infile(FILENAME);
  int h = 0, c = 0, i = 0;
  int cen = 0, cor = 0, edg = 0, mid = 0, pnt = 0;
  while (infile){
    string s;
    if (!getline( infile, s )) break;
    istringstream ss( s );

    while (ss){
      string s;
      if (!getline( ss, s, ',' )) break;
      if (s == "") continue;
      double height = strToDouble(s);
      j[FILENAME]["Channeltype"].push_back(0);
      j[FILENAME]["Location"].push_back(ModuleHeight[h]);
      j[FILENAME]["Cycle"].push_back(ModuleCycle[c]);
      j[FILENAME]["Height"].push_back(height);
      if(i == 0){
        j[FILENAME]["Pad"].push_back(PadCorner[cor]);
        j[FILENAME]["Channel"].push_back(ChanCorner[cor]);
        j[FILENAME]["Which"].push_back(cor);
        cor++;
      } else if(i == 1 || i == 2) {
        j[FILENAME]["Pad"].push_back(PadEdge[edg]);
        j[FILENAME]["Channel"].push_back(ChanEdge[edg]);
        j[FILENAME]["Which"].push_back(edg);
        edg++;
      } else if(i == 3){
        j[FILENAME]["Pad"].push_back(PadCorner[cor]);
        j[FILENAME]["Channel"].push_back(ChanCorner[cor]);
        j[FILENAME]["Which"].push_back(cor);
        cor++;
      } else if(i == 4 || i == 5){
        j[FILENAME]["Pad"].push_back(PadMiddle[mid]);
        j[FILENAME]["Channel"].push_back(ChanMiddle[mid]);
        j[FILENAME]["Which"].push_back(mid);
         mid++;
      } else if(i == 6){
        j[FILENAME]["Pad"].push_back(PadPoint[pnt]);
        j[FILENAME]["Channel"].push_back(ChanPoint[pnt]);
        j[FILENAME]["Which"].push_back(pnt);
         pnt++;
      } else if(i == 7){
        j[FILENAME]["Pad"].push_back(PadCenter[cen]);
        j[FILENAME]["Channel"].push_back(ChanCenter[cen]);
        j[FILENAME]["Which"].push_back(cen);
         cen++;
      }
      
    }

    i++;
    if(i == 1) h++;
    if(i == 3) h--;
    if(i == 4) h+=2;
    if(i == 6 || i == 7) h++;
    if(i == 8) {c++; i = 0; h = 0; cen = 0; cor = 0; edg = 0; mid = 0; pnt = 0;}
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

int getIndex(vector<string> v, string K){

    for(int i = 0; i != v.size(); i++){
      if(v[i] == K) return i;
    }

    return -1;
}

void ThermalCyclePlotter(std::string indir_ = "testDir", string filename = "", string suffix = "module805", TString label = "Module X"){
  gROOT->SetBatch(1);

  TString suffix_ = TString(suffix);
  std::string indir = indir_ + "/*";
  std::vector<std::string> files = readFileTotal(indir);
  for(auto s : files) cout << s << endl;
  std::vector<double> val_up_total, val_down_total;
  json j, jtot, j_orig;
  std::ofstream jout;
  jout.open(indir_ + "/" + suffix_ + ".json");
  for(unsigned i = 0; i != files.size(); i++){
    TString name = TString(files[i]);
    if(!name.Contains(filename)) continue;
    j = readFile(files[i]);
    for (json::iterator it = j.begin(); it != j.end(); ++it) {
      jtot[suffix] = it.value();
    }
  }
  jout << jtot.dump(3);
  jout.close();
}

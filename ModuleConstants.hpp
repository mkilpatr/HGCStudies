// ModuleConstants.hpp
namespace constants {
  TString outputDir = "";
  TString eosDir = "/eos/uscms/store/user/mkilpatr/13TeV/";
  TString localDir = "";
  const bool debug = false;
  const TString whichOverlap = "bas_kap_stack_hist";

  //Define Constant Widths
  const double pcb_w_const = 166.80, sensor_w_const = 166.57, kapton_w_const = 167.19, baseplate_w_const = 166.94;

  const double pcb_a_const_semi = 96.7799, sensor_a_const_semi = pcb_a_const_semi + (sensor_w_const - pcb_w_const), kapton_a_const_semi = pcb_a_const_semi + (kapton_w_const - pcb_w_const), baseplate_a_const_semi = pcb_a_const_semi + (baseplate_w_const - pcb_w_const);

  double baseplate_center = 3.05, baseplate_err = 0.05;
  double baseplate_centerX = 0., baseplate_centerY = 0.;
  double baseplate_pinhole_err = 0.050, baseplate_pinhole_theta = 2*TMath::Pi();
  double cassette_err = 0.050; 
  double pcb_w = pcb_w_const, sensor_w = sensor_w_const, kapton_w = kapton_w_const, baseplate_w = baseplate_w_const;
  double pcb_w_err = 0.150, sensor_w_err = 0.02, kapton_w_err = 0.05, baseplate_w_err = 0.05;
  double PcbToBas_shift_r = 0.04, PcbToBas_shift_theta = 2*TMath::Pi();
  double sensor_shift_x= 0.500, sensor_shift_y = 0.050, sensor_shift_theta = 0.050; //theta is radians
  double baseValue_total = sensor_w;
  double base_kap_min = 0.100;
  double gold_min = 0.25, gold_size = 1.500, backside_min = 0.250;
  double senToKap_min = 0.0;
  pair<double, double> backside_x_err = make_pair(0.2133, 0.1769);
  pair<double, double> backside_y_err = make_pair(0.250, 0.000);
  
  int max = 300000;
  double weight = 30000./double(max);

  int nbins = 650;
  double step = 0.002;
  double width_new = 1.4;
  double width_backY = 4.0;
  double axis = width_new/2 + nbins*step + 0.1;

  int nPeaks = 6;
  //     if(geo == "Five" || geo == "Semi") nPeaks = 3;
  //else if(geo == "Half" || geo == "Three") nPeaks = 2;

  //TString whichGroup = "Kaptonplus"; //group by different Kapton widths
  //TString whichGroup = "PCBplus"; //group by different PCB widths
  //TString whichGroup = "Base", whichGroupSpecific = "Baseplus000"; //group by different Baseplate widths
  //TString whichGroup = "Peak", whichGroupSpecific = "Peak1"; //group by different Peak
  TString whichGroup = "Center", whichGroupSpecific = "CenterXplus000"; //group by different center XY locations

  //TString whichGroupSpecific = "Kaptonplus"; //group by different Kapton widths
  //TString whichGroupSpecific = "PCBplus"; //group by different PCB widths

  //vector<string> Geometry = {"Full", "Five", "Semi", "Half", "Three"};
  vector<string> Geometry = {"Full"};
  vector<TString> whichComp = {"sen_kap_stack_hist", "sen_bas_stack_hist", "bas_kap_stack_hist", "bas_kap_sen_stack_hist", "sen_pcb_stack_hist", "pcb_bas_stack_hist", "pcb_kap_stack_hist", "kap_pcb_hist", "sen_pcb_hist", "sen_pcb_kap_x_hist", "sen_pcb_kap_y_hist"};
  //vector<TString> whichComp = {"bas_kap_stack_hist", "bas_kap_sen_stack_hist", "sen_pcb_stack_hist", "pcb_bas_stack_hist"};
  vector<TString> whichPlot = whichComp;//{"bas_kap_stack_hist", "bas_kap_sen_stack_hist", "sen_pcb_stack_hist", "pcb_bas_stack_hist"};
  vector<string> Dist = {
               "Gaussian_Baseplus000_CenterXplus000_CenterYplus000_midSensor",  "Gaussian_Baseplus000_CenterXplus000_CenterYplus000_newSensor",
               "Gaussian_Baseplus000_CenterXplus050_CenterYplus050_midSensor",  "Gaussian_Baseplus000_CenterXplus050_CenterYplus050_newSensor",
               "Gaussian_Baseplus000_CenterXplus100_CenterYplus030_midSensor",  "Gaussian_Baseplus000_CenterXplus100_CenterYplus030_newSensor",
               "Gaussian_Baseplus000_CenterXplus100_CenterYplus070_midSensor",  "Gaussian_Baseplus000_CenterXplus100_CenterYplus070_newSensor",
               "Gaussian_Baseplus000_CenterXplus150_CenterYplus050_midSensor",  "Gaussian_Baseplus000_CenterXplus150_CenterYplus050_newSensor",
               "Gaussian_Baseplus025_CenterXplus000_CenterYplus000_midSensor",  "Gaussian_Baseplus025_CenterXplus000_CenterYplus000_newSensor",
               "Gaussian_Baseplus025_CenterXplus050_CenterYplus050_midSensor",  "Gaussian_Baseplus025_CenterXplus050_CenterYplus050_newSensor",
               "Gaussian_Baseplus025_CenterXplus100_CenterYplus030_midSensor",  "Gaussian_Baseplus025_CenterXplus100_CenterYplus030_newSensor",
               "Gaussian_Baseplus025_CenterXplus100_CenterYplus070_midSensor",  "Gaussian_Baseplus025_CenterXplus100_CenterYplus070_newSensor",
               "Gaussian_Baseplus025_CenterXplus150_CenterYplus050_midSensor",  "Gaussian_Baseplus025_CenterXplus150_CenterYplus050_newSensor",
               "Gaussian_Baseplus050_CenterXplus000_CenterYplus000_midSensor",  "Gaussian_Baseplus050_CenterXplus000_CenterYplus000_newSensor",
               "Gaussian_Baseplus050_CenterXplus050_CenterYplus050_midSensor",  "Gaussian_Baseplus050_CenterXplus050_CenterYplus050_newSensor",
               "Gaussian_Baseplus050_CenterXplus100_CenterYplus030_midSensor",  "Gaussian_Baseplus050_CenterXplus100_CenterYplus030_newSensor",
               "Gaussian_Baseplus050_CenterXplus100_CenterYplus070_midSensor",  "Gaussian_Baseplus050_CenterXplus100_CenterYplus070_newSensor",
               "Gaussian_Baseplus050_CenterXplus150_CenterYplus050_midSensor",  "Gaussian_Baseplus050_CenterXplus150_CenterYplus050_newSensor",
                           };
  
  const vector<TString> Order = {
               "Gaussian_Baseplus000_CenterXplus000_CenterYplus000_newSensor","Gaussian_Baseplus000_CenterXplus000_CenterYplus000_midSensor",  
               "Gaussian_Baseplus000_CenterXplus050_CenterYplus050_newSensor","Gaussian_Baseplus000_CenterXplus050_CenterYplus050_midSensor",  
               "Gaussian_Baseplus000_CenterXplus100_CenterYplus030_newSensor","Gaussian_Baseplus000_CenterXplus100_CenterYplus030_midSensor",  
               "Gaussian_Baseplus000_CenterXplus100_CenterYplus070_newSensor","Gaussian_Baseplus000_CenterXplus100_CenterYplus070_midSensor",  
               "Gaussian_Baseplus000_CenterXplus150_CenterYplus050_newSensor","Gaussian_Baseplus000_CenterXplus150_CenterYplus050_midSensor",  
               "Gaussian_Baseplus025_CenterXplus000_CenterYplus000_newSensor","Gaussian_Baseplus025_CenterXplus000_CenterYplus000_midSensor",  
               "Gaussian_Baseplus025_CenterXplus050_CenterYplus050_newSensor","Gaussian_Baseplus025_CenterXplus050_CenterYplus050_midSensor",  
               "Gaussian_Baseplus025_CenterXplus100_CenterYplus030_newSensor","Gaussian_Baseplus025_CenterXplus100_CenterYplus030_midSensor",  
               "Gaussian_Baseplus025_CenterXplus100_CenterYplus070_newSensor","Gaussian_Baseplus025_CenterXplus100_CenterYplus070_midSensor",  
               "Gaussian_Baseplus025_CenterXplus150_CenterYplus050_newSensor","Gaussian_Baseplus025_CenterXplus150_CenterYplus050_midSensor",  
               "Gaussian_Baseplus050_CenterXplus000_CenterYplus000_newSensor","Gaussian_Baseplus050_CenterXplus000_CenterYplus000_midSensor",  
               "Gaussian_Baseplus050_CenterXplus050_CenterYplus050_newSensor","Gaussian_Baseplus050_CenterXplus050_CenterYplus050_midSensor",  
               "Gaussian_Baseplus050_CenterXplus100_CenterYplus030_newSensor","Gaussian_Baseplus050_CenterXplus100_CenterYplus030_midSensor",  
               "Gaussian_Baseplus050_CenterXplus100_CenterYplus070_newSensor","Gaussian_Baseplus050_CenterXplus100_CenterYplus070_midSensor",  
               "Gaussian_Baseplus050_CenterXplus150_CenterYplus050_newSensor","Gaussian_Baseplus050_CenterXplus150_CenterYplus050_midSensor",  
                           };
  
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
    {"bas_kap_sen_stack_hist",  "Base to Kap w/ Sen"},
    {"sen_kap_stack_hist",  "Sen to Kap"},
    {"sen_pcb_stack_hist",  "Sen to PCB"},
    {"allOverlaps_stack_hist",  "All Overlaps"},
  };

  map<TString, string> compMap {
    {"kap_pcb_hist",            "ShieldBond"},
    {"sen_pcb_hist",            "GuardBond"},
    {"sen_pcb_kap_x_hist",      "BacksideBondX"},
    {"sen_pcb_kap_y_hist",      "BacksideBondY"},
    {"pcb_bas_stack_hist",      "PCBToBas"},
    {"sen_bas_stack_hist",      "SenToBas"},
    {"pcb_kap_stack_hist",      "PCBToKap"},
    {"bas_kap_stack_hist",      "BasToKap"},
    {"bas_kap_sen_stack_hist",  "BasToKapWSen"},
    {"sen_kap_stack_hist",      "SenToKap"},
    {"sen_pcb_stack_hist",      "SenToPCB"},
  };

  map<TString, Color_t> colors {
    {"minAllOverlaps_stack_hist",       kBlue },
    {"maxAllOverlaps_stack_hist",       kBlue },
    {"kap_pcb_hist",       kBlue },
    {"sen_pcb_hist",       kGreen + 3},
    {"sen_pcb_kap_x_hist", kRed},
    {"sen_pcb_kap_y_hist", kCyan},
    {"pcb_bas_stack_hist", kViolet-4},
    {"sen_bas_stack_hist", kSpring-9},
    {"pcb_kap_stack_hist", kMagenta},
    {"bas_kap_stack_hist", kAzure},
    {"bas_kap_sen_stack_hist", kAzure},
    {"sen_kap_stack_hist", kBlue },
    {"sen_pcb_stack_hist", kGreen + 3 },
  };

  map<TString, pair<double, double>> axisX {
    {"minAllOverlaps_stack_hist",       {}},
    {"maxAllOverlaps_stack_hist",       {}},
    {"kap_pcb_hist",       		{0., 1.}},
    {"sen_pcb_hist",       		{0., 1.}},
    {"sen_pcb_kap_x_hist", 		{0.5, 1.5}},
    {"sen_pcb_kap_y_hist", 		{1.5, 2.5}},
    {"pcb_bas_stack_hist", 		{-0.3, 0.2}},
    {"sen_bas_stack_hist", 		{-0.2, 0.2}},
    {"pcb_kap_stack_hist", 		{-0.3, 0.2}},
    {"bas_kap_stack_hist", 		{-0.3, 0.1}},
    {"bas_kap_sen_stack_hist", 		{-0.3, 0.3}},
    {"sen_kap_stack_hist", 		{-0.3, 0.3}},
    {"sen_pcb_stack_hist", 		{-0.2, 0.2}},
  };

  const std::map<TString, TString> latexMap{
    {"Custom", ""},
    {"Distribution", ""},
    {"PCBDist", ""},
    {"KaptonMultiDist", ""},
    {"TotalBadModules", ""},
    {"Sensor", ""},
    {"forward", ""},
    {"backward", ""},
    {"CenterXplus000", R"(center $= (0)"},
    {"CenterXplus050", R"(center $= (50)"},
    {"CenterXplus100", R"(center $= (100)"},
    {"CenterXplus150", R"(center $= (150)"},
    {"CenterYplus000", R"( 0) \mu m$)"},
    {"CenterYplus050", R"( 50) \mu m$)"},
    {"CenterYplus030", R"( 30) \mu m$)"},
    {"CenterYplus050", R"( 50) \mu m$)"},
    {"CenterYplus070", R"( 70) \mu m$)"},
    {"otherCenter", R"(center $= + 29 \mu m$)"},
    {"PCBplus000", R"(Nominal Hexaboard width)"},
    {"Baseplus000", R"(Nominal Baseplate width)"},
    {"PCBminus25", R"(PCB - 25 $\mu m$)"},
    {"PCBplus25", R"(Nominal Hexaboard width + 25 $\mu m$)"},
    {"PCBminus50", R"(PCB - 50 $\mu m$)"},
    {"PCBplus050", R"(Nominal Hexaboard width + 50 $\mu m$)"},
    {"PCBminus75", R"(PCB - 75 $\mu m$)"},
    {"PCBplus75", R"(Nominal Hexaboard width + 75 $\mu m$)"},
    {"PCBplus100", R"(Nominal Hexaboard width + 100 $\mu m$)"},
    {"PCBplus125", R"(Nominal Hexaboard width + 125 $\mu m$)"},
    {"PCBplus150", R"(Nominal Hexaboard width + 150 $\mu m$)"},
    {"PCBplus175", R"(Nominal Hexaboard width + 175 $\mu m$)"},
    {"PCBplus200", R"(Nominal Hexaboard width + 200 $\mu m$)"},
    {"PCBplus225", R"(Nominal Hexaboard width + 225 $\mu m$)"},
    {"PCBplus250", R"(Nominal Hexaboard width + 250 $\mu m$)"},
    {"PCBplus275", R"(Nominal Hexaboard width + 275 $\mu m$)"},
    {"PCBplus300", R"(Nominal Hexaboard width + 300 $\mu m$)"},
    {"Baseplus025", R"(Nominal Baseplate width + 25 $\mu m$)"},
    {"Baseplus050", R"(Nominal Baseplate width + 50 $\mu m$)"},
    {"Kaptonplus000", R"(Nominal)"},
    {"Kaptonminus25", R"(Kapton - 25 $\mu m$)"},
    {"Kaptonplus025", R"(Nominal Kapton width + 25 $\mu m$)"},
    {"Kaptonminus50", R"(Kapton - 50 $\mu m$)"},
    {"Kaptonplus050", R"(Nominal Kapton width + 50 $\mu m$)"},
    {"Kaptonminus75", R"(Kapton - 75 $\mu m$)"},
    {"Kaptonplus075", R"(Nominal Kapton width + 75 $\mu m$)"},
    {"Kaptonminus150", R"(Kapton - 150 $\mu m$)"},
    {"Kaptonminus170", R"(Kapton - 170 $\mu m$)"},
    {"Kaptonminus75", R"(Kapton - 75 $\mu m$)"},
    {"Kaptonplus100", R"(Nominal Kapton width + 100 $\mu m$)"},
    {"Kaptonplus125", R"(Nominal Kapton width + 125 $\mu m$)"},
    {"Kaptonplus150", R"(Nominal Kapton width + 150 $\mu m$)"},
    {"Kaptonplus160", R"(Nominal Kapton width + 160 $\mu m$)"},
    {"Kaptonplus170", R"(Nominal Kapton width + 170 $\mu m$)"},
    {"Kaptonplus180", R"(Nominal Kapton width + 180 $\mu m$)"},
    {"Kaptonplus190", R"(Nominal Kapton width + 190 $\mu m$)"},
    {"Kaptonplus200", R"(Nominal Kapton width + 200 $\mu m$)"},
    {"Kaptonplus210", R"(Nominal Kapton width + 210 $\mu m$)"},
    {"Kaptonplus220", R"(Nominal Kapton width + 220 $\mu m$)"},
    {"Kaptonplus230", R"(Nominal Kapton width + 230 $\mu m$)"},
    {"Kaptonplus240", R"(Nominal Kapton width + 240 $\mu m$)"},
    {"Kaptonplus0250", R"(Nominal Kapton width + 250 $\mu m$)"},
    {"Kaptonplus260", R"(Nominal Kapton width + 260 $\mu m$)"},
    {"Kaptonplus270", R"(Nominal Kapton width + 270 $\mu m$)"},
    {"Kaptonplus280", R"(Nominal Kapton width + 280 $\mu m$)"},
    {"Kaptonplus290", R"(Nominal Kapton width + 290 $\mu m$)"},
    {"Kaptonplus300", R"(Nominal Kapton width + 300 $\mu m$)"},
    {"Kaptonplus400", R"(Nominal Kapton width + 400 $\mu m$)"},
    {"Kaptonplus500", R"(Nominal Kapton width + 500 $\mu m$)"},
    {"senTokap000", R"(\\ Incomplete coverage = 0 $\mu m$)"},
    {"senTokap100", R"(\\ Incomplete coverage = 100 $\mu m$)"},
    {"senTokap150", R"(\\ Incomplete coverage = 150 $\mu m$)"},
    {"senTokap185", R"(\\ Incomplete coverage = 185 $\mu m$)"},
    {"Nominal", "Nominal"},
    {"Gaussian", "Gaussian"},
    {"Flat", "Flat"},
    {"CustomGaus", "Custom Gaussian"},
    {"CustomFlat", "Custom Flat"},
    {"newSensor", R"([50 $\mu m$])"},
    {"midSensor", R"([100 $\mu m$])"},
    {"oldSensor", R"([200 $\mu m$])"},
    {"newSplitSensor", R"(}} & [50 $\mu m$])"},
    {"midSplitSensor", R"( & [100 $\mu m$])"},
    {"oldSplitSensor", R"(}} & [200 $\mu m$])"},
    {"min", ""},
    {"max", ""},
    {"Fit", ""},
    {"Worst", ""},
    {"integrate", ""},
    {"Peak1", "Peak 1"},
    {"Peak2", "Peak 2"},
    {"Peak3", "Peak 3"},
    {"Peak4", "Peak 4"},
    {"Peak5", "Peak 5"},
    {"Peak6", "Peak 6"},
  };
  
  const std::map<TString, TString> plotMap{
    {"Custom", ""},
    {"Distribution", ""},
    {"PCBDist", ""},
    {"KaptonMultiDist", ""},
    {"TotalBadModules", ""},
    {"Sensor", ""},
    {"forward", ""},
    {"backward", ""},
    {"ShieldBond", "Shield Bond"},
    {"GuardBond",  "Guard Bond"},
    {"BacksideBondX", "Backside Bond X"},
    {"BacksideBondY", "Backside Bond Y"},
    {"SenToBas", "Sen to Base"},
    {"BasToKapMax", "Max Base to Kap"},
    {"BasToKapWSen", "Base to Kap with Sen"},
    {"BasToKap", "Base to Kap"},
    {"SenToPCB", "PCB to Sen"},
    {"PCBToBas", "Base to PCB"},
    {"PCBToKap", "Kap to PCB"},
    {"SenToKap", "Kap to Sen"},
    {"CenterXplus000", R"(center = (0)"},
    {"CenterXplus050", R"(center = (50)"},
    {"CenterXplus100", R"(center = (100)"},
    {"CenterXplus150", R"(center = (150)"},
    {"CenterYplus000", R"( 0) #mum)"},
    {"CenterYplus050", R"( 50) #mum)"},
    {"CenterYplus030", R"( 30) #mum)"},
    {"CenterYplus050", R"( 50) #mum)"},
    {"CenterYplus070", R"( 70) #mum)"},
    {"otherCenter", R"(center = + 29 #mum)"},
    {"PCBplus000", R"(Nominal Hexaboard width)"},
    {"Baseplus000", R"(Nominal Baseplate width)"},
    {"PCBminus25", R"(PCB - 25 #mum)"},
    {"PCBplus25", R"(Nominal Hexaboard width + 25 #mum)"},
    {"PCBminus50", R"(PCB - 50 #mum)"},
    {"PCBplus050", R"(Nominal Hexaboard width + 50 #mum)"},
    {"PCBminus75", R"(PCB - 75 #mum)"},
    {"PCBplus75", R"(Nominal Hexaboard width + 75 #mum)"},
    {"PCBplus100", R"(Nominal Hexaboard width + 100 #mum)"},
    {"PCBplus125", R"(Nominal Hexaboard width + 125 #mum)"},
    {"PCBplus150", R"(Nominal Hexaboard width + 150 #mum)"},
    {"PCBplus175", R"(Nominal Hexaboard width + 175 #mum)"},
    {"PCBplus200", R"(Nominal Hexaboard width + 200 #mum)"},
    {"PCBplus225", R"(Nominal Hexaboard width + 225 #mum)"},
    {"PCBplus250", R"(Nominal Hexaboard width + 250 #mum)"},
    {"PCBplus275", R"(Nominal Hexaboard width + 275 #mum)"},
    {"PCBplus300", R"(Nominal Hexaboard width + 300 #mum)"},
    {"Baseplus025", R"(Nominal Baseplate width + 25 #mum)"},
    {"Baseplus050", R"(Nominal Baseplate width + 50 #mum)"},
    {"Kaptonplus000", R"(Nominal)"},
    {"Kaptonminus25", R"(Kapton - 25 #mum)"},
    {"Kaptonplus025", R"(Nominal Kapton width + 25 #mum)"},
    {"Kaptonminus50", R"(Kapton - 50 #mum)"},
    {"Kaptonplus050", R"(Nominal Kapton width + 50 #mum)"},
    {"Kaptonminus75", R"(Kapton - 75 #mum)"},
    {"Kaptonplus075", R"(Nominal Kapton width + 75 #mum)"},
    {"Kaptonminus150", R"(Kapton - 150 #mum)"},
    {"Kaptonminus170", R"(Kapton - 170 #mum)"},
    {"Kaptonminus75", R"(Kapton - 75 #mum)"},
    {"Kaptonplus100", R"(Nominal Kapton width + 100 #mum)"},
    {"Kaptonplus125", R"(Nominal Kapton width + 125 #mum)"},
    {"Kaptonplus150", R"(Nominal Kapton width + 150 #mum)"},
    {"Kaptonplus160", R"(Nominal Kapton width + 160 #mum)"},
    {"Kaptonplus170", R"(Nominal Kapton width + 170 #mum)"},
    {"Kaptonplus180", R"(Nominal Kapton width + 180 #mum)"},
    {"Kaptonplus190", R"(Nominal Kapton width + 190 #mum)"},
    {"Kaptonplus200", R"(Nominal Kapton width + 200 #mum)"},
    {"Kaptonplus210", R"(Nominal Kapton width + 210 #mum)"},
    {"Kaptonplus220", R"(Nominal Kapton width + 220 #mum)"},
    {"Kaptonplus230", R"(Nominal Kapton width + 230 #mum)"},
    {"Kaptonplus240", R"(Nominal Kapton width + 240 #mum)"},
    {"Kaptonplus0250", R"(Nominal Kapton width + 250 #mum)"},
    {"Kaptonplus260", R"(Nominal Kapton width + 260 #mum)"},
    {"Kaptonplus270", R"(Nominal Kapton width + 270 #mum)"},
    {"Kaptonplus280", R"(Nominal Kapton width + 280 #mum)"},
    {"Kaptonplus290", R"(Nominal Kapton width + 290 #mum)"},
    {"Kaptonplus300", R"(Nominal Kapton width + 300 #mum)"},
    {"Kaptonplus400", R"(Nominal Kapton width + 400 #mum)"},
    {"Kaptonplus500", R"(Nominal Kapton width + 500 #mum)"},
    {"senTokap000", R"(Incomplete coverage = 0 #mum)"},
    {"senTokap100", R"(Incomplete coverage = 100 #mum)"},
    {"senTokap150", R"(Incomplete coverage = 150 #mum)"},
    {"senTokap185", R"(Incomplete coverage = 185 #mum)"},
    {"Nominal", ""},
    {"Gaussian", "Gaussian"},
    {"Flat", "Flat"},
    {"CustomGaus", "Custom Gaussian"},
    {"CustomFlat", "Custom Flat"},
    {"newSensor", R"(Sen-to-kap placement = 50 #mum)"},
    {"midSensor", R"(Sen-to-kap placement = 100 #mum)"},
    {"oldSensor", R"(Sen-to-kap placement = 200 #mum)"},
    {"min", ""},
    {"max", ""},
    {"Fit", "Fit"},
    {"Worst", "Worst"},
    {"integrate", ""},
    {"Peak1", "Peak 1"},
    {"Peak2", "Peak 2"},
    {"Peak3", "Peak 3"},
    {"Peak4", "Peak 4"},
    {"Peak5", "Peak 5"},
    {"Peak6", "Peak 6"},
  };


}

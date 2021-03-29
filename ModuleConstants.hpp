// ModuleConstants.hpp
namespace constants {
  TString outputDir = "";
  const TString eosDir = "root://cmseos.fnal.gov//eos/uscms/store/user/mkilpatr/13TeV/ModuleTolerances_complete_300K_25Mar21/";
  const TString localDir = "ModuleTolerances_complete_300K_25Mar21_small";
  const bool debug = false;
  const TString whichOverlap = "sen_pcb_kap";

  //Define Constant Widths
  const double pcb_w_const = 166.64, sensor_w_const = 166.57, kapton_w_const = 166.94, baseplate_w_const = 166.94;

  double baseplate_center = 3.05, baseplate_err = 0.05;
  double baseplate_pinhole_err = 0.100, baseplate_pinhole_theta = 2*TMath::Pi();
  double pcb_w = pcb_w_const, sensor_w = sensor_w_const, kapton_w = kapton_w_const, baseplate_w = baseplate_w_const;
  double pcb_w_err = 0.1, sensor_w_err = 0.02, kapton_w_err = 0.05, baseplate_w_err = 0.05;
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

  //vector<string> Geometry = {"Full", "Five", "Semi", "Half", "Three"};
  vector<string> Geometry = {"Full"};
  vector<TString> whichComp = {"sen_kap_stack_hist", "sen_bas_stack_hist", "bas_kap_stack_hist", "sen_pcb_stack_hist", "pcb_bas_stack_hist", "pcb_kap_stack_hist", "kap_pcb_hist", "sen_pcb_hist", "sen_pcb_kap_x_hist", "sen_pcb_kap_y_hist"};
  vector<string> Dist = {
                           "Gaussian_Kaptonminus0_senTokap185_midSensor",              "Gaussian_Kaptonminus0_senTokap185_newSensor",
               "Gaussian_PCBplus100_Kaptonplus200_senTokap185_midSensor",  "Gaussian_PCBplus100_Kaptonplus100_senTokap185_newSensor",
               "Gaussian_PCBplus100_Kaptonplus300_senTokap185_midSensor",  "Gaussian_PCBplus100_Kaptonplus150_senTokap185_newSensor",
               "Gaussian_PCBplus100_Kaptonplus400_senTokap185_midSensor",  "Gaussian_PCBplus100_Kaptonplus200_senTokap185_newSensor",
                           };
  
  const vector<TString> Order = {
               "Gaussian_PCBplus100_Kaptonplus200_senTokap185_midSensor",  "Gaussian_PCBplus100_Kaptonplus100_senTokap185_newSensor",
               "Gaussian_PCBplus100_Kaptonplus300_senTokap185_midSensor",  "Gaussian_PCBplus100_Kaptonplus150_senTokap185_newSensor",
               "Gaussian_PCBplus100_Kaptonplus400_senTokap185_midSensor",  "Gaussian_PCBplus100_Kaptonplus200_senTokap185_newSensor",
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
    {"sen_kap_stack_hist", kBlue },
    {"sen_pcb_stack_hist", kGreen + 3 },
  };

  const std::map<TString, TString> latexMap{
    {"Custom", ""},
    {"Distribution", ""},
    {"PCBDist", ""},
    {"KaptonMultiDist", ""},
    {"TotalBadModules", ""},
    {"Sensor", ""},
    {"otherCenter", R"(center $= + 29 \mu m$)"},
    {"PCBminus0", R"(PCB - 0 $\mu m$)"},
    {"PCBminus25", R"(PCB - 12.5 $\mu m$)"},
    {"PCBplus25", R"(PCB + 12.5 $\mu m$)"},
    {"PCBminus50", R"(PCB - 25 $\mu m$)"},
    {"PCBplus50", R"(PCB + 25 $\mu m$)"},
    {"PCBminus75", R"(PCB - 37.5 $\mu m$)"},
    {"PCBplus75", R"(PCB + 37.5 $\mu m$)"},
    {"PCBplus100", R"(PCB + 50 $\mu m$)"},
    {"PCBplus125", R"(PCB + 62.5 $\mu m$)"},
    {"PCBplus150", R"(PCB + 75 $\mu m$)"},
    {"PCBplus175", R"(PCB + 87.5 $\mu m$)"},
    {"PCBplus200", R"(PCB + 100 $\mu m$)"},
    {"PCBplus225", R"(PCB + 112.5 $\mu m$)"},
    {"PCBplus250", R"(PCB + 125 $\mu m$)"},
    {"PCBplus275", R"(PCB + 137.5 $\mu m$)"},
    {"PCBplus300", R"(PCB + 150 $\mu m$)"},
    {"Kaptonminus0", R"(Nominal)"},
    {"Kaptonminus25", R"(Kapton - 25 $\mu m$)"},
    {"Kaptonplus25", R"(Nominal Kapton width + 25 $\mu m$)"},
    {"Kaptonminus50", R"(Kapton - 50 $\mu m$)"},
    {"Kaptonplus50", R"(Nominal Kapton width + 50 $\mu m$)"},
    {"Kaptonminus150", R"(Kapton - 150 $\mu m$)"},
    {"Kaptonminus170", R"(Kapton - 170 $\mu m$)"},
    {"Kaptonminus75", R"(Kapton - 75 $\mu m$)"},
    {"Kaptonplus75", R"(Nominal Kapton width + 75 $\mu m$)"},
    {"Kaptonplus100", R"(Nominal Kapton width + 100 $\mu m$)"},
    {"Kaptonplus125", R"(Nominal Kapton width + 125 $\mu m$)"},
    {"Kaptonplus150", R"(Nominal Kapton width + 150 $\mu m$)"},
    {"Kaptonplus175", R"(Nominal Kapton width + 175 $\mu m$)"},
    {"Kaptonplus200", R"(Nominal Kapton width + 200 $\mu m$)"},
    {"Kaptonplus225", R"(Nominal Kapton width + 225 $\mu m$)"},
    {"Kaptonplus300", R"(Nominal Kapton width + 300 $\mu m$)"},
    {"Kaptonplus400", R"(Nominal Kapton width + 400 $\mu m$)"},
    {"senTokap000", R"(\\Incomplete coverage = 0 $\mu m$}})"},
    {"senTokap100", R"(\\Incomplete coverage = 100 $\mu m$}})"},
    {"senTokap150", R"(\\Incomplete coverage = 150 $\mu m$}})"},
    {"senTokap185", R"(\\Incomplete coverage = 185 $\mu m$}})"},
    {"Nominal", "Nominal"},
    {"Gaussian", "Gaussian"},
    {"Flat", "Flat"},
    {"CustomGaus", "Custom Gaussian"},
    {"CustomFlat", "Custom Flat"},
    {"newSensor", R"([50 $\mu m$])"},
    {"midSensor", R"([100 $\mu m$])"},
    {"oldSensor", R"([200 $\mu m$])"},
    {"newSplitSensor", R"( & [50 $\mu m$])"},
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
  };
  
  const std::map<TString, TString> plotMap{
    {"Custom", ""},
    {"Distribution", ""},
    {"PCBDist", ""},
    {"KaptonMultiDist", ""},
    {"TotalBadModules", ""},
    {"Sensor", ""},
    {"ShieldBond", "Shield Bond"},
    {"GuardBond",  "Guard Bond"},
    {"BacksideBondX", "Backside Bond X"},
    {"BacksideBondY", "Backside Bond Y"},
    {"SenToBas", "Sen to Base"},
    {"BasToKapMax", "Max Base to Kap"},
    {"BasToKap", "Base to Kap"},
    {"SenToPCB", "PCB to Sen"},
    {"PCBToBas", "Base to PCB"},
    {"PCBToKap", "Kap to PCB"},
    {"SenToKap", "Kap to Sen"},
    {"otherCenter", R"(center $= + 29 #mum)"},
    {"PCBminus0", R"(PCB - 0 #mum)"},
    {"PCBminus25", R"(PCB - 12.5 #mum)"},
    {"PCBplus25", R"(PCB + 12.5 #mum)"},
    {"PCBminus50", R"(PCB - 25 #mum)"},
    {"PCBplus50", R"(PCB + 25 #mum)"},
    {"PCBminus75", R"(PCB - 37.5 #mum)"},
    {"PCBplus75", R"(PCB + 37.5 #mum)"},
    {"PCBplus100", R"(PCB + 50 #mum)"},
    {"PCBplus125", R"(PCB + 62.5 #mum)"},
    {"PCBplus150", R"(PCB + 75 #mum)"},
    {"PCBplus175", R"(PCB + 87.5 #mum)"},
    {"PCBplus200", R"(PCB + 100 #mum)"},
    {"PCBplus225", R"(PCB + 112.5 #mum)"},
    {"PCBplus250", R"(PCB + 125 #mum)"},
    {"PCBplus275", R"(PCB + 137.5 #mum)"},
    {"PCBplus300", R"(PCB + 150 #mum)"},
    {"Kaptonminus0", R"(Nominal)"},
    {"Kaptonminus25", R"(Kapton - 25 #mum)"},
    {"Kaptonplus25", R"(Nominal Kapton width + 25 #mum)"},
    {"Kaptonminus50", R"(Kapton - 50 #mum)"},
    {"Kaptonplus50", R"(Nominal Kapton width + 50 #mum)"},
    {"Kaptonminus150", R"(Kapton - 150 #mum)"},
    {"Kaptonminus170", R"(Kapton - 170 #mum)"},
    {"Kaptonminus75", R"(Kapton - 75 #mum)"},
    {"Kaptonplus75", R"(Nominal Kapton width + 75 #mum)"},
    {"Kaptonplus100", R"(Nominal Kapton width + 100 #mum)"},
    {"Kaptonplus125", R"(Nominal Kapton width + 125 #mum)"},
    {"Kaptonplus150", R"(Nominal Kapton width + 150 #mum)"},
    {"Kaptonplus175", R"(Nominal Kapton width + 175 #mum)"},
    {"Kaptonplus200", R"(Nominal Kapton width + 200 #mum)"},
    {"Kaptonplus225", R"(Nominal Kapton width + 225 #mum)"},
    {"Kaptonplus300", R"(Nominal Kapton width + 300 #mum)"},
    {"Kaptonplus400", R"(Nominal Kapton width + 400 #mum)"},
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
    {"Peak1", ""},
    {"Peak2", ""},
    {"Peak3", ""},
  };


}

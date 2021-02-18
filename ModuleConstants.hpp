// ModuleConstants.hpp
namespace constants {
  TString outputDir = "";
  const TString eosDir = "root://cmseos.fnal.gov//eos/uscms/store/user/mkilpatr/13TeV/ModuleTolerances_complete_100K_021821/";
  const TString baseDir = "Gaussian_Kaptonminus0_oldSensor";
  const TString localDir = "ModuleTolerances_complete_100K_021821";
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
  pair<double, double> backside_x_err = make_pair(0.2133, 0.1769);
  pair<double, double> backside_y_err = make_pair(0.250, 0.000);
  
  int max = 100000;
  double weight = 30000./double(max);

  int nbins = 650;
  double step = 0.002;
  double width_new = 1.4;
  double width_backY = 4.0;
  double axis = width_new/2 + nbins*step + 0.1;

  //vector<string> Geometry = {"Full", "Five", "Semi", "Half", "Three"};
  vector<string> Geometry = {"Full"};
  vector<TString> whichComp = {"sen_kap_stack_hist", "sen_bas_stack_hist", "bas_kap_stack_hist", "sen_pcb_stack_hist", "pcb_bas_stack_hist", "pcb_kap_stack_hist", "kap_pcb_hist", "sen_pcb_hist", "sen_pcb_kap_x_hist", "sen_pcb_kap_y_hist"};
  vector<string> Dist = {"Gaussian_Kaptonminus0_oldSensor"};
  //vector<string> Dist = {"Gaussian_Kaptonminus0_oldSensor", "Flat_Kaptonminus0_oldSensor", "Gaussian_Kaptonminus0_newSensor", "Gaussian_Kaptonminus0_midSensor", "CustomGaus_Kaptonminus0_oldSensor", "CustomFlat_Kaptonminus0_oldSensor", "Gaussian_Kaptonplus175_oldSensor", "Gaussian_Kaptonplus200_oldSensor", "Gaussian_Kaptonplus225_oldSensor", "Gaussian_Kaptonplus100_oldSensor", "Gaussian_Kaptonplus125_oldSensor", "Gaussian_Kaptonplus150_oldSensor", "Gaussian_Kaptonplus25_oldSensor", "Gaussian_Kaptonplus50_oldSensor", "Gaussian_Kaptonplus75_oldSensor", "Gaussian_Kaptonminus25_oldSensor", "Gaussian_Kaptonminus50_oldSensor", "Gaussian_Kaptonminus75_oldSensor", "Gaussian_Kaptonminus150_oldSensor", "Gaussian_Kaptonminus170_oldSensor", "Gaussian_Kaptonplus175_newSensor", "Gaussian_Kaptonplus200_newSensor", "Gaussian_Kaptonplus225_newSensor", "Gaussian_Kaptonplus100_newSensor", "Gaussian_Kaptonplus125_newSensor", "Gaussian_Kaptonplus150_newSensor", "Gaussian_Kaptonplus25_newSensor", "Gaussian_Kaptonplus50_newSensor", "Gaussian_Kaptonplus75_newSensor", "Gaussian_Kaptonminus25_newSensor", "Gaussian_Kaptonminus50_newSensor", "Gaussian_Kaptonminus75_newSensor", "Gaussian_Kaptonminus150_newSensor", "Gaussian_Kaptonminus170_newSensor", "Gaussian_Kaptonplus175_midSensor", "Gaussian_Kaptonplus200_midSensor", "Gaussian_Kaptonplus225_midSensor", "Gaussian_Kaptonplus100_midSensor", "Gaussian_Kaptonplus125_midSensor", "Gaussian_Kaptonplus150_midSensor", "Gaussian_Kaptonplus25_midSensor", "Gaussian_Kaptonplus50_midSensor", "Gaussian_Kaptonplus75_midSensor", "Gaussian_Kaptonminus25_midSensor", "Gaussian_Kaptonminus50_midSensor", "Gaussian_Kaptonminus75_midSensor", "Gaussian_Kaptonminus150_midSensor", "Gaussian_Kaptonminus170_midSensor", "Gaussian_PCBplus25_oldSensor", "Gaussian_PCBplus50_oldSensor", "Gaussian_PCBplus75_oldSensor", "Gaussian_PCBplus25_newSensor", "Gaussian_PCBplus50_newSensor", "Gaussian_PCBplus75_newSensor", "Gaussian_PCBplus25_midSensor", "Gaussian_PCBplus50_midSensor", "Gaussian_PCBplus75_midSensor", "Gaussian_PCBminus0_oldSensor", "Gaussian_PCBminus0_newSensor", "Gaussian_PCBminus0_midSensor", "Gaussian_PCBminus25_oldSensor", "Gaussian_PCBminus50_oldSensor", "Gaussian_PCBminus75_oldSensor", "Gaussian_PCBminus25_newSensor", "Gaussian_PCBminus50_newSensor", "Gaussian_PCBminus75_newSensor", "Gaussian_PCBminus25_midSensor", "Gaussian_PCBminus50_midSensor", "Gaussian_PCBminus75_midSensor", "Gaussian_PCBplus125_Kaptonplus175_oldSensor", "Gaussian_PCBplus125_Kaptonplus200_oldSensor", "Gaussian_PCBplus125_Kaptonplus225_oldSensor", "Gaussian_PCBplus100_Kaptonplus175_oldSensor", "Gaussian_PCBplus100_Kaptonplus200_oldSensor", "Gaussian_PCBplus100_Kaptonplus225_oldSensor", "Gaussian_PCBplus75_Kaptonplus175_oldSensor", "Gaussian_PCBplus75_Kaptonplus200_oldSensor", "Gaussian_PCBplus75_Kaptonplus225_oldSensor", "Gaussian_PCBplus50_Kaptonplus175_oldSensor", "Gaussian_PCBplus50_Kaptonplus200_oldSensor", "Gaussian_PCBplus50_Kaptonplus225_oldSensor", "Gaussian_PCBplus25_Kaptonplus175_oldSensor", "Gaussian_PCBplus25_Kaptonplus200_oldSensor", "Gaussian_PCBplus25_Kaptonplus225_oldSensor", "Gaussian_PCBplus125_Kaptonplus175_midSensor", "Gaussian_PCBplus125_Kaptonplus200_midSensor", "Gaussian_PCBplus125_Kaptonplus225_midSensor", "Gaussian_PCBplus100_Kaptonplus175_midSensor", "Gaussian_PCBplus100_Kaptonplus200_midSensor", "Gaussian_PCBplus100_Kaptonplus225_midSensor", "Gaussian_PCBplus75_Kaptonplus175_midSensor", "Gaussian_PCBplus75_Kaptonplus200_midSensor", "Gaussian_PCBplus75_Kaptonplus225_midSensor", "Gaussian_PCBplus50_Kaptonplus175_midSensor", "Gaussian_PCBplus50_Kaptonplus200_midSensor", "Gaussian_PCBplus50_Kaptonplus225_midSensor", "Gaussian_PCBplus25_Kaptonplus175_midSensor", "Gaussian_PCBplus25_Kaptonplus200_midSensor", "Gaussian_PCBplus25_Kaptonplus225_midSensor", "Gaussian_PCBplus125_Kaptonplus175_newSensor", "Gaussian_PCBplus125_Kaptonplus200_newSensor", "Gaussian_PCBplus125_Kaptonplus225_newSensor", "Gaussian_PCBplus100_Kaptonplus175_newSensor", "Gaussian_PCBplus100_Kaptonplus200_newSensor", "Gaussian_PCBplus100_Kaptonplus225_newSensor", "Gaussian_PCBplus75_Kaptonplus175_newSensor", "Gaussian_PCBplus75_Kaptonplus200_newSensor", "Gaussian_PCBplus75_Kaptonplus225_newSensor", "Gaussian_PCBplus50_Kaptonplus175_newSensor", "Gaussian_PCBplus50_Kaptonplus200_newSensor", "Gaussian_PCBplus50_Kaptonplus225_newSensor", "Gaussian_PCBplus25_Kaptonplus175_newSensor", "Gaussian_PCBplus25_Kaptonplus200_newSensor", "Gaussian_PCBplus25_Kaptonplus225_newSensor"};

  const vector<TString> Order = {
  			 "Gaussian_Kaptonminus170_oldSensor",            "Gaussian_Kaptonminus170_midSensor",            "Gaussian_Kaptonminus170_newSensor",  
  			 "Gaussian_Kaptonminus150_oldSensor",            "Gaussian_Kaptonminus150_midSensor",            "Gaussian_Kaptonminus150_newSensor",  
  			  "Gaussian_Kaptonminus75_oldSensor",             "Gaussian_Kaptonminus75_midSensor",             "Gaussian_Kaptonminus75_newSensor",  
  			  "Gaussian_Kaptonminus50_oldSensor",             "Gaussian_Kaptonminus50_midSensor",             "Gaussian_Kaptonminus50_newSensor",  
  			  "Gaussian_Kaptonminus25_oldSensor",             "Gaussian_Kaptonminus25_midSensor",             "Gaussian_Kaptonminus25_newSensor",  
  			   "Gaussian_Kaptonminus0_oldSensor",              "Gaussian_Kaptonminus0_midSensor",              "Gaussian_Kaptonminus0_newSensor", 
  			   "Gaussian_Kaptonplus25_oldSensor",              "Gaussian_Kaptonplus25_midSensor",              "Gaussian_Kaptonplus25_newSensor",  
  			   "Gaussian_Kaptonplus50_oldSensor",              "Gaussian_Kaptonplus50_midSensor",              "Gaussian_Kaptonplus50_newSensor",  
  			   "Gaussian_Kaptonplus75_oldSensor",              "Gaussian_Kaptonplus75_midSensor",              "Gaussian_Kaptonplus75_newSensor",
  			  "Gaussian_Kaptonplus100_oldSensor",             "Gaussian_Kaptonplus100_midSensor",             "Gaussian_Kaptonplus100_newSensor",
  			  "Gaussian_Kaptonplus125_oldSensor",             "Gaussian_Kaptonplus125_midSensor",             "Gaussian_Kaptonplus125_newSensor",
  			  "Gaussian_Kaptonplus150_oldSensor",             "Gaussian_Kaptonplus150_midSensor",             "Gaussian_Kaptonplus150_newSensor",
  			  "Gaussian_Kaptonplus175_oldSensor",             "Gaussian_Kaptonplus175_midSensor",             "Gaussian_Kaptonplus175_newSensor",
  			  "Gaussian_Kaptonplus200_oldSensor",             "Gaussian_Kaptonplus200_midSensor",             "Gaussian_Kaptonplus200_newSensor",
  			  "Gaussian_Kaptonplus225_oldSensor",             "Gaussian_Kaptonplus225_midSensor",             "Gaussian_Kaptonplus225_newSensor",
  			     "Gaussian_PCBminus75_oldSensor",                "Gaussian_PCBminus75_midSensor",                "Gaussian_PCBminus75_newSensor",  
  			     "Gaussian_PCBminus50_oldSensor",                "Gaussian_PCBminus50_midSensor",                "Gaussian_PCBminus50_newSensor",  
  			     "Gaussian_PCBminus25_oldSensor",                "Gaussian_PCBminus25_midSensor",                "Gaussian_PCBminus25_newSensor",  
  			      "Gaussian_PCBminus0_oldSensor",                 "Gaussian_PCBminus0_midSensor",                 "Gaussian_PCBminus0_newSensor", 
  			      "Gaussian_PCBplus25_oldSensor",                 "Gaussian_PCBplus25_midSensor",                 "Gaussian_PCBplus25_newSensor",  
  			      "Gaussian_PCBplus50_oldSensor",                 "Gaussian_PCBplus50_midSensor",                 "Gaussian_PCBplus50_newSensor",  
  			      "Gaussian_PCBplus75_oldSensor",                 "Gaussian_PCBplus75_midSensor",                 "Gaussian_PCBplus75_newSensor",
                "Gaussian_PCBplus25_Kaptonplus175_oldSensor",   "Gaussian_PCBplus25_Kaptonplus175_midSensor",   "Gaussian_PCBplus25_Kaptonplus175_newSensor",
                "Gaussian_PCBplus25_Kaptonplus200_oldSensor",   "Gaussian_PCBplus25_Kaptonplus200_midSensor",   "Gaussian_PCBplus25_Kaptonplus200_newSensor",
                "Gaussian_PCBplus25_Kaptonplus225_oldSensor",   "Gaussian_PCBplus25_Kaptonplus225_midSensor",   "Gaussian_PCBplus25_Kaptonplus225_newSensor",
                "Gaussian_PCBplus50_Kaptonplus175_oldSensor",   "Gaussian_PCBplus50_Kaptonplus175_midSensor",   "Gaussian_PCBplus50_Kaptonplus175_newSensor",
                "Gaussian_PCBplus50_Kaptonplus200_oldSensor",   "Gaussian_PCBplus50_Kaptonplus200_midSensor",   "Gaussian_PCBplus50_Kaptonplus200_newSensor",
                "Gaussian_PCBplus50_Kaptonplus225_oldSensor",   "Gaussian_PCBplus50_Kaptonplus225_midSensor",   "Gaussian_PCBplus50_Kaptonplus225_newSensor",
                "Gaussian_PCBplus75_Kaptonplus175_oldSensor",   "Gaussian_PCBplus75_Kaptonplus175_midSensor",   "Gaussian_PCBplus75_Kaptonplus175_newSensor",
                "Gaussian_PCBplus75_Kaptonplus200_oldSensor",   "Gaussian_PCBplus75_Kaptonplus200_midSensor",   "Gaussian_PCBplus75_Kaptonplus200_newSensor",
                "Gaussian_PCBplus75_Kaptonplus225_oldSensor",   "Gaussian_PCBplus75_Kaptonplus225_midSensor",   "Gaussian_PCBplus75_Kaptonplus225_newSensor",
               "Gaussian_PCBplus100_Kaptonplus175_oldSensor",  "Gaussian_PCBplus100_Kaptonplus175_midSensor",  "Gaussian_PCBplus100_Kaptonplus175_newSensor",
               "Gaussian_PCBplus100_Kaptonplus200_oldSensor",  "Gaussian_PCBplus100_Kaptonplus200_midSensor",  "Gaussian_PCBplus100_Kaptonplus200_newSensor",
               "Gaussian_PCBplus100_Kaptonplus225_oldSensor",  "Gaussian_PCBplus100_Kaptonplus225_midSensor",  "Gaussian_PCBplus100_Kaptonplus225_newSensor",
               "Gaussian_PCBplus125_Kaptonplus175_oldSensor",  "Gaussian_PCBplus125_Kaptonplus175_midSensor",  "Gaussian_PCBplus125_Kaptonplus175_newSensor",
               "Gaussian_PCBplus125_Kaptonplus200_oldSensor",  "Gaussian_PCBplus125_Kaptonplus200_midSensor",  "Gaussian_PCBplus125_Kaptonplus200_newSensor",
               "Gaussian_PCBplus125_Kaptonplus225_oldSensor",  "Gaussian_PCBplus125_Kaptonplus225_midSensor",  "Gaussian_PCBplus125_Kaptonplus225_newSensor",
  		               "Flat_Kaptonminus0_oldSensor",            "CustomGaus_Kaptonminus0_oldSensor",            "CustomFlat_Kaptonminus0_oldSensor", 
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
    {"PCBminus25", R"(PCB - 25 $\mu m$)"},
    {"PCBplus25", R"(PCB + 25 $\mu m$)"},
    {"PCBminus50", R"(PCB - 50 $\mu m$)"},
    {"PCBplus50", R"(PCB + 50 $\mu m$)"},
    {"PCBminus75", R"(PCB - 75 $\mu m$)"},
    {"PCBplus75", R"(PCB + 75 $\mu m$)"},
    {"PCBplus100", R"(PCB + 100 $\mu m$)"},
    {"PCBplus125", R"(PCB + 125 $\mu m$)"},
    {"Kaptonminus0", R"(Kapton - 0 $\mu m$)"},
    {"Kaptonminus25", R"(Kapton - 25 $\mu m$)"},
    {"Kaptonplus25", R"(Kapton + 25 $\mu m$)"},
    {"Kaptonminus50", R"(Kapton - 50 $\mu m$)"},
    {"Kaptonplus50", R"(Kapton + 50 $\mu m$)"},
    {"Kaptonminus150", R"(Kapton - 150 $\mu m$)"},
    {"Kaptonminus170", R"(Kapton - 170 $\mu m$)"},
    {"Kaptonminus75", R"(Kapton - 75 $\mu m$)"},
    {"Kaptonplus75", R"(Kapton + 75 $\mu m$)"},
    {"Kaptonplus100", R"(Kapton + 100 $\mu m$)"},
    {"Kaptonplus125", R"(Kapton + 125 $\mu m$)"},
    {"Kaptonplus150", R"(Kapton + 150 $\mu m$)"},
    {"Kaptonplus175", R"(Kapton + 175 $\mu m$)"},
    {"Kaptonplus200", R"(Kapton + 200 $\mu m$)"},
    {"Kaptonplus225", R"(Kapton + 225 $\mu m$)"},
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
    {"oldSplitSensor", R"( & [200 $\mu m$])"},
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
    {"otherCenter", R"(center $= + 29 #mum)"},
    {"PCBminus0", R"(PCB - 0 #mum)"},
    {"PCBminus25", R"(PCB - 25 #mum)"},
    {"PCBplus25", R"(PCB + 25 #mum)"},
    {"PCBminus50", R"(PCB - 50 #mum)"},
    {"PCBplus50", R"(PCB + 50 #mum)"},
    {"PCBminus75", R"(PCB - 75 #mum)"},
    {"PCBplus75", R"(PCB + 75 #mum)"},
    {"PCBplus100", R"(PCB + 100 #mum)"},
    {"PCBplus125", R"(PCB + 125 #mum)"},
    {"Kaptonminus0", R"(Kapton - 0 #mum)"},
    {"Kaptonminus25", R"(Kapton - 25 #mum)"},
    {"Kaptonplus25", R"(Kapton + 25 #mum)"},
    {"Kaptonminus50", R"(Kapton - 50 #mum)"},
    {"Kaptonplus50", R"(Kapton + 50 #mum)"},
    {"Kaptonminus150", R"(Kapton - 150 #mum)"},
    {"Kaptonminus170", R"(Kapton - 170 #mum)"},
    {"Kaptonminus75", R"(Kapton - 75 #mum)"},
    {"Kaptonplus75", R"(Kapton + 75 #mum)"},
    {"Kaptonplus100", R"(Kapton + 100 #mum)"},
    {"Kaptonplus125", R"(Kapton + 125 #mum)"},
    {"Kaptonplus150", R"(Kapton + 150 #mum)"},
    {"Kaptonplus175", R"(Kapton + 175 #mum)"},
    {"Kaptonplus200", R"(Kapton + 200 #mum)"},
    {"Kaptonplus225", R"(Kapton + 225 #mum)"},
    {"Nominal", ""},
    {"Gaussian", "Gaussian"},
    {"Flat", "Flat"},
    {"CustomGaus", "Custom Gaussian"},
    {"CustomFlat", "Custom Flat"},
    {"newSensor", R"([50 #mum])"},
    {"midSensor", R"([100 #mum])"},
    {"oldSensor", R"([200 #mum])"},
    {"Fit", "Fit"},
    {"Worst", "Worst"},
    {"integrate", ""},
    {"Peak1", "Peak 1"},
    {"Peak2", "Peak 2"},
    {"Peak3", "Peak 3"},
  };


}

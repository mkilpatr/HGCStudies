// ModuleConstants.hpp
namespace constants {
  TString outputDir = "";
  const TString baseDir = "moduleTolerances_complete_020821_TestAll";
  const bool debug = false;
  const TString whichOverlap = "sen_pcb_kap";

  //Define Constant Widths
  const double pcb_w_const = 166.64, sensor_w_const = 166.57, kapton_w_const = 166.94, baseplate_w_const = 166.94;

  double baseplate_center = 3.05, baseplate_err = 0.05;
  double baseplate_pinhole_err = 0.100, baseplate_pinhole_theta = 2*TMath::Pi();
  double pcb_w = pcb_w_const, sensor_w = sensor_w_const, kapton_w = kapton_w_const, baseplate_w = baseplate_w_const;
  double pcb_w_err = 0.1, sensor_w_err = 0.02, kapton_w_err = 0.05, baseplate_w_err = 0.05;
  double PcbToSen_shift_r = 0.04, PcbToSen_shift_theta = 2*TMath::Pi();
  double sensor_shift_x= 0.500, sensor_shift_y = 0.050, sensor_shift_theta = 0.050; //theta is radians
  double baseValue_total = sensor_w;
  double base_kap_min = 0.100;
  double gold_min = 0.25, gold_size = 1.500, backside_min = 0.250;
  pair<double, double> backside_x_err = make_pair(0.2133, 0.1769);
  pair<double, double> backside_y_err = make_pair(0.250, 0.000);
  
  int max = 3000;

  int nbins = 650;
  double step = 0.002;
  double width_new = 1.4;
  double width_backY = 4.0;
  double axis = width_new/2 + nbins*step + 0.1;

  //vector<string> Geometry = {"Full", "Five", "Semi", "Half", "Three"};
  vector<string> Dist = {"Gaussian_Kaptonminus0_oldSensor", "Flat", "Gaussian_Kaptonminus0_newSensor", "Gaussian_Kaptonminus0_midSensor",
                         "CustomGaus", "CustomFlat",
                         "Gaussian_Kaptonplus25_oldSensor", "Gaussian_Kaptonplus50_oldSensor", "Gaussian_Kaptonplus75_oldSensor",
                         "Gaussian_Kaptonminus25_oldSensor", "Gaussian_Kaptonminus50_oldSensor", "Gaussian_Kaptonminus75_oldSensor", "Gaussian_Kaptonminus150_oldSensor",  "Gaussian_Kaptonminus170_oldSensor",
                         "Gaussian_Kaptonplus25_newSensor", "Gaussian_Kaptonplus50_newSensor", "Gaussian_Kaptonplus75_newSensor",
                         "Gaussian_Kaptonminus25_newSensor", "Gaussian_Kaptonminus50_newSensor", "Gaussian_Kaptonminus75_newSensor", "Gaussian_Kaptonminus150_newSensor",  "Gaussian_Kaptonminus170_newSensor",
                         "Gaussian_Kaptonplus25_midSensor", "Gaussian_Kaptonplus50_midSensor", "Gaussian_Kaptonplus75_midSensor",
                         "Gaussian_Kaptonminus25_midSensor", "Gaussian_Kaptonminus50_midSensor", "Gaussian_Kaptonminus75_midSensor", "Gaussian_Kaptonminus150_midSensor",  "Gaussian_Kaptonminus170_midSensor",
                         "Gaussian_PCBplus25_oldSensor", "Gaussian_PCBplus50_oldSensor", "Gaussian_PCBplus75_oldSensor",
                         "Gaussian_PCBplus25_newSensor", "Gaussian_PCBplus50_newSensor", "Gaussian_PCBplus75_newSensor",
                         "Gaussian_PCBplus25_midSensor", "Gaussian_PCBplus50_midSensor", "Gaussian_PCBplus75_midSensor",
                         "Gaussian_PCBminus0_oldSensor", "Gaussian_PCBminus0_newSensor", "Gaussian_PCBminus0_midSensor",
                         "Gaussian_PCBminus25_oldSensor", "Gaussian_PCBminus50_oldSensor", "Gaussian_PCBminus75_oldSensor",
                         "Gaussian_PCBminus25_newSensor", "Gaussian_PCBminus50_newSensor", "Gaussian_PCBminus75_newSensor",
                         "Gaussian_PCBminus25_midSensor", "Gaussian_PCBminus50_midSensor", "Gaussian_PCBminus75_midSensor",
                        }; 
  //vector<string> Dist = {"Gaussian_Kaptonminus0_oldSensor"};
  vector<string> Geometry = {"Full"};

  const vector<TString> Order = {
  			 "Gaussian_Kaptonminus170_newSensor", "Gaussian_Kaptonminus170_midSensor", "Gaussian_Kaptonminus170_oldSensor",  
  			 "Gaussian_Kaptonminus150_newSensor", "Gaussian_Kaptonminus150_midSensor", "Gaussian_Kaptonminus150_oldSensor",  
  			 "Gaussian_Kaptonminus75_newSensor", "Gaussian_Kaptonminus75_midSensor", "Gaussian_Kaptonminus75_oldSensor",  
  			 "Gaussian_Kaptonminus50_newSensor", "Gaussian_Kaptonminus50_midSensor", "Gaussian_Kaptonminus50_oldSensor",  
  			 "Gaussian_Kaptonminus25_newSensor", "Gaussian_Kaptonminus25_midSensor", "Gaussian_Kaptonminus25_oldSensor",  
  			 "Gaussian_Kaptonminus0_newSensor", "Gaussian_Kaptonminus0_midSensor", "Gaussian_Kaptonminus0_oldSensor", 
  			 "Gaussian_Kaptonplus25_newSensor", "Gaussian_Kaptonplus25_midSensor", "Gaussian_Kaptonplus25_oldSensor",  
  			 "Gaussian_Kaptonplus50_newSensor", "Gaussian_Kaptonplus50_midSensor", "Gaussian_Kaptonplus50_oldSensor",  
  			 "Gaussian_Kaptonplus75_newSensor", "Gaussian_Kaptonplus75_midSensor", "Gaussian_Kaptonplus75_oldSensor",
  			 "Gaussian_PCBminus75_newSensor", "Gaussian_PCBminus75_midSensor", "Gaussian_PCBminus75_oldSensor",  
  			 "Gaussian_PCBminus50_newSensor", "Gaussian_PCBminus50_midSensor", "Gaussian_PCBminus50_oldSensor",  
  			 "Gaussian_PCBminus25_newSensor", "Gaussian_PCBminus25_midSensor", "Gaussian_PCBminus25_oldSensor",  
  			 "Gaussian_PCBminus0_newSensor", "Gaussian_PCBminus0_midSensor", "Gaussian_PCBminus0_oldSensor", 
  			 "Gaussian_PCBplus25_newSensor", "Gaussian_PCBplus25_midSensor", "Gaussian_PCBplus25_oldSensor",  
  			 "Gaussian_PCBplus50_newSensor", "Gaussian_PCBplus50_midSensor", "Gaussian_PCBplus50_oldSensor",  
  			 "Gaussian_PCBplus75_newSensor", "Gaussian_PCBplus75_midSensor", "Gaussian_PCBplus75_oldSensor",
  			 "CustomGaus", "CustomFlat", 
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
    {"Flat", "Flat"},
    {"CustomGaus", "Custom Gaussian"},
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


}

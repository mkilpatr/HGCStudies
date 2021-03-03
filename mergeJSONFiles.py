import os
import json

# json file witll_BkgPred_v2e bkg predictions and signal yields
def write_json(data, filename='ModuleTolerances_complete_100K_030321_small/TestFull/Full_BadmoduleTolerances.json'):
    with open(filename,'w') as f:
    	json.dump(data, f, indent=2, ensure_ascii=False)

eosDir = "/eos/uscms/store/user/mkilpatr/13TeV/ModuleTolerances_complete_100K_030321/"
Dist = ["Gaussian_Kaptonminus0_oldSensor", "Gaussian_Kaptonminus0_newSensor", "Gaussian_Kaptonminus0_midSensor",
        "Gaussian_PCBplus75_Kaptonplus175_oldSensor",   "Gaussian_PCBplus75_Kaptonplus175_midSensor",   "Gaussian_PCBplus75_Kaptonplus175_newSensor",
        "Gaussian_PCBplus75_Kaptonplus200_oldSensor",   "Gaussian_PCBplus75_Kaptonplus200_midSensor",   "Gaussian_PCBplus75_Kaptonplus200_newSensor",
        "Gaussian_PCBplus75_Kaptonplus225_oldSensor",   "Gaussian_PCBplus75_Kaptonplus225_midSensor",   "Gaussian_PCBplus75_Kaptonplus225_newSensor",
       "Gaussian_PCBplus100_Kaptonplus175_oldSensor",  "Gaussian_PCBplus100_Kaptonplus175_midSensor",  "Gaussian_PCBplus100_Kaptonplus175_newSensor",
       "Gaussian_PCBplus100_Kaptonplus200_oldSensor",  "Gaussian_PCBplus100_Kaptonplus200_midSensor",  "Gaussian_PCBplus100_Kaptonplus200_newSensor",
       "Gaussian_PCBplus100_Kaptonplus225_oldSensor",  "Gaussian_PCBplus100_Kaptonplus225_midSensor",  "Gaussian_PCBplus100_Kaptonplus225_newSensor",
       "Gaussian_PCBplus125_Kaptonplus175_oldSensor",  "Gaussian_PCBplus125_Kaptonplus175_midSensor",  "Gaussian_PCBplus125_Kaptonplus175_newSensor",
       "Gaussian_PCBplus125_Kaptonplus200_oldSensor",  "Gaussian_PCBplus125_Kaptonplus200_midSensor",  "Gaussian_PCBplus125_Kaptonplus200_newSensor",
       "Gaussian_PCBplus125_Kaptonplus225_oldSensor",  "Gaussian_PCBplus125_Kaptonplus225_midSensor",  "Gaussian_PCBplus125_Kaptonplus225_newSensor",
       "Gaussian_PCBplus150_Kaptonplus175_oldSensor",  "Gaussian_PCBplus150_Kaptonplus175_midSensor",  "Gaussian_PCBplus150_Kaptonplus175_newSensor",
       "Gaussian_PCBplus150_Kaptonplus200_oldSensor",  "Gaussian_PCBplus150_Kaptonplus200_midSensor",  "Gaussian_PCBplus150_Kaptonplus200_newSensor",
       "Gaussian_PCBplus150_Kaptonplus225_oldSensor",  "Gaussian_PCBplus150_Kaptonplus225_midSensor",  "Gaussian_PCBplus150_Kaptonplus225_newSensor",
       "Gaussian_PCBplus175_Kaptonplus175_oldSensor",  "Gaussian_PCBplus175_Kaptonplus175_midSensor",  "Gaussian_PCBplus175_Kaptonplus175_newSensor",
       "Gaussian_PCBplus175_Kaptonplus200_oldSensor",  "Gaussian_PCBplus175_Kaptonplus200_midSensor",  "Gaussian_PCBplus175_Kaptonplus200_newSensor",
       "Gaussian_PCBplus175_Kaptonplus225_oldSensor",  "Gaussian_PCBplus175_Kaptonplus225_midSensor",  "Gaussian_PCBplus175_Kaptonplus225_newSensor",
       "Gaussian_PCBplus200_Kaptonplus175_oldSensor",  "Gaussian_PCBplus200_Kaptonplus175_midSensor",  "Gaussian_PCBplus200_Kaptonplus175_newSensor",
       "Gaussian_PCBplus200_Kaptonplus200_oldSensor",  "Gaussian_PCBplus200_Kaptonplus200_midSensor",  "Gaussian_PCBplus200_Kaptonplus200_newSensor",
       "Gaussian_PCBplus200_Kaptonplus225_oldSensor",  "Gaussian_PCBplus200_Kaptonplus225_midSensor",  "Gaussian_PCBplus200_Kaptonplus225_newSensor",
       "Gaussian_PCBplus225_Kaptonplus175_oldSensor",  "Gaussian_PCBplus225_Kaptonplus175_midSensor",  "Gaussian_PCBplus225_Kaptonplus175_newSensor",
       "Gaussian_PCBplus225_Kaptonplus200_oldSensor",  "Gaussian_PCBplus225_Kaptonplus200_midSensor",  "Gaussian_PCBplus225_Kaptonplus200_newSensor",
       "Gaussian_PCBplus225_Kaptonplus225_oldSensor",  "Gaussian_PCBplus225_Kaptonplus225_midSensor",  "Gaussian_PCBplus225_Kaptonplus225_newSensor",
       "Gaussian_PCBplus250_Kaptonplus175_oldSensor",  "Gaussian_PCBplus250_Kaptonplus175_midSensor",  "Gaussian_PCBplus250_Kaptonplus175_newSensor",
       "Gaussian_PCBplus250_Kaptonplus200_oldSensor",  "Gaussian_PCBplus250_Kaptonplus200_midSensor",  "Gaussian_PCBplus250_Kaptonplus200_newSensor",
       "Gaussian_PCBplus250_Kaptonplus225_oldSensor",  "Gaussian_PCBplus250_Kaptonplus225_midSensor",  "Gaussian_PCBplus250_Kaptonplus225_newSensor",
       "Gaussian_PCBplus275_Kaptonplus175_oldSensor",  "Gaussian_PCBplus275_Kaptonplus175_midSensor",  "Gaussian_PCBplus275_Kaptonplus175_newSensor",
       "Gaussian_PCBplus275_Kaptonplus200_oldSensor",  "Gaussian_PCBplus275_Kaptonplus200_midSensor",  "Gaussian_PCBplus275_Kaptonplus200_newSensor",
       "Gaussian_PCBplus275_Kaptonplus225_oldSensor",  "Gaussian_PCBplus275_Kaptonplus225_midSensor",  "Gaussian_PCBplus275_Kaptonplus225_newSensor",
       "Gaussian_PCBplus300_Kaptonplus175_oldSensor",  "Gaussian_PCBplus300_Kaptonplus175_midSensor",  "Gaussian_PCBplus300_Kaptonplus175_newSensor",
       "Gaussian_PCBplus300_Kaptonplus200_oldSensor",  "Gaussian_PCBplus300_Kaptonplus200_midSensor",  "Gaussian_PCBplus300_Kaptonplus200_newSensor",
       "Gaussian_PCBplus300_Kaptonplus225_oldSensor",  "Gaussian_PCBplus300_Kaptonplus225_midSensor",  "Gaussian_PCBplus300_Kaptonplus225_newSensor"]

with open("/eos/uscms/store/user/mkilpatr/13TeV/ModuleTolerances_complete_100K_030321/Gaussian_Kaptonminus0_oldSensor_Full/Full_BadmoduleTolerances.json", "a+") as new:
    newData = json.load(new)
    newData['Worst']["Bad Components"] = {}
    newData['Worst']["Bad Overlaps"] = {}

for d in Dist:
    with open(eosDir+d+"_Full/Full_BadmoduleTolerances.json", "r") as lepcr:
        lep_insert = json.load(lepcr)
        newData['Worst']["Bad Overlaps"][d] = {}
        newData['Worst']["Bad Components"].update(lep_insert['Worst']["Bad Components"])
        newData['Worst']["Bad Overlaps"][d].update(lep_insert['Worst']["Bad Overlaps"][d])
write_json(newData)


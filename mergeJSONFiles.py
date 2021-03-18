import os
import json
import argparse

parser = argparse.ArgumentParser(
        description='Produce or print limits based on existing datacards')
parser.add_argument("-d", "--dir", dest="baseDir", default='ModuleTolerances_complete_100K',
                         help="Search region to use when running the maximum likelihood fit. [Default: Search region with 183]")
args = parser.parse_args()


# json file witll_BkgPred_v2e bkg predictions and signal yields
def write_json(data, filename=args.baseDir + '_small/TestFull/Full_BadmoduleTolerances.json'):
    with open(filename,'w') as f:
    	json.dump(data, f, indent=2, ensure_ascii=False)

eosDir = "/eos/uscms/store/user/mkilpatr/13TeV/" + args.baseDir + "/"
Dist = [
            "Gaussian_Kaptonminus0_senTokap150_oldSensor",              "Gaussian_Kaptonminus0_senTokap150_midSensor",              "Gaussian_Kaptonminus0_senTokap150_newSensor",
"Gaussian_PCBplus100_Kaptonplus100_senTokap150_oldSensor",  "Gaussian_PCBplus100_Kaptonplus100_senTokap150_midSensor",  "Gaussian_PCBplus100_Kaptonplus100_senTokap150_newSensor",
"Gaussian_PCBplus100_Kaptonplus125_senTokap150_oldSensor",  "Gaussian_PCBplus100_Kaptonplus125_senTokap150_midSensor",  "Gaussian_PCBplus100_Kaptonplus125_senTokap150_newSensor",
"Gaussian_PCBplus100_Kaptonplus150_senTokap150_oldSensor",  "Gaussian_PCBplus100_Kaptonplus150_senTokap150_midSensor",  "Gaussian_PCBplus100_Kaptonplus150_senTokap150_newSensor",
"Gaussian_PCBplus100_Kaptonplus175_senTokap150_oldSensor",  "Gaussian_PCBplus100_Kaptonplus175_senTokap150_midSensor",  "Gaussian_PCBplus100_Kaptonplus175_senTokap150_newSensor",
"Gaussian_PCBplus100_Kaptonplus200_senTokap150_oldSensor",  "Gaussian_PCBplus100_Kaptonplus200_senTokap150_midSensor",  "Gaussian_PCBplus100_Kaptonplus200_senTokap150_newSensor",
"Gaussian_PCBplus100_Kaptonplus225_senTokap150_oldSensor",  "Gaussian_PCBplus100_Kaptonplus225_senTokap150_midSensor",  "Gaussian_PCBplus100_Kaptonplus225_senTokap150_newSensor",
            "Gaussian_Kaptonminus0_senTokap185_oldSensor",              "Gaussian_Kaptonminus0_senTokap185_midSensor",              "Gaussian_Kaptonminus0_senTokap185_newSensor",
"Gaussian_PCBplus100_Kaptonplus100_senTokap185_oldSensor",  "Gaussian_PCBplus100_Kaptonplus100_senTokap185_midSensor",  "Gaussian_PCBplus100_Kaptonplus100_senTokap185_newSensor",
"Gaussian_PCBplus100_Kaptonplus125_senTokap185_oldSensor",  "Gaussian_PCBplus100_Kaptonplus125_senTokap185_midSensor",  "Gaussian_PCBplus100_Kaptonplus125_senTokap185_newSensor",
"Gaussian_PCBplus100_Kaptonplus150_senTokap185_oldSensor",  "Gaussian_PCBplus100_Kaptonplus150_senTokap185_midSensor",  "Gaussian_PCBplus100_Kaptonplus150_senTokap185_newSensor",
"Gaussian_PCBplus100_Kaptonplus175_senTokap185_oldSensor",  "Gaussian_PCBplus100_Kaptonplus175_senTokap185_midSensor",  "Gaussian_PCBplus100_Kaptonplus175_senTokap185_newSensor",
"Gaussian_PCBplus100_Kaptonplus200_senTokap185_oldSensor",  "Gaussian_PCBplus100_Kaptonplus200_senTokap185_midSensor",  "Gaussian_PCBplus100_Kaptonplus200_senTokap185_newSensor",
"Gaussian_PCBplus100_Kaptonplus225_senTokap185_oldSensor",  "Gaussian_PCBplus100_Kaptonplus225_senTokap185_midSensor",  "Gaussian_PCBplus100_Kaptonplus225_senTokap185_newSensor",
]

os.makedirs(os.path.join(args.baseDir + "_small", "TestFull"))

with open("/eos/uscms/store/user/mkilpatr/13TeV/" + args.baseDir + "/Gaussian_Kaptonminus0_senTokap150_oldSensor_Full/Full_BadmoduleTolerances.json", "a+") as new:
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


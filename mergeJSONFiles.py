import os
import json
import argparse
import re

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
Dist = os.listdir(eosDir)

os.makedirs(os.path.join(args.baseDir + "_small", "TestFull"))

with open("/eos/uscms/store/user/mkilpatr/13TeV/" + args.baseDir + "/" + Dist[0] + "/Full_BadmoduleTolerances.json", "a+") as new:
    newData = json.load(new)
    newData['Worst']["Bad Components"] = {}
    newData['Worst']["Bad Overlaps"] = {}

for d in Dist:
    with open(eosDir+d+"/Full_BadmoduleTolerances.json", "r") as lepcr:
        lep_insert = json.load(lepcr)
        loc = d.replace('_Full', '')
        newData['Worst']["Bad Overlaps"][loc] = {}
        newData['Worst']["Bad Components"].update(lep_insert['Worst']["Bad Components"])
        newData['Worst']["Bad Overlaps"][loc].update(lep_insert['Worst']["Bad Overlaps"][loc])
write_json(newData)


import os
import io
import json
import argparse
import re
import pandas as pd

parser = argparse.ArgumentParser(
        description='Produce or print limits based on existing datacards')
parser.add_argument("-d", "--dir", dest="baseDir", default='.',
                         help="Base directory where files are located")
parser.add_argument("-t", "--type", dest="type", default="",
                         help="Determine JSON file format depending on test type")
parser.add_argument("-k", "--kind", dest="kind", default="Unconstrained warp",
                         help="Determine JSON file format depending on test type")
parser.add_argument("-l", "--hb_type", dest="hb_type", default="LD",
                         help="hexaboard types [LD, HD]")
args = parser.parse_args()

Dist = os.listdir(args.baseDir)

ModuleHeight = ["Corner", "Edges", "Middle", "Point", "Center"]
ModuleCycle = ["OGP", "0 Cycles", "1 Cycle", "2 Cycles", "3 Cycles", "Remount", "4 Cycles", "5 Cycles", "6 Cycles"] #Module 812
if args.type == "Module805": ModuleCycle = ["0 Cycles", "1 Cycle", "2 Cycles", "3 Cycles", "4 Cycles", "5 Cycles", "10 Cycles", "After Fix", "11 Cycles", "15 Cycles"] #Module 805
elif args.type == "Module801": ModuleCycle = ["0 Cycles", "1 Cycle", "2 Cycles", "3 Cycles", "4 Cycles", "5 Cycles", "10 Cycles", "After Fix", "11 Cycles", "15 Cycles"] #Module 801
elif args.type == "Deformation": ModuleCycle = ["OGP"] #Deformation studies
elif args.type == "Pandas": ModuleCycle = ["OGP"] #Deformation studies
PadCenter = [88]
PadCorner = [1, 18, 111, 197, 190, 81, 8, 95, 189, 191, 96, 9]
PadEdge = [3, 39, 140, 195, 169, 52, 5, 65, 168, 193, 141, 29]
PadMiddle = [22, 49, 124, 175, 158, 68, 24, 78, 152, 173, 129, 42]
PadPoint = [45, 76, 122, 148, 131, 71]
if args.hb_type == "HD":
    PadCenter = [228]
    PadPoint = [100, 164, 281, 339, 272, 155]

if args.type == "Pandas": 
    if args.hb_type == "LD":
        PadCorner = [198, 190, 180, 96, 81, 9, 1, 8, 18, 95, 111, 189]
    elif args.hb_type == "HD":
        PadCorner = [444, 432, 418, 217, 193, 13, 1, 12, 25, 216, 240, 431]

def readFile(FILENAME):
    file = open(args.baseDir + "/" + FILENAME, 'r')
    label = FILENAME.replace(".csv", "")
    j = {label: {"Cycle": [], "Height": [], "Pad": []}}
    h = c = i = 0
    cen = cor = edg = mid = pnt = 0
    while (True):
        next_line = file.readline()
        if not next_line: break
        ss = next_line.split(",")
    
        for s in ss:
            if (s == ""): continue
            height = float(s)
            j[label]["Cycle"].append(ModuleCycle[c])
            j[label]["Height"].append(height)
            if(i == 0):
                j[label]["Pad"].append(PadCorner[cor])
                cor+=1
            elif(i == 1 or i == 2):
                j[label]["Pad"].append(PadEdge[edg])
                edg+=1
            elif(i == 3):
                j[label]["Pad"].append(PadCorner[cor])
                cor+=1
            elif(i == 4 or i == 5):
                j[label]["Pad"].append(PadMiddle[mid])
                mid+=1
            elif(i == 6):
                j[label]["Pad"].append(PadPoint[pnt])
                pnt+=1
            elif(i == 7):
                j[label]["Pad"].append(PadCenter[cen])
                cen+=1
            
        i+=1
        if(i == 1): h+=1
        if(i == 3): h-=1
        if(i == 4): h+=2
        if(i == 6 or i == 7): h+=1
        if(i == 8):
            c+=1 
            h = i = 0
            cen = cor = edg = mid = pnt = 0
    return j

indexCorner = [3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20]
indexEdge = [4, 7, 10, 13, 16, 19]
indexPoint = [21, 22, 23, 24, 25, 26]
OGPPadEdge = [4, 51, 155, 194, 156, 40]
if args.type == "Pandas": 
    indexCorner = [0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17]
    indexEdge = [1, 4, 7, 10, 13, 16]
    indexPoint = [18, 19, 20, 21, 22, 23]
    if args.hb_type == "LD":
        OGPPadEdge = [194, 156, 40, 4, 51, 155]
    elif args.hb_type == "HD":
        OGPPadEdge = [438, 330, 91, 7, 108, 348]
  

def readFileOGP(FILENAME):
  file = open(args.baseDir + "/" + FILENAME, 'r')
  label = FILENAME.replace(".txt", "")
  j = {label: {"Cycle": [], "Height": [], "Pad": []}}
  h = 0
  i = 3
  cen = cor = edg = mid = pnt = 0
  height = 0.
  while (True):
      next_line = file.readline()
      if not next_line: break
      if("Point" not in next_line): continue

      s = re.search(r'[0-9].[0-9]+', next_line)
      height = float(s.group())

      if(bool(s)):
        j[label]["Cycle"].append(ModuleCycle[0])
        j[label]["Height"].append(height)
        if(i in indexCorner):
          j[label]["Pad"].append(PadCorner[cor])
          cor+=1
        elif(i in indexEdge):
          j[label]["Pad"].append(OGPPadEdge[edg])
          edg+=1
        elif(i in indexPoint):
          j[label]["Pad"].append(PadPoint[pnt])
          pnt+=1
        elif(i == 27):
          j[label]["Pad"].append(PadCenter[cen])
          cen+=1

      i+=1
      if(i in indexCorner): h = 0
      elif(i in indexEdge): h = 1
      elif(i in indexPoint): h = 3
      elif(i == 27): h = 4

  return j

def readPandas(filename, type = "Unconstrained warp"):
  df = pd.read_csv(args.baseDir + "/" + filename)
 
  label = filename.replace(".csv", "") + "_" + type.replace(" ", "_")
 
  j = {label: {"Cycle": [], "Height": [], "Pad": []}}
  h = cen = cor = edg = mid = pnt = 0
  height = 0.

  for i in range(0,df[type].size):
      height = df[type][i]
      point = int(df["#"][i]) - 1
      if pd.isna(height) or str(height) == " ": height = float(-999.)

      j[label]["Cycle"].append(ModuleCycle[0])
      j[label]["Height"].append(float(height))
      if(i in indexCorner):
        j[label]["Pad"].append(PadCorner[cor])
        cor+=1
      elif(i in indexEdge):
        j[label]["Pad"].append(OGPPadEdge[edg])
        edg+=1
      elif(i in indexPoint):
        j[label]["Pad"].append(PadPoint[pnt])
        pnt+=1
      elif(i == 24):
        j[label]["Pad"].append(PadCenter[cen])
        cen+=1
        break

      if(point in indexCorner): h = 0
      elif(point in indexEdge): h = 1
      elif(point in indexPoint): h = 3
      elif(point == 24): h = 4

  return j
  

for d in Dist:
    if os.path.isdir(args.baseDir + "/" + d): continue
    if "txt" in d or "csv" in d: 
        print(d)
        name = d.replace(".txt", "").replace(".csv", "")
        label = args.kind.replace(" ", "_")
        if args.type == "Pandas": name += "_" + label
        fout = open(args.baseDir + "/" + name + ".json", 'w')
        j = []
        if args.type == "Deformation": j = readFileOGP(d)
        elif args.type == "Pandas": j = readPandas(d, args.kind)
        else: j = readFile(d)

        if args.type == "Deformation" or args.type == "Pandas": fout.write(json.dumps(j, sort_keys=False, indent=1))
        else:                         fout.write(json.dumps(j, sort_keys=False))

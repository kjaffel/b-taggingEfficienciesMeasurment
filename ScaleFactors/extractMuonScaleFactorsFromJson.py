#! /bin/env python

import math
import argparse
import json
import re

parser = argparse.ArgumentParser()
parser.add_argument('file', help='Json file containing muon scale factors')
parser.add_argument('-s', '--suffix', help='Suffix to append at the end of the output filename', required=True)
parser.add_argument('-m', '--mode', choices=['id', 'iso'], help='Retrieve either ID or isolation SF')

args = parser.parse_args()

with open(args.file) as _f:
    data = json.load(_f)

if args.mode == 'id':
    num = "TightID"
    den = "TrackerMuons"
if args.mode == 'iso':
    num = 'TightRelIso'
    den = 'TightIDandIPCut'

entry = "NUM_{}_DEN_{}".format(num, den)
scaleFactors = data[entry]["abseta_pt"]

eta_binning = set()
pt_binning = set()

def getBoundaries(var, string):
    start = re.search(var + r':\[(.*),.*\]', string)
    stop = re.search(var + r':\[.*,(.*)\]', string)
    return (start.group(1), stop.group(1))

for key in scaleFactors:
    bounds = getBoundaries("abseta", key)
    eta_binning.add(bounds[0])
    eta_binning.add(bounds[1])
eta_binning = sorted(list(eta_binning), key=lambda x: float(x))

for key in next(iter(scaleFactors.values())):
    bounds = getBoundaries("pt", key)
    pt_binning.add(bounds[0])
    pt_binning.add(bounds[1])
pt_binning = sorted(list(pt_binning), key=lambda x: float(x))

json_content = {'dimension': 2, 'variables': ['AbsEta', 'Pt'], 'binning': {'x': [float(x) for x in eta_binning], 'y': [float(x) for x in pt_binning]}, 'data': [], 'error_type': 'absolute'}
json_content_data = json_content['data']

for i in range(0, len(eta_binning) - 1):
    eta_data = {'bin': [eta_binning[i], eta_binning[i + 1]], 'values': []}
    eta_entry = scaleFactors["abseta:[{},{}]".format(eta_binning[i], eta_binning[i+1])]
    
    for j in range(0, len(pt_binning) - 1):
        pt_entry = eta_entry["pt:[{},{}]".format(pt_binning[j], pt_binning[j+1])]
        value = pt_entry["value"]
        stat = pt_entry["stat"]
        syst = pt_entry["syst"]
        error = math.sqrt(stat**2 + syst**2)
        pt_data = {'bin': [pt_binning[j], pt_binning[j + 1]], 'value': value, 'error_low': error, 'error_high': error}
        eta_data['values'].append(pt_data)

    json_content_data.append(eta_data)

# Save JSON file
filename = 'Muon_%s.json' % args.suffix
with open(filename, 'w') as j:
    json.dump(json_content, j, indent=2)

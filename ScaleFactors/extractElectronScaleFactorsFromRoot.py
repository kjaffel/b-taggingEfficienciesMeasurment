#! /bin/env python

import math
import argparse
import ROOT

parser = argparse.ArgumentParser()
parser.add_argument('file', help='ROOT file containing electron scale factors')
parser.add_argument('-s', '--suffix', help='Suffix to append at the end of the output filename', required=True)

args = parser.parse_args()

# Moriond17: last pt bin is the same as the previous one, but with 100% error and should be ignored
IGNORE_LAST_PT_BIN = False

f = ROOT.TFile.Open(args.file)

for key in f.GetListOfKeys():
    wp = key.GetName()

    if not 'SF' in wp:
        continue

    if not '2D' in wp:
        continue

    h = f.Get(wp)

    # Get binning
    eta_binning = []
    for i in range(1, h.GetXaxis().GetNbins() + 1):
        if len(eta_binning) == 0:
            eta_binning.append(h.GetXaxis().GetBinLowEdge(i))
            eta_binning.append(h.GetXaxis().GetBinUpEdge(i))
        else:
            eta_binning.append(h.GetXaxis().GetBinUpEdge(i))

    pt_binning = []
    for i in range(1, h.GetYaxis().GetNbins() + 1):
        if len(pt_binning) == 0:
            pt_binning.append(h.GetYaxis().GetBinLowEdge(i))
            pt_binning.append(h.GetYaxis().GetBinUpEdge(i))
        else:
            pt_binning.append(h.GetYaxis().GetBinUpEdge(i))

    if IGNORE_LAST_PT_BIN and len(pt_binning) > 2:
        pt_binning.pop()

    eta = 'ClusEta' if eta_binning[0] < 0 else 'AbsClusEta'
    json_content = {'dimension': 2, 'variables': [eta, 'Pt'], 'binning': {'x': eta_binning, 'y': pt_binning}, 'data': [], 'error_type': 'absolute'}
    json_content_data = json_content['data']

    for i in range(0, len(eta_binning) - 1):
        eta_data = {'bin': [eta_binning[i], eta_binning[i + 1]], 'values': []}
        mean_eta = (eta_binning[i] + eta_binning[i + 1]) / 2.
        for j in range(0, len(pt_binning) - 1):
            mean_pt = (pt_binning[j] + pt_binning[j + 1]) / 2.
            bin = h.FindBin(mean_eta, mean_pt)
            error = h.GetBinError(bin)
            pt_data = {'bin': [pt_binning[j], pt_binning[j + 1]], 'value': h.GetBinContent(bin), 'error_low': error, 'error_high': error}

            eta_data['values'].append(pt_data)

        json_content_data.append(eta_data)

    # Save JSON file
    filename = 'Electron_%s_%s.json' % (wp, args.suffix)
    with open(filename, 'w') as j:
        import json
        json.dump(json_content, j, indent=2)

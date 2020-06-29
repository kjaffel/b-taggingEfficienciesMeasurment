#! /usr/bin/env python
import sys, os, json
import copy
import datetime
import subprocess
import numpy as np
import math
import glob
import re
import argparse
import ROOT
from ROOT import TCanvas, TPad, TLine, TH1F


#inputs: TTBar files with path
#open file as r+
#read it
#get all the rho_steps without "__", i.e. the nominal ones
#clone them


def produceUpDownHistograms(nominal, uncertainty):
    
    #Get number of bins
    n_bins = nominal.GetNcells()

    ups = []
    downs = []

    #i=0: underflow bin
    #i=n_bins+1: overflow bin
    for i in range(1, n_bins-1):
        up = nominal.Clone()
        up.SetDirectory(ROOT.nullptr)
        up.Reset()

        down = nominal.Clone()
        down.SetDirectory(ROOT.nullptr)
        down.Reset()
        
        for j in range(1, n_bins-1):
            if j==i:
                c = nominal.GetBinContent(i)
                cup = c*(1+float(uncertainty))
                cdown = c*(1-float(uncertainty))
                up.SetBinContent(i, cup)
                down.SetBinContent(i, cdown)
            else:
                up.SetBinContent(j, nominal.GetBinContent(j))
                down.SetBinContent(j, nominal.GetBinContent(j))
            
        up.SetName(nominal.GetName()+"__ttbarNormBin{}up".format(i))
        down.SetName(nominal.GetName()+"__ttbarNormBin{}down".format(i))

        up.SetDirectory(0)
        down.SetDirectory(0)

        ups.append(up)
        downs.append(down)

        del up
        del down

    return (ups, downs)


def CloneHistogram(nominal):

    #We want to add this uncertainty only for ttbar. If the file is some other MC,
    #create the up and down histograms with no variation.

    #Get number of bins
    n_bins = nominal.GetNcells()

    ups = []
    downs = []

    #i=0: underflow bin
    #i=n_bins+1: overflow bin
    for i in range(1, n_bins-1):
        up = nominal.Clone()
        up.SetDirectory(ROOT.nullptr)

        down = nominal.Clone()
        down.SetDirectory(ROOT.nullptr)
        
        for j in range(1, n_bins-1):
            up.SetBinContent(j, nominal.GetBinContent(j))
            down.SetBinContent(j, nominal.GetBinContent(j))
            
        up.SetName(nominal.GetName()+"__ttbarNormBin{}up".format(i))
        down.SetName(nominal.GetName()+"__ttbarNormBin{}down".format(i))

        up.SetDirectory(0)
        down.SetDirectory(0)

        ups.append(up)
        downs.append(down)

        del up
        del down

    return (ups, downs)


def UpdateFileFromFileRegex(file, regex, veto=None):
    """
    Return all histograms found in a file whose name matches a regexp.

    Arguments:

    file: Path to the considered file
    regex: Regexp to match the histogram names
    veto: Also a regexp. If specified, will not consider histograms matching the veto.

    Returns: Dictionary with (key, value) = (name, histogram)
    """

    _files = set()
    myRe = re.compile(regex)
    try:
        myReVeto = re.compile(veto)
    except TypeError:
        myReVeto = None

    foundHistos = {}

    nHistos=0

    r_file = ROOT.TFile.Open(file, "update")
    _files.add(r_file)
    if not r_file or not r_file.IsOpen():
        raise Exception("Could not open file {}".format(file))

    print (file)
    content = r_file.GetListOfKeys()
    print (type(content))

    for key in content:
        name = key.GetName()
        print ("name: ", name)

        if myRe.match(name) is not None:
            if myReVeto is not None:
                if myReVeto.match(name) is not None:
                    continue
            
            if "__" in name:
                continue
            
            item = key.ReadObj()
            
            if not item.InheritsFrom("TH1"):
                continue
            
            item.SetDirectory(0)
            foundHistos[name] = item

            [ups, downs] = produceUpDownHistograms(item, 0.053)

            for up in ups:
                up.Write()
            for down in downs:
                down.Write()

            del ups
            del downs

    r_file.Close()

    return r_file


if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Modify the ttbar histograms for add a bin by bin uncertainty of 5.3%')
    parser.add_argument('-d', '--directory', required=True, help='Name of the folder that contains the histograms')
    
    args = parser.parse_args()

    rootDir = args.directory
    slurmDir = os.path.join(rootDir, "2017Results/passMETcut_tofinalSel_removeSL_triggerSfs_primarySLsatasetsstillON/results/")   
    print ( slurmDir)
    for filename in os.listdir(slurmDir):
        #no systematics for data
        if not filename.startswith("SingleElectron") and not filename.startswith("SingleMuon") and not filename.startswith("MuonEG") and not filename.startswith("EGamma"): # last for 2018 
            #only for ttbar
            if filename.startswith("TTToFullLeptonic") or filename.startswith("TTToSemiLeptonic") or filename.startswith("TTToHadronic"):
                UpdateFileFromFileRegex(slurmDir+filename, "muel_Only2_bJet2_METcut_phi_DeepCSVT")

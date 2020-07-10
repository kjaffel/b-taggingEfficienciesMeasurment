import os
import sys

btagPath = os.path.dirname(__file__)
if btagPath not in sys.path:
        sys.path.append(btagPath)

import utils
from utils import *
import HistogramTools as HT
import logging
logger = logging.getLogger("Btag Plotter")
from bamboo.analysismodules import NanoAODHistoModule
from bamboo import treefunctions as op
from bamboo.plots import Plot, SummedPlot
from bamboo.plots import EquidistantBinning as EqBin
from bamboo.plots import VariableBinning as VarBin


def MakeBtagEfficienciesMaps(jets, bjets, categories, era):

    bFlavJets = op.select(jets, lambda j: j.hadronFlavour == 5)
    cFlavJets = op.select(jets, lambda j: j.hadronFlavour == 4)
    lFlavJets = op.select(jets, lambda j: j.hadronFlavour == 0)
    
    cutOnJetsLen = {
                    'Only2': op.rng_len(jets) ==2,
                    'atleast2': op.rng_len(jets) >1
                }

    btaggingWPs = {
            "DeepCSV":{ # era: (loose, medium, tight)
                        "2016": (0.2217, 0.6321, 0.8953), 
                        "2017":(0.1522, 0.4941, 0.8001), 
                        "2018":(0.1241, 0.4184, 0.7527) },
            "DeepFlavour":{
                        "2016":(0.0614, 0.3093, 0.7221), 
                        "2017":(0.0521, 0.3033, 0.7489), 
                        "2018": (0.0494, 0.2770, 0.7264) }
            }
        
    # look at specific slice in pt/|eta|
    jetSel = lambda j: op.AND(op.in_range(30., j.pt, 200.), op.in_range(-1., j.eta, 1.))
    selBJets = op.select(bFlavJets, jetSel)
    selLightJets = op.select(lFlavJets, jetSel)
        
    selBJetsMaxDR = op.select(selBJets, lambda jet: op.NOT(op.rng_any(jets, lambda j: op.AND(jet != j, op.deltaR(jet.p4, j.p4) < 0.4))))
    selBJetsMinDR = op.select(selBJets, lambda jet: op.rng_any(jets, lambda j: op.AND(jet != j, op.deltaR(jet.p4, j.p4) < 0.4)))
    selLightJetsMaxDR = op.select(selLightJets, lambda jet: op.NOT(op.rng_any(jets, lambda j: op.AND(jet != j, op.deltaR(jet.p4, j.p4) < 0.4))))
    selLightJetsMinDR = op.select(selLightJets, lambda jet: op.rng_any(jets, lambda j: op.AND(jet != j, op.deltaR(jet.p4, j.p4) < 0.4)))
    
    jetPairs = op.combine(jets, N=2)
    bJetPairs = op.combine(bFlavJets, N=2)
    lightJetPairs = op.combine(lFlavJets, N=2)
    minJetDR = op.rng_min(jetPairs, lambda pair: op.deltaR(pair[0].p4, pair[1].p4))
    minBJetDR = op.rng_min(bJetPairs, lambda pair: op.deltaR(pair[0].p4, pair[1].p4))
    minLightJetDR = op.rng_min(lightJetPairs, lambda pair: op.deltaR(pair[0].p4, pair[1].p4))
    
    for nbrj, jcut in cutOnJetsLen.items():
        for channel, (dilepton, catSel) in categories.items():
            if channel == "elmu":
                MuonElectronTwoJetsSel = catSel.refine(f"{nbrj}Jets_{channel}Sel", cut=[ jcut])
            elif channel == "muel":
                ElectronMuonTwoJetsSel = catSel.refine(f"{nbrj}Jets_{channel}Sel", cut=[ jcut])
    
        plots = []
        for flav, flavJets in zip(['b', 'c', 'light'], [bFlavJets, cFlavJets, lFlavJets]):
        # b tagging efficiencies as a function of flavour/pt/|eta|
            binning = (VarBin([30,50,70,100,140,200,300,600,1000]), EqBin(5, 0, 2.5))
    
            pt = op.map(flavJets, lambda j: j.pt)
            eta = op.map(flavJets, lambda j: op.abs(j.eta))
            plots.append(Plot.make2D(f"elmu_{nbrj}j_jet_pt_eta_{flav}", (pt, eta), ElectronMuonTwoJetsSel, binning))
            plots.append(Plot.make2D(f"muel_{nbrj}j_jet_pt_eta_{flav}", (pt, eta), MuonElectronTwoJetsSel, binning))
            plots.append(SummedPlot(f"pair_lept_OSOF_{nbrj}j_jet_pt_eta_{flav}", plots[-2:-1]))
            
            for tagger in btaggingWPs.keys():
                for (wp, discr) in zip(["L", "M", "T"], btaggingWPs[tagger][era]):
                    if tagger == "DeepFlavour":
                        selJets = op.select(flavJets, lambda j: j.btagDeepFlavB >= discr)
                        lambda_ = lambda j: j.btagDeepFlavB
                        suffix = "deepFlav"
                    elif tagger =="DeepCSV":
                        selJets = op.select(flavJets, lambda j: j.btagDeepB >= discr)
                        lambda_ = lambda j: j.btagDeepB
                        suffix = "deepB"

                    pt = op.map(selJets, lambda j: j.pt)
                    eta = op.map(selJets, lambda j: op.abs(j.eta))
                    plots.append(Plot.make2D(f"elmu_{nbrj}j_jet_pt_eta_{flav}_{tagger}{wp}", (pt, eta), ElectronMuonTwoJetsSel, binning))
                    plots.append(Plot.make2D(f"muel_{nbrj}j_jet_pt_eta_{flav}_{tagger}{wp}", (pt, eta), MuonElectronTwoJetsSel, binning))
                    plots.append(SummedPlot(f"pair_lept_OSOF_{nbrj}j_jet_pt_eta_{flav}_{tagger}{wp}", plots[-2:-1]))

        for vari, sel in zip(["elmu", "muel"], [ElectronMuonTwoJetsSel, MuonElectronTwoJetsSel]):
            plots.append(Plot.make1D(vari + f"_{nbrj}_minJetDR", minJetDR, sel, EqBin(40, 0.4, 2.), plotopts=utils.getOpts(vari)))
            plots.append(Plot.make1D(vari + f"_{nbrj}_minBJetDR", minBJetDR, sel, EqBin(40, 0.4, 2.), plotopts=utils.getOpts(vari)))
            plots.append(Plot.make1D(vari + f"_{nbrj}_minLightJetDR", minLightJetDR, sel, EqBin(40, 0.4, 2.), plotopts=utils.getOpts(vari)))

        #for tagger in btaggingWPs.keys():
        #    lambda_ = (lambda j: j.btagDeepFlavB if tagger== "DeepFlavour" else (lambda j: j.btagDeepB))
        #    for channel, sel in zip(['muel', 'elmu'], [MuonElectronTwoJetsSel,ElectronMuonTwoJetsSel]):
        #        plots.append(Plot.make1D(f"{channel}_2j_minDR_bJet_{tagger}", op.map(selBJetsMinDR, lambda_), ElectronMuonTwoJetsSel, EqBin(30, 0., 1.)))
        #        plots.append(Plot.make1D(f"{channel}_2j_minDR_lightJet_{tagger}", op.map(selLightJetsMinDR, lambda_), ElectronMuonTwoJetsSel, EqBin(30, 0., 1.)))
        #        plots.append(Plot.make1D(f"{channel}_2j_maxDR_bJet_{tagger}", op.map(selBJetsMaxDR, lambda_), ElectronMuonTwoJetsSel, EqBin(30, 0., 1.)))
        #        plots.append(Plot.make1D(f"{channel}_2j_maxDR_lightJet_{tagger}", op.map(selLightJetsMaxDR, lambda_), ElectronMuonTwoJetsSel, EqBin(30, 0., 1.)))

    return plots

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
    ## Get list of plots (taken from bamboo.HistogramsModule)
        if not self.plotList:
            tup, smpName, smpCfg = self.getATree()
            tree, noSel, backend, runAndLS = self.prepareTree(tup, sample=smpName, sampleCfg=smpCfg)
            self.plotList = self.definePlots(tree, noSel, sample=smpName, sampleCfg=smpCfg)

        # merge processes using their cross sections
        utils.normalizeAndMergeSamples(self.plotList, self.readCounters, config, resultsdir, os.path.join(resultsdir, "mergedProcesses.root"))
        def getRatio(f, n, d, suffix):
            num = f.Get(n)
            den = f.Get(d)
            ratio = num.Clone(num.GetName() + suffix)
            ratio.Divide(den)
            return ratio

        for proc in list(config["samples"].keys()) + ["mergedProcesses"]:
            tf = HT.openFileAndGet(os.path.join(resultsdir, proc +".root"), "update")

            # compute ratio histogram needed to correct the jet multiplicity
            #getRatio(tf, "1lep_nJets", "1lep_shape_nJets", "_sfCorr").Write()

            # compute efficiencies (divide histo after cut by total histo)
            for flav in ["b", "c", "light"]:
                for wp in ["L", "M", "T"]:
                    for nbrj in ["Only2", "atleast2"]:
                        getRatio(tf, f"pair_lept_OSOF_{nbrj}j_jet_pt_eta_{flav}_wp{wp}", f"pair_lept_OSOF_{nbrj}j_jet_pt_eta_{flav}", "_eff").Write()         
            
        tf.Close()

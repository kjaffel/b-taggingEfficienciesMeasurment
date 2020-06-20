import os
import sys

sys.path.append('/home/ucl/cp3/kjaffel/bamboodev/b-taggingEfficienciesMeasurment')
import utils
from ControlPLots import safeget
import ControlPLots as cp
import HistogramTools as HT

from bamboo.analysismodules import NanoAODHistoModule
from bamboo import treefunctions as op
from bamboo.plots import Plot, SummedPlot
from bamboo.plots import EquidistantBinning as EqBin
from bamboo.plots import VariableBinning as VarBin


def MakeBtagEfficienciesPlots(self, jets, bjets, categories):

    bFlavJets = op.select(jets, lambda j: j.hadronFlavour == 5)
    cFlavJets = op.select(jets, lambda j: j.hadronFlavour == 4)
    lFlavJets = op.select(jets, lambda j: j.hadronFlavour == 0)

    bJets_ = safeget(bjets, "DeepFlavour", "M")
    
    for channel, (dilepton, catSel) in categories.items():
        if channel=="elmu":
            MuonElectronTwoJetsSel = catSel.refine("twoJet{0}Sel_".format("ElMu"), cut=[ op.rng_len(jets) > 1 ])
        else:
            ElectronMuonTwoJetsSel = catSel.refine("twoJet{0}Sel_".format("MuEl"), cut=[ op.rng_len(jets) > 1 ])

    plots = []
    for flav, flavJets in zip(['b', 'c', 'light'], [bFlavJets, cFlavJets, lFlavJets]):
    # b tagging efficiencies as a function of flavour/pt/|eta|
        binning = (VarBin([30, 40, 60, 80, 100, 200, 350, 1000]), EqBin(5, 0, 2.5))

        pt = op.map(flavJets, lambda j: j.pt)
        eta = op.map(flavJets, lambda j: op.abs(j.eta))
        plots.append(Plot.make2D(f"1el1mu_2j_jet_pt_eta_{flav}", (pt, eta), ElectronMuonTwoJetsSel, binning))
        plots.append(Plot.make2D(f"1mu1el_2j_jet_pt_eta_{flav}", (pt, eta), MuonElectronTwoJetsSel, binning))
        plots.append(SummedPlot(f"pair_lept_OSSF_2j_jet_pt_eta_{flav}", plots[-2:-1]))

        # I am trying for 2018 data if it works --> i should propagate this to # era and to both tagger !
        for (wp, deepThr) in zip(["L", "M", "T"], [0.0494, 0.277, 0.7264]):
        
            selJets = op.select(flavJets, lambda j: j.btagDeepFlavB >= deepThr)
            pt = op.map(selJets, lambda j: j.pt)
            eta = op.map(selJets, lambda j: op.abs(j.eta))
            plots.append(Plot.make2D(f"1el1mu_2j_jet_pt_eta_{flav}_wp{wp}", (pt, eta), ElectronMuonTwoJetsSel, binning))
            plots.append(Plot.make2D(f"1mu1el_2j_jet_pt_eta_{flav}_wp{wp}", (pt, eta), MuonElectronTwoJetsSel, binning))
            plots.append(SummedPlot(f"pair_lept_OSSF_2j_jet_pt_eta_{flav}_wp{wp}", plots[-2:-1]))

    jetPairs = op.combine(jets, N=2)
    minJetDR = op.rng_min(jetPairs, lambda pair: op.deltaR(pair[0].p4, pair[1].p4))
    bJetPairs = op.combine(bFlavJets, N=2)
    minBJetDR = op.rng_min(bJetPairs, lambda pair: op.deltaR(pair[0].p4, pair[1].p4))
    lightJetPairs = op.combine(lFlavJets, N=2)
    minLightJetDR = op.rng_min(lightJetPairs, lambda pair: op.deltaR(pair[0].p4, pair[1].p4))

    for vari, sel in zip(["1el1mu", "1mu1el"], [ElectronMuonTwoJetsSel, MuonElectronTwoJetsSel]):
        plots.append(Plot.make1D(vari +"_minJetDR", minJetDR, sel, EqBin(40, 0.4, 2.), plotopts=utils.getOpts(vari)))
        plots.append(Plot.make1D(vari + "_minBJetDR", minBJetDR, sel, EqBin(40, 0.4, 2.), plotopts=utils.getOpts(vari)))
        plots.append(Plot.make1D(vari + "_minLightJetDR", minLightJetDR, sel, EqBin(40, 0.4, 2.), plotopts=utils.getOpts(vari)))

        plots += cp.makeJetPlots(sel, jets, vari, maxJet=2, binScaling=2)
        plots += cp.makeBJetPlots(sel, jets, vari)
        plots.append(Plot.make1D(vari + f"_nDeepFlavM", op.rng_len(bJets_), sel, EqBin(6, 0, 6), title="Number of deepFlavourM b jets", plotopts=utils.getOpts(vari)))


    # look at specific slice in pt/|eta|
    jetSel = lambda j: op.AND(op.in_range(80., j.pt, 200.), op.in_range(-1., j.eta, 1.))
    selBJets = op.select(bFlavJets, jetSel)
    selLightJets = op.select(lFlavJets, jetSel)

    # deepFlavour shape for specific jet multiplicities
    for i in range(0, 2):
        nJetSel = catSel.refine(f"muon_eq{i}jets", cut=op.rng_len(jets)==i)
        plots.append(Plot.make1D(f"1mu_eq{i}j_bJet_deepFlav", op.map(selBJets, lambda j: j.btagDeepFlavB), nJetSel, EqBin(30, 0., 1.)))
        plots.append(Plot.make1D(f"1mu_eq{i}j_lightJet_deepFlav", op.map(selLightJets, lambda j: j.btagDeepFlavB), nJetSel, EqBin(30, 0., 1.)))

    # b and light jets which have - or not - another jet within DeltaR < 0.6
    selBJetsMaxDR = op.select(selBJets, lambda jet: op.NOT(op.rng_any(jets, lambda j: op.AND(jet != j, op.deltaR(jet.p4, j.p4) < 0.6))))
    selBJetsMinDR = op.select(selBJets, lambda jet: op.rng_any(jets, lambda j: op.AND(jet != j, op.deltaR(jet.p4, j.p4) < 0.6)))
    selLightJetsMaxDR = op.select(selLightJets, lambda jet: op.NOT(op.rng_any(jets, lambda j: op.AND(jet != j, op.deltaR(jet.p4, j.p4) < 0.6))))
    selLightJetsMinDR = op.select(selLightJets, lambda jet: op.rng_any(jets, lambda j: op.AND(jet != j, op.deltaR(jet.p4, j.p4) < 0.6)))
    # deepFlavour shape for the latter
    plots.append(Plot.make1D(f"1el1mu_2j_minDR_bJet_deepFlav", op.map(selBJetsMinDR, lambda j: j.btagDeepFlavB), ElectronMuonTwoJetsSel, EqBin(30, 0., 1.)))
    plots.append(Plot.make1D(f"1el1mu_2j_minDR_lightJet_deepFlav", op.map(selLightJetsMinDR, lambda j: j.btagDeepFlavB), ElectronMuonTwoJetsSel, EqBin(30, 0., 1.)))
    plots.append(Plot.make1D(f"1el1mu_2j_maxDR_bJet_deepFlav", op.map(selBJetsMaxDR, lambda j: j.btagDeepFlavB), ElectronMuonTwoJetsSel, EqBin(30, 0., 1.)))
    plots.append(Plot.make1D(f"1el1mu_2j_maxDR_lightJet_deepFlav", op.map(selLightJetsMaxDR, lambda j: j.btagDeepFlavB), ElectronMuonTwoJetsSel, EqBin(30, 0., 1.)))

    return plots

    #def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
    ## Get list of plots (taken from bamboo.HistogramsModule)
    #    if not self.plotList:
    #        tup, smpName, smpCfg = self.getATree()
    #        tree, noSel, backend, runAndLS = self.prepareTree(tup, sample=smpName, sampleCfg=smpCfg)
    #        self.plotList = self.definePlots(tree, noSel, sample=smpName, sampleCfg=smpCfg)

    #    # merge processes using their cross sections
    #    utils.normalizeAndMergeSamples(self.plotList, self.readCounters, config, resultsdir, os.path.join(resultsdir, "mergedProcesses.root"))
    #    def getRatio(f, n, d, suffix):
    #        num = f.Get(n)
    #        den = f.Get(d)
    #        ratio = num.Clone(num.GetName() + suffix)
    #        ratio.Divide(den)
    #        return ratio

    #    for proc in list(config["samples"].keys()) + ["mergedProcesses"]:
    #        tf = HT.openFileAndGet(os.path.join(resultsdir, proc +".root"), "update")

    #        # compute ratio histogram needed to correct the jet multiplicity
    #        getRatio(tf, "1lep_nJets", "1lep_shape_nJets", "_sfCorr").Write()

    #        # compute efficiencies (divide histo after cut by total histo)
    #        for flav in ["b", "c", "light"]:
    #            for wp in ["L", "M", "T"]:
    #                getRatio(tf, f"pair_lept_OSSF_2j_jet_pt_eta_{flav}_wp{wp}", f"pair_lept_OSSF_2j_jet_pt_eta_{flav}", "_eff").Write()         
    #        
    #        tf.Close()

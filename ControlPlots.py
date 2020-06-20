import os
import sys

btagPath = os.path.dirname(__file__)
if btagPath not in sys.path:
    sys.path.append(btagPath)

from bambooToOls import Plot
from bamboo.plots import SummedPlot
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

import utils
from utils import *

def makeJetPlots(sel, jets, suffix, uname):
    maxJet=2
    binScaling=1
    plots = []
    for i in range(maxJet):
        plots.append(Plot.make1D(f"{uname}_{suffix}_jet{i+1}_pt", jets[i].pt, sel,
                    EqBin(60 // binScaling, 30., 730. - max(2, i) * 100), title=f"{utils.getCounter(i+1)} jet p_{{T}} (GeV)",
                    plotopts=utils.getOpts(uname, **{"log-y": True})))
        plots.append(Plot.make1D(f"{uname}_{suffix}_jet{i+1}_eta", jets[i].eta, sel,
                    EqBin(50 // binScaling, -2.4, 2.4), title=f"{utils.getCounter(i+1)} jet eta",
                    plotopts=utils.getOpts(uname, **{"log-y": True})))
        plots.append(Plot.make1D(f"{uname}_{suffix}_jet{i+1}_phi", jets[i].phi, sel,
                    EqBin(50 // binScaling, -3.1416, 3.1416), title=f"{utils.getCounter(i+1)} jet #phi", 
                    plotopts=utils.getOpts(uname, **{"log-y": True})))

    return plots
    
def makeBJetPlots(sel, bjets, wp, suffix, uname):
    
    plots = []
    binScaling=1
    for key in sel.keys():
        tagger=key.replace(wp, "")
        bjets_ = safeget(bjets, tagger, wp)
        for i in range(2):
            plots.append(Plot.make1D(f"{uname}_{suffix}_bJet{i+1}_pT_{key}".format(key=key), 
                        bjets_[i].pt,
                        safeget(sel, key), 
                        EqBin(60 // binScaling, 30., 730. - max(2, i) * 100), 
                        title=f"{utils.getCounter(i+1)} bJet pT [GeV]", 
                        plotopts=utils.getOpts(uname, **{"log-y": True})))
                
            plots.append(Plot.make1D(f"{uname}_{suffix}_bJet{i+1}_eta_{key}".format(key=key), 
                        bjets_[i].eta,
                        safeget(sel, key), 
                        EqBin(50 // binScaling, -2.4, 2.4), 
                        title=f"{utils.getCounter(i+1)} bJet eta", 
                        plotopts=utils.getOpts(uname, **{"log-y": True})))
            
            plots.append(Plot.make1D(f"{uname}_{suffix}_bJet{i+1}_phi_{key}".format(key=key), 
                        bjets_[i].phi, 
                        safeget(sel, key),
                        EqBin(50 // binScaling, -3.1416, 3.1416), 
                        title=f"{utils.getCounter(i+1)} bJet phi", 
                        plotopts=utils.getOpts(uname, **{"log-y": True})))
    return plots

def makeDiscriminatorPlots(sel, bjets, wp, btaggingWPs, suffix, uname, era):
    plots = []
    for key in sel.keys():
            
        tagger=key.replace(wp, "")
        bjets_ = safeget(bjets, tagger, wp)
        idx = ( 0 if wp=="L" else ( 1 if wp=="M" else 2))
        bin0 = btaggingWPs[tagger][era][idx]
        for i in range(2):
            if tagger =="DeepFlavour":
                plots.append(Plot.make1D(f"{uname}_{suffix}_jet{i+1}_discr_deepFlav{wp}", 
                                bjets_[i].btagDeepFlavB,
                                safeget(sel, "DeepFlavour{0}".format(wp)),
                                EqBin(60, bin0, 1.), 
                                title="DeepFlavourBDisc {0}".format(wp), 
                                plotopts=utils.getOpts(uname, **{"log-y": True})))
            else:
                plots.append(Plot.make1D(f"{uname}_{suffix}_jet{i+1}_discr_deepCSV{wp}", 
                                bjets_[i].btagDeepB,
                                safeget(sel, "DeepCSV{0}".format(wp)),
                                EqBin(60, bin0, 1.), 
                                title= "DeepCSVBDisc {0}".format(wp), 
                                plotopts=utils.getOpts(uname, **{"log-y": True})))
        # Let's do L, M, T discr in one plot
    return plots

def makeLeptonPlots(sel, leptons, suffix, uname):
    binScaling=1
    plots = []
    if "mu" in uname:
        flav = "Muon"
    if "el" in uname:
        flav = "Electron"
    for i in range(2):
        leptonptCut = (25. if i == 0 else( 10. if uname[-1]=='u' else( 15. if uname[-1]=='l' else(0.))))
        plots.append(Plot.make1D(f"{uname}_%s_lep{i+1}_pt"%suffix, leptons[i].pt, sel,
                        EqBin(60 // binScaling, leptonptCut, 530.), title="%s p_{T} [GeV]" % flav,
                        plotopts=utils.getOpts(uname, **{"log-y": True})))
        plots.append(Plot.make1D(f"{uname}_%s_lep{i+1}_eta"%suffix, leptons[i].eta, sel,
                        EqBin(50 // binScaling, -2.4, 2.4), title="%s eta" % flav,
                        plotopts=utils.getOpts(uname, **{"log-y": True})))
        plots.append(Plot.make1D(f"{uname}_%s_lep{i+1}_phi"%suffix, leptons[i].phi, sel,
                        EqBin(50 // binScaling, -3.1416, 3.1416), title="%s #phi" % flav,
                        plotopts=utils.getOpts(uname, **{"log-y": True})))
    return plots
        
def makeJetmultiplictyPlots(sel, jets, suffix, uname):
    binScaling=1
    plots=[]
    plots.append(Plot.make1D(f"{uname}_%s_Jet_mulmtiplicity"%suffix, op.rng_len(jets), sel,
                    EqBin(10, 0., 10.), title="Jet mulmtiplicity",
                    plotopts=utils.getOpts(uname, **{"log-y": True})))
    return plots
def makePrimaryANDSecondaryVerticesPlots(sel, uname):
    binScaling=1
    plots=[]
    sv_mass=op.map(t.SV, lambda sv: sv.mass)
    sv_eta=op.map(t.SV, lambda sv: sv.eta)
    sv_phi=op.map(t.SV, lambda sv: sv.phi)
    sv_pt=op.map(t.SV, lambda sv: sv.pt)

    plots.append(Plot.make1D(f"{uname}_number_primary_reconstructed_vertices", 
                    t.PV.npvsGood, sel, 
                    EqBin(50 // binScaling, 0., 60.), title="reconstructed vertices",
                    plotopts=utils.getOpts(uname, **{"log-y": True})))
    plots.append(Plot.make1D(f"{uname}_secondary_vertices_mass", 
                    sv_mass, sel, 
                    EqBin(50 // binScaling, 0., 450.), title="SV mass",
                    plotopts=utils.getOpts(uname, **{"log-y": True})))
    plots.append(Plot.make1D(f"{uname}_secondary_vertices_eta", 
                    sv_eta, sel, 
                    EqBin(50 // binScaling,-2.4, 2.4), title="SV eta",
                    plotopts=utils.getOpts(uname, **{"log-y": True})))
    plots.append(Plot.make1D(f"{uname}_secondary_vertices_phi", 
                    sv_phi, sel, 
                    EqBin(50 // binScaling, -3.1416, 3.1416), title="SV #phi",
                    plotopts=utils.getOpts(uname, **{"log-y": True})))
    plots.append(Plot.make1D(f"{uname}_secondary_vertices_pt", 
                    sv_pt, sel, 
                    EqBin(50 // binScaling, 0., 450.), title="SV p_{T} [GeV]",
                    plotopts=utils.getOpts(uname, **{"log-y": True})))

    return plots

def makeMETPlots(sel, leptons, met, corrMET, uname):
    binScaling=1
    plots = []

    plots.append(Plot.make1D(f"{uname}_MET_pt", met.pt, sel,
                EqBin(60 // binScaling, 0., 600.), title="MET p_{T} [GeV]",
                plotopts=utils.getOpts(uname, **{"log-y": True})))
    plots.append(Plot.make1D(f"{uname}_MET_phi", met.phi, sel,
                EqBin(60 // binScaling, -3.1416, 3.1416), title="MET #phi",
                plotopts=utils.getOpts(uname, **{"log-y": True})))
    for i in range(2):
        plots.append(Plot.make1D(f"{uname}_MET_lep{i+1}_deltaPhi",
                op.Phi_mpi_pi(leptons[i].phi - met.phi), sel, EqBin(60 // binScaling, -3.1416, 3.1416),
                title="#Delta #phi (lepton, MET)", plotopts=utils.getOpts(uname, **{"log-y": True})))
    
        MT = op.sqrt( 2. * met.pt * leptons[i].p4.Pt() * (1. - op.cos(op.Phi_mpi_pi(met.phi - leptons[i].p4.Phi()))) )
        plots.append(Plot.make1D(f"{uname}_MET_MT_lep{i+1}", MT, sel,
                EqBin(60 // binScaling, 0., 600.), title="Lepton M_{T} [GeV]",
                plotopts=utils.getOpts(uname, **{"log-y": True})))

    return plots


def MakeEXTRAMETPlots(sel, corrmet, met, suffix, uname):
    plots = []
    binScaling=1

    plots.append(Plot.make1D(f"{uname}_{suffix}_met_pT",
                met.pt, sel, 
                EqBin(60 // binScaling, 0., 600.), title="MET p_{T} [GeV]",
                plotopts=utils.getOpts(uname, **{"log-y": True})))
                        
    plots.append(Plot.make1D(f"{uname}_{suffix}_xycorrmet_pT",
                corrmet.pt, sel,
                EqBin(60 // binScaling, 0., 600.), title="corrMET p_{T} [GeV]",
                plotopts=utils.getOpts(uname, **{"log-y": True})))

    return plots


def makedeltaRPlots( self, sel, jets, leptons, suffix, uname):
    plots = []
    plots.append(Plot.make1D(f"{uname}_{suufix}_jet1jet2_deltaR", op.deltaR(jets[0].p4, jets[1].p4),
                sel, EqBin(50, 0., 8.), title="#Delta R (leading jet, sub-leading jet)",
                plotopts=utils.getOpts(uname, **{"log-y": False})))
    plots.append(Plot.make1D(f"{uname}_{suffix}_jet1lep1_deltaR", op.deltaR(jets[0].p4, leptons[0].p4),
                sel, EqBin(50, 0., 8.), title="#Delta R (leading jet, leading lepton)",
                plotopts=utils.getOpts(uname, **{"log-y": False})))
    plots.append(Plot.make1D(f"{uname}_{suffix}_jet1lep2_deltaR", op.deltaR(jets[0].p4, leptons[1].p4),
                sel, EqBin(50, 0., 8.), title="#Delta R (leading jet, sub-leading lepton)",
                plotopts=utils.getOpts(uname, **{"log-y": False})))
    plots.append(Plot.make1D(f"{uname}_{suffix}_jet2lep1_deltaR", op.deltaR(jets[1].p4, leptons[0].p4),
                sel, EqBin(50, 0., 8.), title="#Delta R (sub-leading jet, leading lepton)",
                plotopts=utils.getOpts(uname, **{"log-y": False})))
    plots.append(Plot.make1D(f"{uname}_{suffix}_jet2lep2_deltaR", op.deltaR(jets[1].p4, leptons[1].p4),
                sel, EqBin(50, 0., 8.), title="#Delta R (leading jet, sub-leading lepton)",
                plotopts=utils.getOpts(uname, **{"log-y": False})))
    plots.append(Plot.make1D(f"{uname}_{suffix}_lep1lep2_deltaR", op.deltaR(leptons[0].p4, leptons[1].p4),
                sel, EqBin(50, 0., 8.), title="#Delta R (leading lepton, sub-leading lepton)",
                plotopts=utils.getOpts(uname, **{"log-y": False})))
    return plots

def makehistosforTTbarEstimation(selections, ll, bjets, wp, met, suffix, uname):
    plots=[]
    for key, sel in selections.items():
        tagger = key.replace(wp, "")
        bjets_ = safeget(bjets, tagger, wp)
        bb = bjets[key.replace(wp, "")][wp]
        plots += [ Plot.make1D(f"{nm}_{uname}_{suffix}_lljj_{tagger}_btag{wp}_mll_and_{met}",
            llbbVar, sel, binning, title=title, plotopts=utils.getOpts(uname))
            for nm, (llbbVar, binning, title) in {
                "jj_M"   : (op.invariant_mass(bb[0].p4+bb[1].p4), EqBin(40, 10., 1000.), "Mjj [GeV]"),
                "lljj_M" : ((ll[0].p4+ll[1].p4+bb[0].p4+bb[1].p4).M(), EqBin(50, 100., 1500.), "Mlljj [GeV]"),
                "ll_M"   : (op.invariant_mass(ll[0].p4+ll[1].p4), EqBin(60, 70., 120.), "Mll [GeV]"),
                "jj_DR"  : (op.deltaR(bb[0].p4, bb[1].p4), EqBin(50, 0., 6.), "jj deltaR"),
                "jj_pt"  : ((bjets_[0].p4 + bjets_[1].p4).Pt(), EqBin(50, 0., 450.), "dijets p_{T} [GeV]"),
                "jet1_pt": (bb[0].pt, EqBin( 50, 20., 500.), "jet1 p_{T} [GeV]"),
                "jet2_pt": (bb[1].pt, EqBin( 50, 20., 300.), "jet2 p_{T} [GeV]"),
                "lep1_pt": (ll[0].pt, EqBin( 50, 20., 400.), "Leading Lepton p_{T} [GeV]"),
                "lep2_pt": (ll[1].pt, EqBin( 50, 10., 200.), "Sub-leading Lepton p_{T} [GeV]"),
                }.items()
            ]
    return plots 

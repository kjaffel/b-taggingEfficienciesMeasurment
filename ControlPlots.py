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
                    plotopts=utils.getOpts(uname, **{"log-y": False})))
        plots.append(Plot.make1D(f"{uname}_{suffix}_jet{i+1}_eta", jets[i].eta, sel,
                    EqBin(50 // binScaling, -2.4, 2.4), title=f"{utils.getCounter(i+1)} jet eta",
                    plotopts=utils.getOpts(uname, **{"log-y": False})))
        plots.append(Plot.make1D(f"{uname}_{suffix}_jet{i+1}_phi", jets[i].phi, sel,
                    EqBin(50 // binScaling, -3.1416, 3.1416), title=f"{utils.getCounter(i+1)} jet #phi", 
                    plotopts=utils.getOpts(uname, **{"log-y": False})))

    jj_p4 = jets[0].p4+jets[1].p4

    plots += [ Plot.make1D("{0}_{1}_jj{nm}".format(uname, suffix, nm=nm), var, sel, binning,
                title=f"di-jets {title}", plotopts=utils.getOpts(uname))
                for nm, (var, binning, title) in {
                    "PT" : (jj_p4.Pt() , EqBin(60 // binScaling, 0., 450.), "P_{T} [GeV]"),
                    "Phi": (jj_p4.Phi(), EqBin(50 // binScaling, -3.1416, 3.1416), "#phi"),
                    "Eta": (jj_p4.Eta(), EqBin(50 // binScaling, -3., 3.), "Eta")
                }.items()
            ]

    return plots
    
def makeBJetPlots(sel, bb, key, channel, cut, suffix):
    plots = []
    binScaling=1
    for i in range(2):
        plots.append(Plot.make1D(f"{channel}_{suffix}_bJet{i+1}_{cut}_pT_{key}",
                    bb[i].pt, sel, EqBin(60 // binScaling, 30., 730. - max(2, i) * 100),
                    title=f"{utils.getCounter(i+1)} bJet pT [GeV]",
                    plotopts=utils.getOpts(channel, **{"log-y": False})))
            
        plots.append(Plot.make1D(f"{channel}_{suffix}_bJet{i+1}_{cut}_eta_{key}",
                    bb[i].eta, sel, EqBin(50 // binScaling, -2.4, 2.4),
                    title=f"{utils.getCounter(i+1)} bJet eta",
                    plotopts=utils.getOpts(channel, **{"log-y": False})))
        
        plots.append(Plot.make1D(f"{channel}_{suffix}_bJet{i+1}_{cut}_phi_{key}",
                    bb[i].phi, sel, EqBin(50 // binScaling, -3.1416, 3.1416),
                    title=f"{utils.getCounter(i+1)} bJet phi",
                    plotopts=utils.getOpts(channel, **{"log-y": False})))

    jj_p4 = bb[0].p4+bb[1].p4

    plots += [ Plot.make1D(f"{channel}_{suffix}{cut}_bb{nm}_{key}",
                var, sel, binning, title=f"di-bjet {title}", plotopts=utils.getOpts(channel))
                for nm, (var, binning, title) in {
                    "PT" : (jj_p4.Pt() , EqBin(60 // binScaling, 0., 450.), "P_{T} [GeV]"),
                    "Phi": (jj_p4.Phi(), EqBin(50 // binScaling, -3.1416, 3.1416), "#phi"),
                    "Eta": (jj_p4.Eta(), EqBin(50 // binScaling, -3., 3.), "Eta")
                }.items()
            ]

    return plots

def makeLeptonPlots(sel, leptons, suffix, uname):
    binScaling=1
    plots = []

    for i in range(2):
        leptonptCut = (25. if i == 0 else( 10. if uname[-1]=='u' else( 15. if uname[-1]=='l' else(0.))))
        plots.append(Plot.make1D(f"{uname}_%s_lep{i+1}_pt"%suffix, leptons[i].pt, sel,
                        EqBin(60 // binScaling, leptonptCut, 530.), title=f"{utils.getCounter(i+1)} Lepton pT [GeV]",
                        plotopts=utils.getOpts(uname, **{"log-y": False})))
        plots.append(Plot.make1D(f"{uname}_%s_lep{i+1}_eta"%suffix, leptons[i].eta, sel,
                        EqBin(50 // binScaling, -2.4, 2.4), title=f"{utils.getCounter(i+1)} Lepton eta",
                        plotopts=utils.getOpts(uname, **{"log-y": False})))
        plots.append(Plot.make1D(f"{uname}_%s_lep{i+1}_phi"%suffix, leptons[i].phi, sel,
                        EqBin(50 // binScaling, -3.1416, 3.1416), title=f"{utils.getCounter(i+1)} Lepton #phi",
                        plotopts=utils.getOpts(uname, **{"log-y": False})))
    
    plots.append(Plot.make1D(f"{uname}_llpT_%s"%(suffix), (leptons[0].p4 + leptons[1].p4).Pt(), sel,
                    EqBin(60 // binScaling, 0., 450.), title= "dilepton P_{T} [GeV]",
                    plotopts=utils.getOpts(uname)))
    plots.append(Plot.make1D(f"{uname}_llphi_%s"%(suffix), (leptons[0].p4 + leptons[1].p4).Phi(), sel,
                    EqBin(50 // binScaling, -3.1416, 3.1416), title= "dilepton #phi ",
                    plotopts=utils.getOpts(uname)))
    plots.append(Plot.make1D(f"{uname}_lleta_%s"%(suffix), (leptons[0].p4 + leptons[1].p4).Eta(), sel,
                    EqBin(50 // binScaling, -2.5, 2.5), title= "dilepton eta",
                    plotopts=utils.getOpts(uname)))

    return plots

def makePrimaryANDSecondaryVerticesPlots(sel, t, uname):
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
                plotopts=utils.getOpts(uname, **{"log-y": False})))
    plots.append(Plot.make1D(f"{uname}_MET_phi", met.phi, sel,
                EqBin(60 // binScaling, -3.1416, 3.1416), title="MET #phi",
                plotopts=utils.getOpts(uname, **{"log-y": False})))
    for i in range(2):
        plots.append(Plot.make1D(f"{uname}_MET_lep{i+1}_deltaPhi",
                op.Phi_mpi_pi(leptons[i].phi - met.phi), sel, EqBin(60 // binScaling, -3.1416, 3.1416),
                title="#Delta #phi (lepton, MET)", plotopts=utils.getOpts(uname, **{"log-y": True})))
    
        MT = op.sqrt( 2. * met.pt * leptons[i].p4.Pt() * (1. - op.cos(op.Phi_mpi_pi(met.phi - leptons[i].p4.Phi()))) )
        plots.append(Plot.make1D(f"{uname}_MET_MT_lep{i+1}", MT, sel,
                EqBin(60 // binScaling, 0., 600.), title="Lepton M_{T} [GeV]",
                plotopts=utils.getOpts(uname, **{"log-y": False})))

    return plots


def MakeEXTRAMETPlots(sel, corrmet, met, suffix, uname):
    plots = []
    binScaling=1

    plots.append(Plot.make1D(f"{uname}_{suffix}_met_pT",
                met.pt, sel, 
                EqBin(60 // binScaling, 0., 600.), title="MET p_{T} [GeV]",
                plotopts=utils.getOpts(uname, **{"log-y": False})))
                        
    plots.append(Plot.make1D(f"{uname}_{suffix}_xycorrmet_pT",
                corrmet.pt, sel,
                EqBin(60 // binScaling, 0., 600.), title="corrMET p_{T} [GeV]",
                plotopts=utils.getOpts(uname, **{"log-y": False})))

    return plots


def makedeltaRPlots(sel, jets, leptons, suffix, uname):
    plots = []
    plots.append(Plot.make1D(f"{uname}_{suffix}_jet1jet2_deltaR", op.deltaR(jets[0].p4, jets[1].p4),
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

def makehistosforTTbarEstimation(selection, ll, bb, uname, channel):
    return [ Plot.make1D(f"{nm}_{channel}_{uname}",
        llbbVar, selection, binning, title=title, plotopts=utils.getOpts(channel))
        for nm, (llbbVar, binning, title) in {
            "jj_M"   : (op.invariant_mass(bb[0].p4+bb[1].p4), EqBin(40, 10., 1000.), "Mjj [GeV]"),
            "lljj_M" : ((ll[0].p4+ll[1].p4+bb[0].p4+bb[1].p4).M(), EqBin(50, 100., 1500.), "Mlljj [GeV]"),
            "ll_M"   : (op.invariant_mass(ll[0].p4+ll[1].p4), EqBin(60, 70., 120.), "Mll [GeV]"),
            "jj_DR"  : (op.deltaR(bb[0].p4, bb[1].p4), EqBin(50, 0., 6.), "jj deltaR"),
            "jj_pt"  : ((bb[0].p4 + b[1].p4).Pt(), EqBin(50, 0., 450.), "dijets p_{T} [GeV]"),
            "jet1_pt": (bb[0].pt, EqBin( 50, 20., 500.), "jet1 p_{T} [GeV]"),
            "jet2_pt": (bb[1].pt, EqBin( 50, 20., 300.), "jet2 p_{T} [GeV]"),
            "lep1_pt": (ll[0].pt, EqBin( 50, 20., 400.), "Leading Lepton p_{T} [GeV]"),
            "lep2_pt": (ll[1].pt, EqBin( 50, 10., 200.), "Sub-leading Lepton p_{T} [GeV]"),
            }.items()
        ]

def make2DMAPS(sel, jets, suffix, uname):
    plots=[]
    plots.append( Plot.make2D(f"{uname}_{suffix}_jetETA_vs_Phi",(op.map(jets, lambda j : j.eta), op.map(jets, lambda j : j.phi)), sel, 
        (EqBin(50, -2.5, 2.5), EqBin(50,  -2.5, 2.5)), title=" Eta vs #phi", 
        plotopts=utils.getOpts(uname, **{"log-y": True})))
    #plots.append( Plot.make2D(f"{uname}_{suffix}_jetPT_vs_ETA",(op.map(jets, lambda j : j.pt), op.map(jets, lambda j : j.eta)), sel, 
    #    (EqBin(50, 30., 600.), EqBin(50,  -2.5, 2.5)), title=" p_{T} vs Eta", 
    #    plotopts=utils.getOpts(uname, **{"log-y": True})))
    lowerPtBinEdges=[30,50,70,100,140,200,300, 600]
    for i in range(len(lowerPtBinEdges)-1):
        ptBin= str(lowerPtBinEdges[i])+"to"+str(lowerPtBinEdges[i+1])
        plots.append( Plot.make2D(f"{uname}_{suffix}_jetPT_vs_ETA_ptBin{ptBin}",(op.map(jets, lambda j : j.pt), op.map(jets, lambda j : j.eta)), sel, 
            (EqBin(50, lowerPtBinEdges[i], lowerPtBinEdges[i+1]), EqBin(50, -2.5, 2.5)), title=" p_{T} vs Eta", 
            plotopts=utils.getOpts(uname, **{"log-y": True})))

    return plots

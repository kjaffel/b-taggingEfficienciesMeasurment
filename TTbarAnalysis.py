from bamboo.analysismodules import NanoAODHistoModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.scalefactors import binningVariables_nano

from bamboo import treefunctions as op
from bamboo.plots import EquidistantBinning as EqB
from bamboo import scalefactors
import logging
logger = logging.getLogger("ttbar Plotter")

from itertools import chain
from functools import partial
import os.path
import collections
from collections import defaultdict
import math
import argparse
import sys
import os

btagPath = os.path.dirname(__file__)
if btagPath not in sys.path:
    sys.path.append(btagPath)

from ControlPlots import *
from utils import *
from BtagEfficiencies import MakeBtagEfficienciesMaps
import bamboo.scalefactors

binningVariables = {
      "Eta"       : lambda obj : obj.eta
    , "ClusEta"   : lambda obj : obj.eta + obj.deltaEtaSC
    , "AbsEta"    : lambda obj : op.abs(obj.eta)
    , "AbsClusEta": lambda obj : op.abs(obj.eta + obj.deltaEtaSC)
    , "Pt"        : lambda obj : obj.pt
    }

def getTriggersAndPrimaryDatasets(year, fullEra, evt, isMC=False):
    if fullEra:
        era = fullEra[0] ## catch things like "C1" and "C2"
    else:
        era = ""
    hlt = evt.HLT
    def _getSel(hltSel):
        if str(hltSel) != hltSel:
            return [ getattr(hlt, sel) for sel in hltSel ]
        else:
            return [ getattr(hlt, hltSel) ]
    def forEra(hltSel, goodEras):
        if isMC or era in goodEras:
            return _getSel(hltSel)
        else:
            return []
    def notForEra(hltSel, badEras):
        if isMC or era not in badEras:
            return _getSel(hltSel)
        else:
            return []
    def fromRun(hltSel, firstRun, theEra, otherwise=True):
        if isMC:
            return _getSel(hltSel)
        elif fullEra == theEra:
            sel = _getSel(hltSel)
            return [ op.AND((evt.run >= firstRun), (op.OR(*sel) if len(sel) > 1 else sel[0])) ]
        elif otherwise:
            return _getSel(hltSel)
        else:
            return []
    if year == "2017":
        return { ## only consider eras B-F
                               # HLT_IsoMu24 is off for a 3.48/fb, HLT_IsoMu24_eta2p1 off for ~9/fb, HLT_TkMu100 not existing for first ~5/fb
            #"SingleMuon"     :([ #hlt.IsoMu24,  # prescaled for approx. 50/pb of Run2017B.  
            #                     hlt.IsoMu24_eta2p1, hlt.IsoMu27, 
            #                     hlt.Mu50 ]+ notForEra(("OldMu100", "TkMu100"), "B") ), # FIXME look if you want to keep these 
            
            #"SingleElectron" :( fromRun("Ele32_WPTight_Gsf_L1DoubleEG", 302026, "C2", otherwise=(era in "DEF"))+ # I need to require the electron to pass the L1 seeds of HLT_Ele35_WPTight_Gsf_v 
            #                   [ hlt.Ele35_WPTight_Gsf]+
            #                   notForEra("Ele115_CaloIdVT_GsfTrkIdT", "B")+[ hlt.Photon200 ]), # single electron (high pt)

            "MuonEG"         :([ hlt.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, 
                                 hlt.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, 
                                 hlt.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ ]+
                               notForEra(("Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", 
                                          "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", 
                                          "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL" #prescaled 
                                          ), "B")
                               # Introduced from menu version 3 on (missing for first ~14/fb)
                               #+fromRun(("Mu27_Ele37_CaloIdL_MW", "Mu37_Ele27_CaloIdL_MW"), 302026, "C2", otherwise=(era in ("D", "E", "F")))
                               ), #FIXME look as well if you should keep them 
            }
    elif year == "2018":
        return {
            "SingleMuon"     : [ hlt.IsoMu24, hlt.IsoMu27, hlt.Mu50, hlt.OldMu100, hlt.TkMu100 ], # OldMu100 and TkMu100 are recommend to recover inefficiencies at high pt (https://indico.cern.ch/event/766895/contributions/3184188/attachments/1739394/2814214/IdTrigEff_HighPtMu_Min_20181023_v2.pdf)
            "EGamma"         : [ hlt.Ele32_WPTight_Gsf, 
                                 hlt.Ele28_eta2p1_WPTight_Gsf_HT150, # Pieter he's not uusing this one : nut it's recommended by EgHLT however it  need additional cut for 2017   
                                 hlt.Ele115_CaloIdVT_GsfTrkIdT, hlt.Photon200],

            "MuonEG"         : [ hlt.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, 
                                 hlt.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
                                 
                                 hlt.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
                                 hlt.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,

                                 hlt.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, 
                                 hlt.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,

                                 hlt.Mu27_Ele37_CaloIdL_MW, 
                                 hlt.Mu37_Ele27_CaloIdL_MW ]
            }

def localizeSF(aPath, era):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "ScaleFactors", era, aPath)

def localize_trigger(aPath, era):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "TriggerEfficienciesStudies", era, aPath)

def localize_PileupJetID(aPath):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "ScaleFactors/run2PileupJetID", aPath)

def localize_HHMoriond17Trigger(aPath):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "TriggerEfficienciesStudies", aPath)

myScaleFactors = {
    "2017": {
            "electron_ID": {"cut_medium": localizeSF("Electron_EGamma_SF2D_2017_Medium_Fall17V2.json", "2017")},
            #"electron_ID": {"cut_tight": localizeSF("Electron_EGamma_SF2D_2017_Tight_Fall17V2.json", "2017")},
            "electron_reco_ptAbove20": localizeSF("Electron_EGamma_SF2D_RECO_2017RunBCDEF_ptG20_fromPOG.json", "2017"),
            "electron_trigger": localizeSF("Electron_ele28_ht150_OR_ele32.json", "2017"),
            "muon_ID": {"cut_medium": localizeSF("Muon_NUM_MediumID_DEN_genTracks_pt_abseta_{uncer}_2017RunBCDEF.json".format(uncer=uncer), "2017") for uncer in ("syst", "stat")},
            "muon_iso":{"iso_tight_id_medium": localizeSF("Muon_NUM_TightRelIso_DEN_MediumID_pt_abseta_{uncer}_2017RunBCDEF.json".format(uncer=uncer), "2017") for uncer in ("syst", "stat")},
            #"muon_ID": {"cut_tight": localizeSF("Muon_NUM_TightID_DEN_genTracks_pt_abseta_{uncer}_2017RunBCDEF.json".format(uncer=uncer), "2017") for uncer in ("syst", "stat")},
            #"muon_iso":{"iso_tight_id_tight": localizeSF("Muon_NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_{uncer}_2017RunBCDEF.json".format(uncer=uncer), "2017") for uncer in ("syst", "stat")},
            "muon_trigger": localizeSF("Muon_IsoMu27_RunBtoF.json", "2017"),
            "elmu_trigger": localizeSF("Triggers_emu2017.json", "2017"),
            
            # https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018#UL_2017_scales_and_smearings
            # https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2
            # https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffsRun2Legacy
            
            "muon_ID_17UL": {"cut_medium": localizeSF("Muon_NUM_MediumID_DEN_genTracks_pt_abseta_{uncer}_2017ULRunBCDEF.json".format(uncer=uncer), "2017UL") for uncer in ("syst", "stat")},
            "muon_iso_17UL":{"iso_tight_id_medium": localizeSF("Muon_NUM_TightRelIso_DEN_MediumID_pt_abseta_{uncer}_2017ULRunBCDEF.json".format(uncer=uncer), "2017UL") for uncer in ("syst", "stat")},
            #"muon_ID_17UL": {"cut_tight": localizeSF("Muon_NUM_TightID_DEN_genTracks_pt_abseta_{uncer}_2017ULRunBCDEF.json".format(uncer=uncer), "2017UL") for uncer in ("syst", "stat")},
            #"muon_iso_17UL":{"iso_tight_id_tight": localizeSF("Muon_NUM_TightRelTkIso_DEN_HighPtIDandIPCut_pair_newTuneP_probe_pt_abseta_{uncer}_2017ULRunBCDEF.json".format(uncer=uncer), "2017UL") for uncer in ("syst", "stat")},
            "electron_reco_ptAbove20_17UL": localizeSF("Electron_EGamma_SF2D_reco_ptAbove20_2017UL.json", "2017UL"),
            "electron_reco_ptBelow20_17UL": localizeSF("Electron_EGamma_SF2D_reco_ptBelow20_2017UL.json", "2017UL"),
            "electron_ID_17UL": {"cut_medium": localizeSF("Electron_EGamma_SF2D_2017UL_Medium.json", "2017UL")},
            #"electron_ID_17UL": {"cut_tight": localizeSF("Electron_EGamma_SF2D_2017UL_Tight.json", "2017UL")},
            
            },
    "2018": {
            
            "electron_ID": {"cut_medium": localizeSF("Electron_EGamma_SF2D_2018_Medium_Fall17V2.json", "2018")},
            "electron_reco": localizeSF("Electron_EGamma_SF2D_RECO_2018_fromPOG.json", "2018"),
            "electron_trigger": localize_trigger("Electron_ele28_ht150_OR_ele32_etaCut.json", "2018"),
            "muon_ID": {"cut_medium": localizeSF("Muon_NUM_MediumID_DEN_TrackerMuons_pt_abseta_{uncer}_2018RunABCD.json".format(uncer=uncer), "2018") for uncer in ("syst", "stat")},
            "muon_iso":{"iso_tight_id_medium": localizeSF("Muon_NUM_TightRelIso_DEN_MediumID_pt_abseta_{uncer}_2018RunABCD.json".format(uncer=uncer), "2018") for uncer in ("syst", "stat")},
            "muon_trigger":[(["Run315264to316360"], localizeSF("Muon_IsoMu24_BeforeMuonHLTUpdate.json", "2018")),
                                (["Run316361to325175"], localizeSF("Muon_IsoMu24_AfterMuonHLTUpdate.json", "2018"))],
            #"muon_trigger": {"afterHLTupdate":localize_trigger("{trig}_PtEtaBins_2018AfterMuonHLTUpdate.json".format(trig=trig),"2018")
            #                for trig in ("IsoMu24_OR_IsoTkMu24","Mu50_OR_OldMu100_OR_TkMu100" )},
            },
              
            "mueleLeg_HHMoriond17_2016" : tuple(localize_HHMoriond17Trigger("{wp}.json".format(wp=wp)) 
                                                for wp in ("Muon_XPathIsoMu23leg", "Muon_XPathIsoMu8leg", "Electron_IsoEle23Leg", "Electron_IsoEle12Leg")),
            "elemuLeg_HHMoriond17_2016" : tuple(localize_HHMoriond17Trigger("{wp}.json".format(wp=wp)) 
                                                for wp in ("Electron_IsoEle23Leg", "Electron_IsoEle12Leg", "Muon_XPathIsoMu23leg", "Muon_XPathIsoMu8leg")),
            "doubleEleLeg_HHMoriond17_2016" : tuple(localize_HHMoriond17Trigger("{wp}.json".format(wp=wp))
                                                for wp in ("Electron_IsoEle23Leg", "Electron_IsoEle12Leg", "Electron_IsoEle23Leg", "Electron_IsoEle12Leg")),
            "doubleMuLeg_HHMoriond17_2016" : tuple(localize_HHMoriond17Trigger("{wp}.json".format(wp=wp))
                                                for wp in ("Muon_DoubleIsoMu17Mu8_IsoMu17leg", "Muon_DoubleIsoMu17TkMu8_IsoMu8legORTkMu8leg", "Muon_DoubleIsoMu17Mu8_IsoMu17leg", "Muon_DoubleIsoMu17TkMu8_IsoMu8legORTkMu8leg")),

        }

def getScaleFactor(objType, key, periods=None, combine=None, additionalVariables=dict(), getFlavour=None, isElectron=False, systName=None):
    return bamboo.scalefactors.get_scalefactor(objType, key, periods=periods, combine=combine,
                                        additionalVariables=additionalVariables,
                                        sfLib=myScaleFactors,
                                        paramDefs=bamboo.scalefactors.binningVariables_nano,
                                        getFlavour=getFlavour,
                                        isElectron=isElectron,
                                        systName=systName)

class makeYieldPlots:
    def __init__(self):
        self.calls = 0
        self.plots = []
    def addYields(self, sel, name, title):
        """
            Make Yield plot and use it also in the latex yield table
            sel     = refine selection
            name    = name of the PDF to be produced
            title   = title that will be used in the LateX yield table
        """
        self.plots.append(Plot.make1D("Yield_"+name,   
                        op.c_int(0),
                        sel,
                        EqB(1, 0., 1.),
                        title = title + " Yield",
                        plotopts = {"for-yields":True, "yields-title":title, 'yields-table-order':self.calls}))
        self.calls += 1
    def returnPlots(self):
        return self.plots

def getL1PreFiringWeight(tree):
    return op.systematic(tree.L1PreFiringWeight_Nom,
                        name="L1PreFiring",
                        up=tree.L1PreFiringWeight_Up,
                        down=tree.L1PreFiringWeight_Dn)

lowerPtBinEdges=[30,50,70,100,140]#,200,300, 600] doesn't worth going further , check pt vs eta 2D maps 
def getPtLabel(iPT):
    ptBinLabel=""
    if(iPT==-1):
        ptBinLabel+="Inclusive"
    else:
        ptBinLabel+=str(lowerPtBinEdges[iPT])+"to"
        if(iPT==len(lowerPtBinEdges)-1):
            ptBinLabel+="Inf"
        else:
            ptBinLabel+=str(lowerPtBinEdges[iPT+1])
    return ptBinLabel

class TTbarDileptonMeasurment(NanoAODHistoModule):
    """
    This class is for b-tagging efficiencies measurments with a tag counting method in ttbar dilepton 
    events at 13TeV 
    samples are for 2018 && 2017 data in NanoAODv5 format 

    """
    def __init__(self, args):
        super(TTbarDileptonMeasurment, self).__init__(args)
        self.plotDefaults = {
                            "y-axis"           : "Events",
                            "log-y"            : "both",
                            "y-axis-show-zero" : True,
                            "save-extensions"  : ["pdf", "png"],
                            "show-ratio"       : True,
                            "sort-by-yields"   : False,
                            }

        self.doSysts = self.args.systematic
        self.doYields = self.args.yields
        self.doHLT = self.args.hlt
    def addArgs(self, parser):
        super(TTbarDileptonMeasurment, self).addArgs(parser)
        parser.add_argument("--backend", type=str, default="dataframe", help="Backend to use, 'dataframe' (default) or 'lazy'")
        parser.add_argument("-s", "--systematic", action="store_true", help="Produce systematic variations")
        parser.add_argument("-y", "--yields", action="store_true", help="Produce Yield latex table at different stage of selection")
        parser.add_argument("-hlt", "--hlt", action="store_true", help="Produce 2D Maps for HLT desicions")

    def prepareTree(self, tree, sample=None, sampleCfg=None):
        era = sampleCfg.get("era") if sampleCfg else None
        isMC = self.isMC(sample)
        metName = "METFixEE2017" if era =="2017" and "_UL17" not in sample else "MET"
        year = sampleCfg.get("era")
        eraInYear = "" if isMC else next(tok for tok in sample.split("_") if tok.startswith(year))[4:]
        ## initializes tree.Jet.calc so should be called first (better: use super() instead)
        # JEC's Recommendation for Full RunII: https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        # JER : -----------------------------: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
       
        from bamboo.treedecorators import NanoAODDescription, nanoRochesterCalc, nanoJetMETCalc, nanoJetMETCalc_METFixEE2017

        tree,noSel,be,lumiArgs = NanoAODHistoModule.prepareTree(self, tree, 
                                                                sample=sample, 
                                                                sampleCfg=sampleCfg, 
                                                                description=NanoAODDescription.get("v5", year=(era if era else "2016"), 
                                                                isMC=isMC, 
                                                                systVariations=[ nanoRochesterCalc, 
                                                                    (nanoJetMETCalc_METFixEE2017 if era == "2017" and "_UL17" not in sample else nanoJetMETCalc) ]),
                                                                lazyBackend   = (self.args.backend == "lazy")) ## will do Jet and MET variations, and the Rochester correction
        triggersPerPrimaryDataset = {}
        jec, smear, jesUncertaintySources = None, None, None

        from bamboo.analysisutils import configureJets, configureType1MET, configureRochesterCorrection
        isNotWorker = (self.args.distributed != "worker") 
        

        if era == "2018":
            configureRochesterCorrection(tree._Muon, os.path.join(os.path.dirname(__file__), "data", "RoccoR2018.txt"), isMC=isMC, backend=be, uName=sample)

            if self.isMC(sample):
                jec="Autumn18_V19_MC"
                smear="Autumn18_V7b_MC"
                jesUncertaintySources=["Total"]

            else:
                if "2018A" in sample:
                    jec="Autumn18_RunA_V19_DATA"
                elif "2018B" in sample:
                    jec="Autumn18_RunB_V19_DATA"
                elif "2018C" in sample:
                    jec="Autumn18_RunC_V19_DATA"
                elif "2018D" in sample:
                    jec="Autumn18_RunD_V19_DATA"

        elif era == "2017":
            
            configureRochesterCorrection(tree._Muon, os.path.join(os.path.dirname(__file__), "data", "RoccoR2017.txt"), isMC=isMC, backend=be, uName=sample)
            
            if self.isMC(sample):
                if "_UL17" in sample:
                    jec= "Summer19UL17_V5_MC"
                    smear= "Summer19UL17_JRV2_MC"
                    jesUncertaintySources=["Total"]
                else:
                    jec="Fall17_17Nov2017_V32_MC"
                    smear="Fall17_V3b_MC"
                    jesUncertaintySources=["Total"]

            else:
                if "_UL17" in sample:
                    if "2017B" in sample:
                        jec = "Summer19UL17_RunB_V5_DATA"
                    elif "2017C" in sample:
                        jec = "Summer19UL17_RunC_V5_DATA"
                    elif "2017D" in sample:
                        jec = "Summer19UL17_RunD_V5_DATA"
                    elif "2017E" in sample:
                        jec = "Summer19UL17_RunE_V5_DATA"
                    elif "2017F" in sample:
                        jec = "Summer19UL17_RunF_V5_DATA"
                else:
                    if "2017B" in sample:
                        jec="Fall17_17Nov2017B_V32_DATA"
                    elif "2017C" in sample:
                        jec="Fall17_17Nov2017C_V32_DATA"
                    elif "2017D" in sample or "2017E" in sample:
                        jec="Fall17_17Nov2017DE_V32_DATA"
                    elif "2017F" in sample:
                        jec="Fall17_17Nov2017F_V32_DATA"

        try:
            configureJets(tree._Jet, "AK4PFchs", jec=jec, smear=smear, jesUncertaintySources=jesUncertaintySources, mayWriteCache=isNotWorker, isMC=isMC, backend=be, uName=sample)
        except Exception as ex:
            logger.exception("Problem while configuring jet correction and variations")
        
        ## Configure MET
        try:
            configureType1MET(getattr(tree, f"_{metName}"), jec=jec, smear=smear, jesUncertaintySources=jesUncertaintySources, mayWriteCache=isNotWorker, isMC=isMC, backend=be, uName=sample)
        except Exception as ex:
            logger.exception("Problem while configuring MET correction and variations")
       
        if self.doHLT:
            if self.isMC(sample):
                if  self.doSysts:
                    logger.info("Adding FSR, ISR and QCD scale variations")
                    noSel = utils.addTheorySystematics(self, tree, noSel, qcdScale=False, PSISR=False, PSFSR=False)
                    getpartOf_selection(noSel)
        else:
            triggersPerPrimaryDataset = getTriggersAndPrimaryDatasets(year, eraInYear, tree, isMC=isMC) 
            if self.isMC(sample):
                noSel = noSel.refine("genWeight", weight=tree.genWeight, cut=op.OR(*chain.from_iterable(triggersPerPrimaryDataset.values())), autoSyst=self.doSysts)
                if self.doSysts:
                    logger.info("Adding FSR, ISR and QCD scale variations")
                    noSel = utils.addTheorySystematics(self, tree, noSel, qcdScale=True, PSISR=False, PSFSR=False)
            else:
                noSel = noSel.refine("withTrig", cut=makeMultiPrimaryDatasetTriggerSelection(sample, triggersPerPrimaryDataset) )
            #cuts = []
            #if era == '2017':
            #    cuts.append(op.OR(
            #            op.AND(HLT.Ele32_WPTight_Gsf_L1DoubleEG, op.rng_any(TrigObj, lambda obj: op.AND(op.deltaR(obj.p4, ele.p4) < 0.1, obj.filterBits & 1024))),
            #                    HLT.Ele28_eta2p1_WPTight_Gsf_HT150))
            #    noSel = noSel.refine("additionalTrig", cut=cuts)
            
            #cuts = []
            #if era == '2017':
            #    cuts.append(op.OR(
            #            op.AND(HLT.Ele32_WPTight_Gsf_L1DoubleEG, op.rng_any(TrigObj, lambda obj: op.AND(obj.id == 11, obj.filterBits & 1024))),
            #                    HLT.Ele35_WPTight_Gsf))
            #    noSel = noSel.refine("additionalTrig", cut=cuts)
        return tree,noSel,be,lumiArgs
    
    def definePlots(self, t, noSel, sample=None, sampleCfg=None):    

        from bamboo.analysisutils import forceDefine
        from bamboo.plots import EquidistantBinning as EqBin
        from bambooToOls import Plot
        from bamboo.plots import SummedPlot
        from bamboo import treefunctions as op
        from METFilter_xyCorr import METFilter, METcorrection
        from bamboo.analysisutils import makePileupWeight
        from systematics import getTriggerSystematcis

        isMC = self.isMC(sample)
        era = sampleCfg.get("era") if sampleCfg else None

        noSel = noSel.refine("passMETFlags", cut=METFilter(t.Flag, era, isMC, sample) )
        puWeightsFile = None
        mcprofile= None
        TriggerFrom_HHMoriond17 =True
        TriggerFrom_tth = False
        yield_object = makeYieldPlots()
        binScaling = 1
        plots = []
        self.eff_plots= []

        if era == "2018":
            suffix= '2018_Autumn18'
            puWeightsFile = os.path.join(os.path.dirname(__file__), "data/PileupFullRunII/", "puweights2018_Autumn18.json")
        elif era == "2017":
            suffix ='2017_Fall17'
            if sample in ['DYJetsToLL_M-10to50-madgraphMLM', 'ST_tW_antitop_5f', 'WZ', 'ZZ'] and '_UL17' not in sample: 
                if "pufile" not in sampleCfg:
                    raise KeyError("Could not find 'pufile' entry for sample %s in the YAML file"%sampleCfg["sample"])
                mcprofile= os.path.join(os.path.dirname(__file__), "data/PileupFullRunII/mcprofile/", "%s_2017.json"%sample)
            else:
                puWeightsFile = os.path.join(os.path.dirname(__file__), "data/PileupFullRunII/", "puweights2017_Fall17.json")

        if self.isMC(sample):
            if mcprofile is not None:
                PUWeight = makePileupWeight(mcprofile, t.Pileup_nTrueInt, variation="Nominal",
                                                   nameHint="puWeight_{0}".format(sample.replace('-','_')))
            else:
                PUWeight = makePileupWeight(puWeightsFile, t.Pileup_nTrueInt, systName="puweights%s"%suffix)
            noSel = noSel.refine("puWeight", weight=PUWeight)


        if self.isMC(sample) and sampleCfg["group"] in ['TTJets_SL', 'TTJets_DL', 'TTJets_AH']:
            # https://indico.cern.ch/event/904971/contributions/3857701/attachments/2036949/3410728/TopPt_20.05.12.pdf
            genTop_pt = op.select(t.GenPart, lambda gp : op.AND((gp.statusFlags & (0x1<<13)), gp.pdgId==6))
            gen_antiTop_pt = op.select(t.GenPart, lambda gp : op.AND((gp.statusFlags & (0x1<<13)), gp.pdgId==-6))
            
            gen_pt = op.map(genTop_pt, lambda top: top.pt)
            forceDefine(gen_pt, noSel)
            plots.append(Plot.make1D("gen_ToppT_noReweighting".format(suffix), gen_pt, noSel, EqB(60 // binScaling, 0., 1200.), title="gen Top p_{T} [GeV]"))

            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#Use_case_3_ttbar_MC_is_used_to_m 
            scalefactor = lambda t : op.exp(-2.02274e-01 + 1.09734e-04*t.pt -1.30088e-07*t.pt**2 + (5.83494e+01/(t.pt+1.96252e+02)))
            top_weight = lambda top, antitop : op.sqrt(scalefactor(top)*scalefactor(antitop))
    
            noSel = noSel.refine("top_reweighting", weight=top_weight(genTop_pt[0], gen_antiTop_pt[0]))
            forceDefine(gen_pt, noSel)
            plots.append(Plot.make1D("gen_ToppT_withReweighting".format(suffix), gen_pt, noSel, EqB(60 // binScaling, 0., 1200.), title="gen Top p_{T} [GeV]"))

        forceDefine(t._Muon.calcProd, noSel)

        sorted_muons = op.sort(t.Muon, lambda mu : -mu.pt)
        muons = op.select(sorted_muons, lambda mu : op.AND(mu.pt > 10., op.abs(mu.eta) < 2.4, mu.mediumId, mu.pfRelIso04_all<0.15, op.abs(mu.sip3d) < 4.))
        
        suffix = ("_17UL" if "_UL17" in sample else (""))
        muonIDSF = getScaleFactor("lepton", (era, "muon_ID%s"%suffix, "cut_medium"), systName="muon_ID")
        muonIsoSF = getScaleFactor("lepton", (era, "muon_iso%s"%suffix, "iso_tight_id_medium"), systName="muon_iso")
        #muonIDSF = getScaleFactor("lepton", (era, "muon_ID%s"%suffix, "cut_tight"), systName="muon_ID")
        #muonIsoSF = getScaleFactor("lepton", (era, "muon_iso%s"%suffix, "iso_tight_id_tight"), systName="muon_iso")
        #muonTriggerSF = getScaleFactor("lepton", (era, "muon_trigger"), combine="weight", systName="muon_trigger")

        sorted_electrons = op.sort(t.Electron, lambda ele : -ele.pt)
        electrons = op.select(sorted_electrons, lambda ele : op.AND(ele.pt > 15., op.abs(ele.eta) < 2.5 , ele.cutBased>=3, op.abs(ele.sip3d) < 4., 
                                                                    op.OR(op.AND(op.abs(ele.dxy) < 0.05, op.abs(ele.dz) < 0.1), 
                                                                          op.AND(op.abs(ele.dxy) < 0.05, op.abs(ele.dz) < 0.2) ))) 
        
        eleIDSF = getScaleFactor("lepton", (era, "electron_ID%s"%suffix, "cut_medium"), isElectron=True, systName="ele_ID")
        #eleIDSF = getScaleFactor("lepton", (era, "electron_ID%s"%suffix, "cut_tight"), isElectron=True, systName="ele_ID")
        eleRecoSF = getScaleFactor("lepton", (era, "electron_reco_ptAbove20%s"%suffix), isElectron=True, systName="ele_reco")
        #eleTriggerSF = getScaleFactor("lepton", (era, "electron_trigger"), isElectron=True, systName="ele_trigger")

        sorted_jets=op.sort(t.Jet, lambda j : -j.pt)
        jetsSel = op.select(sorted_jets, lambda j : op.AND(j.pt > 30., op.abs(j.eta)< 2.5, j.jetId &4))
                                                            #,op.AND( op.in_range(30., j.pt, 50.), j.puId==7 ))) #j.jetId &2)) old cut 
        
        # exclude from the jetsSel any jet that happens to include within its reconstruction cone a muon or an electron.
        jets= op.select(jetsSel, lambda j : op.AND(op.NOT(op.rng_any(electrons, lambda ele : op.deltaR(j.p4, ele.p4) < 0.4 )), 
                                                   op.NOT(op.rng_any(muons, lambda mu : op.deltaR(j.p4, mu.p4) < 0.4 ))))

        #jetpuId = getScaleFactor("jet", ("JetId_InHighPileup_{0}_94X".format(era).replace("94X", "102X" if era=="2018" else "94X"), "puId_T"), systName="JetpuID_T")

        t_bDisc = {
            "DeepFlavour" : lambda j : j.btagDeepFlavB,
            "DeepCSV"     : lambda j : j.btagDeepB
            }
        
        bTaggingWPs_era = {
            "2016" : {
                "DeepCSV"     : { "L":0.2217, "M":0.6321, "T":0.8953 },
                "DeepFlavour" : { "L":0.0614, "M":0.3093, "T":0.7221 },
                },
            "2017" : {
                "DeepCSV"     : { "L":0.1522, "M":0.4941, "T":0.8001 },
                "DeepFlavour" : { "L":0.0521, "M":0.3033, "T":0.7489 },
                },
            "2018" : {
                "DeepCSV"     : { "L":0.1241, "M":0.4184, "T":0.7527 },
                "DeepFlavour" : { "L":0.0494, "M":0.2770, "T":0.7264 },
                }
            }
        #DoubleB_2017MistagSF = {
        #        }

        wp_text = {"L": "loose", "M": "medium", "T": "tight"}
        
        MET = t.METFixEE2017 if era == "2017" and "UL17" not in sample else t.MET
        corrMET=METcorrection(MET, t.PV, sample, era, self.isMC(sample))
        
        L1Prefiring = 1.
        HLTZvtxSF = 1.
        if era in ["2016", "2017"]:
            if sample in ['TTToFullLeptonic_UL17', 'TTToSemiLeptonic_UL17', 'WW_UL17', 'WZ_UL17', 'ZZ_UL17']:
                L1Prefiring = 1.
            else:
                L1Prefiring = getL1PreFiringWeight(t)
        if era =='2017':
            HLTZvtxSF = op.systematic(op.c_float(0.991), name='HLT_Zvtx_eff', up=op.c_float(0.992), down=op.c_float(0.990)) 
        
        #The contribution from Z + jets events is reduced by applying a veto around the Z boson mass when the two leptons have the same flavour, ( | Mll -MZ | > 10 GeV).
        osdilep = lambda l1,l2 : op.AND(l1.charge != l2.charge, op.invariant_mass(l1.p4, l2.p4) >12. )
        Zboson_Vetocut = lambda l1, l2 : op.AND(osdilep(l1, l2), op.NOT(op.in_range(75., op.invariant_mass(l1.p4, l2.p4), 115.)))
        osdilep_Z = lambda l1,l2 : op.AND(l1.charge != l2.charge, op.in_range(75., op.invariant_mass(l1.p4, l2.p4), 115.))

        osll = {
                # combined cat
                #"elmu": op.combine((electrons, muons), pred=lambda ele,mu : osdilep(ele,mu)),

                "elmu": op.combine((electrons, muons), pred=lambda ele,mu : op.AND(osdilep(ele,mu), ele.pt > mu.pt, op.NOT(op.in_range(75., op.invariant_mass(ele.p4, mu.p4), 115.)))),
                "muel": op.combine((muons, electrons), pred=lambda mu,ele : op.AND(osdilep(mu,ele), mu.pt > ele.pt, op.NOT(op.in_range(75., op.invariant_mass(mu.p4, ele.p4), 115.)))),

                #"mumu" : op.combine(muons, N=2, pred=osdilep_Z),
                #"elel" : op.combine(electrons, N=2, pred=osdilep_Z)
                #"mumu_SR" : op.combine(muons, N=2, pred=Zboson_Vetocut),
                #"mumu_CR" : op.combine(muons, N=2, pred=Zboson_Vetocut),
                #"elel_SR" : op.combine(electrons, N=2, pred=osdilep_Z),
                #"elel_CR" : op.combine(electrons, N=2, pred=osdilep_Z),
                }
        
        if TriggerFrom_HHMoriond17:
            elemuTrigSF = getScaleFactor("dilepton", ("elemuLeg_HHMoriond17_2016"), systName="elmutrig")
            mueleTrigSF = getScaleFactor("dilepton", ("mueleLeg_HHMoriond17_2016"), systName="mueltrig")
            doubleMuTrigSF = getScaleFactor("dilepton", ("doubleMuLeg_HHMoriond17_2016"), systName="mumutrig")
            doubleEleTrigSF = getScaleFactor("dilepton", ("doubleEleLeg_HHMoriond17_2016"), systName="eleltrig")
        elif TriggerFrom_tth:
            elemuTrigSF = getTriggerSystematcis(era, osll.get('elmu')[0], 'elmu', "tth")
            mueleTrigSF = getTriggerSystematcis(era, osll.get('muel')[0], 'muel', "tth")
            # FIXME get SAME FLAVOUR TRIGGER SFS
        else:
            elemuTrigSF = getScaleFactor("lepton", (era, "elmu_trigger"), systName="elmutrig")
            mueleTrigSF = getScaleFactor("lepton", (era, "elmu_trigger"), systName="mueltrig")
            # FIXME get SAME FLAVOUR TRIGGER SFS
        if TriggerFrom_HHMoriond17:
            osllSFs = {
                    "elmu" : (lambda lep :  [ HLTZvtxSF, L1Prefiring, muonIDSF(lep[1]), muonIsoSF(lep[1]) , eleRecoSF(lep[0]), eleIDSF(lep[0]), elemuTrigSF(lep)]),# eleTriggerSF(lep[0]), muonTriggerSF(lep[1])]),
                    "muel" : (lambda lep :  [ HLTZvtxSF, L1Prefiring, muonIDSF(lep[0]), muonIsoSF(lep[0]) , eleRecoSF(lep[1]), eleIDSF(lep[1]), mueleTrigSF(lep)]),# eleTriggerSF(lep[1]), muonTriggerSF(lep[0])])
                    "mumu" : (lambda lep :  [ L1Prefiring, muonIDSF(lep[0]), muonIsoSF(lep[0]) , muonIDSF(lep[1]), muonIsoSF(lep[1]),  doubleMuTrigSF(lep)]),# eleTriggerSF(lep[1]), muonTriggerSF(lep[0])])
                    "elel" : (lambda lep :  [ HLTZvtxSF, L1Prefiring, eleRecoSF(lep[0]), eleIDSF(lep[0]), eleRecoSF(lep[1]), eleIDSF(lep[1]), doubleEleTrigSF(lep)]),# eleTriggerSF(lep[1]), muonTriggerSF(lep[0])])
                }
        else:
            osllSFs = {
                    "elmu" : (lambda lep :  [ HLTZvtxSF, L1Prefiring, muonIDSF(lep[1]), muonIsoSF(lep[1]) , eleRecoSF(lep[0]), eleIDSF(lep[0]), elemuTrigSF(lep[0]), elemuTrigSF(lep[1])]),# eleTriggerSF(lep[0]), muonTriggerSF(lep[1]) ]),
                    "muel" : (lambda lep :  [ HLTZvtxSF, L1Prefiring, muonIDSF(lep[0]), muonIsoSF(lep[0]) , eleRecoSF(lep[1]), eleIDSF(lep[1]), mueleTrigSF(lep[0]), mueleTrigSF(lep[1])]),# eleTriggerSF(lep[1]), muonTriggerSF(lep[0])])
                    "mumu" : (lambda lep :  [ L1Prefiring, muonIDSF(lep[0]), muonIsoSF(lep[0]) , muonIDSF(lep[1]), muonIsoSF(lep[1]),  doubleMuTrigSF(lep[0]), doubleMuTrigSF(lep[1])]),# eleTriggerSF(lep[1]), muonTriggerSF(lep[0])])
                    "elel" : (lambda lep :  [ HLTZvtxSF, L1Prefiring, eleRecoSF(lep[0]), eleIDSF(lep[0]), eleRecoSF(lep[1]), eleIDSF(lep[1]), doubleEleTrigSF(lep[0]), doubleEleTrigSF(lep[1])]),# eleTriggerSF(lep[1]), muonTriggerSF(lep[0])])
                }

        osllcatrng = lambda catrng : op.AND(op.rng_len(catrng) > 0, catrng[0][0].pt > 25.)
        hasOSLL = noSel.refine("hasOSLL", cut=op.OR(*( osllcatrng(rng) for rng in osll.values())))
        
        forceDefine(t._Jet.calcProd, hasOSLL)
        forceDefine(getattr(t, "_{0}".format("METFixEE2017" if era =="2017" and "UL17" not in sample else "MET")).calcProd, hasOSLL)
        
        categories = dict((channel, (leadpair[0],
            hasOSLL.refine("{0}_EventSelection".format(channel), cut=osllcatrng(leadpair), 
                                                                 weight=(osllSFs[channel if channel.find('_')==-1 else (channel.split('_')[0])](leadpair[0]) if isMC else None))
                                                                        )) for channel, leadpair in osll.items())

        make_PrimaryANDSecondaryVerticesPlots = True
        make_LeptonsPlots_LepSel = True #*
        make_LeptonsPlots_LeptonPlusJetsSel = True
        make_JetsPlots = True #*
        make_METPlots = False
        make_bJetsPlots = True #*
        make_btagScore = True #*
        make_TwoTagCountPlots = False
        make_ttbarEstimationPlots = False #*
        make_efficienciesmaps = False
        make_deltaRPlots = False
        make_2Dmaps = False
        make_significance = False
        #==================================================================================================================================
        #==================================================================================================================================
        if make_efficienciesmaps and isMC:
            plots.extend(MakeBtagEfficienciesMaps(jets, categories, era))

        optstex = {
                "elel" : 'e^+e^-',
                "mumu" : '\mu^+\mu^-',
                "muel" : '\mu^+e^-',
                "elmu" : 'e^+\mu^-'
                }
            
        jLenCuts = {
                "Only2"    : lambda rng : op.rng_len(rng) == 2,
                "atleast2" : lambda rng : op.rng_len(rng) >= 2
                }

        plotsToSum_twoTagCount = defaultdict(list)

        for channel, (leptons, cat) in categories.items():

            if make_significance:
                if '_' not in channel:
                    pass
                else:
                    for nbrj, njcut in jLenCuts.items():
                        TwoLepTwoJets = cat.refine(f"{nbrj}Jets{channel}Sel_toOptimizeMETcut", cut=njcut(jets))
                        for tagger,workingPoints in bTaggingWPs_era[era].items():
                            jets_by_desc_discri = op.sort(jets, partial((lambda j,fun=None : -fun(j)), fun=t_bDisc[tagger]))
                            for wp,thresh in workingPoints.items():
                                Twobjets = op.select(jets_by_desc_discri, partial((lambda j,fun=None,val=None : fun(j) >= val), fun=t_bDisc[tagger], val=thresh))
                                TwoLepTwoBJets = TwoLepTwoJets.refine(f"TwoLeptons{nbrj}BJets_{key}_{channel}_toOptimizeMETcut", cut=njcut(Twobjets))
                                for suffix, Object in {"MET": MET, "corrMET": corrMET}.items():
                                    plots.append(Plot.make1D(f"{channel}_{nbrj}j_{key}_{suffix}_pt", Object.pt, TwoLepTwoBJets,
                                                            EqBin(60 // binScaling, 0., 600.), title="%s p_{T} [GeV]"%suffix,
                                                            plotopts=utils.getOpts(channel.split('_')[0], **{"log-y": False})))
                                        
            if channel.find('_')!=-1:
                continue
            if self.doHLT:
                SingleLeptonsTriggers = {
                    "SingleElectron": [ t.Ele32_WPTight_Gsf_L1DoubleEG, t.HLT.Ele35_WPTight_Gsf, t.HLT.Ele115_CaloIdVT_GsfTrkIdT, t.HLT.Photon200],
                    "SingleMuon":     [ t.HLT.IsoMu24_eta2p1, t.HLT.IsoMu27, t.HLT.Mu50, t.HLT.OldMu100, t.HLT.TkMu100]
                    } 
                DoubleLeptonsTriggers = {
                    "MuonEG": [ t.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, t.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, t.HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
                                t.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, t.HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, t.HLT.Mu27_Ele37_CaloIdL_MW, t.HLT.Mu37_Ele27_CaloIdL_MW]
                    }

                passSingleLepton = Plot.make1D("OSll_HLTSingleLeptonDecision_{0}".format(channel), 
                                                op.OR(*chain.from_iterable(SingleLeptonsTriggers.values())), 
                                                cat, EqB(2, 0., 2.), title="Single lepton trigger decision", plotopts=utils.getOpts(channel))
                passDiLepton = Plot.make1D("OSll_HLTDiLeptonDecision_{0}".format(channel), 
                                                op.OR(*chain.from_iterable(DoubleLeptonsTriggers.values())), 
                                                cat, EqB(2, 0., 2.), title="Dilepton trigger decision", plotopts=utils.getOpts(channel))
                passboth = Plot.make1D("OSll_HLTDiLeptonANDSiLeptonDecision_{0}".format(channel), 
                                                op.OR(*chain.from_iterable(DoubleLeptonsTriggers.values()), *chain.from_iterable(SingleLeptonsTriggers.values())), 
                                                cat, EqB(2, 0., 2.), title="Dilepton || Single Lepton trigger decision", plotopts=utils.getOpts(channel))
                
                Triggerscuts = {
                    "SingleLeptonsTriggers": [t.Ele32_WPTight_Gsf_L1DoubleEG, t.HLT.Ele35_WPTight_Gsf, t.HLT.Ele115_CaloIdVT_GsfTrkIdT, t.HLT.Photon200, 
                                              t.HLT.IsoMu24_eta2p1, t.HLT.IsoMu27, t.HLT.Mu50, t.HLT.OldMu100, t.HLT.TkMu100],
                    "DoubleLeptonsTriggers": [ t.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, t.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, t.HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
                                                t.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, t.HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, t.HLT.Mu27_Ele37_CaloIdL_MW, t.HLT.Mu37_Ele27_CaloIdL_MW]
                            }
                plots.append(passSingleLepton)
                plots.append(passDiLepton)
                plots.append(passboth)
                pt_bin_lead = VariableBinning([25., 30., 40., 50., 60., 80., 100., 120., 150., 200.])
                pt_bin_trail = VariableBinning(([10., 20.] if channel != "elmu" else [15., 20.])+[30., 40., 50., 60., 80., 100., 120., 150., 200.])
                denom_lead = Plot.make1D("{0}_leadpt_DEN".format(channel), leptons[0].Pt, cat, pt_bin_lead)
                denom_trail = Plot.make1D("{0}_trailpt_DEN".format(channel), leptons[1].Pt, cat, pt_bin_trail)
                denom_2D = Plot.make2D("{0}_pt_DEN".format(channel), (leptons[0].Pt, leptons[1].Pt), cat, (pt_bin_lead, pt_bin_trail))
                plots += [ denom_lead, denom_trail, denom_2D ]
                for trigName, trigSel in trigCuts.items():
                    pref = "{0}_{1}_trig".format(channel, trigName)
                    catOn_HLT = cat.refine("dilep{0}".format(pref), cut=trigSel)
                    num_lead = Plot.make1D("{0}_lead_NUM".format(pref), leptons[0].Pt, catOn_HLT, pt_bin_lead)
                    num_trail = Plot.make1D("{0}_trail_NUM".format(pref), leptons[1].Pt, catOn_HLT, pt_bin_trail)
                    num_2D = Plot.make2D("{0}_pt_NUM".format(pref), (leptons[0].Pt, leptons[1].Pt), catOn_HLT, (pt_bin_lead, pt_bin_trail))
                    plots += [ num_lead, num_trail, num_2D ]
                    self.eff_plots += [ (num_lead, denom_lead), (num_trail, denom_trail), (num_2D, denom_2D) ]
            
            plots.append(Plot.make1D(f"{channel}_NOcutOnJetsLen_Jet_multiplicity",
                op.rng_len(jets), cat,
                EqBin(10, 0., 10.), title="Jet multiplicity",
                plotopts=utils.getOpts(channel, **{"log-y": True})))
            if self.doYields:
                yield_object.addYields(cat,"hasOs%s"%channel,"OS leptons + M_{ll} cut (channel : %s)"%optstex[channel])
            if make_LeptonsPlots_LepSel:
                plots.extend(makeLeptonPlots(cat, leptons, '_2lepOSOFSel_', channel))

            for nbrj, njcut in jLenCuts.items():
                sel2l2j = cat.refine(f"{nbrj}Jets{channel}Sel", cut=njcut(jets))
                if self.doYields:
                    yield_object.addYields(sel2l2j,"2Lep%s_%sJets"%(nbrj, channel),"2Lep(OS) + %s Jets (channel : %s)"%(nbrj, optstex[channel]))

                if nbrj == 'atleast2':
                    plots.append(Plot.make1D(f"{channel}_{nbrj}Jets_Jet_multiplicity",
                        op.rng_len(jets), sel2l2j,
                        EqBin(10, 0., 10.), title="Jet multiplicity",
                        plotopts=utils.getOpts(channel, **{"log-y": True})))
                
                if make_LeptonsPlots_LeptonPlusJetsSel:
                    plots.extend(makeLeptonPlots(sel2l2j, leptons, '_2lep{0}jeSel_'.format(nbrj), channel))
                if make_JetsPlots:
                    plots.extend(makeJetPlots(sel2l2j, jets, nbrj, channel))
                if make_deltaRPlots:
                    plots.extend(makedeltaRPlots(sel2l2j, jets, leptons, nbrj, channel))
                if make_2Dmaps:
                    plots.extend(make2DMAPS(sel2l2j, jets, nbrj, channel))
                if make_PrimaryANDSecondaryVerticesPlots:
                    plots.extend(makePrimaryANDSecondaryVerticesPlots(sel2l2j, t, nbrj, channel))

                if make_bJetsPlots or make_btagScore or make_ttbarEstimationPlots:
                    for tagger,workingPoints in bTaggingWPs_era[era].items():
                        jets_by_desc_discri = op.sort(jets, partial((lambda j,fun=None : -fun(j)), fun=t_bDisc[tagger]))
                        for wp,thresh in workingPoints.items():
                            key = f"{tagger}{wp}"
                            print(f"Btagging: Era= {era}, Tagger={tagger}, Pass_{wp_text[wp]}_working_point={thresh}")
                            print(f"btag_{era}_{('102X' if era == '2018' else '94X')}", f"{tagger}_{wp_text[wp]}")
                            bjets_pass = op.select(jets_by_desc_discri, partial((lambda j,fun=None,val=None : fun(j) >= val), fun=t_bDisc[tagger], val=thresh))
                            sel2l2b = sel2l2j.refine(f"TwoLeptons{nbrj}BJets_{key}_{channel}", cut=njcut(bjets_pass))
                            sel2l2b_noMET = sel2l2b.refine(f"TwoLeptons{nbrj}Bjets_{key}_{channel}_plusmetcut", cut=[ corrMET.pt < 80. ])
                            if make_METPlots:     
                                plots.extend(makeMETPlots(sel2l2b, leptons, MET, corrMET, nbrj, key, channel))
                            for metcut, sel in {"METcut": sel2l2b_noMET, "noMETcut": sel2l2b}.items():
                                if self.doYields:
                                    yield_object.addYields(sel, f"{channel}_{nbrj}bJets_{tagger}{wp}","2Lep(OS) + {nbrj} bJets {tagger}{wp} + {metcut} (channel : {optstex[channel]})")
                                if make_bJetsPlots:
                                    plots.append(Plot.make1D(f"{channel}_{metcut}_{nbrj}j_{tagger}{wp}_Jet_mulmtiplicity",
                                        op.rng_len(bjets_pass), sel,
                                        EqBin(10, 0., 10.), title="Jet mulmtiplicity",
                                        plotopts=utils.getOpts(channel, **{"log-y": True})))
                                    plots += makeBJetPlots(sel, bjets_pass, key, channel, metcut, nbrj)
                                if make_btagScore:
                                    for i in range(2):
                                        plots.append(Plot.make1D(f"{channel}_{nbrj}j_jet{i+1:d}_{metcut}_discr_{key}",
                                            t_bDisc[tagger](bjets_pass[i]), sel,
                                            EqBin(60, thresh, 1.), title=f"{tagger}Disc {wp}",
                                            plotopts=utils.getOpts(channel)))
                            if make_ttbarEstimationPlots:
                                sel2l2b_MET = sel2l2b.refine(f"TwoLeptons{nbrj}Bjets_{key}_{channel}_highMETtail", cut=[ corrMET.pt > 80. ])
                                for metcut, sel in {"METcut": sel2l2b_noMET, "HighMET": sel2l2b_MET}.items():
                                    plots += makehistosforTTbarEstimation(sel, leptons, bjets_pass, "{nbrj}_jets_{tagger}_btag{wp}_mll_and_{met}", channel)

                ## two-tag count method
                sel2l2j_METCut = sel2l2j.refine(f"{nbrj}Jets{channel}selection_METcut", cut=[ corrMET.pt < 80.])
                if make_TwoTagCountPlots:
                    for iPT in range(-1, len(lowerPtBinEdges)):
                        ptBin = getPtLabel(iPT)
                        if ptBin == "Inclusive":
                            sel_twotag_pt = sel2l2j_METCut
                        else:
                            if "Inf"in ptBin:
                                minpt = 140.
                                maxpt = 200.
                            else:
                                minpt = lowerPtBinEdges[iPT]
                                maxpt = lowerPtBinEdges[iPT+1]
                            sel_twotag_pt = sel2l2j_METCut.refine(f"{sel2l2j_METCut.name}_2Jets_{ptBin}ptBin_cut{nbrj}", cut=[ op.in_range(minpt, jets[0].pt, maxpt), op.in_range(minpt, jets[1].pt, maxpt) ])
                        if isMC:
                            nTrueB = op.rng_count(jets[:2], lambda j : j.hadronFlavour == 5)
                            nTrueC = op.rng_count(jets[:2], lambda j : j.hadronFlavour == 4)
                            trueBC_cat = op.multiSwitch((nTrueB != 0, nTrueB+1), (nTrueC != 0, op.c_int(1)), op.c_int(0))
                        else:
                            nTrueB = op.c_int(0)
                            trueBC_cat = op.c_int(0)
                        sel_twotag_pt_twoTrueB = sel_twotag_pt.refine(f"{sel_twotag_pt.name}_twoTrueB", cut=(nTrueB == 2))
                        for tagger,workingPoints in bTaggingWPs_era[era].items():
                            for wp,thresh in workingPoints.items():
                                key = f"{tagger}{wp}"
                                nTaggedB = op.rng_count(jets[:2], partial((lambda j,fun=None,val=None : fun(j) >= val), fun=t_bDisc[tagger], val=thresh))
                                ## now two plots: two true b-jets -> n pass; two tagged: n true
                                sel_twotag_pt_twoTaggedB = sel_twotag_pt.refine(f"{sel_twotag_pt.name}_twoTaggedB_{key}", cut=(op.rng_len(bjets_pass) == 2))
                                pNTagged_TwoTrue = Plot.make1D(f"{channel}_Events_2l{nbrj}j_nbtagged_TruthFlav_2Beauty_{key}_{ptBin}",
                                        nTaggedB, sel_twotag_pt_twoTrueB,
                                        EqBin(3, 0., 3.), title="N b-tags for events with two true b's",
                                        plotopts=utils.getOpts(channel, **{"log-y": True}))
                                pNTrueB_TwoTagged = Plot.make1D(f"{channel}_Events_2l{nbrj}j_2btagged_TruthFlav_nBeauty_{key}_{ptBin}",
                                        nTrueB, sel_twotag_pt_twoTaggedB,
                                        EqBin(3, 0., 3.), title="N true b's for events with two btags",
                                        plotopts=utils.getOpts(channel, **{"log-y": True}))
                                pNnotTrueB_TwoTagged = Plot.make1D(f"{channel}_Events_2l{nbrj}j_2btagged_TruthFlav_notBeauty_{key}_{ptBin}",
                                        trueBC_cat, sel_twotag_pt_twoTaggedB,
                                        EqBin(4, 0., 4.), title="N b-tags (cross MC truth)",
                                        plotopts=utils.getOpts(channel, **{"log-y": True}))
                                pNTagged_CrossTruth = Plot.make1D(f"{channel}_Events_2l{nbrj}j_nBtagged_cross_TruthFlav_{key}_{ptBin}",
                                        nTaggedB*4+trueBC_cat, sel_twotag_pt_twoTaggedB,
                                        EqBin(12, 0., 12.), title="N b-tags (cross MC truth)",
                                        plotopts=utils.getOpts(channel, **{"log-y": True}))
                                plots += [ pNTagged_TwoTrue, pNTrueB_TwoTagged, pNnotTrueB_TwoTagged, pNTagged_CrossTruth]

                                plotsToSum_twoTagCount[(channel, nbrj, iPT, key, "TwoTrueB")].append(pNTagged_TwoTrue)
                                plotsToSum_twoTagCount[(channel, nbrj, iPT, key, "TwoTagged")].append(pNTrueB_TwoTagged)
                                plotsToSum_twoTagCount[(channel, nbrj, iPT, key, "TwoTaggedNotTrueB")].append(pNnotTrueB_TwoTagged)
                                plotsToSum_twoTagCount[(channel, nbrj, iPT, key, "NTaggedCrossTruth")].append(pNTagged_CrossTruth)
                                
                                #eff_b= op.sqrt(((nTrueB if not isMC else None) - trueBC_cat)/nTaggedB)
                                #btag_score = op.map(bjets_pass, t_bDisc[tagger]) 
                                #plots.append(Plot.make2D(f"{channel}_2l{nbrj}j_btagScore_vs_effb_{key}_{ptBin}",
                                #            (btag_score, eff_b), sel_twotag_pt_twoTaggedB,
                                #            (EqBin(50, thresh, 1.), EqBin(50, 0., 1.)), title=f"eff_b .vs. {key} discr {ptBin} [GeV]",
                                #            plotopts=utils.getOpts(channel, **{"log-y": True})))
        ## make summedPlots (over lepton flavour channels) of plotsToSum_twoTagCount
        for pkey, t2plots in plotsToSum_twoTagCount.items():
            channel, nbrj, iPT, key, pName = pkey
            ptBin = getPtLabel(iPT)
            if 'muel' in channel or 'elmu' in channel:
                plots.append(SummedPlot(f"{pName}_2lOSOF{nbrj}_{key}_{ptBin}", t2plots))
            else:
                plots.append(SummedPlot(f"{pName}_2lOSSF{nbrj}_{key}_{ptBin}", t2plots))

        if self.doYields:
            plots.extend(yield_object.returnPlots())

        return plots

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        # run plotIt as defined in HistogramsModule - this will also ensure that self.plotList is present
        super(TTbarDileptonMeasurment, self).postProcess(taskList, config, workdir, resultsdir)

        from bamboo.plots import CutFlowReport, DerivedPlot
        import bambooToOls

        plotstoNormalized = []
        for plots in self.plotList:
            if 'btagged_TruthFlav_' in plots.name:
                plotstoNormalized.append(plots)
        if not os.path.isdir(os.path.join(resultsdir, "normalizedFor2TagCount")):
            os.makedirs(os.path.join(resultsdir, "normalizedFor2TagCount"))
        normalizeAndMergeSamples(plotstoNormalized, self.readCounters, config, resultsdir, os.path.join(resultsdir, "normalizedFor2TagCount"))


        plotList_2D = [ ap for ap in self.plotList if ( isinstance(ap, Plot) or isinstance(ap, DerivedPlot) ) and len(ap.binnings) == 2 ]
        logger.debug("Found {0:d} plots to save".format(len(plotList_2D)))

        from bamboo.analysisutils import loadPlotIt
        p_config, samples, plots_2D, systematics, legend = loadPlotIt(config, plotList_2D, eras=self.args.eras, workdir=workdir, resultsdir=resultsdir, readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes, plotDefaults=self.plotDefaults)
        from plotit.plotit import Stack
        from bamboo.root import gbl
        for plot in plots_2D:
            logger.debug(f"Saving plot {plot.name}")
            obsStack = Stack(smp.getHist(plot) for smp in samples if smp.cfg.type == "DATA")
            expStack = Stack(smp.getHist(plot) for smp in samples if smp.cfg.type == "MC")
            cv = gbl.TCanvas(f"c{plot.name}")
            cv.Divide(2)
            cv.cd(1)
            expStack.obj.Draw("COLZ0")
            cv.cd(2)
            obsStack.obj.Draw("COLZ0")
            cv.Update()
            cv.SaveAs(os.path.join(resultsdir, f"{plot.name}.png"))

        Topsamples = ['TTJets_SL', 'TTJets_DL', 'TTJets_AH']
        for plotName in ['gen_ToppT_noReweighting', 'gen_ToppT_withReweighting']:
            #for smp in Topsamples:
            print (samples)
            expStack = Stack(smp.getHist(plotName)for smp in samples if smp.cfg.type == "MC")
            cv = gbl.TCanvas(f"c{plotName}")
            cv.cd()
            expStack.obj.Draw("COLZ0")
            cv.SaveAs(os.path.join(resultsdir, f"{plotName}.png"))

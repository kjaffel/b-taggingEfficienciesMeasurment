from bamboo.analysismodules import NanoAODHistoModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.scalefactors import binningVariables_nano

from bamboo import treefunctions as op
from bamboo import scalefactors
import logging
logger = logging.getLogger("ttbar Plotter")

from itertools import chain
import os.path
import collections
import math
import argparse
import sys
import os

btagPath = os.path.dirname(__file__)
if btagPath not in sys.path:
    sys.path.append(btagPath)

from ControlPlots import *
from utils import *
import bamboo.scalefactors

binningVariables = {
      "Eta"       : lambda obj : obj.eta
    , "ClusEta"   : lambda obj : obj.eta + obj.deltaEtaSC
    , "AbsEta"    : lambda obj : op.abs(obj.eta)
    , "AbsClusEta": lambda obj : op.abs(obj.eta + obj.deltaEtaSC)
    , "Pt"        : lambda obj : obj.pt
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
            "electron_reco": localizeSF("Electron_EGamma_SF2D_RECO_2017RunBCDEF_ptG20_fromPOG.json", "2017"),
            "electron_trigger": localizeSF("Electron_ele28_ht150_OR_ele32.json", "2017"),
            "muon_ID": {"cut_medium": localizeSF("Muon_NUM_MediumID_DEN_genTracks_pt_abseta_{uncer}_2017RunBCDEF.json".format(uncer=uncer), "2017") for uncer in ("syst", "stat")},
            "muon_iso":{"iso_tight_id_medium": localizeSF("Muon_NUM_TightRelIso_DEN_MediumID_pt_abseta_{uncer}_2017RunBCDEF.json".format(uncer=uncer), "2017") for uncer in ("syst", "stat")},
            "muon_trigger": localizeSF("Muon_IsoMu27_RunBtoF.json", "2017"),
            },
    "2018": {
            #***********  leptons ID , ISO , Scales Factors **************************
            "electron_ID": {"cut_medium": localizeSF("Electron_EGamma_SF2D_2018_Medium_Fall17V2.json", "2018")},
            "electron_reco": localizeSF("Electron_EGamma_SF2D_RECO_2018_fromPOG.json", "2018"),
            "electron_trigger": localize_trigger("Electron_ele28_ht150_OR_ele32.json", "2018"),
            "muon_ID": {"cut_medium": localizeSF("Muon_NUM_MediumID_DEN_TrackerMuons_pt_abseta_{uncer}_2018RunABCD.json".format(uncer=uncer), "2018") for uncer in ("syst", "stat")},
            "muon_iso":{"iso_tight_id_medium": localizeSF("Muon_NUM_TightRelIso_DEN_MediumID_pt_abseta_{uncer}_2018RunABCD.json".format(uncer=uncer), "2018") for uncer in ("syst", "stat")},
            "muon_trigger" :tuple(localize_trigger("{trig}_PtEtaBins_2018AfterMuonHLTUpdate.json".format(trig=trig),"2018")
                             for trig in ("IsoMu24_OR_IsoTkMu24","Mu50_OR_OldMu100_OR_TkMu100" )),
            },
            # double leptons trigger 
            "mueleLeg_HHMoriond17_2016" : tuple(localize_HHMoriond17Trigger("{wp}.json".format(wp=wp)) 
                                                for wp in ("Muon_XPathIsoMu23leg", "Muon_XPathIsoMu8leg", "Electron_IsoEle23Leg", "Electron_IsoEle12Leg")),
            "elemuLeg_HHMoriond17_2016" : tuple(localize_HHMoriond17Trigger("{wp}.json".format(wp=wp)) 
                                                for wp in ("Electron_IsoEle23Leg", "Electron_IsoEle12Leg", "Muon_XPathIsoMu23leg", "Muon_XPathIsoMu8leg")),
        }

def getScaleFactor(objType, key, periods=None, combine=None, additionalVariables=dict(), getFlavour=None, isElectron=False, systName=None):
    return bamboo.scalefactors.get_scalefactor(objType, key, periods=periods, combine=combine,
                                        additionalVariables=additionalVariables,
                                        sfLib=myScaleFactors,
                                        paramDefs=bamboo.scalefactors.binningVariables_nano,
                                        getFlavour=getFlavour,
                                        isElectron=isElectron,
                                        systName=systName)

def getL1PreFiringWeight(tree):
    return op.systematic(tree.L1PreFiringWeight_Nom,
                        name="L1PreFiring",
                        up=tree.L1PreFiringWeight_Up,
                        down=tree.L1PreFiringWeight_Dn)

lowerPtBinEdges=[30,50,70,100,140,200,300, 600]
def ReturnPtLabel(iPT):
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


def ReturnVarAtIndex(tagger, wp, bjets, idx):
    bjets_=safeget(bjets, tagger, wp)
    if tagger=="DeepCSV":
        op.sort(bjets_, lambda j: - j.btagDeepB)
        return bjets_[idx].btagDeepB
    elif tagger=="DeepFlavour":
        op.sort(bjets_, lambda j: - j.btagDeepFlavB)
        return bjets_[idx].btagDeepFlavB
    else:
        raise RuntimeError("Something went wrong in returning {0} discriminator !".format(tagger))

def TwoTagCount(uname, key, OsOflep, jets, rawbjets, jnotbtagged, _2l2nobjets, bjets, _2l2jsel, _2l2bjsel, WP, nbrj, isMC):
    from bambooToOls import Plot
    from bamboo.plots import SummedPlot
    from bamboo.plots import EquidistantBinning as EqBin
    from bamboo import treefunctions as op
    plots=[]
    for iPT in range(-1,len(lowerPtBinEdges)):
        ptBin= ReturnPtLabel(iPT)
        histosName = uname+"_{0}_twoTags_".format(nbrj) + key + "_"+ ptBin
        if ptBin=="Inclusive": 
            _2l2jsel_ptBin  = _2l2jsel
            _2l2bjsel_ptBin = _2l2bjsel
            _2l2nobjets_ptBin = _2l2nobjets
        else:
            if ptBin =="600toInf":
                minpt = 600.
                maxpt = 1000.
            else:
                minpt = lowerPtBinEdges[iPT]
                maxpt = lowerPtBinEdges[iPT+1]
                            
            _2l2bjsel_ptBin = ((_2l2bjsel.get(key)).refine("{0}_2bJets_{1}_{2}ptBin_cut{3}".format(uname, key, ptBin, nbrj), cut=[op.in_range(minpt, bjets[0].pt, maxpt), op.in_range(minpt, bjets[1].pt, maxpt)] )) 
            #_2l2bjsel_ptBin = ((_2l2bjsel.get(key)).refine("{0}_2bJets_{1}_{2}ptBin".format(uname, key, ptBin), cut=[op.in_range(minpt, jets[0].pt, maxpt), op.in_range(minpt, jets[1].pt, maxpt)] ))
            _2l2nobjets_ptBin = ((_2l2nobjets.get(key)).refine("{0}_2Jets_NOTbtagged_{1}_{2}ptBin_cut{3}".format(uname, key, ptBin, nbrj), cut=[op.in_range(minpt, jnotbtagged[0].pt, maxpt), op.in_range(minpt, jnotbtagged[1].pt, maxpt)] ))
            _2l2jsel_ptBin = (_2l2jsel_ptBin.refine("{0}_2Jets_{1}_{2}ptBin_cut{3}".format(uname, key, ptBin, nbrj), cut=[op.in_range(minpt, jets[0].pt, maxpt), op.in_range(minpt, jets[1].pt, maxpt)] )) 
            
            #tagger=key.replace(WP, "")
            #passing_bjets=safeget(rawbjets, tagger, WP)
                
            #_2l2bjsel_ptBin = {
            #    key :  _2l2jsel_ptBin.refine("{0}_2bJets_{1}_{2}ptBin".format(uname, key, ptBin),
            #                                                        cut=[ op.rng_len(passing_bjets) ==2 ]),
            #                }
        
        logger.info("Start filling histos ...")
        
        isLight  = lambda abshf : op.OR(abshf == 21, op.in_range(1, abshf, 3))
        isCharm  = lambda abshf : abshf == 4
        isBeauty = lambda abshf : abshf == 5
        TwoTagCrossFlavour = {
                    "2b"    : lambda f1, f2 : op.AND(isBeauty(f1), isBeauty(f2)),
                    "1b_1c" : lambda f1, f2 : op.OR(op.AND(isBeauty(f1), isCharm(f2)), op.AND(isBeauty(f2), isCharm(f1))),
                    "1b_1l" : lambda f1, f2 : op.OR(op.AND(isBeauty(f1), isLight(f2)), op.AND(isBeauty(f2), isLight(f1))),
                    "1c_1l" : lambda f1, f2 : op.OR(op.AND(isCharm(f1), isLight(f2)), op.AND(isCharm(f2), isLight(f1))),
                    "2c"    : lambda f1, f2 : op.AND(isCharm(f1), isCharm(f2)),
                    "2l"    : lambda f1, f2 : op.AND(isLight(f1), isLight(f2)),
                    }
        sel = (_2l2bjsel_ptBin.get(key) if ptBin == "Inclusive" else (_2l2bjsel_ptBin))
        nobjets_sel = (_2l2nobjets_ptBin.get(key) if ptBin == "Inclusive" else (_2l2nobjets_ptBin))
        #sel = _2l2bjsel_ptBin.get(key)
        twoTagCrossFlavour_categories = dict((catName, sel.refine("{0}_TwoTagCross_{1}_{2}_{3}_cut{4}".format(uname, catName, key, ptBin, nbrj), cut=catSel(op.abs(bjets[0].hadronFlavour), op.abs(bjets[1].hadronFlavour))if isMC else None )) for catName, catSel in TwoTagCrossFlavour.items())
            
        TwoBeautyWith012PassBtagWP = dict((catName, _2l2jsel_ptBin.refine("{0}_TwoBeauty_{1}_{2}_{3}_cut{4}".format(uname, catName, key, ptBin, nbrj), cut=catSel(op.abs(jets[0].hadronFlavour), op.abs(jets[1].hadronFlavour))if isMC else None )) for catName, catSel in TwoTagCrossFlavour.items())
        
        twobeauty_failsbtag = dict((catName, nobjets_sel.refine("{0}_2lep_2jets_failsbtagged_truthflav_{1}_{2}_{3}_cut{4}".format(uname, catName, key, ptBin, nbrj), cut=catSel(op.abs(jnotbtagged[0].hadronFlavour), op.abs(jnotbtagged[1].hadronFlavour))if isMC else None )) for catName, catSel in TwoTagCrossFlavour.items())

        # other mc (falv at least 1b) - 2btagged
        Twobtagged_Truth_BeautyCharm = Plot.make1D("{0}_Events_2l{1}j_2btagged_TruthFlav_BeautyCharm_{2}_{3}".format(uname, nbrj, key, ptBin),
                                                    op.c_int(0)if isMC else op.c_int(-1), twoTagCrossFlavour_categories.get("1b_1c"),
                                                    EqBin(4, -1., 3.),
                                                    title="N Fake %s b-tags (mc truth)"%WP,
                                                    plotopts=utils.getOpts(uname, **{"log-y": True}))

        Twobtagged_Truth_BeautyLight = Plot.make1D("{0}_Events_2l{1}j_2btagged_TruthFlav_BeautyLight_{2}_{3}".format(uname, nbrj, key, ptBin),
                                                    op.c_int(0)if isMC else op.c_int(-1), twoTagCrossFlavour_categories.get("1b_1l"),
                                                    EqBin(4, -1., 3.),
                                                    title="N Fake %s b-tags (mc truth)"%WP,
                                                    plotopts=utils.getOpts(uname, **{"log-y": True}))
        
        Twobtagged_Truth_CharmLight = Plot.make1D("{0}_Events_2l{1}j_2btagged_TruthFlav_CharmLight_{2}_{3}".format(uname, nbrj, key, ptBin),
                                                    op.c_int(0)if isMC else op.c_int(-1), twoTagCrossFlavour_categories.get("1c_1l"),
                                                    EqBin(4, -1., 3.),
                                                    title="N Fake %s b-tags (mc truth)"%WP,
                                                    plotopts=utils.getOpts(uname, **{"log-y": True}))

        Twobtagged_Truth_2Charm = Plot.make1D("{0}_Events_2l{1}j_2btagged_TruthFlav_2Charm_{2}_{3}".format(uname, nbrj, key, ptBin),
                                                    op.c_int(0)if isMC else op.c_int(-1), twoTagCrossFlavour_categories.get("2c"),
                                                    EqBin(4, -1., 3.),
                                                    title="N Fake %s b-tags (mc truth)"%WP,
                                                    plotopts=utils.getOpts(uname, **{"log-y": True}))
        
        Twobtagged_Truth_2Light = Plot.make1D("{0}_Events_2l{1}j_2btagged_TruthFlav_2Light_{2}_{3}".format(uname, nbrj, key, ptBin),
                                                    op.c_int(0)if isMC else op.c_int(-1), twoTagCrossFlavour_categories.get("2l"),
                                                    EqBin(4, -1., 3.),
                                                    title="N Fake %s b-tags (mc truth)"%WP,
                                                    plotopts=utils.getOpts(uname, **{"log-y": True}))
        # 2bmc truth flav 2b with (0 or 1 btagged ): means fails btagging requirements 
        NOTbtagged_Truth_2Beauty = Plot.make1D("{0}_Events_2l{1}j_NOTbtagged_TruthFlav_2Beauty_{2}_{3}".format(uname, nbrj, key, ptBin),
                                                    op.c_int(1)if isMC else op.c_int(-1), twobeauty_failsbtag.get("2b"),
                                                    EqBin(4, -1., 3.),
                                                    title="N Fake %s b-tags (mc truth)"%WP,
                                                    plotopts=utils.getOpts(uname, **{"log-y": True}))
        # 2bmc truth flav 2b - 2btagged 
        Twobtagged_Truth_2Beauty = Plot.make1D("{0}_Events_2l{1}j_2btagged_TruthFlav_2Beauty_{2}_{3}".format(uname, nbrj, key, ptBin),
                                                    op.c_int(2)if isMC else op.c_int(-1), twoTagCrossFlavour_categories.get("2b"),
                                                    EqBin(4, -1., 3.),
                                                    title="N Fake %s b-tags (mc truth)"%WP,
                                                    plotopts=utils.getOpts(uname, **{"log-y": True}))
        
        #plots.append( Twobtagged_Truth_BeautyCharm)
        #plots.append( Twobtagged_Truth_BeautyLight)
        #plots.append( Twobtagged_Truth_CharmLight)
        #plots.append( Twobtagged_Truth_2Charm)
        #plots.append( Twobtagged_Truth_2Light)

        plots.append(SummedPlot("{0}_{1}_twoTags_{2}_{3}".format(uname, nbrj, key, ptBin),
                    [Twobtagged_Truth_BeautyCharm, Twobtagged_Truth_BeautyLight, Twobtagged_Truth_CharmLight, Twobtagged_Truth_CharmLight, Twobtagged_Truth_2Light, 
                     Twobtagged_Truth_2Beauty, NOTbtagged_Truth_2Beauty],
                    title="N Fake %s b-tags (mc truth)"%WP,
                    plotopts=utils.getOpts(uname, **{"log-y": True})))
    return plots 


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
    def addArgs(self, parser):
        super(TTbarDileptonMeasurment, self).addArgs(parser)
        parser.add_argument("--backend", type=str, default="dataframe", help="Backend to use, 'dataframe' (default) or 'lazy'")
        parser.add_argument("-s", "--systematic", action="store_true", help="Produce systematic variations")

    def prepareTree(self, tree, sample=None, sampleCfg=None):
        era = sampleCfg.get("era") if sampleCfg else None
        isMC = self.isMC(sample)
        metName = "METFixEE2017" if era == "2017" else "MET"
        ## initializes tree.Jet.calc so should be called first (better: use super() instead)
        # JEC's Recommendation for Full RunII: https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        # JER : -----------------------------: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
       
        from bamboo.treedecorators import NanoAODDescription, nanoRochesterCalc, nanoJetMETCalc, nanoJetMETCalc_METFixEE2017

        tree,noSel,be,lumiArgs = NanoAODHistoModule.prepareTree(self, tree, 
                                                                    sample=sample, 
                                                                    sampleCfg=sampleCfg, 
                                                                    description=NanoAODDescription.get("v5", year=(era if era else "2016"), 
                                                                    isMC=isMC, 
                                                                    systVariations=[ nanoRochesterCalc, (nanoJetMETCalc_METFixEE2017 if era == "2017" else nanoJetMETCalc) ]),
                                                                    lazyBackend   = (self.args.backend == "lazy")) ## will do Jet and MET variations, and the Rochester correction
        triggersPerPrimaryDataset = {}
        jec, smear, jesUncertaintySources = None, None, None

        from bamboo.analysisutils import configureJets, configureType1MET, configureRochesterCorrection
        isNotWorker = (self.args.distributed != "worker") 
        

        if era == "2018":
            configureRochesterCorrection(tree._Muon, os.path.join(os.path.dirname(__file__), "data", "RoccoR2018.txt"), isMC=isMC, backend=be, uName=sample)

            triggersPerPrimaryDataset = {
                "MuonEG"     : [ tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
                                 tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
                                 
                                 tree.HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
                                 tree.HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,

                                 tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
                                 tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
                                 
                                 tree.HLT.Mu27_Ele37_CaloIdL_MW, 
                                 tree.HLT.Mu37_Ele27_CaloIdL_MW],
                "EGamma":[ tree.HLT.Ele32_WPTight_Gsf  ],
                # OldMu100 and TkMu100 are recommend to recover inefficiencies at high pt 
                # here: (https://indico.cern.ch/event/766895/contributions/3184188/attachments/1739394/2814214/IdTrigEff_HighPtMu_Min_20181023_v2.pdf)
                "SingleMuon": [ tree.HLT.IsoMu24, 
                                tree.HLT.IsoMu27, 
                                tree.HLT.Mu50, 
                                tree.HLT.OldMu100, 
                                tree.HLT.TkMu100 ], 
                }

            if self.isMC(sample):
                jec="Autumn18_V8_MC"
                smear="Autumn18_V1_MC"
                jesUncertaintySources=["Total"]

            else:
                if "2018A" in sample:
                    jec="Autumn18_RunA_V8_DATA"

                elif "2018B" in sample:
                    jec="Autumn18_RunB_V8_DATA"

                elif "2018C" in sample:
                    jec="Autumn18_RunC_V8_DATA"

                elif "2018D" in sample:
                    jec="Autumn18_RunD_V8_DATA"

        elif era == "2017":
            
            configureRochesterCorrection(tree._Muon, os.path.join(os.path.dirname(__file__), "data", "RoccoR2017.txt"), isMC=isMC, backend=be, uName=sample)
            
            # https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2017
            triggersPerPrimaryDataset = {
                "MuonEG"     : [ tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
                                 tree.HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
                                 tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ ],
                
                "SingleElectron": [ tree.HLT.Ele35_WPTight_Gsf,
                                    tree.HLT.Ele28_eta2p1_WPTight_Gsf_HT150 ],
                "SingleMuon" :    [ tree.HLT.IsoMu27,
                                    tree.HLT.Mu50   ],

                
            }
            
            if "2017B" not in sample:
             ## all are removed for 2017 era B
                triggersPerPrimaryDataset["MuonEG"] += [ 
                        tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
                        tree.HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
                        tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL ]

            if self.isMC(sample):
                jec="Fall17_17Nov2017_V32_MC"
                smear="Fall17_V3_MC"
                jesUncertaintySources=["Total"]

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
        
        
        if self.isMC(sample):
            noSel = noSel.refine("genWeight", weight=tree.genWeight, cut=op.OR(*chain.from_iterable(triggersPerPrimaryDataset.values())), autoSyst=self.doSysts)
            if self.doSysts:
                logger.info("Adding FSR, ISR and QCD scale variations")
                noSel = utils.addTheorySystematics(self, tree, noSel)
        else:
            noSel = noSel.refine("withTrig", cut=makeMultiPrimaryDatasetTriggerSelection(sample, triggersPerPrimaryDataset) )
            
        return tree,noSel,be,lumiArgs
    
    def definePlots(self, t, noSel, sample=None, sampleCfg=None):    

        from bamboo.analysisutils import forceDefine
        from bamboo.plots import EquidistantBinning as EqBin
        from bambooToOls import Plot
        from bamboo.plots import SummedPlot
        from bamboo import treefunctions as op
        from METFilter_xyCorr import METFilter, METcorrection
        from bamboo.analysisutils import makePileupWeight

        isMC = self.isMC(sample)
        era = sampleCfg.get("era") if sampleCfg else None
        noSel = noSel.refine("passMETFlags", cut=METFilter(t.Flag, era, isMC) )
        puWeightsFile = None
        mcprofile= None
        
        if era == "2018":
            suffix= '2018_Autumn18'
            puWeightsFile = os.path.join(os.path.dirname(__file__), "data/PileupFullRunII/", "puweights2018_Autumn18.json")
        elif era == "2017":
            suffix ='2017_Fall17'
            if sample in ['DYJetsToLL_M-10to50-madgraphMLM', 'ST_tW_antitop_5f', 'WZ', 'ZZ']: 
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

        plots = []
        forceDefine(t._Muon.calcProd, noSel)

        sorted_muons = op.sort(t.Muon, lambda mu : -mu.pt)
        muons = op.select(sorted_muons, lambda mu : op.AND(mu.pt > 10., op.abs(mu.eta) < 2.4, mu.mediumId, mu.pfRelIso04_all<0.15, op.abs(mu.sip3d) < 4.))
                     
        muonIDSF = getScaleFactor("lepton", (era, "muon_ID", "cut_medium"), systName="muon_ID")
        muonIsoSF = getScaleFactor("lepton", (era, "muon_iso", "iso_tight_id_medium"), systName="muon_iso")
        #     if combPrefix == "":
        #UnboundLocalError: local variable 'combPrefix' referenced before assignment
        #FIXME muonTriggerSF = getScaleFactor("lepton", (era, "muon_trigger"), systName="muon_trigger")

        sorted_electrons = op.sort(t.Electron, lambda ele : -ele.pt)
        electrons = op.select(sorted_electrons, lambda ele : op.AND(ele.pt > 15., op.abs(ele.eta) < 2.5 , ele.cutBased>=3, op.abs(ele.sip3d) < 4., 
                                                                    op.OR(op.AND(op.abs(ele.dxy) < 0.05, op.abs(ele.dz) < 0.1), 
                                                                          op.AND(op.abs(ele.dxy) < 0.05, op.abs(ele.dz) < 0.2) ))) 
        
        eleIDSF = getScaleFactor("lepton", (era, "electron_ID", "cut_medium"), isElectron=True, systName="ele_ID")
        eleRecoSF = getScaleFactor("lepton", (era, "electron_reco"), isElectron=True, systName="ele_reco")
        #FIXME eleTriggerSF = getScaleFactor("lepton", (era, "electron_trigger"), isElectron=True, systName="ele_trigger")

        sorted_jets=op.sort(t.Jet, lambda j : -j.pt)
        jetsSel = op.select(sorted_jets, lambda j : op.AND(j.pt > 30., op.abs(j.eta)< 2.5, j.jetId &4))
                                                            #,op.AND( op.in_range(30., j.pt, 50.), j.puId==7 ))) #j.jetId &2)) old cut 
        
        # exclude from the jetsSel any jet that happens to include within its reconstruction cone a muon or an electron.
        jets= op.select(jetsSel, lambda j : op.AND(op.NOT(op.rng_any(electrons, lambda ele : op.deltaR(j.p4, ele.p4) < 0.4 )), 
                                                   op.NOT(op.rng_any(muons, lambda mu : op.deltaR(j.p4, mu.p4) < 0.4 ))))

        #jetpuId = getScaleFactor("jet", ("JetId_InHighPileup_{0}_94X".format(era).replace("94X", "102X" if era=="2018" else "94X"), "puId_T"), systName="JetpuID_T")
        # order jets 
        cleanedJetsByDeepFlav = op.sort(jets, lambda j: -j.btagDeepFlavB)
        cleanedJetsByDeepB = op.sort(jets, lambda j: -j.btagDeepB)
        
        btaggingWPs = {
                "DeepCSV":{ # era: (loose, medium, tight)
                            #"2016": (0.2217, 0.6321, 0.8953), 
                            "2017":(0.1522, 0.4941, 0.8001), 
                            "2018":(0.1241, 0.4184, 0.7527) },
                "DeepFlavour":{
                            #"2016":(0.0614, 0.3093, 0.7221), 
                            "2017":(0.0521, 0.3033, 0.7489), 
                            "2018": (0.0494, 0.2770, 0.7264) }
                }
        
        bjets = {}
        failings = {}
        # bjets ={ "DeepFlavour": {"L": jets pass loose  , "M":  jets pass medium  , "T":jets pass tight    }     
        #           "DeepCSV":    {"L":    ---           , "M":         ---        , "T":   ----            }
        #        }
        WorkingPoints = ["L", "M", "T"]
        for tagger  in btaggingWPs.keys():
            
            bJets_deepflavour ={}
            failing_deepflavour = {}
            bJets_deepcsv ={}
            failing_deepcsv = {}
            
            for wp in sorted(WorkingPoints):
                
                suffix = ("loose" if wp=='L' else ("medium" if wp=='M' else "tight"))
                idx = ( 0 if wp=="L" else ( 1 if wp=="M" else 2))
                
                if tagger=="DeepFlavour":
                    
                    print ("Btagging: Era= {0}, Tagger={1}, Pass_{2}_working_point={3}".format(era, tagger, suffix, btaggingWPs[tagger][era][idx] ))
                    print ("btag_{0}_94X".format(era).replace("94X", "102X" if era=="2018" else "94X"), "{0}_{1}".format('DeepJet', suffix))
                    
                    bJets_deepflavour[wp] = op.select(cleanedJetsByDeepFlav, lambda j : j.btagDeepFlavB >= btaggingWPs[tagger][era][idx] )
                    failing_deepflavour[wp] = op.select(cleanedJetsByDeepFlav, lambda j : j.btagDeepFlavB < btaggingWPs[tagger][era][idx] )
                    Jet_DeepFlavourBDisc = { "BTagDiscri": lambda j : j.btagDeepFlavB }
                    
                    bjets[tagger]=bJets_deepflavour
                    failings[tagger] = failing_deepflavour
                    
                else:
                    print ("Btagging: Era= {0}, Tagger={1}, Pass_{2}_working_point={3}".format(era, tagger, suffix, btaggingWPs[tagger][era][idx] ))
                    print ("btag_{0}_94X".format(era).replace("94X", "102X" if era=="2018" else "94X"), "{0}_{1}".format('DeepJet', suffix))
                    
                    bJets_deepcsv[wp] = op.select(cleanedJetsByDeepB, lambda j : j.btagDeepB >= btaggingWPs[tagger][era][idx] )   
                    failing_deepcsv[wp] = op.select(cleanedJetsByDeepB, lambda j : j.btagDeepB < btaggingWPs[tagger][era][idx] )
                    Jet_DeepCSVBDis = { "BTagDiscri": lambda j : j.btagDeepB }
                    
                    bjets[tagger]=bJets_deepcsv
                    failings[tagger] = failing_deepcsv
        

        MET = t.MET if era != "2017" else t.METFixEE2017
        corrMET=METcorrection(MET,t.PV,sample,era,self.isMC(sample))
        
        L1Prefiring = 1.
        HLTZvtxSF = 1.
        if era in ["2016", "2017"]:
            L1Prefiring = getL1PreFiringWeight(t)
        if era =='2017':
            HLTZvtxSF = op.systematic(op.c_float(0.991), name='HLT_Zvtx_eff', up=op.c_float(0.992), down=op.c_float(0.990)) 
        
        #The contribution from Z + jets events is reduced by applying a veto around the Z boson mass when the two leptons have the same flavour, ( | Mll -MZ | > 10 GeV).
        osdilep = lambda l1,l2 : op.AND(l1.charge != l2.charge, op.invariant_mass(l1.p4, l2.p4) >12. )
        Zboson_Vetocut = lambda l1, l2 : op.AND(l1.charge != l2.charge, op.abs(op.invariant_mass(l1.p4, l2.p4) - 90.) > 10. )
        osdilep_Z = lambda l1,l2 : op.AND(l1.charge != l2.charge, op.in_range(70., op.invariant_mass(l1.p4, l2.p4), 120.))

        osll = {
                "elmu": op.combine((electrons, muons), pred=lambda ele,mu : op.AND(osdilep(ele,mu), ele.pt > mu.pt)),
                "muel": op.combine((muons, electrons), pred=lambda mu,ele : op.AND(osdilep(mu,ele), mu.pt > ele.pt)),

                #"mumu" : op.combine(muons, N=2, pred= osdilep_Z),
                #"elel" : op.combine(electrons, N=2, pred=osdilep_Z)
             }
        
        elemuTrigSF = getScaleFactor("dilepton", ("elemuLeg_HHMoriond17_2016"), systName="elmutrig")
        mueleTrigSF = getScaleFactor("dilepton", ("mueleLeg_HHMoriond17_2016"), systName="mueltrig")
        #doubleMuTrigSF = getScaleFactor("dilepton", ("doubleMuLeg_HHMoriond17_2016"), systName="mumutrig")
        #doubleEleTrigSF = getScaleFactor("dilepton", ("doubleEleLeg_HHMoriond17_2016"), systName="eleltrig")

        osllSFs = {
                "elmu" : (lambda lep :  [ HLTZvtxSF, L1Prefiring, muonIDSF(lep[1]), muonIsoSF(lep[1]) , eleRecoSF(lep[0]), eleIDSF(lep[0]), elemuTrigSF(lep)]),#, eleTriggerSF(lep[0]), muonTriggerSF(lep[1]) ]),
                "muel" : (lambda lep :  [ HLTZvtxSF, L1Prefiring, muonIDSF(lep[0]), muonIsoSF(lep[0]) , eleRecoSF(lep[1]), eleIDSF(lep[1]), mueleTrigSF(lep)])#,  eleTriggerSF(lep[1]), muonTriggerSF(lep[0])])
                
                #"mumu" : (lambda ll : [ muonIDSF(ll[0]), muonIDSF(ll[1]), muonIsoSF(ll[0]), muonIsoSF(ll[1]), doubleMuTrigSF(ll)]),
                #"elel" : (lambda ll : [ eleIDSF(ll[0]), eleIDSF(ll[1]), eleRecoSF(ll[0]), eleRecoSF(ll[1]), doubleEleTrigSF(ll)])
            }

        osllcatrng = lambda catrng : op.AND(op.rng_len(catrng) > 0, catrng[0][0].pt > 25.)
        hasOSLL = noSel.refine("hasOSLL", cut=op.OR(*( osllcatrng(rng) for rng in osll.values())))
        
        forceDefine(t._Jet.calcProd, hasOSLL)
        forceDefine(getattr(t, "_{0}".format("MET" if era != "2017" else "METFixEE2017")).calcProd, hasOSLL)


        categories = dict((channel, (leadpair[0],
                            hasOSLL.refine("{0}_EventSelection".format(channel), cut=osllcatrng(leadpair), weight=(osllSFs[channel](leadpair[0]) if isMC else None))
                            )) for channel, leadpair in osll.items())

        make_PrimaryANDSecondaryVerticesPlots = False
        make_LeptonsPlots_LepSel = True
        make_LeptonsPlots_LeptonPlusJetsSel = True
        make_JetsPlots = True
        make_METPlots = False
        make_bJetsPlots = True
        make_btagScore = True
        make_TwoTagCountPlots = True
        make_ttbarEstimationPlots =True

        for channel, (leptons, cat) in categories.items():

            if make_PrimaryANDSecondaryVerticesPlots:
                plots.extend(makePrimaryANDSecondaryVerticesPlots(cat, channel))
            if make_LeptonsPlots_LepSel:
                plots.extend(makeLeptonPlots(cat, leptons, '_2lepOSOFSel_', channel))

            plots.extend(makeJetmultiplictyPlots(cat, jets, '_NOcutOnJetsLen_', channel))
            
            cutOnJetsLen = {
                    #'Only2': op.rng_len(jets) ==2,
                    'atleast2': op.rng_len(jets) >1
                    }
            
            for nbrj, jcut in cutOnJetsLen.items():
            
                TwoLeptonsTwoJets=cat.refine("{0}Jets{1}Sel".format( nbrj, channel), cut=[ jcut])
                plots.extend(makeJetmultiplictyPlots(TwoLeptonsTwoJets, jets, '{0}Jets_'.format(nbrj), channel))
           
                if make_LeptonsPlots_LeptonPlusJetsSel:
                    plots.extend(makeLeptonPlots(TwoLeptonsTwoJets, leptons, '_2lep{0}jeSel_'.format(nbrj), channel))
                if make_JetsPlots:
                    plots.extend(makeJetPlots(TwoLeptonsTwoJets, jets, nbrj, channel))
                if make_METPlots:     
                    plots.extend(makeMETPlots(TwoLeptonsTwoJets, leptons, MET, corrMET, nbrj, channel))
                    plots.extend(MakeEXTRAMETPlots(TwoLeptonsTwoJets, corrMET, MET, nbrj, channel))
            
                for WP in WorkingPoints:
                    bJets_PassdeepflavourWP=safeget(bjets, "DeepFlavour", WP)
                    bJets_PassdeepcsvWP=safeget(bjets, "DeepCSV", WP)
                    
                    bJets_NOTPassdeepflavourWP=safeget(failings, "DeepFlavour", WP)
                    bJets_NOTPassdeepcsvWP=safeget(failings, "DeepCSV", WP)
                
                    cutOnbJetsLen_PassWP = {
                            'Only2': {
                                'DeepFlavour':op.rng_len(bJets_PassdeepflavourWP) == 2,
                                'DeepCSV': op.rng_len(bJets_PassdeepcsvWP) ==2
                                },
                            'atleast2':{ 
                                'DeepFlavour':op.rng_len(bJets_PassdeepflavourWP) >1,
                                'DeepCSV':op.rng_len(bJets_PassdeepcsvWP) >1
                                }
                            }
                    cutOnbJetsLen_FailWP = {
                            'Only2': {
                                'DeepFlavour':op.rng_len(bJets_NOTPassdeepflavourWP) ==2,
                                'DeepCSV':op.rng_len(bJets_NOTPassdeepcsvWP) ==2
                                },
                            'atleast2':{
                                'DeepFlavour':op.rng_len(bJets_NOTPassdeepflavourWP) >1,
                                'DeepCSV':op.rng_len(bJets_NOTPassdeepcsvWP) >1
                                }
                            }

                    TwoLeptonsTwoBJets = {
                        "DeepFlavour{0}".format(WP) :  TwoLeptonsTwoJets.refine("TwoLeptons{0}BJets_DeepFlavour{1}_{2}".format(nbrj, WP, channel),
                                                                            cut=[ cutOnbJetsLen_PassWP[nbrj]['DeepFlavour']]),
                        "DeepCSV{0}".format(WP)     :  TwoLeptonsTwoJets.refine("TwoLeptonsT{0}BJets_DeepCSV{1}_{2}".format(nbrj, WP, channel), 
                                                                            cut=[ cutOnbJetsLen_PassWP[nbrj]['DeepCSV']])
                                        }
                    
                    TwoLeptonsTwoFailingBJets = {
                        "DeepFlavour{0}".format(WP) :  TwoLeptonsTwoJets.refine("TwoLeptons{0}BJets_fail_DeepFlavour{1}_{2}".format(nbrj, WP, channel),
                                                                            cut=[ cutOnbJetsLen_FailWP[nbrj]['DeepFlavour']]),
                        "DeepCSV{0}".format(WP)     :  TwoLeptonsTwoJets.refine("TwoLeptons{0}BJets_fail_DeepCSV{1}_{2}".format(nbrj, WP, channel), 
                                                                            cut=[ cutOnbJetsLen_FailWP[nbrj]['DeepCSV']])
                                        }

                    if make_bJetsPlots:
                        plots.extend(makeBJetPlots(TwoLeptonsTwoBJets, bjets, WP, nbrj, channel))
                    if make_btagScore:
                        plots.extend(makeDiscriminatorPlots(TwoLeptonsTwoBJets, bjets, WP, btaggingWPs, nbrj, channel, era))
            
                    TwoLeptonsTwoBjets_METcut = dict((key, selNoMET.refine("TwoLeptons{0}Bjets_{1}_{2}_plusmetcut".format(nbrj, key, channel), 
                                                            cut=[ corrMET.pt < 80. ])) for key, selNoMET in TwoLeptonsTwoBJets.items()) 
                    TwoLeptonsTwoBjets_Inverted_METcut = dict((key, selNoMET.refine("TwoLeptons{0}Bjets_{1}_{2}_highMETtail".format(nbrj, key, channel), 
                                                                    cut=[ corrMET.pt > 80. ])) for key, selNoMET in TwoLeptonsTwoBJets.items()) 
                    
                    if make_ttbarEstimationPlots:
                        for sel, metcut in zip ([TwoLeptonsTwoBjets_METcut, TwoLeptonsTwoBjets_Inverted_METcut] , [ 'METcut', 'HighMET']):
                            plots.extend(makehistosforTTbarEstimation(sel, leptons, bjets, WP, metcut, nbrj, channel))        
                    
                    if make_TwoTagCountPlots:
                        # TwoTag Count Method --> start filling 2tag crossFlavour histograms
                        for key in TwoLeptonsTwoBJets.keys():
                            tagger=key.replace(WP, "")
                            bjets_=safeget(bjets, tagger, WP)
                            failings_=safeget(failings, tagger, WP)
                
                            plots +=(TwoTagCount(channel, key, cat, jets, bjets, failings_, TwoLeptonsTwoFailingBJets, bjets_, TwoLeptonsTwoJets, TwoLeptonsTwoBJets, WP, nbrj, isMC))

        return plots

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        # run plotIt as defined in HistogramsModule - this will also ensure that self.plotList is present
        super(TTbarDileptonMeasurment, self).postProcess(taskList, config, workdir, resultsdir)

        from bamboo.plots import CutFlowReport, DerivedPlot
        import bambooToOls

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

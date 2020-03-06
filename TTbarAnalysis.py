from bamboo.analysismodules import NanoAODHistoModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.scalefactors import binningVariables_nano

from bamboo import treefunctions as op
from bamboo import scalefactors

from bamboo.logging import getLogger
logger = getLogger(__name__)


from itertools import chain
import os.path
import collections
import math
import argparse
import sys
import os

sys.path.append('/home/ucl/cp3/kjaffel/bamboodev/b-taggingEfficienciesMeasurment')
import utils

binningVariables = {
      "Eta"       : lambda obj : obj.eta
    , "ClusEta"   : lambda obj : obj.eta + obj.deltaEtaSC
    , "AbsEta"    : lambda obj : op.abs(obj.eta)
    , "AbsClusEta": lambda obj : op.abs(obj.eta + obj.deltaEtaSC)
    , "Pt"        : lambda obj : obj.pt
    }

def localizeSF(aPath, era="2018"):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "ScaleFactors", era, aPath)

def localize_trigger(aPath):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "TriggerEfficienciesStudies", aPath)

myScaleFactors = {
    "2018": {
            #***********  leptons ID , ISO , Scales Factors **************************
            "electron_ID": {"cut_medium": localizeSF("Electron_EGamma_SF2D_2018_Medium_Fall17V2.json", "2018")},
            "electron_reco": localizeSF("Electron_EGamma_SF2D_Reco.json", "2018"),
            "muon_ID": {"cut_medium": localizeSF("Muon_NUM_MediumID_DEN_TrackerMuons_pt_abseta_{uncer}_2018RunABCD.json".format(uncer=uncer), "2018") for uncer in ("syst", "stat")},
            "muon_iso":{"iso_tight_id_medium": localizeSF("Muon_NUM_TightRelIso_DEN_MediumID_pt_abseta_{uncer}_2018RunABCD.json".format(uncer=uncer), "2018") for uncer in ("syst", "stat")},
            },
            #************ Trigger Scales Factors **********************************
            # single leptons trigger  
            "electron_trigger": localize_trigger("Electron_ele28_ht150_OR_ele32_etaCut.json"),
            #"electron_trigger": localize_trigger("Electron_ele28_ht150_OR_ele32.json", "2018"),
            "muon_trigger" :tuple(localize_trigger("{trig}_PtEtaBins_2018AfterMuonHLTUpdate.json".format(trig=trig))
                             for trig in ("IsoMu24_OR_IsoTkMu24","Mu50_OR_OldMu100_OR_TkMu100" )),
            #"muon_trigger": [(["Run315264to316360"], localize_trigger("Muon_IsoMu24_BeforeMuonHLTUpdate.json")),
                             #(["Run316361to325175"], localize_trigger("Muon_IsoMu24_AfterMuonHLTUpdate.json"))],
            # FIXME these are old versions should be updated ! can be passed to all eras
            # double leptons trigger 
            "mueleLeg_HHMoriond17_2016" : tuple(localize_trigger("{wp}.json".format(wp=wp)) 
                                          for wp in ("Muon_XPathIsoMu23leg", "Muon_XPathIsoMu8leg", "Electron_IsoEle23Leg", "Electron_IsoEle12Leg")),
            "elemuLeg_HHMoriond17_2016" : tuple(localize_trigger("{wp}.json".format(wp=wp)) 
                                          for wp in ("Electron_IsoEle23Leg", "Electron_IsoEle12Leg", "Muon_XPathIsoMu23leg", "Muon_XPathIsoMu8leg")),
    }

    # fill in some defaults: myScalefactors and bamboo.scalefactors.binningVariables_nano
def getScaleFactor(objType, key, periods=None, combine=None, isElectron=False, systName=None):
    return scalefactors.get_scalefactor(objType, key, combine=combine,
                                        sfLib=myScaleFactors, paramDefs=scalefactors.binningVariables_nano,
                                        isElectron=isElectron, systName=systName)

def safeget(dct, *keys):
    for key in keys:
        try:
            dct = dct[key]
        except KeyError:
            return None
    return dct


def METFilter(flags, era):
    # from https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    if era == '2018':
        cuts = [
                flags.goodVertices,
                flags.globalSuperTightHalo2016Filter, # not tested need to be careful
                flags.HBHENoiseFilter,
                flags.HBHENoiseIsoFilter,
                flags.EcalDeadCellTriggerPrimitiveFilter,
                flags.BadPFMuonFilter,
                flags.ecalBadCalibFilterV2 ]
    return cuts

class METcorrection(object):
    # https://lathomas.web.cern.ch/lathomas/METStuff/XYCorrections/XYMETCorrection.h
    def __init__(self,rawMET,pv,sample,era,isMC):
        if era == "2018":
            if isMC:
                xcorr = (-0.296713,  0.141506)
                ycorr = (-0.115685, -0.0128193)
            else:
                if '2018A' in sample:
                    xcorr= (-0.362865,  1.94505)
                    ycorr= (-0.0709085, 0.307365)
                elif'2018B' in sample:
                    xcorr = (-0.492083, 2.93552)
                    ycorr = (-0.17874,  0.786844)
                elif '2018C' in sample:
                    xcorr = (-0.521349, 1.44544)
                    ycorr = (-0.118956, 1.96434)
                else:
                    xcorr = (-0.531151,  1.37568)
                    ycorr = (-0.0884639, 1.57089)
            
            METxcorr=xcorr[0] *pv.npvs+xcorr[1]
            METycorr=ycorr[0] *pv.npvs+ycorr[1]
            
            corrMETx=rawMET.pt*op.cos(rawMET.phi) +METxcorr
            corrMETy=rawMET.pt*op.sin(rawMET.phi) +METycorr
                
            self.pt=op.sqrt(corrMETx**2 +corrMETy**2)
            atan=op.atan(corrMETy/corrMETx)
            self.phi=op.multiSwitch((corrMETx> 0,atan),(corrMETy> 0,atan+math.pi),atan-math.pi)

class TTbarDileptonMeasurment(NanoAODHistoModule):

    def __init__(self, args):
        super(TTbarDileptonMeasurment, self).__init__(args)
        self.plotDefaults = {
                            "y-axis"           : "Events",
                            "log-y"            : "both",
                            "y-axis-show-zero" : True,
                            "save-extensions"  : ["pdf"],
                            "show-ratio"       : True,
                            "sort-by-yields"   : False,
                            }

    def prepareTree(self, tree, sample=None, sampleCfg=None):
        era = sampleCfg.get("era") if sampleCfg else None
        isMC = self.isMC(sample)
        metName = "METFixEE2017" if era == "2017" else "MET"
        ## initializes tree.Jet.calc so should be called first (better: use super() instead)
        # JEC's Recommendation for Full RunII: https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        # JER : -----------------------------: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
        
        tree,noSel,be,lumiArgs = NanoAODHistoModule.prepareTree(self, tree, sample=sample, sampleCfg=sampleCfg, calcToAdd=["nJet", metName, "nMuon"])
        triggersPerPrimaryDataset = {}
        jec, smear, jesUncertaintySources = None, None, None

        from bamboo.analysisutils import configureJets, configureType1MET, configureRochesterCorrection
        isNotWorker = (self.args.distributed != "worker") 
        

        if era == "2018":
            configureRochesterCorrection(tree._Muon, os.path.join(os.path.dirname(__file__), "data", "RoccoR2018.txt"), isMC=isMC, backend=be, uName=sample)


            ## TODO add single muon single el trigger
            triggersPerPrimaryDataset = {
                "MuonEG"     : [ tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
                                 tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
                                 
                                 tree.HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
                                 tree.HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,

                                 tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
                                 tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
                                 
                                 tree.HLT.Mu27_Ele37_CaloIdL_MW, 
                                 tree.HLT.Mu37_Ele27_CaloIdL_MW],
                "SingleElectron":[ tree.HLT.Ele32_WPTight_Gsf  ],
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
            noSel = noSel.refine("genWeight", weight=tree.genWeight, cut=op.OR(*chain.from_iterable(triggersPerPrimaryDataset.values())))
        else:
            noSel = noSel.refine("withTrig", cut=makeMultiPrimaryDatasetTriggerSelection(sample, triggersPerPrimaryDataset) )
            
        return tree,noSel,be,lumiArgs
    
    def definePlots(self, t, noSel, sample=None, sampleCfg=None):    

        from bamboo.analysisutils import forceDefine
        from bamboo.plots import EquidistantBinning as EqBin
        from bamboo.plots import Plot, SummedPlot
        from bamboo import treefunctions as op

        era = sampleCfg.get("era") if sampleCfg else None
        noSel = noSel.refine("passMETFlags", cut=METFilter(t.Flag, era) )
        puWeightsFile = None
        
        if era == "2018":
            puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "puweights2018.json")
        
        if self.isMC(sample) and puWeightsFile is not None:
            from bamboo.analysisutils import makePileupWeight
            noSel = noSel.refine("puWeight", weight=makePileupWeight(puWeightsFile, t.Pileup_nTrueInt, systName="pileup"))
        
        isMC = self.isMC(sample)
        plots = []
        forceDefine(t._Muon.calcProd, noSel)

        sorted_muons = op.sort(t.Muon, lambda mu : -mu.pt)
        muons = op.select(sorted_muons, lambda mu : op.AND(mu.pt > 20., op.abs(mu.eta) < 2.4, mu.mediumId, mu.pfRelIso04_all<0.15, op.abs(mu.sip3d) < 4.))
                     
        muonIDSF = getScaleFactor("lepton", (era, "muon_ID", "cut_medium"), systName="muon_ID")
        muonIsoSF = getScaleFactor("lepton", (era, "muon_iso", "iso_tight_id_medium"), systName="muon_iso")
        # single muon trigger
        #muonTriggerSF = getScaleFactor("lepton", (era, "muon_trigger"), systName="muon_trigger")

        sorted_electrons = op.sort(t.Electron, lambda ele : -ele.pt)
        electrons = op.select(sorted_electrons, lambda ele : op.AND(ele.pt > 20., op.abs(ele.eta) < 2.5 , ele.cutBased>=3, op.abs(ele.sip3d) < 4., op.OR(op.AND(op.abs(ele.dxy) < 0.05, op.abs(ele.dz) < 0.1), op.AND(op.abs(ele.dxy) < 0.05, op.abs(ele.dz) < 0.2) ))) 
        
        eleIDSF = getScaleFactor("lepton", (era, "electron_ID", "cut_medium"), isElectron=True, systName="ele_ID")
        eleRecoSF = getScaleFactor("lepton", (era, "electron_reco"), isElectron=True, systName="ele_reco")
        # single electron trigger 
        #eleTriggerSF = getScaleFactor("lepton", (era, "electron_trigger"), isElectron=True, systName="ele_trigger")

        sorted_jets=op.sort(t.Jet, lambda j : -j.pt)
        jetsSel = op.select(sorted_jets, lambda j : op.AND(j.pt > 30., op.abs(j.eta)< 2.5, j.jetId &2))   
        
        # exclude from the jetsSel any jet that happens to include within its reconstruction cone a muon or an electron.
        jets= op.select(jetsSel, lambda j : op.AND(op.NOT(op.rng_any(electrons, lambda ele : op.deltaR(j.p4, ele.p4) < 0.3 )), op.NOT(op.rng_any(muons, lambda mu : op.deltaR(j.p4, mu.p4) < 0.3 ))))
        
        # order jets 
        cleanedJetsByDeepFlav = op.sort(jets, lambda j: -j.btagDeepFlavB)
        cleanedJetsByDeepB = op.sort(jets, lambda j: -j.btagDeepB)
        
        btaggingWPs = {
                "DeepCSV":{ # era: (loose, medium, tight)
                            #"2016": (0.2217, 0.6321, 0.8953), 
                            #"2017":(0.1522, 0.4941, 0.8001), 
                            "2018":(0.1241, 0.4184, 0.7527)
                            },
                "DeepFlavour":{
                            #"2016":(0.0614, 0.3093, 0.7221), 
                            #"2017":(0.0521, 0.3033, 0.7489), 
                            "2018": (0.0494, 0.2770, 0.7264)
                            }
                }
        
        bjets = {}
        # bjets ={ "DeepFlavour": {"L": jets pass loose  , "M":  jets pass medium  , "T":jets pass tight    }     
        #           "DeepCSV":    {"L":    ---           , "M":         ---        , "T":   ----            }
        #        }
        WorkingPoints = ["L", "M", "T"]
        for tagger  in btaggingWPs.keys():
            
            bJets_deepflavour ={}
            bJets_deepcsv ={}
            
            for wp in sorted(WorkingPoints):
                
                suffix = ("loose" if wp=='L' else ("medium" if wp=='M' else "tight"))
                idx = ( 0 if wp=="L" else ( 1 if wp=="M" else 2))
                
                if tagger=="DeepFlavour":
                    
                    print ("Btagging: Era= {0}, Tagger={1}, Pass_{2}_working_point={3}".format(era, tagger, suffix, btaggingWPs[tagger][era][idx] ))
                    print ("btag_{0}_94X".format(era).replace("94X", "102X" if era=="2018" else "94X"), "{0}_{1}".format('DeepJet', suffix))
                    
                    bJets_deepflavour[wp] = op.select(cleanedJetsByDeepFlav, lambda j : j.btagDeepFlavB >= btaggingWPs[tagger][era][idx] )
                    Jet_DeepFlavourBDisc = { "BTagDiscri": lambda j : j.btagDeepFlavB }
                    
                    bjets[tagger]=bJets_deepflavour
                    
                else:
                    print ("Btagging: Era= {0}, Tagger={1}, Pass_{2}_working_point={3}".format(era, tagger, suffix, btaggingWPs[tagger][era][idx] ))
                    print (suffix, era)
                    print ("btag_{0}_94X".format(era).replace("94X", "102X" if era=="2018" else "94X"), "{0}_{1}".format('DeepJet', suffix))
                    
                    bJets_deepcsv[wp] = op.select(cleanedJetsByDeepB, lambda j : j.btagDeepB >= btaggingWPs[tagger][era][idx] )   
                    Jet_DeepCSVBDis = { "BTagDiscri": lambda j : j.btagDeepB }
                    
                    bjets[tagger]=bJets_deepcsv
        

        MET = t.MET if era != "2017" else t.METFixEE2017
        corrMET=METcorrection(MET,t.PV,sample,era,self.isMC(sample))
        
        #The contribution from Z + jets events is reduced by applying a veto around the Z boson mass when the two leptons have the same flavour, ( | Mll -MZ | > 10 GeV).
        osdilep = lambda l1,l2 : op.AND(l1.charge != l2.charge, op.invariant_mass(l1.p4, l2.p4) >12. )
        Zboson_Vetocut = lambda l1, l2 : op.abs(op.invariant_mass(l1.p4, l2.p4) - 90.) > 10.
        
        osll = {
                "elmu": op.combine((electrons, muons), pred=lambda ele,mu : op.AND(osdilep(ele,mu), ele.pt > mu.pt)),
                "muel": op.combine((muons, electrons), pred=lambda mu,ele : op.AND(osdilep(mu,ele), mu.pt > ele.pt))
             }
        
        elemuTrigSF = getScaleFactor("dilepton", ("elemuLeg_HHMoriond17_2016"), systName="elmutrig")
        mueleTrigSF = getScaleFactor("dilepton", ("mueleLeg_HHMoriond17_2016"), systName="mueltrig")
       
        osllSFs = {
                "elmu" : (lambda lep :  [ muonIDSF(lep[1]), muonIsoSF(lep[1]) , eleRecoSF(lep[0]), eleIDSF(lep[0]), elemuTrigSF(lep)]), #, eleTriggerSF(lep[0]), muonTriggerSF(lep[1]) ]),
                "muel" : (lambda lep :  [ muonIDSF(lep[0]), muonIsoSF(lep[0]) , eleRecoSF(lep[1]), eleIDSF(lep[1]), mueleTrigSF(lep)])  #  eleTriggerSF(lep[1]), muonTriggerSF(lep[0]) ]),
            }

        osllcatrng = lambda catrng : op.rng_len(catrng) > 0
        hasOSLL = noSel.refine("hasOSLL", cut=op.OR(*( osllcatrng(rng) for rng in osll.values())))
        
        forceDefine(t._Jet.calcProd, hasOSLL)
        forceDefine(getattr(t, "_{0}".format("MET" if era != "2017" else "METFixEE2017")).calcProd, hasOSLL)


        categories = dict((channel, (leadpair[0],
                            hasOSLL.refine("{0}_EventSelection".format(channel), cut=osllcatrng(leadpair), weight=(osllSFs[channel](leadpair[0]) if isMC else None))
                            )) for channel, leadpair in osll.items())

        def makeJetPlots(sel, jets, uname):
            maxJet=2
            binScaling=1
            plots = []
            for i in range(maxJet):
                plots.append(Plot.make1D(f"{uname}_jet{i+1}_pt", jets[i].pt, sel,
                            EqBin(60 // binScaling, 30., 730. - max(2, i) * 100), title=f"{utils.getCounter(i+1)} jet p_{{T}} (GeV)",
                            plotopts=utils.getOpts(uname)))
                plots.append(Plot.make1D(f"{uname}_jet{i+1}_eta", jets[i].eta, sel,
                            EqBin(50 // binScaling, -2.4, 2.4), title=f"{utils.getCounter(i+1)} jet eta",
                            plotopts=utils.getOpts(uname, **{"log-y": False})))
                plots.append(Plot.make1D(f"{uname}_jet{i+1}_phi", jets[i].phi, sel,
                            EqBin(50 // binScaling, -3.1416, 3.1416), title=f"{utils.getCounter(i+1)} jet phi", 
                            plotopts=utils.getOpts(uname, **{"log-y": False})))

            return plots
        
        def makeBJetPlots(sel, bjets, wp, uname):
            
            plots = []
            binScaling=1
            for key in sel.keys():
                tagger=key.replace(wp, "")
                bjets_ = safeget(bjets, tagger, wp)
                for i in range(2):
                    plots.append(Plot.make1D(f"{uname}_bJet{i+1}_pT_{key}".format(key=key), 
                                bjets_[i].pt,
                                safeget(sel, key), 
                                EqBin(60 // binScaling, 30., 730. - max(2, i) * 100), 
                                title=f"{utils.getCounter(i+1)} bJet pT [GeV]", 
                                plotopts=utils.getOpts(uname)))
                    
                    plots.append(Plot.make1D(f"{uname}_bJet{i+1}_eta_{key}".format(key=key), 
                                bjets_[i].eta,
                                safeget(sel, key), 
                                EqBin(50 // binScaling, -2.4, 2.4), 
                                title=f"{utils.getCounter(i+1)} bJet eta", 
                                plotopts=utils.getOpts(uname)))
            
                    plots.append(Plot.make1D(f"{uname}_bJet{i+1}_phi_{key}".format(key=key), 
                                bjets_[i].phi, 
                                safeget(sel, key),
                                EqBin(50 // binScaling, -3.1416, 3.1416), 
                                title=f"{utils.getCounter(i+1)} bJet phi", 
                                plotopts=utils.getOpts(uname)))
            return plots

        def makeDiscriminatorPlots(sel, bjets, wp, uname):
            plots = []
            for key in sel.keys():
                 
                tagger=key.replace(wp, "")
                bjets_ = safeget(bjets, tagger, wp)
                idx = ( 0 if wp=="L" else ( 1 if wp=="M" else 2))
                bin0 = btaggingWPs[tagger][era][idx]
                for i in range(2):
                    if tagger =="DeepFlavour":
                        plots.append(Plot.make1D(f"{uname}_jet{i+1}_discr_deepFlav{wp}", 
                                        bjets_[i].btagDeepFlavB,
                                        safeget(sel, "DeepFlavour{0}".format(wp)),
                                        EqBin(60, bin0, 1.), 
                                        title="DeepFlavourBDisc {0}".format(wp), 
                                        plotopts=utils.getOpts(uname.lower())))
                    else:
                        plots.append(Plot.make1D(f"{uname}_jet{i+1}_discr_deepCSV{wp}", 
                                        bjets_[i].btagDeepB,
                                        safeget(sel, "DeepCSV{0}".format(wp)),
                                        EqBin(60, bin0, 1.), 
                                        title= "DeepCSVBDisc {0}".format(wp), 
                                        plotopts=utils.getOpts(uname.lower())))
    
            return plots

        def makeLeptonPlots(sel, leptons, uname):
            binScaling=1
            plots = []
            if "mu" in uname:
                flav = "Muon"
            if "el" in uname:
                flav = "Electron"
            for i in range(2):
                plots.append(Plot.make1D(f"{uname}_lep{i+1}_pt", leptons[i].pt, sel,
                                EqBin(60 // binScaling, 30., 530.), title="%s p_{T} [GeV]" % flav,
                                plotopts=utils.getOpts(uname)))
                plots.append(Plot.make1D(f"{uname}_lep{i+1}_eta", leptons[i].eta, sel,
                                EqBin(50 // binScaling, -2.4, 2.4), title="%s eta" % flav,
                                plotopts=utils.getOpts(uname, **{"log-y": False})))
                plots.append(Plot.make1D(f"{uname}_lep{i+1}_phi", leptons[i].phi, sel,
                                EqBin(50 // binScaling, -3.1416, 3.1416), title="%s phi" % flav,
                                plotopts=utils.getOpts(uname, **{"log-y": False})))
            return plots
        def makeJetmultiplictyPlots(sel, jets, uname):
            binScaling=1
            plots=[]
            plots.append(Plot.make1D(f"{uname}_Jet_mulmtiplicity", op.rng_len(jets), sel,
                            EqBin(10, 2., 8.), title="Jet mulmtiplicity",
                            plotopts=utils.getOpts(uname)))
            plots.append(Plot.make1D(f"{uname}_nVX", 
                            t.PV.npvs, sel, 
                            EqBin(50 // binScaling, 0., 60.), title="reconstructed vertices",
                            plotopts=utils.getOpts(uname)))
            # FIXME
            #for i in range(2, 5):
            #    event=sel.refine("{i}JetSel".format(i), cut=[ op.rng_len(jets) ==i ])
            #    plots.append(Plot.make1D(f"{uname}_EventSelection", 
            #                    event, event,
            #                    EqBin(50 // binScaling, 0., i), title="Event selection",
            #                    plotopts=utils.getOpts(uname)))
            return plots
        
        def makeMETPlots(sel, leptons, met, corrMET, uname):
            binScaling=1
            plots = []

            plots.append(Plot.make1D(f"{uname}_MET_pt", met.pt, sel,
                        EqBin(60 // binScaling, 0., 600.), title="MET p_{T} [GeV]",
                        plotopts=utils.getOpts(uname)))
            plots.append(Plot.make1D(f"{uname}_MET_phi", met.phi, sel,
                        EqBin(60 // binScaling, -3.1416, 3.1416), title="MET #phi",
                        plotopts=utils.getOpts(uname, **{"log-y": False})))
            for i in range(2):
                plots.append(Plot.make1D(f"{uname}_MET_lep{i+1}_deltaPhi",
                            op.Phi_mpi_pi(leptons[i].phi - met.phi), sel, EqBin(60 // binScaling, -3.1416, 3.1416),
                            title="#Delta #phi (lepton, MET)", plotopts=utils.getOpts(uname, **{"log-y": False})))
                
                MT = op.sqrt( 2. * met.pt * leptons[i].p4.Pt() * (1. - op.cos(op.Phi_mpi_pi(met.phi - leptons[i].p4.Phi()))) )
                plots.append(Plot.make1D(f"{uname}_MET_MT_lep{i+1}", MT, sel,
                            EqBin(60 // binScaling, 0., 600.), title="Lepton M_{T} [GeV]",
                            plotopts=utils.getOpts(uname)))

            return plots


        def MakeEXTRAMETPlots(sel, corrmet, met, uname):
            plots = []
            binScaling=1

            plots.append(Plot.make1D(f"{uname}_met_pT",
                        met.pt, 
                        sel, 
                        EqBin(60 // binScaling, 0., 600.), title="MET p_{T} [GeV]",
                        plotopts=utils.getOpts(uname)))
                        
            plots.append(Plot.make1D(f"{uname}_xycorrmet_pT",
                        corrmet.pt, 
                        sel,
                        EqBin(60 // binScaling, 0., 600.), title="corrMET p_{T} [GeV]",
                        plotopts=utils.getOpts(uname)))

            return plots

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


        def TwoTagCount(uname, histosName, key, jets, bjets, _2l2jsel, _2l2bjsel, iPt, ptBin):
            plots=[]
            
            # check if both jets are falling in the same bins befor filling the histos
            #if(iPt > -1): # -1 is inclusive--> always keep -1
            #    pt1 = bjets[0].pt
            #    pt2 = bjets[1].pt

            #    if(lowerPtBinEdges[iPt] > max(pt1,pt2)):
            #        return
            #    if(ptBin != len(lowerPtBinEdges)):
            #        if(lowerPtBinEdges[iPt+1] < min(pt1,pt2)):
            #            return
            nTrueB = (op.rng_count(jets, lambda j : j.hadronFlavour == 5) if isMC else op.c_int(-1))
            plots.append(Plot.make1D("{0}_Events_2l2j_truthFlav_2b_{1}_{2}".format(uname, key, ptBin), 
                                        nTrueB, 
                                        _2l2jsel,                                        
                                        EqBin(2, 0., 2.),
                                        title="N b-tags (mc truth)",
                                        plotopts=utils.getOpts(uname)))

            plots.append(Plot.make1D("{0}_Events_2l2j_2btagged_truthFlav_1atleast_light_or_charm_{1}_{2}".format(uname, key, ptBin),
                                        nTrueB < 2,
                                        _2l2bjsel.get(key),
                                        EqBin(2, 0., 2.),
                                        title="N b-tags (mc truth)",
                                        plotopts=utils.getOpts(uname)))
            
            ########## 
            if isMC:
                flavour1=(bjets[0].hadronFlavour)
                flavour2=(bjets[1].hadronFlavour)
    
                TwoTagCrossFlavour= {
                        # both are b-hadron (truth flav): but at last 1 is btagged --> it can be  0 passing or 1 passing or 2 passing my btagging wp for specific tagger 
                        "2b"   : [flavour1 == 5, flavour2 == 5],
                        # these catgories  : both jet are btagged --> but one at least is light or charm hadron flav
                        "1b_1c": [flavour1==5,  flavour2==4],
                        "1b_1l": [flavour1==5, (flavour2==21 or (flavour2>0 and flavour2<4))],
                        "1c_1l": [flavour1==4, (flavour2==21 or (flavour2>0 and flavour2<4))],
                        "2c"   : [flavour1==4,  flavour2==4],
                        "2l"   : [(flavour1==21 or (flavour1>0 and flavour1<4)), (flavour2==21 or (flavour2>0 and flavour2<4))]
                        }
            
            nTrueB = (op.rng_count(jets, lambda j : j.hadronFlavour == 5) if isMC else op.c_int(-1))
            tagger=key.replace(WP, "")
            if tagger == "DeepFlavour":
                bjetsfortagger= bJets_PassdeepflavourWP
            else:
                bjetsfortagger= bJets_PassdeepcsvWP

            nbtagged = {"0jet_btagged": op.rng_len(bjetsfortagger) ==0 ,
                        "1jet_btagged": op.rng_len(bjetsfortagger) ==1 ,
                        "2jet_btagged": op.rng_len(bjetsfortagger) ==2 ,
                    }

            #nTrueBFlav = (op.rng_count(bjets, lambda j : j.hadronFlavour == 5) if isMC else op.c_int(-1))
            #plots.append(Plot.make1D("{0}_Events_2l2j_truthFlav_2b_{1}_{2}_ver1".format(uname, key, ptBin), 
            #                            nTrueBFlav, 
            #                            _2l2jsel.refine("2l2j_truthFlav_2b_{0}btagged".format(),cut=op.OR(*nbtagged.values())),
            #                            EqBin(2, 0., 2.),
            #                            title="N b-tags (mc truth)",
            #                            plotopts=utils.getOpts(uname)))

            #nFakeJet_2btagged = (op.rng_count(bjets, lambda j : [j[0].hadronFlavour,j[1].hadronFlavour] == op.OR(*TwoTagCrossFlavour.values())) if isMC else op.c_int(-1))
            #plots.append(Plot.make1D("{0}_Events_2l2j_2btagged_truthFlav_1atleast_light_or_charm_{1}_{2}_ver1".format(uname, key, ptBin),
            #                            nFakeJet_2btagged,
            #                            _2l2bjsel.get(key),
            #                            EqBin(2, 0., 2.),
            #                            title="N b-tags (mc truth)",
            #                            plotopts=utils.getOpts(uname)))
            
            return plots 

    # ---- Ask for plots  --- 
        for channel, (leptons, cat) in categories.items():

            # before cutting on the jets len 
            plots.extend(makeJetmultiplictyPlots(cat, jets, channel))
            
            TwoLeptonsTwoJets=cat.refine("twoJet{0}Sel".format(channel), cut=[ op.rng_len(jets) ==2 ])
            plots.extend(makeLeptonPlots(TwoLeptonsTwoJets, leptons, channel))
            plots.extend(makeJetPlots(TwoLeptonsTwoJets, jets, channel))
            
            plots.extend(makeMETPlots(TwoLeptonsTwoJets, leptons, MET, corrMET, channel))
            plots.extend(MakeEXTRAMETPlots(TwoLeptonsTwoJets, corrMET, MET, channel))
            
            for WP in WorkingPoints:
                bJets_PassdeepflavourWP=safeget(bjets, "DeepFlavour", WP)
                bJets_PassdeepcsvWP=safeget(bjets, "DeepCSV", WP)
                
                TwoLeptonsTwoBjets = {
                    "DeepFlavour{0}".format(WP) :  TwoLeptonsTwoJets.refine("TwoLeptonsTwoBjets_DeepFlavour{0}_{1}".format(WP, channel),
                                                                        cut=[ op.rng_len(bJets_PassdeepflavourWP) ==2 ]),
                    "DeepCSV{0}".format(WP)     :  TwoLeptonsTwoJets.refine("TwoLeptonsTwoBjets_DeepCSV{0}_{1}".format(WP, channel), 
                                                                        cut=[ op.rng_len(bJets_PassdeepcsvWP) ==2 ])
                                            }

                plots.extend(makeBJetPlots(TwoLeptonsTwoBjets, bjets, WP, channel))
                plots.extend(makeDiscriminatorPlots(TwoLeptonsTwoBjets, bjets, WP, channel))
                
                # TwoTag Count Method --> start filling 2tag crossFlavour histograms
                for key in TwoLeptonsTwoBjets.keys():
                    tagger=key.replace(WP, "")
                    bjets_=safeget(bjets, tagger, WP)
                    for iPT in range(-1,len(lowerPtBinEdges)):
                        histosName = channel+"_only2_twoTags_"+ key + "_"+ ReturnPtLabel(iPT)
                        TwoTagCount(channel, histosName, key, jets, bjets_, TwoLeptonsTwoJets, TwoLeptonsTwoBjets, iPT, ReturnPtLabel(iPT))

        return plots

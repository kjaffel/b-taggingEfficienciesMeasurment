import ROOT
import sys, os
import math
import csv
import glob
import argparse
import yaml
import logging
logger = logging.getLogger("BTAG Plotter")

lowerPtBinEdges=[30,50,70,100,140]#,200,300,600] 
def get_options():
    parser = argparse.ArgumentParser(description='Create TwoTag Efficiencies maps ')
    parser.add_argument('-i', '--input', action='store', dest='root_path', type=str, default='/home/ucl/cp3/kjaffel/bamboodev/b-taggingEfficienciesMeasurment/newversion/2017UL/controlplots.v9._/results/normalizedFor2TagCount', help='histFactory input path')
    parser.add_argument('-c', '--config', type=str, help='Input YML file')
    parser.add_argument('-e', '--era',   action='store', dest='era',       type=str, default= '2017',   help='you need to pass your era ')
    parser.add_argument('-t', '--tagCount',action='store', dest='tagCount',type=str, default= 'atleast2',   help='What do you want as cut on your jets multiplicty: **atleast2** or **Only2** ')
    parser.add_argument('-cat', '--cat',action='store', dest='cat', type=str, default= 'elmu',   help='choose channel: muel, elmu')
    options = parser.parse_args()
    return options

def getLuminosity(era):
    if era =='2016':
        lumi = 35921.875594646
    elif era =='2017':
        lumi = 41529.152060112
    else:
        lumi = 59740.565201546
    return lumi

def getLuminosityUncertainty(era):
    if era == '2016':
        uncer = 1.025
    elif era == '2017':
        uncer = 1.023
    elif era == '2018':
        uncer=  1.025
    return uncer

def mistagrate_uncer(baseName, era):
    tag=baseName.split("_")
    var=(tag[len(tag)-2] if (baseName.find("up")!=-1 or baseName.find("down")!=-1) else tag[len(tag)-1])
    wp= ("L" if var.find("L")!=-1 else ("M" if var.find("M")!=-1 else "T"))
    pt=var.split(wp)
    ptLabel=var.split(wp)

    if var.find("Inclusive")!=-1:
        x=0
    elif var.find("Inf")!=-1:
        pt=ptLabel[1].split("to")
        x=(int(pt[0]))/2
    else:
        pt=ptLabel[1].split("to")
        ptmin=pt[0]
        ptmax=pt[1]
        x=(int(ptmin)+int(ptmax))/2

    if baseName.find("DeepCSV")!=-1:
        uncer=( (0.0559259+1.96455e-05*x+-3.60571e-08*x*x) if wp=="L" else ( (0.122811+0.000162564*x+-1.66422e-07*x*x) if wp=="M" else (0.207667+0.000165081*x+-1.45307e-07*x*x)))
    else:
        uncer=( (0.0631436+7.87525e-05*x+-1.10405e-07*x*x) if wp=="L" else ( (0.142253+0.000227323*x+-2.71704e-07*x*x) if wp=="M" else (0.223324+0.000231341*x+-2.34362e-07*x*x)))
    return uncer
    
def ReturnPtLabel(iPT):
    ptBinLabel=""
    if(iPT==-1):
        ptBinLabel="Inclusive"
    else:
        ptBinLabel+=str(lowerPtBinEdges[iPT])+"to"
        if(iPT==len(lowerPtBinEdges)-1):
            ptBinLabel+="Inf"
        else:
            ptBinLabel+=str(lowerPtBinEdges[iPT+1])
    return ptBinLabel

def DumpResults(ptLabel, tag, scalefactor, measurementType, line_count, ultralegacy):
    wp= (0 if tag.find("L")==len(tag)-1 else (1 if tag.find("M")==len(tag)-1 else 2))
    ptMin=(ptLabel.split("to"))[0]
    ptMax=(ptLabel.split("to"))[1]
    if ptMax==str("Inf"):
        ptMax=ptMax.replace(str("Inf"),str(200))
    
    tagger=("DeepCSV" if "DeepCSV" in tag else ("DeepFlavour"))
    GlobalTag = ("RunIIAutumn18" if options.era =='2018' else ("RunIISummer19UL17NanoAOD" if options.era=='2017' and ultralegacy else ("RunIIFall17NanoAOD")))
    if measurementType != "Central":
        csv_file = '{0}_{1}_SFs_syst_stat.csv'.format(tagger, GlobalTag)
    else:
        csv_file = '{0}_{1}_SFs_stat.csv'.format(tagger, GlobalTag)
    mode= ('a' if line_count ==0 else('w'))
    with open(csv_file, mode=mode) as f:
        btag_params = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        if line_count == 0:
            btag_params.writerow([tagger,'OperatingPoint', 'measurementType', 'sysType', 'jetFlavor', 'etaMin','etaMax', 'ptMin','ptMax', 'discrMin','discrMax', 'SFs'])
        btag_params.writerow([wp, "tagcount", measurementType, 0, -2.5, 2.5, ptMin, ptMax, 0, 1, scalefactor])
    return csv_file
                    
def main():
    global options
    options = get_options()
    print ( 'era:', options.era, 'lumi:', getLuminosity(options.era))
    options.luminosity= getLuminosity(options.era)
    era = options.era
    root_path =options.root_path
    tagCount = options.tagCount
    cat = options.cat
    DrawPlot_foreachptBin = False
    ultralegacy = True

    def ComputeTwoTagEfficieny(basename, root_path, era, combined):
        Truthmc2bFlavour= ROOT.TH1F("Truth_plotter", "Events", 3, 0., 3.)
        Data2btag= ROOT.TH1F("data_plotter", "Events", 3, 0., 3.)
        Fakemc2btag= ROOT.TH1F("Fake_plotter", "Events", 4, 0., 4.)
       
        histo_data = ROOT.TH1F("data_plotter", "Events", 12, 0., 12.)
        histo_mc = ROOT.TH1F("mc_plotter", "Events", 12, 0., 12.)
        tot={}
        statError={}
        tot["2bmc"]=0   # all the events in which both jets are 2b jet truth flav   with 0 or 1 b-tagged 
        tot["data"]=0   ## all data
        tot["othermc"]=0 
        tot["2bmc_2btag"]=0  ##  number of events in which both jets are b-jets and both are b-tagged
        tot["othermc_2btag"]=0 ##  combinations of different flavour jets (1b_1c, 1b_1l, 1l_1c, ,2l, 2c) where both are btagged 
        tot["data_2btag"]=0   ## data where we have 2b jets && both passing b-tagging wp 
        print(basename)
        for filename in glob.glob(os.path.join(options.root_path, '*.root')):
            f = ROOT.TFile(filename)
            split_filename = filename.split('/')
            smp = str(split_filename[-1])
            if combined:
                if '__skeleton__' in smp :
                    pass
                elif 'MuonEG' in smp:
                    data = f.Get(basename)
                    data.SetDirectory(0)
                    histo_data.Add(data)
                else:
                    mc = f.Get(basename)
                    mc.SetDirectory(0)
                    histo_mc.Add(mc)

            else:
                histName1= basename.replace("nbtagged", "2btagged")    # 2btagged , 2Truth b 
                histName2= basename.replace("nBeauty", "2Beauty")      # 0, 1 , or 2 btagged , 2Truth b
                histName3= basename.replace("nBeauty", "notBeauty").replace("nbtagged", "2btagged")     # 2btagged , Truth both not b 
                if '__skeleton__' in smp :
                    pass
                elif 'MuonEG' in smp:
                    if 'up' not in histName1 and 'down' not in histName1:
                        datahisto = f.Get(histName1)
                        datahisto.SetDirectory(0)
                        Data2btag.Add(datahisto)
                else:
                    mcFake= f.Get(histName3)
                    mcFake.SetDirectory(0)
                    Fakemc2btag.Add(mcFake)
                    
                    mcTruth=f.Get(histName2)
                    mcTruth.SetDirectory(0)
                    Truthmc2bFlavour.Add(mcTruth)
        if combined:
            tot["2bmc"]=tot["2bmc"]+histo_mc.GetBinContent(4)
            tot["2bmc"]=tot["2bmc"]+histo_mc.GetBinContent(8)
            
            tot["2bmc_2btag"]=tot["2bmc_2btag"]+histo_mc.GetBinContent(12)
            for iBin in range(9, 12):
                tot["othermc_2btag"]=tot["othermc_2btag"]+histo_mc.GetBinContent(iBin)  
            tot["data_2btag"]=tot["data_2btag"]+histo_data.GetBinContent(9)
            for name in tot:
                if  name.find("mc")!=-1:
                    tot[name]=tot[name]/histo_mc.Integral()
                else:
                    tot[name]=tot[name]/histo_data.Integral()

        else:
            print ( "hist_Truthmc2bFlavour.Integral():",  Truthmc2bFlavour.Integral())
            print ( "hist_Fakemc2btag.Integral():",  Fakemc2btag.Integral())
            print ( "hist_data2btag.Integral():", Data2btag.Integral() )
            
            tot["2bmc_2btag"]=tot["2bmc_2btag"]+Truthmc2bFlavour.GetBinContent(3)
            for i in range(1, 4):
                tot["2bmc"]= tot["2bmc"] + Truthmc2bFlavour.GetBinContent(i)
            for i in range(1, 4):
                tot["othermc_2btag"]=tot["othermc_2btag"]+Fakemc2btag.GetBinContent(i)
            tot["data_2btag"]=tot["data_2btag"]+Data2btag.GetBinContent(1)
    
            print ( "data_2btag:", tot["data_2btag"], "2bmc_2btag:", tot["2bmc_2btag"], "othermc_2btag:", tot["othermc_2btag"], "2bmc:", tot["2bmc"] )
            
            tot["data_2btag"]=tot["data_2btag"]/Data2btag.Integral()
            tot["2bmc"]= tot["2bmc"]/Truthmc2bFlavour.Integral()
            tot["othermc_2btag"]=tot["othermc_2btag"]/Fakemc2btag.Integral()
            tot["2bmc_2btag"] = tot["2bmc_2btag"]/Truthmc2bFlavour.Integral()

        print ( "data_2btag:", tot["data_2btag"], "2bmc_2btag:", tot["2bmc_2btag"], "othermc_2btag:", tot["othermc_2btag"], "2bmc:", tot["2bmc"] )
        try:
            #eff_b = math.sqrt((tot["data_2btag"]-tot["othermc_2btag"]) / (tot["2bmc"] +tot["2bmc_2btag"]))
            eff=math.sqrt((tot["data_2btag"]-tot["othermc_2btag"])/(tot["2bmc"]))
            effmc=math.sqrt(tot["2bmc_2btag"]/(tot["2bmc"]))
            
            #effUp=math.sqrt((tot["data_2btag"]-tot["othermc_2btag"]*(1-mistagrate(basename, era)))/tot["2bmc"])
            #effDown=math.sqrt((tot["data_2btag"]-tot["othermc_2btag"]*(1+mistagrate(basename, era)))/tot["2bmc"])
        
            # 50 % conservative uncertainty 
            effUp=math.sqrt((tot["data_2btag"]-tot["othermc_2btag"]*1.5)/tot["2bmc"])
            effDown=math.sqrt((tot["data_2btag"]-tot["othermc_2btag"]*0.5)/tot["2bmc"])
           
            suffix= basename.replace("%s_"%cat,"")
            if combined:
                statError[suffix +"data_2btag"] = math.sqrt( tot["data_2btag"] * (1 - tot["data_2btag"]) / histo_data.GetEffectiveEntries())
                statError[suffix +"othermc_2btag"] = math.sqrt( tot["othermc_2btag"] * (1 - tot["othermc_2btag"]) / histo_mc.GetEffectiveEntries())
                statError[suffix +"2bmc"] = math.sqrt( tot["2bmc"] * (1 - tot["2bmc"]) /  histo_mc.GetEffectiveEntries())
            else:
                statError[suffix +"data_2btag"] = math.sqrt( tot["data_2btag"] * (1 - tot["data_2btag"]) / Data2btag.GetEffectiveEntries())
                statError[suffix +"othermc_2btag"] = math.sqrt( tot["othermc_2btag"] * (1 - tot["othermc_2btag"]) / Fakemc2btag.GetEffectiveEntries())
                statError[suffix +"2bmc"] = math.sqrt( tot["2bmc"] * (1 - tot["2bmc"]) /  Truthmc2bFlavour.GetEffectiveEntries())
            statError[suffix +"eff"] = math.sqrt(math.pow(eff*statError[suffix +"2bmc"],2) + (math.pow(statError[suffix+"data_2btag"],2)+math.pow(statError[suffix+"othermc_2btag"],2))/ math.pow(eff,2))/(2*tot["2bmc"])
            logger.info("basename:", basename,  "effdata:",eff, "effmc:", effmc , "effUp:", effUp, "effDown:", effDown)
            return [tot["2bmc_2btag"],tot["2bmc"],tot["data_2btag"],tot["othermc_2btag"],eff,effmc,eff/effmc,effUp/effmc,effDown/effmc, statError[suffix+"eff"]]
        except:
            logger.exception(f"unable to get Efficency from baseName :{basename} **  Eff ~ 0  (will return 1. the default value ) Keep In Mind : This exception can be avoided if you exclude histograms where in ptBin the Integral is almost None " )
       
    
    tags=["DeepCSVL","DeepCSVM","DeepCSVT", "DeepFlavourL","DeepFlavourM","DeepFlavourT"]
    tags.sort()
    systematics = [""]#, "PSFSR", "PSISR", "qcdScale"]
    #systematics = ["", 'jer', 'jesTotal', 'muid', 'muiso', 'elid', 'ele_reco', 'pu', 'elmutrig', 'L1PreFiring', 'HLTZvtxSF']
            
    #with open(options.config) as _f:
    #    plotConfig = yaml.load(_f, Loader=yaml.FullLoader)

    #systematics.append(plotConfig["systematics"])
    tgraphs={}
    results={}

    can=ROOT.TCanvas("can","",800,800)
    can.SetGridy()
    can.SetGridx()
    outputFile=ROOT.TFile("twoTagSFGraphs.root","RECREATE")
    
    for syst in systematics:
        print (syst)
        for key in tags:
            tag = f"{tagCount}_twoTags_{key}"
            for ptBin in range(-1,len(lowerPtBinEdges)):
                ptLabel=ReturnPtLabel(ptBin)
                analysiskeys= (tag,ptLabel,syst)
                basename= f"{cat}_Events_2l{tagCount}j_nbtagged_TruthFlav_nBeauty_{key}_{ptLabel}"
                #basename= f"NTaggedCrossTruth_2lOSSF{tagCount}_{key}_{ptLabel}"
                print (80*"_")
                print ("analysiskeys:", cat, analysiskeys)
                if syst =="":
                    results[analysiskeys]=ComputeTwoTagEfficieny(basename, root_path, era, combined=False)
                else:
                    resultsup=ComputeTwoTagEfficieny(f"{basename}__{syst}up", root_path, era, combined=False)
                    resultsdown=ComputeTwoTagEfficieny(f"{basename}__{syst}down", root_path, era, combined=False)
                    statErrorup=ComputeTwoTagEfficieny(f"{basename}__{syst}up", root_path, era, combined=False)
                    statErrordown=ComputeTwoTagEfficieny(f"{basename}__{syst}down", root_path, era, combined=False)
                    #I just set the numbers to 1 because you don't really need them for systematics unless you're debugging
                    # Instead you should get sth like :
                    # return [tot["2bmc_2btag"],tot["2bmc"],tot["data_2btag"],tot["othermc_2btag"],eff,effmc,eff/effmc,effUp/effmc,effDown/effmc, statError[suffix+"eff"]]
                    results[analysiskeys]=[1.,1.,1.,1.,1.,1.,(resultsup[6] if resultsup else 1.), (resultsdown[6] if resultsdown else 1.)]
                    print (results[analysiskeys])
    multiG={}
    multiG2={}
    for index, key in enumerate(tags):
        tag = f"{tagCount}_twoTags_{key}"
        statuncertainty2up=[]
        statuncertainty2down=[]
        systuncertainty2up=[]
        systuncertainty2down=[]
    
        it=0
        if index==0:
            tgraphs["wp"]=ROOT.TGraphAsymmErrors()
            tgraphs["wp"].SetTitle(";Jet P_{T} [GeV];Scale Factor")
        for syst in systematics:
            iGraph=0
            graphKey=tag+syst
            tgraphs[graphKey]=ROOT.TGraphAsymmErrors()
            tgraphs[graphKey].SetTitle(graphKey+";Jet P_{T} [GeV];Scale Factor")
            if syst==systematics[0]:
                tgraphs[tag+"totalError"]=ROOT.TGraphAsymmErrors()
                tgraphs[tag+"systematics"]=ROOT.TGraphAsymmErrors()
                tgraphs[tag+"statError"]=ROOT.TGraphErrors()
                tgraphs[tag+"eff_b"]=ROOT.TGraphErrors()
            for ptBin in range(0,len(lowerPtBinEdges)):
                ptLabel=ReturnPtLabel(ptBin)
                # nominal must always be run for this to work
                if ptBin>-1:
                    pt1=lowerPtBinEdges[ptBin]
                    pt2=200
                    if pt1!=lowerPtBinEdges[-1]:
                        pt2=lowerPtBinEdges[ptBin+1]
                try:
                    if syst=="":
                        print (tag, ptLabel, results[(tag,ptLabel,"")])
                        nominalSF=results[(tag,ptLabel,"")][6]
                        nominalSFUp=results[(tag,ptLabel,"")][7]
                        nominalSFDown=results[(tag,ptLabel,"")][8]
                        statError=results[(tag,ptLabel,"")][9]
                        eff_b=results[(tag,ptLabel,"")][4]
                            
                        tgraphs[tag+"systematics"].SetPoint(iGraph,(pt1+pt2)/2.,nominalSF)
                        tgraphs[tag+"systematics"].SetPointEXhigh(iGraph,(-pt1+pt2)/2.)
                        tgraphs[tag+"systematics"].SetPointEXlow(iGraph,(-pt1+pt2)/2.)
        
                        tgraphs[tag+"statError"].SetPoint(iGraph,(pt1+pt2)/2, nominalSF)
                        tgraphs[tag+"statError"].SetPointError(iGraph,(pt2-pt1)/2, statError)
    
                        tgraphs[tag+"eff_b"].SetPoint(iGraph,(pt1+pt2)/2, eff_b)
                        tgraphs[tag+"eff_b"].SetPointError(iGraph,(pt2-pt1)/2, statError)
                        
                        statuncertainty2up.append(math.pow(nominalSFUp-nominalSF,2))
                        statuncertainty2down.append(math.pow(-nominalSFDown+nominalSF,2))
    
                        systuncertainty2up.append(0)
                        systuncertainty2down.append(0)
        
                        tgraphs[graphKey].SetPoint(iGraph,(pt1+pt2)/2.,nominalSF)
                        tgraphs[graphKey].SetPointEXhigh(iGraph, (-pt1+pt2)/2.)
                        tgraphs[graphKey].SetPointEXlow( iGraph, (pt2-pt1)/2.)
                        tgraphs[graphKey].SetPointEYhigh(iGraph, nominalSFUp-nominalSF)
                        tgraphs[graphKey].SetPointEYlow( iGraph, -nominalSFDown+nominalSF)
                    
                        tgraphs["wp"].SetPoint(it,(pt1+pt2)/2.,nominalSF)
                        tgraphs["wp"].SetPointEXhigh(it, (-pt1+pt2)/2.)
                        tgraphs["wp"].SetPointEXlow( it, (pt2-pt1)/2.)
                        tgraphs["wp"].SetPointEYhigh(it, nominalSFUp-nominalSF)
                        tgraphs["wp"].SetPointEYlow( it, -nominalSFDown+nominalSF)
                        print ("--------------------------------------------------------------------------")
                        print (("Nominal:{0}".format(nominalSF)))
                        print (("Up     :{0}".format(nominalSFUp)))
                        print (("Down   :{0}".format(nominalSFDown)))
                        iGraph=iGraph+1
                        if index>=3:
                            it=it+1
                        DumpResults(ptLabel, tag, nominalSF, "Central", it, ultralegacy)
                    else:
                        systupSF=results[(tag,ptLabel,syst)][6]
                        systdownSF=results[(tag,ptLabel,syst)][7]
                        tgraphs[tag+"totalError"].SetPoint(iGraph,(pt1+pt2)/2.,nominalSF)
                        tgraphs[tag+"totalError"].SetPointEXhigh(iGraph,(-pt1+pt2)/2.)
                        tgraphs[tag+"totalError"].SetPointEXlow(iGraph,(-pt1+pt2)/2.)
    
                        if (systupSF-systdownSF >0.):
                            tgraphs[graphKey].SetPointEYhigh(iGraph,systupSF-nominalSF)
                            tgraphs[graphKey].SetPointEYlow( iGraph,-systdownSF+nominalSF)
        
                            statuncertainty2up[ptBin] += math.pow(systupSF-nominalSF,2)
                            statuncertainty2down[ptBin] += math.pow(-systdownSF+nominalSF,2)
                            systuncertainty2up[ptBin] += math.pow(systupSF-nominalSF,2)
                            systuncertainty2down[ptBin] += math.pow(-systdownSF+nominalSF,2)
                        else:
                            tgraphs[graphKey].SetPointEYlow( iGraph,-systupSF+nominalSF)
                            tgraphs[graphKey].SetPointEYhigh(iGraph,systdownSF-nominalSF)
                    
                            statuncertainty2up[ptBin] += math.pow(-systupSF+nominalSF,2)
                            statuncertainty2down[ptBin] += math.pow(systdownSF-nominalSF,2)
                            systuncertainty2up[ptBin] += math.pow(-systupSF+nominalSF,2)
                            systuncertainty2down[ptBin] += math.pow(systdownSF-nominalSF,2)
    
                        tgraphs[graphKey].SetPoint(iGraph,(pt1+pt2)/2.,nominalSF)
                        tgraphs[graphKey].SetPointEXhigh(iGraph,(-pt1+pt2)/2.)
                        tgraphs[graphKey].SetPointEXlow(iGraph,(pt2-pt1)/2.)
                        iGraph=iGraph+1
    
                        print ("-------------------------------------------------------------------------------")
                        print(("up   : {0}".format(systupSF)))
                        print(("down    : {0}".format(systdownSF)))
                        print(("up syst -Nominal:{0}".format(systupSF-nominalSF)))
                        print(("Nominal- down syst :{0}".format(-systdownSF+nominalSF)))
    
                        if math.copysign(1,systupSF-nominalSF) != math.copysign(1,nominalSF-systdownSF):
                            loger.info("up and down vary from nominal in the same direction for {graphKey} : Up, Nonimal, Down {systupSF} {nominalSF} {systdownSF}")                 
    
                except:
                    logger.exception("skip plotting the default value SF=1 ! ")
                if DrawPlot_foreachptBin:
                    tgraphs[graphKey].SetName(graphKey)
                    tgraphs[graphKey].Draw("APE")
                    can.Update()
                    can.SaveAs(graphKey+".png")
                    outputFile.cd()
                    tgraphs[graphKey].Write()
                
                    tgraphs["wp"].SetName("wp")
                    tgraphs["wp"].Draw("AP")
                    can.Update()
                    can.SaveAs("wp"+".png")
                    outputFile.cd()
                    tgraphs["wp"].Write()
        
        try:
            allSF=tgraphs[tag+"totalError"].GetY()
            for iGraph in range(0,len(lowerPtBinEdges)):
                tgraphs[tag+"totalError"].SetPointEYhigh(iGraph, math.sqrt(statuncertainty2up[iGraph]))
                tgraphs[tag+"totalError"].SetPointEYlow(iGraph, math.sqrt(statuncertainty2down[iGraph]))
                
                SFsUp=allSF[iGraph]+ math.sqrt(statuncertainty2up[iGraph])
                SFsDown=allSF[iGraph]-math.sqrt(statuncertainty2down[iGraph])

                DumpResults(ptLabel, tag, SFsUp, "Up", it, ultralegacy)
                DumpResults(ptLabel, tag, SFsDown, "Down", it, ultralegacy)
    
            otherSF=tgraphs[tag+"systematics"].GetY()
            for iGraph in range(0,len(lowerPtBinEdges)):
                tgraphs[tag+"systematics"].SetPointEYhigh(iGraph, math.sqrt(systuncertainty2up[iGraph]))
                tgraphs[tag+"systematics"].SetPointEYlow(iGraph, math.sqrt(systuncertainty2down[iGraph]))
                
                SFsUp=otherSF[iGraph]+ math.sqrt(systuncertainty2up[iGraph])
                SFsDown=otherSF[iGraph]-math.sqrt(systuncertainty2down[iGraph])
                
                DumpResults(ptLabel, tag, SFsUp, "Up", it, ultralegacy)
                DumpResults(ptLabel, tag, SFsDown, "Down", it, ultralegacy)
        except:
            logger.exception('scale factor not computed ; will skipp ')
        # save the statError , othersystematics and totalError once per tag
        
        optstex = {
                "elel" : 'e^+e^-',
                "mumu" : '\mu^+\mu^-',
                "muel" : '\mu^+e^-',
                "elmu" : 'e^+\mu^-'
                }

        legend = ROOT.TLegend(0.4,0.72,0.89,0.89)
        #legend.SetFillStyle(0)
        #legend.SetBorderSize(0)
        legend.SetTextSize(0.020)
        legend.SetTextFont(352)
        legend.SetTextAlign(11)
        legend.SetHeader("#splitline{CMS Preliminary __ 2017UL Data 41.53 fb^{-1} (13 TeV)}{#splitline{%s}{%s channel, 2 AK4 jets (pT > 30 GeV)}}"%(key, cat), "C");
        
        tgraphs[tag+"eff_b"].SetName(tag+"eff_b")
        tgraphs[tag+"eff_b"].SetTitle(tag)
        tgraphs[tag+"eff_b"].SetFillColor(2)
        tgraphs[tag+"eff_b"].SetFillStyle(3001)
        tgraphs[tag+"eff_b"].Draw("a2")
        tgraphs[tag+"eff_b"].Draw("p")
        #legend.AddEntry(tgraphs[tag+"eff_b"] ,"efficienciy.","lep" )
        #legend.Draw()
        tgraphs[tag+"eff_b"].GetXaxis().SetTitle( 'p_{T} [GeV]' )
        tgraphs[tag+"eff_b"].GetYaxis().SetTitle( 'Efficiency' )
        tgraphs[tag+"eff_b"].SetMaximum(1.1)
        tgraphs[tag+"eff_b"].SetMinimum(0.7)
        can.SaveAs(tag+"_eff_b"+".png")
        outputFile.cd()
        tgraphs[tag+"eff_b"].Write()
        
        tgraphs[tag+"statError"].SetName(tag+"statError")
        tgraphs[tag+"statError"].SetTitle(tag)
        tgraphs[tag+"statError"].SetFillColor(2)
        tgraphs[tag+"statError"].SetFillStyle(3001)
        tgraphs[tag+"statError"].Draw("a2")
        tgraphs[tag+"statError"].Draw("p")
        legend.AddEntry(tgraphs[tag+"statError"] ,"stat.","lep" )
        legend.Draw("same")
        tgraphs[tag+"statError"].GetXaxis().SetTitle( 'p_{T} [GeV]' )
        tgraphs[tag+"statError"].GetYaxis().SetTitle( 'Scale Factor' )
        tgraphs[tag+"statError"].SetMaximum(1.4)
        tgraphs[tag+"statError"].SetMinimum(0.7)
        can.SaveAs(tag+"_statError"+".png")
        outputFile.cd()
        tgraphs[tag+"statError"].Write()
   
        
        if len(systematics) !=1:
            tgraphs[tag+"totalError"].SetName(tag+"totalError")
            tgraphs[tag+"totalError"].SetTitle(tag)
            tgraphs[tag+"totalError"].SetFillColor(4)
            tgraphs[tag+"totalError"].SetFillStyle(3001)
            tgraphs[tag+"totalError"].Draw("a2")
            tgraphs[tag+"totalError"].Draw("p")
            legend.AddEntry(tgraphs[tag+"totalError"] ,"totalError.","f" )
            legend.Draw("same")
            tgraphs[tag+"totalError"].GetXaxis().SetTitle( 'p_{T} [GeV]' )
            tgraphs[tag+"totalError"].GetYaxis().SetTitle( 'Scale Factor' )
            tgraphs[tag+"totalError"].SetMaximum(1.4)
            tgraphs[tag+"totalError"].SetMinimum(0.8)
            can.Update()
            can.SaveAs(tag+"_totalError"+".png")
            outputFile.cd()
            tgraphs[tag+"totalError"].Write()
            
            tgraphs[tag+"systematics"].SetName(tag+"systematics")
            tgraphs[tag+"systematics"].SetTitle(tag)
            tgraphs[tag+"systematics"].SetFillColor(7)
            tgraphs[tag+"systematics"].SetFillStyle(3001)
            tgraphs[tag+"systematics"].Draw("a2")
            tgraphs[tag+"systematics"].Draw("p")
            legend.AddEntry(tgraphs[tag+"systematics"] ,"syst.","f" )
            legend.Draw("same")
            tgraphs[tag+"systematics"].GetXaxis().SetTitle( 'p_{T} [GeV]' )
            tgraphs[tag+"systematics"].GetYaxis().SetTitle( 'Scale Factor' )
            tgraphs[tag+"systematics"].SetMaximum(1.4)
            tgraphs[tag+"systematics"].SetMinimum(0.7)
            can.Update()
            can.SaveAs(tag+"_systematics"+".png")
            outputFile.cd()
            tgraphs[tag+"systematics"].Write()
            
            # all together 
            tgraphs[tag].SetFillColor(3)          
            tgraphs[tag+"statError"].SetFillColor(2) 
            tgraphs[tag+"totalError"].SetFillColor(4)
            tgraphs[tag+"systematics"].SetFillColor(7)
            
            tgraphs[tag].SetFillStyle(3001)          
            tgraphs[tag+"statError"].SetFillStyle(3001)
            tgraphs[tag+"totalError"].SetFillStyle(3001)
            tgraphs[tag+"systematics"].SetFillStyle(3001)
            
            multiG[tag]=ROOT.TMultiGraph()
            multiG[tag].SetMaximum(1.4)
            multiG[tag].SetMinimum(0.8)
            multiG[tag].SetTitle(tag)
            multiG[tag].Add(tgraphs[tag+"totalError"] ,"pe")
            multiG[tag].Add(tgraphs[tag+"systematics"],"pe")
            multiG[tag].Add(tgraphs[tag+""]           ,"pe")
            multiG[tag].Add(tgraphs[tag+"statError"]  ,"pe")
        
            if key==tags[0]:    
                legend.AddEntry(tgraphs[tag+"totalError"],"Total uncer.","f" )
                legend.AddEntry(tgraphs[tag+"systematics"],"Sys. uncer.","f" )
                legend.AddEntry(tgraphs[tag+""]          ,"Non-bb bkgs sub","f" )
                legend.AddEntry(tgraphs[tag+"statError"] ,"Stat.","f" )
                multiG[tag].Draw("a4")
                multiG[tag].Draw("p")
                legend.Draw("same")
                can.Update()
                can.SaveAs(tag+"_summary"+".png")
                outputFile.cd()        
        
    outputFile.Close()
    sortedResults=results.keys()
    
    for result in sortedResults:
        if result[2]!="":
            print (result[0]+result[1]+result[2]+"up", ":",results[result][6] , "//", result[0]+result[1]+result[2]+"down", ":",results[result][7])
    
if __name__ == '__main__':
    main()


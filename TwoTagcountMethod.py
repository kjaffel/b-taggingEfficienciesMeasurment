import ROOT
import sys
import math
import csv

#ROOT.gROOT.ProcessLine('.L BTagCalibrationStandalone.cpp+')
#calib = ROOT.BTagCalibration("MyAlgo")

fileName=sys.argv[1]
print 'successfully opening file', fileName
tFile=ROOT.TFile.Open(fileName)
lowerPtBinEdges=[20,30,50,70,100,140,200,300,600] 

def mistagrate_uncer(baseName):
    tag=baseName.split("_")
    var=(tag[len(tag)-2] if (baseName.find("Hi")!=-1 or baseName.find("Lo")!=-1) else tag[len(tag)-1])
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
    
    if baseName.find("deepCSV")!=-1:
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

def ComputeTwoTagEfficieny(basename):
    for sample in samples:
        histName=basename+"/"+basename+"_"+sample
        print histName
        hist=tFile.Get(histName)
        if sample==samples[0]: 
            histBkg=hist.Clone("histBkg")
        else:
            histBkg.Add(hist)

    tot={}
    statError={}
    tot["2bmc"]=0   # all the events in which both jets are 2b jet flav   with 0 & 1 passing 
    tot["othermc"]=0 # 2 jets (bX,cX, ll) 1 at least is passing
    tot["data"]=0   ## all data
    tot["2bmc_2btag"]=0  ##  number of events in which both jets are b-jets and both are passing 
    tot["othermc_2btag"]=0 ##  2jets (bX,CX,ll) both are btagged 
    tot["data_2btag"]=0   ## data where we have 2b jets flav && both passing (btagged)

    for iBin in range(15):
        if iBin%4==1:
            tot["2bmc"]=tot["2bmc"]+histBkg.GetBinContent(iBin)
            if iBin>9:
                tot["2bmc_2btag"]=tot["2bmc_2btag"]+histBkg.GetBinContent(iBin)
        else:
            tot["othermc"]=tot["othermc"]+histBkg.GetBinContent(iBin)
            if iBin>9: 
                tot["othermc_2btag"]=tot["othermc_2btag"]+histBkg.GetBinContent(iBin)
    twoTagDataHist=basename+"/"+basename
    hist=tFile.Get(twoTagDataHist)
    print basename 
    for iBin in range(15):
        print iBin,hist.GetBinContent(iBin),histBkg.GetBinContent(iBin)
        tot["data"]=tot["data"]+hist.GetBinContent(iBin)
        if iBin>9:
            tot["data_2btag"]=tot["data_2btag"]+hist.GetBinContent(iBin)

    if (histBkg.Integral()>0):
        hist.Integral()/histBkg.Integral()
    else:
        pass

    for name in tot:
        if  name.find("mc")!=-1:
            tot[name]=tot[name]/histBkg.Integral()
        else:
            tot[name]=tot[name]/hist.Integral()
    try:
        eff=math.sqrt((tot["data_2btag"]-tot["othermc_2btag"])/(tot["2bmc"]))
        effmc=math.sqrt(tot["2bmc_2btag"]/(tot["2bmc"]))
        
        #effUp=math.sqrt((tot["data_2btag"]-tot["othermc_2btag"]*(1-mistagrate_uncer(basename)))/tot["2bmc"])
        #effDown=math.sqrt((tot["data_2btag"]-tot["othermc_2btag"]*(1+mistagrate_uncer(basename)))/tot["2bmc"])
        
        # 50 % conservative uncertainty 
        effUp=math.sqrt((tot["data_2btag"]-tot["othermc_2btag"]*1.5)/tot["2bmc"])
        effDown=math.sqrt((tot["data_2btag"]-tot["othermc_2btag"]*0.5)/tot["2bmc"])
        
        statError[basename.replace("emu_","")+"data_2btag"] = math.sqrt( tot["data_2btag"] * (1 - tot["data_2btag"]) / hist.Integral())
        statError[basename.replace("emu_","")+"othermc_2btag"] = math.sqrt( tot["othermc_2btag"] * (1 - tot["othermc_2btag"]) / histBkg.GetEffectiveEntries())
        statError[basename.replace("emu_","")+"2bmc"] = math.sqrt( tot["2bmc"] * (1 - tot["2bmc"]) /  histBkg.GetEffectiveEntries())
        statError[basename.replace("emu_","")+"eff"] = math.sqrt(math.pow(eff*statError[basename.replace("emu_","")+"2bmc"],2)  + (math.pow(statError[basename.replace("emu_","")+"data_2btag"],2)+math.pow(statError[basename.replace("emu_","")+"othermc_2btag"],2)  ) / math.pow(eff,2)    )/(2*tot["2bmc"])
        print [tot["2bmc_2btag"],tot["2bmc"],tot["data_2btag"],tot["othermc_2btag"],eff,effmc,eff/effmc,effUp/effmc,effDown/effmc, statError[basename.replace("emu_","")+"eff"]]
        print "basename:", basename,  "eff:",eff, "effmc:", effmc , "effUp:", effUp, "effDown:", effDown
        return [tot["2bmc_2btag"],tot["2bmc"],tot["data_2btag"],tot["othermc_2btag"],eff,effmc,eff/effmc,effUp/effmc,effDown/effmc, statError[basename.replace("emu_","")+"eff"]]
    except:
        print "ERROR: ERROR in computing efficiencies from",basename

#tags=["twoTags_deepCSVL","only2_twoTags_deepCSVL","twoTags_deepCSVM","only2_twoTags_deepCSVM","twoTags_deepCSVT","only2_twoTags_deepCSVT","twoTags_deepFlavourL","only2_twoTags_deepFlavourL","twoTags_deepFlavourM","only2_twoTags_deepFlavourM","twoTags_deepFlavourT","only2_twoTags_deepFlavourT"]
tags=["only2_twoTags_deepCSVL","only2_twoTags_deepCSVM","only2_twoTags_deepCSVT","only2_twoTags_deepFlavourL","only2_twoTags_deepFlavourM","only2_twoTags_deepFlavourT"]
tags.sort()

samples=["t#bar{t} SL","t#bar{t} DL","t#bar{t} AH","ZZ","WZ","WW","tW","DY","W2Jets", "W3Jets", "W4Jets"]
systs=["","_fsrCon","_fsrDef","_fsrRed","_isrCon","_isrDef","_isrRed","_qcdScale"]

hiORlo=["Hi","Lo"]

tgraphs={}
results={}

can=ROOT.TCanvas("can","",800,800)
outputFile=ROOT.TFile("twoTagSFGraphs.root","RECREATE")

taggers=['DeepFlavour','DeepCSV']
#for tagger in taggers:
#    with open('{tagger}_94XSF_Moriond19_2018DATA_syst_stat.csv'.format(tagger=tagger), mode='a') as csv_file:
#        btag_params = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
#        btag_params.writerow([tagger,'OperatingPoint', 'measurementType', 'sysType', 'jetFlavor', 'etaMin','etaMax', 'ptMin','ptMax', 'discrMin','discrMax', 'SFs'])

#compute and save all the results
for syst in systs:
    for tag in tags:
        for ptBin in range(-1,len(lowerPtBinEdges)):
            ptLabel=ReturnPtLabel(ptBin)
            analysiskeys= (tag,ptLabel,syst)

            print "---------------------------------------------------------------------------------------"
            print "analysiskeys:" , analysiskeys
            basename="emu_"+tag+ptLabel+syst
            if syst=="":
                results[analysiskeys]=ComputeTwoTagEfficieny(basename)
            else:
                print basename+"Hi"
                resultsHi=ComputeTwoTagEfficieny(basename+"Hi")
                resultsLo=ComputeTwoTagEfficieny(basename+"Lo")
                statErrorHi=ComputeTwoTagEfficieny(basename+"Hi")
                statErrorLo=ComputeTwoTagEfficieny(basename+"Lo")
                #I just set the numbers to 1 because you don't really need them for systematics unless you're debugging
                #return[tot["2bmc_2btag"],tot["2bmc"],tot["data_2btag"],tot["othermc_2btag"],eff,effmc,eff/effmc,effUp/effmc,effDown/effmc]
                results[analysiskeys]=[1.,1.,1.,1.,1.,1.,(resultsHi[6] if resultsHi else 1.), (resultsLo[6] if resultsLo else 1.)]
                print results[analysiskeys]
#make the tgraphs for all systematics
multiG={}
multiG2={}
for index, tag in enumerate(tags):
    uncertainty2Hi=[]
    uncertainty2Lo=[]
    otheruncertainty2Hi=[]
    otheruncertainty2Lo=[]
    
    it=0
    if index==0:
        tgraphs["wp"]=ROOT.TGraphAsymmErrors()
        tgraphs["wp"].SetTitle(";Jet P_{T};Scale Factor")
    for syst in systs:
        iGraph=0
        graphKey=tag+syst
        tgraphs[graphKey]=ROOT.TGraphAsymmErrors()
        tgraphs[graphKey].SetTitle(graphKey+";Jet P_{T};Scale Factor")
        if syst==systs[0]:
            tgraphs[tag+"totalError"]=ROOT.TGraphAsymmErrors()
            tgraphs[tag+"othersysts"]=ROOT.TGraphAsymmErrors()
            tgraphs[tag+"statError"]=ROOT.TGraphErrors()
        #for ptBin in range(-1,len(lowerPtBinEdges)): -1 for Inclusive
        for ptBin in range(0,len(lowerPtBinEdges)):
            ptLabel=ReturnPtLabel(ptBin)
            # nominal must always be run for this to work
            if ptBin>-1:
                pt1=lowerPtBinEdges[ptBin]
                pt2=1000
                if pt1!=lowerPtBinEdges[-1]:
                    pt2=lowerPtBinEdges[ptBin+1]
            try:
                nominalSF=results[(tag,ptLabel,"")][6]
                nominalSFUp=results[(tag,ptLabel,"")][7]
                nominalSFDown=results[(tag,ptLabel,"")][8]
                statError=results[(tag,ptLabel,"")][9]
                if syst=="":
                    tgraphs[tag+"othersysts"].SetPoint(iGraph,(pt1+pt2)/2.,nominalSF)
                    tgraphs[tag+"othersysts"].SetPointEXhigh(iGraph,(-pt1+pt2)/2.)
                    tgraphs[tag+"othersysts"].SetPointEXlow(iGraph,(-pt1+pt2)/2.)
    
                    tgraphs[tag+"statError"].SetPoint(iGraph,(pt1+pt2)/2, nominalSF)
                    tgraphs[tag+"statError"].SetPointError(iGraph,(pt2-pt1)/2, statError)

                    uncertainty2Hi.append(math.pow(nominalSFUp-nominalSF,2))
                    uncertainty2Lo.append(math.pow(-nominalSFDown+nominalSF,2))
    
                    otheruncertainty2Hi.append(0)
                    otheruncertainty2Lo.append(0)
    
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
                    print "--------------------------------------------------------------------------"
                    print ("Nominal:{0}".format(nominalSF))
                    print ("Up     :{0}".format(nominalSFUp))
                    print ("Down   :{0}".format(nominalSFDown))
                    iGraph=iGraph+1
                    if index>=3:
                        it=it+1

                    wp= (0 if tag.find("L")==len(tag)-1 else (1 if tag.find("M")==len(tag)-1 else 2))
                    meas="twotag"
                    var="Central"
                    flav=0
                    etaMin=-2.4
                    etaMax=2.4
                    ptMin=(ptLabel.split("to"))[0]
                    ptMax=(ptLabel.split("to"))[1]
                    if ptMax==str("Inf"):
                        ptMax=ptMax.replace(str("Inf"),str(1000))
                    discrMin=0
                    discrMax=1
                    SFs=nominalSF
                
                    if index >= 3:
                        tagger="DeepFlavour"
                    else:
                        tagger="DeepCSV"

                    with open('{tagger}_94XSF_Moriond19_2018DATA_syst_stat.csv'.format(tagger=tagger), mode='a') as csv_file:
                        btag_params = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                        btag_params.writerow([wp, meas, var, flav, etaMin, etaMax, ptMin,ptMax ,discrMin,discrMax, SFs])
                
                else:
                    systHiSF=results[(tag,ptLabel,syst)][6]
                    systLoSF=results[(tag,ptLabel,syst)][7]
                    tgraphs[tag+"totalError"].SetPoint(iGraph,(pt1+pt2)/2.,nominalSF)
                    tgraphs[tag+"totalError"].SetPointEXhigh(iGraph,(-pt1+pt2)/2.)
                    tgraphs[tag+"totalError"].SetPointEXlow(iGraph,(-pt1+pt2)/2.)

                    if (systHiSF-systLoSF >0.):
                        tgraphs[graphKey].SetPointEYhigh(iGraph,systHiSF-nominalSF)
                        tgraphs[graphKey].SetPointEYlow( iGraph,-systLoSF+nominalSF)
    
                        
                        uncertainty2Hi[ptBin] += math.pow(systHiSF-nominalSF,2)
                        uncertainty2Lo[ptBin] += math.pow(-systLoSF+nominalSF,2)
                        otheruncertainty2Hi[ptBin] += math.pow(systHiSF-nominalSF,2)
                        otheruncertainty2Lo[ptBin] += math.pow(-systLoSF+nominalSF,2)
                    else:
                        tgraphs[graphKey].SetPointEYlow( iGraph,-systHiSF+nominalSF)
                        tgraphs[graphKey].SetPointEYhigh(iGraph,systLoSF-nominalSF)
                    
                        uncertainty2Hi[ptBin] += math.pow(-systHiSF+nominalSF,2)
                        uncertainty2Lo[ptBin] += math.pow(systLoSF-nominalSF,2)
                        otheruncertainty2Hi[ptBin] += math.pow(-systHiSF+nominalSF,2)
                        otheruncertainty2Lo[ptBin] += math.pow(systLoSF-nominalSF,2)

                    tgraphs[graphKey].SetPoint(iGraph,(pt1+pt2)/2.,nominalSF)
                    tgraphs[graphKey].SetPointEXhigh(iGraph,(-pt1+pt2)/2.)
                    tgraphs[graphKey].SetPointEXlow(iGraph,(pt2-pt1)/2.)
                    iGraph=iGraph+1

                    print "-------------------------------------------------------------------------------"
                    print("High   : {0}".format(systHiSF))
                    print("Low    : {0}".format(systLoSF))
                    print("High syst -Nominal:{0}".format(systHiSF-nominalSF))
                    print("Nominal- Low syst :{0}".format(-systLoSF+nominalSF))

                    if math.copysign(1,systHiSF-nominalSF) != math.copysign(1,nominalSF-systLoSF):
                        print "Hi and Lo vary from nominal in the same direction",graphKey,"hi, nom, lo",systHiSF,nominalSF,systLoSF                 

            except:
                print "skip plotting the default value SF=1 ! "
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
            tgraphs[tag+"totalError"].SetPointEYhigh(iGraph, math.sqrt(uncertainty2Hi[iGraph]))
            tgraphs[tag+"totalError"].SetPointEYlow(iGraph, math.sqrt(uncertainty2Lo[iGraph]))
            wp= (0 if tag.find("L")==len(tag)-1 else (1 if tag.find("M")==len(tag)-1 else 2))
            meas="twotag"
            SFsUp=allSF[iGraph]+ math.sqrt(uncertainty2Hi[iGraph])
            SFsDown=allSF[iGraph]-math.sqrt(uncertainty2Lo[iGraph])
            flav=0
            etaMin=-2.4
            etaMax=2.4
            ptLabel=ReturnPtLabel(iGraph)
            ptMin=(ptLabel.split("to"))[0]
            ptMax=(ptLabel.split("to"))[1]
            if ptMax==str("Inf"):
                ptMax=ptMax.replace(str("Inf"),str(1000))
                discrMin=0
                discrMax=1
        
            if index >= 3:
                tagger="DeepFlavour"
            else:
                tagger="DeepCSV"

            with open('{tagger}_94XSF_Moriond19_2018DATA_syst_stat.csv'.format(tagger=tagger), mode='a') as csv_file:
                btag_params = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                btag_params.writerow([wp, meas, "Up", flav, etaMin, etaMax, ptMin,ptMax ,discrMin,discrMax, SFsUp])
                btag_params.writerow([wp, meas, "Down", flav, etaMin, etaMax, ptMin,ptMax ,discrMin,discrMax, SFsDown])

        otherSF=tgraphs[tag+"othersysts"].GetY()
        for iGraph in range(0,len(lowerPtBinEdges)):
            tgraphs[tag+"othersysts"].SetPointEYhigh(iGraph, math.sqrt(otheruncertainty2Hi[iGraph]))
            tgraphs[tag+"othersysts"].SetPointEYlow(iGraph, math.sqrt(otheruncertainty2Lo[iGraph]))
            wp= (0 if tag.find("L")==len(tag)-1 else (1 if tag.find("M")==len(tag)-1 else 2))
            meas="twotag"
            SFsUp=otherSF[iGraph]+ math.sqrt(otheruncertainty2Hi[iGraph])
            SFsDown=otherSF[iGraph]-math.sqrt(otheruncertainty2Lo[iGraph])
            flav=0
            etaMin=-2.4
            etaMax=2.4
            ptLabel=ReturnPtLabel(iGraph)
            ptMin=(ptLabel.split("to"))[0]
            ptMax=(ptLabel.split("to"))[1]
            if ptMax==str("Inf"):
                ptMax=ptMax.replace(str("Inf"),str(1000))
                discrMin=0
                discrMax=1

            if index >= 3:
                tagger="DeepFlavour"
            else:
                tagger="DeepCSV"
    
            with open('{tagger}_94XSF_Moriond19_2018DATA_othersyst.csv'.format(tagger=tagger), mode='a') as csv_file:
                csv_file.truncate()
                btag_params = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                btag_params.writerow([wp, meas, "Up", flav, etaMin, etaMax, ptMin,ptMax ,discrMin,discrMax, SFsUp])
                btag_params.writerow([wp, meas, "Down", flav, etaMin, etaMax, ptMin,ptMax ,discrMin,discrMax, SFsDown])

    except:
        print ' ** Warning ***** unable to dump data in {tagger}_94XSF_Moriond19_2018DATA_othersyst.csv'.format(tagger=tagger)
    # save the statError , othersystematics and totalError once per tag
    tgraphs[tag+"statError"].SetName(tag+"statError")
    tgraphs[tag+"statError"].Draw("APE")
    tgraphs[tag+"statError"].GetXaxis().SetTitle( 'pT(GeV)' )
    tgraphs[tag+"statError"].GetYaxis().SetTitle( 'Scale Factor' )
    tgraphs[tag+"statError"].GetYaxis().SetLimits( 0,7)
    can.Update()
    can.SaveAs(tag+"_statError"+".png")
    outputFile.cd()
    tgraphs[tag+"statError"].Write()

    tgraphs[tag+"totalError"].SetName(tag+"totalError")
    tgraphs[tag+"totalError"].Draw("APE")
    tgraphs[tag+"totalError"].GetXaxis().SetTitle( 'pT(GeV)' )
    tgraphs[tag+"totalError"].GetYaxis().SetTitle( 'Scale Factor' )
    tgraphs[tag+"totalError"].GetYaxis().SetLimits( 0,7)
    can.Update()
    can.SaveAs(tag+"_totalError"+".png")
    outputFile.cd()
    tgraphs[tag+"totalError"].Write()
    
    tgraphs[tag+"othersysts"].SetName(tag+"othersysts")
    tgraphs[tag+"othersysts"].Draw("APE")
    tgraphs[tag+"othersysts"].GetXaxis().SetTitle( 'pT(GeV)' )
    tgraphs[tag+"othersysts"].GetYaxis().SetTitle( 'Scale Factor' )
    tgraphs[tag+"othersysts"].GetYaxis().SetLimits( 0,7)
    can.Update()
    can.SaveAs(tag+"_othersysts"+".png")
    outputFile.cd()
    tgraphs[tag+"othersysts"].Write()

    #tgraphs[tag+"totalError"].SetFillColor(6)
    #tgraphs[tag+"totalError"].SetFillStyle(3005)
    # set Marker style
    tgraphs[tag+"totalError"].SetMarkerStyle(21)
    tgraphs[tag+"othersysts"].SetMarkerStyle(22)
    tgraphs[tag].SetMarkerStyle(23)          
    tgraphs[tag+"statError"].SetMarkerStyle(24) 
    # set the line color 
    tgraphs[tag+"totalError"].SetLineColor(1)
    tgraphs[tag+"othersysts"].SetLineColor(633)
    tgraphs[tag].SetLineColor(417)          
    tgraphs[tag+"statError"].SetLineColor(601) 
    # set the line width
    tgraphs[tag+"totalError"].SetLineWidth(5)
    tgraphs[tag+"othersysts"].SetLineWidth(3)
    tgraphs[tag].SetLineWidth(2)
    tgraphs[tag+"statError"].SetLineWidth(1)

    multiG[tag]=ROOT.TMultiGraph()
    multiG[tag].SetMaximum(10)
    multiG[tag].SetMinimum(0)
    multiG[tag].Add(tgraphs[tag+"totalError"],"pe")
    multiG[tag].Add(tgraphs[tag+"othersysts"],"pe")
    multiG[tag].Add(tgraphs[tag+""]          ,"pe")
    multiG[tag].Add(tgraphs[tag+"statError"] ,"pe")

    if tag==tags[0]:    
        leg=ROOT.TLegend(0.5,0.7,0.9,0.9)
        leg.AddEntry(tgraphs[tag+"totalError"],"Total unc.","l" )
        leg.AddEntry(tgraphs[tag+"othersysts"],"systs:isr/fsr & qcdScale","l" )
        leg.AddEntry(tgraphs[tag+""]          ,"Non-bb bkg sub","l" )
        leg.AddEntry(tgraphs[tag+"statError"] ,"Statistics","l" )
    multiG[tag].Draw("a4")
    leg.Draw("same")
    #multiG[tag].GetYaxis().SetLimits( 0, 7)
    #multiG[tag].GetXaxis().SetTitle( 'pT(GeV)' )
    #multiG[tag].GetYaxis().SetTitle( 'Scale Factor' )
    can.Update()
    can.SaveAs(tag+"_summary"+".png")
    outputFile.cd()        

outputFile.Close()
sortedResults=results.keys()
sortedResults.sort()

for result in sortedResults:
    if result[2]!="":
        print result[0]+result[1]+result[2]+"Hi", ":",results[result][6] , "//", result[0]+result[1]+result[2]+"Lo", ":",results[result][7]

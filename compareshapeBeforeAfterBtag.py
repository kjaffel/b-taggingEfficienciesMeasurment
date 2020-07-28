import os
import sys
import glob
import argparse
import copy
import yaml
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

class StudyBtagRegions():
    def __init__(self,variable,variable_name,output,plot_data,plot_MC):
        self.outname = "Plots_2017UL/"+output
        self.era     = '2017'
        self.variable = variable
        self.variable_name = variable_name
        self.data_hist_dict = {}
        self.MC_hist_dict = {}
        self.plot_data = plot_data
        self.plot_MC = plot_MC

        self.extractAllInformation()
        self.processHistograms()
        self.saveHistograms()

    def extractAllInformation(self):
        self.info_dict = {}
        for (directory,legend),histogram in self.variable.items():
            directory = os.path.abspath(directory)
            self.info_dict[legend] = self.extractDirInformation(directory,histogram)
    
    def extractDirInformation(self,directory,histogram):
        # Get YAML info #
        yaml_dict = self.loadYaml(os.path.join(directory,'plots.yml'),histogram)

        # Loop over file to recover histograms #
        hist_dict = {}
        path_file = os.path.abspath(os.path.join(directory, 'results','*root'))
        for f in glob.glob(path_file):
            name = os.path.basename(f)
            #print ("Looking at %s"%name)
            # Check if in YAML #
            if name not in yaml_dict['samples'].keys():
                print ("[WARNING] File %s will be ignored"%(name))
                continue
            
            h = self.getHistogram(f,histogram)
            hist_dict[name] =h 

        #print ('hist_dict', hist_dict)
        yaml_dict.update({'histograms':hist_dict})

        return yaml_dict


    def loadYaml(self,yaml_path,histogram):
        # Parse YAML #
        with open(yaml_path,"r") as handle:
            full_dict = yaml.load(handle,Loader=yaml.FullLoader)
        # Get Lumi per era #
        lumi_dict = full_dict["configuration"]["luminosity"]
        
        # Get plot options #
        opt_to_keep = ['x-axis','y-axis']
        try:
            options_dict = {k:full_dict["plots"][histogram][k] for k in full_dict["plots"][histogram].keys() & opt_to_keep}
        except KeyError:
            print ("Could not find hist %s in YAML %s, will proceed without"%(histogram,yaml_path))
            options_dict = {k:'' for k in opt_to_keep}

        # Get data per sample #
        sample_dict = {}
        info_to_keep = ['cross-section','generated-events','group','type','era']
        for sample,data in full_dict['files'].items():
            sample_dict[sample] = {k:data[k] for k in data.keys() & info_to_keep}

        return {'luminosity':lumi_dict,'plot_options':options_dict,'samples':sample_dict}
            
    def getHistogram(self,rootfile,histname):
        f = ROOT.TFile(rootfile)
        #print (rootfile)
        if f.GetListOfKeys().Contains(histname):
            h = copy.deepcopy(f.Get(histname))
        else:
            print ("Could not find hist %s in %s"%(histname,rootfile))
            h = None
        f.Close()
        return h

    def processHistograms(self):
        colors = [221, 213, 2, ROOT.kRed, 9, 69, 90, 208] 
        import collections 
        
        for i,(key, val) in enumerate(self.info_dict.items()):
            print ('key:', key)
            #samp_dict=collections.defaultdict(dict)
            #for smp, val2 in val['samples'].items():
                
            #    if 'MuMu' in self.variable_name:
            #        if not str(smp).startswith("DoubleMuon") and not str(smp).startswith("DY"):
            #            continue
            #    elif 'ElEl' in self.variable_name:
            #        if not str(smp).startswith("DoubleEGamma") and not str(smp).startswith("DY"):
            #            continue
            #    else:
            #        if not str(smp).startswith("DoubleEGamma") and not str(smp).startswith("DoubleMuon") and not str(smp).startswith("DY"):
            #            continue
            #    samp_dict[smp]= val2
            
            data_hist = None
            MC_hist = None
            hist_dict = val['histograms'] 
            lumi_dict = val['luminosity'] 
            plot_dict = val['plot_options'] 
            samp_dict = val['samples']
            
            #print ( samp_dict)
            for sample, data_dict in samp_dict.items():
                h = hist_dict[sample]
                if h is None:
                    continue
                if data_hist is None:
                    data_hist = ROOT.TH1F(self.variable_name,self.variable_name,h.GetNbinsX(),h.GetBinLowEdge(1),h.GetBinLowEdge(h.GetNbinsX())+h.GetBinWidth(h.GetNbinsX()))
                if MC_hist is None:
                    MC_hist = ROOT.TH1F(self.variable_name,self.variable_name,h.GetNbinsX(),h.GetBinLowEdge(1),h.GetBinLowEdge(h.GetNbinsX())+h.GetBinWidth(h.GetNbinsX()))

                if data_dict['type'] == 'data':
                    data_hist.Add(h)
                elif data_dict['type'] == 'mc':
                    # Let's normalize to uni just to compare the shape
                    #factor = data_dict["cross-section"]*lumi_dict[self.era]/data_dict["generated-events"]
                    #MC_hist.Add(h,factor)
                    MC_hist.Add(h)
                        

                    
            data_hist.SetLineWidth(2)
            #data_hist.SetFillColorAlpha(colors[i], 0.35)
            data_hist.SetLineColor(colors[i])
            data_hist.GetXaxis().SetTitle(plot_dict['x-axis'])
            data_hist.GetYaxis().SetTitle(plot_dict['y-axis'])
            data_hist.SetTitle(self.variable_name)
            MC_hist.SetLineWidth(2)
            if self.plot_MC and self.plot_data:
                i *= 2
                MC_hist.SetLineColor(colors[i+1])
                print( "color", colors[i+1])
            else:
                MC_hist.SetLineColor(colors[i])
            MC_hist.GetXaxis().SetTitle(plot_dict['x-axis'])
            MC_hist.GetYaxis().SetTitle(plot_dict['y-axis'])
            MC_hist.SetTitle(self.variable_name)
            # Normalize to unity #
            if data_hist.Integral() != 0:
                data_hist.Scale(1/data_hist.Integral())
            if MC_hist.Integral() != 0:
                MC_hist.Scale(1/MC_hist.Integral())
            # Save #
            self.data_hist_dict[key] = data_hist
            self.MC_hist_dict[key] = MC_hist

    def saveHistograms(self):
        num_plots = len(self.data_hist_dict.keys())*self.plot_data + len(self.MC_hist_dict.keys())*self.plot_MC
        plot_ratio = num_plots==2
        #if num_plots>=4:
        #    legend = ROOT.TLegend(0.2,0.8,0.89,0.89)
        #    legend.SetNColumns(2)
        #    legend.SetTextSize(0.012)
        #else:
        legend = ROOT.TLegend(0.5,0.8,0.89,0.89)
        legend.SetTextSize(0.015)
        legend.SetHeader("N- BJets","C")

        # Canvas and pad #
        C = ROOT.TCanvas("c1", "c1", 600, 600)
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.0, 1, 1.0)
        pad1.SetBottomMargin(0.15)
        pad1.SetLeftMargin(0.15)
        pad1.SetRightMargin(0.1)
        if plot_ratio:
            pad1.SetBottomMargin(0.32)
        #pad1.SetGridx()
        #pad1.SetGridy()
        pad1.Draw()
        pad1.cd()

        # Get Max values #
        max_data = max([h.GetMaximum() for h in self.data_hist_dict.values()])
        max_MC = max([h.GetMaximum() for h in self.MC_hist_dict.values()])
        if self.plot_data and self.plot_MC:
            amax = max(max_data,max_MC)
        elif self.plot_data:
            amax = max_data
        elif self.plot_MC:
            amax = max_MC

        opt = "hist"
        if self.plot_data:
            print( self.data_hist_dict.keys())
            getOnekey = list(self.data_hist_dict)[0]
            print ( getOnekey)
            hist_data = self.data_hist_dict[getOnekey]
            hist_data.SetMaximum(amax*1.2)
            hist_data.SetMinimum(0.)
            hist_data.Draw(opt)
            if opt.find("same") == -1:
                opt += " same"
            legend.AddEntry(hist_data, "Data",'l')
        
        for key in self.MC_hist_dict.keys():
            if self.plot_MC:
                hist_MC = self.MC_hist_dict[key]
                hist_MC.SetMaximum(amax*1.2)
                hist_MC.SetMinimum(0.)
                hist_MC.Draw(opt)
                if opt.find("same") == -1:
                    opt += " same"
            if self.plot_MC:
                legend.AddEntry(hist_MC,"%s "%key,'l')

        legend.Draw()
        if plot_ratio:
            keys = list(self.data_hist_dict.keys())
            if self.plot_data and not self.plot_MC:
                hist1 = self.data_hist_dict[keys[0]]
                hist2 = self.data_hist_dict[keys[1]]
            elif self.plot_MC and not self.plot_data:
                hist1 = self.MC_hist_dict[keys[0]]
                hist2 = self.MC_hist_dict[keys[1]]
            elif self.plot_data and self.plot_MC:
                hist1 = self.data_hist_dict[keys[0]]
                hist2 = self.MC_hist_dict[keys[0]]
            ratio = hist1.Clone()
            ratio.Sumw2()
            ratio.Divide(hist2)

            # Redraw axis to avoid clipping 0
            hist1.GetXaxis().SetLabelSize(0.)
            hist1.GetXaxis().SetTitle('')
            
            pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.3)
            pad2.SetTopMargin(0)
            pad2.SetBottomMargin(0.4)
            pad2.SetLeftMargin(0.15)
            pad1.SetRightMargin(0.1)
            #pad2.SetGridx()
            pad2.SetGridy()
            pad2.Draw()
            pad2.cd()

            ratio.SetLineColor(ROOT.kBlack)
            ratio.SetMinimum(0.6)
            ratio.SetMaximum(1.4)
            ratio.SetStats(0)
            ratio.Draw("ep")

            ratio.SetTitle("")
            ratio.GetYaxis().SetTitle("Data/Mc")
            ratio.GetYaxis().SetNdivisions(505)
            ratio.GetYaxis().SetTitleSize(20)
            ratio.GetYaxis().SetTitleFont(43)
            ratio.GetYaxis().SetTitleOffset(1.8)
            ratio.GetYaxis().SetLabelFont(43)
            ratio.GetYaxis().SetLabelSize(15)
            
            ratio.GetXaxis().SetNdivisions(510)
            ratio.GetXaxis().SetTitleSize(20)
            ratio.GetXaxis().SetTitleFont(43)
            ratio.GetXaxis().SetTitleOffset(4.)
            ratio.GetXaxis().SetLabelFont(43)
            ratio.GetXaxis().SetLabelSize(15)

        C.Print(self.outname+".pdf")
        C.Print(self.outname+".png")
            


if __name__ == "__main__":
    
    path = '/home/ucl/cp3/kjaffel/bamboodev/b-taggingEfficienciesMeasurment/newversion/2017UL/controlplots.v9._/'

    for var_ToPlots in ['Jet_mulmtiplicity']:
        for cat in ['elmu', 'muel']:
            
            variable_name = '{0} channel : {1}'.format(cat, var_ToPlots)
            if var_ToPlots == 'Jet_mulmtiplicity':
                no_btag = {(path,'no btag'):  '{0}_atleast2Jets_Jet_mulmtiplicity'.format(cat)}
                for tagger in  ['DeepCSV', 'DeepFlavour']:
                    loOse_btag = {(path,'%s Loose'%tagger): '{0}_noMETcut_atleast2j_{1}L_{2}'.format(cat, tagger, var_ToPlots)}
                    medium_btag = {(path,'%s Meidum'%tagger): '{0}_noMETcut_atleast2j_{1}M_{2}'.format(cat, tagger, var_ToPlots)}
                    tight_btag = {(path,'%s Tight'%tagger): '{0}_noMETcut_atleast2j_{1}T_{2}'.format(cat, tagger, var_ToPlots)}
                    
                    instance = StudyBtagRegions({**no_btag, **loOse_btag, **medium_btag, **tight_btag},
                                            variable_name,"{0}_comparebtag_{1}_{2}".format(cat, tagger, var_ToPlots),
                                            plot_data=True, plot_MC=True)

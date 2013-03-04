'''

Base class to do WH plotting.

Author: Evan K. Friis, UW

Takes as input a set of ROOT [files] with analysis histgrams, and the corresponding
lumicalc.sum [lumifiles] that hve the effective lumi for each sample.

If [blind] is true, data in the p1p2p3 region will not be plotted.

'''

import rootpy.plotting.views as views
import rootpy.plotting as plotting
#from rootpy.tree import TreeChain
import rootpy.io
from FinalStateAnalysis.PlotTools.Plotter import Plotter
from FinalStateAnalysis.PlotTools.PoissonView import PoissonView
from FinalStateAnalysis.PlotTools.HistToTGRaphErrors import HistToTGRaphErrors, HistStackToTGRaphErrors
from FinalStateAnalysis.PlotTools.InflateErrorView import InflateErrorView
from FinalStateAnalysis.MetaData.data_styles import data_styles
#from RecoLuminosity.LumiDB import argparse''' 
import sys
import os
import glob
import pprint
import ROOT
import fnmatch
import math

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)

def quad(*xs):
    return math.sqrt(sum(x*x for x in xs))


class TauEffPlotter(Plotter):
    def __init__(self, channel):
        self.samples = [ 'Zjets*', 'WplusJets*', 'TTplusJets*']#'WH*',
        mcSamples    = [ 'Zjets*', 'WplusJets*', 'TTplusJets*']
        self.samples += ['data_SingleMu*']
        self.jobid = os.environ['jobid']
        self.channel = channel
        self.period = '7TeV' if '7TeV' in jobid else '8TeV'
        self.sqrts = 7 if '7TeV' in jobid else 8
        files = []
        lumifiles = []
        for x in self.samples:
            files += glob.glob('results/%s/TauEffZ%s/%s.root' % (jobid, channel, x))
            lumifiles += glob.glob('inputs/%s/%s.lumicalc.sum' % (jobid, x))
        self.outputdir    = 'results/%s/plots/%s' % (jobid, channel.lower() )
        self.base_out_dir = self.outputdir
        #pprint.pprint(files)
        super(TauEffPlotter, self).__init__(files, lumifiles, self.outputdir, None)
        self.mc_samples = mcSamples

    def set_subdir(self, folder):
        self.outputdir = '/'.join([self.base_out_dir, folder])

    def plot_mc_vs_data(self, folder, variable, rebin=1, xaxis='', leftside=True, xrange=None):
        #print folder
        super(TauEffPlotter, self).plot_mc_vs_data(folder, variable, rebin, xaxis, leftside, xrange)
        self.add_cms_blurb(self.sqrts)

    def get_w_mt_mc_scale(self, iso_name):
        mt_h_name= 'mMtToMET' if self.channel == 'MT' else 'm1MtToMET'
        wjets_mc = self.get_view('WplusJets*')
        wjets_hi = wjets_mc.Get('%s/ss/HiMT/bjets0/%s' % (iso_name,mt_h_name))
        wjets_lo = wjets_mc.Get('%s/ss/LoMT/bjets0/%s' % (iso_name,mt_h_name))
        return wjets_lo.Integral() / wjets_hi.Integral()

    def plot_with_estimate(self, folder, variable, rebin=1, xaxis='', leftside=True, xrange=None, show_error=True):
        iso_name  = folder.split('/')[0]
        wjets_mc  = self.rebin_view(self.get_view('WplusJets*') , rebin)
        zjets_mc  = self.rebin_view(self.get_view('Zjets*')     , rebin)
        data      = self.rebin_view(self.get_view('data')       , rebin)
        ttbar     = self.rebin_view(self.get_view('TTplusJets*'), rebin)
        ttbar     = views.SubdirectoryView(ttbar, folder)
        neg_zjets = views.ScaleView(zjets_mc, -1)
        zjets_mc  = views.SubdirectoryView(zjets_mc, folder)

        #makes WJets estimation view
        wjet_reg  = '%s/ss/HiMT/bjets0/' % iso_name
        #print iso_name
        nzj_at_w  = views.SubdirectoryView(neg_zjets, wjet_reg)
        data_at_w = views.SubdirectoryView(data, wjet_reg)
        w_scale   = self.get_w_mt_mc_scale(iso_name)
        wjet_est  = views.SumView(data_at_w, nzj_at_w)
        wjet_est  = views.ScaleView(wjet_est, w_scale)
        wjet_est  = views.StyleView(wjet_est, **data_styles['WplusJets*'])
        wjet_est  = views.TitleView(wjet_est, 'WplusJets')
        
        #makes QCD Estimation view
        qcd_reg   = '%s/ss/LoMT/bjets0/' % iso_name
        nzj_qcd   = views.SubdirectoryView(neg_zjets, qcd_reg)
        data_qcd  = views.SubdirectoryView(data, qcd_reg)
        nwjet_est = views.ScaleView(wjet_est, -1)
        qcd_est   = views.SumView(data_qcd, nzj_qcd)
        qcd_est   = views.SumView(data_qcd, nwjet_est)
        qcd_est   = views.StyleView(qcd_est, **data_styles['QCD*'])
        qcd_est   = views.TitleView(qcd_est, 'QCD')
        
        stack     = views.StackView(zjets_mc,wjet_est,ttbar, qcd_est).Get(variable)
        stack.Draw()
        stack.GetHistogram().GetXaxis().SetTitle(xaxis)
        if xrange:
            stack.GetXaxis().SetRange(xrange[0], xrange[1])
            stack.Draw()
        self.keep.append(stack)

        #print os.path.join(folder,variable)
        data_h = data.Get(os.path.join(folder,variable))
        data_h.Draw('same')
        self.keep.append(data_h)
        # Make sure we can see everything
        if data_h.GetMaximum() > stack.GetMaximum():
            stack.SetMaximum(1.2*data_h.GetMaximum())
            
        if show_error:
            stack_sum = sum(stack.GetHists())
            stack_sum.SetMarkerSize(0)
            stack_sum.SetFillColor(1)
            stack_sum.SetFillStyle(3013)
            stack_sum.legendstyle = 'f'
            self.keep.append(stack_sum)
            stack_sum.Draw('pe2,same')
            
        # Add legend
        self.add_legend([data_h, stack], leftside, entries=5)



    def write_shapes(self, variable, rebin, outdir, unblinded=False):
        ''' Write final shape histos for [variable] into a TDirectory [outputdir] '''
        sig_view = self.make_signal_views(rebin, unblinded)
        outdir.cd()
        wz = sig_view['wz'].Get(variable)
        zz = sig_view['zz'].Get(variable)
        obs = sig_view['data'].Get(variable)
        Zjets = sig_view['Zjets'].Get(variable)

        wz.SetName('WZ')
        zz.SetName('ZZ')
        obs.SetName('data_obs')
        Zjets.SetName('Zjets')

        #print sig_view.keys()
        for mass in [110, 120, 130, 140]:
            vh = sig_view['vh%i' % mass].Get(variable)
            vh.SetName('ZH_htt%i' % mass)
            vh.Write()
            ww = sig_view['vh%i_hww' % mass].Get(variable)
            ww.SetName('ZH_hww%i' % mass)
            ww.Write()

        wz.Write()
        zz.Write()
        obs.Write()
        Zjets.Write()



        
if __name__ <> "__main__":
    sys.exit(0)
jobid = os.environ['jobid']

leptons = {
    'mt' : {
        1     : { 'name' : '#mu',      'label' : 'm'},
        2     : { 'name' : '#tau_{h}', 'label' : 't'},
        'ids' : ['LooseIso', 'MediumIso', 'TightIso', 'LooseMVAIso', 'MediumMVAIso', 'TightMVAIso',],
        },
    'mm' : {
        1     : { 'name' : '#mu_{1}', 'label' : 'm1'},
        2     : { 'name' : '#mu_{2}', 'label' : 'm2'},
        'ids' : ['h2Tau'],
        },
    }

channels = ['mm', 'mt']

for channel in channels:
    leps    = leptons[channel]
    name1   = leps[1]['name']
    label1  = leps[1]['label']
    name2   = leps[2]['name']
    label2  = leps[2]['label']
    toPlot  = { ("%sPt"       % label1)           : { 'xaxis' : 'p_{%s T} (GeV)'    % name1        , 'rebin' : 10},
                ("%sPt"       % label2)           : { 'xaxis' : 'p_{%s T} (GeV)'    % name2        , 'rebin' : 10},
                ("%sAbsEta"   % label1)           : { 'xaxis' : '|#eta|_{%s} (GeV)' % name1        , 'rebin' : 10},
                ("%sAbsEta"   % label2)           : { 'xaxis' : '|#eta|_{%s} (GeV)' % name2        , 'rebin' : 10},
                ("%s_%s_Mass" % (label1,label2) ) : { 'xaxis' : 'M_{%s%s} (GeV)'    % (name2,name2), 'rebin' : 5 },
                }

    plotter = TauEffPlotter(channel.upper())
    #pprint.pprint(plotter.views)

    #Make Ztt WJets and TTbar region plots
    for iso in leps['ids']:
        for region in ['LoMT','HiMT']:
            for bjet in (['bjets0', 'bjets1', 'bjets2'] if region == 'HiMT' else ['bjets0']):
                for var, kwargs in toPlot.iteritems():
                    folder = '/'.join([iso,'os',region,bjet])
                    plotter.set_subdir(folder)
                    plotter.plot_mc_vs_data(folder, var, **kwargs)
                    plotter.save('mc_vs_data_os_%s_%s_%s_%s' % (iso, region, bjet, var) )
                    if 'os/LoMT/bjets0' in folder:
                        plotter.plot_with_estimate(folder, var, **kwargs)
                        plotter.save('final_data_os_%s_%s_%s_%s' % (iso, region, bjet, var) )
                
    #Make QCD region plots
    for iso in leps['ids']:
        region = 'LoMT'
        bjet   = 'bjets0'
        for var, kwargs in toPlot.iteritems():
            folder = '/'.join([iso,'ss',region,bjet])
            plotter.set_subdir('/'.join([iso,'qcd']))
            plotter.plot_mc_vs_data(folder, var, **kwargs)
            plotter.save('mc_vs_data_os_%s_%s_%s_%s' % (iso, region, bjet, var) )

    #Make WJets region plots
    for iso in leps['ids']:
        region = 'HiMT'
        bjet   = 'bjets0'
        for var, kwargs in toPlot.iteritems():
            folder = '/'.join([iso,'ss',region,bjet])
            plotter.set_subdir('/'.join([iso,'wjets']))
            plotter.plot_mc_vs_data(folder, var, **kwargs)
            plotter.save('mc_vs_data_os_%s_%s_%s_%s' % (iso, region, bjet, var) )


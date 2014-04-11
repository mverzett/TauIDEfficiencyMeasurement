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
import FinalStateAnalysis.Utilities.prettyjson as prettyjson
from TauEffPlotterBase import TauEffPlotterBase, remove_name_entry
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


class TauEffPlotterMM(TauEffPlotterBase):
    def __init__(self):
        super(TauEffPlotterMM, self).__init__('MM')
        self.sample_mapping['Zjets_M50'] = 'zmm'

    def get_qcd_estimation(self):
        qcd_reg   = 'h2Tau/ss'

        #MC ZJets view
        neg_zjet = views.ScaleView(self.get_view('Zjets*'),-1)
        neg_zjet = views.SubdirectoryView(neg_zjet,qcd_reg)

        #DATA WJets estimation
        neg_wjet = views.ScaleView(self.get_view('WplusJets*'),-1)
        neg_wjet = views.SubdirectoryView(neg_wjet,qcd_reg)
        
        #Data view
        data_view= views.SubdirectoryView(self.get_view('data'),qcd_reg)

        qcd_est  = views.SumView(data_view, neg_zjet)
        qcd_est  = views.SumView(qcd_est  , neg_wjet)
        qcd_est  = views.StyleView(qcd_est, **remove_name_entry(data_styles['QCD*']))
        qcd_est  = views.TitleView(qcd_est, 'QCD')
        return qcd_est

    def make_folder_views(self, folder, rebin):
        iso_name  = folder.split('/')[0]
        zjets_mc  = self.get_view_dir('Zjets*'     , rebin, folder)
        data      = self.get_view_dir('data'       , rebin, folder)
        ttbar     = self.get_view_dir('TTplusJets*', rebin, folder)
        wjet_mc   = self.get_view_dir('WplusJets*' , rebin, folder)
        diboson   = views.TitleView(
            views.StyleView(
                views.SumView(
                    self.get_view_dir('WZ*' , rebin, folder),
                    self.get_view_dir('WW*' , rebin, folder),
                    self.get_view_dir('ZZ*' , rebin, folder)
                    ),
                **remove_name_entry(data_styles['WZ*'])
            ),
            'diboson'
            )
        #makes QCD Estimation view
        qcd_est   = self.rebin_view(self.get_qcd_estimation(), rebin)
        return {
            'Z_jets'   : zjets_mc,
            'ttbar'    : ttbar,
            'diboson'  : diboson,
            'WJets'    : wjet_mc,
            'QCD'      : qcd_est,
            'data'     : data,            
            }


    def plot_with_estimate(self, folder, variable, rebin=1, xaxis='', leftside=True,
                           xrange=None, show_error=True,  logscale=False, yrange=None):
        folder_views = self.make_folder_views(folder, rebin)
        zjets_mc = folder_views['Z_jets' ]
        ttbar    = folder_views['ttbar'  ]
        diboson  = folder_views['diboson']
        wjet_mc  = folder_views['WJets'  ]
        qcd_est  = folder_views['QCD'    ]
        data     = folder_views['data'   ]       

        stack     = views.StackView(qcd_est, wjet_mc, ttbar, diboson, zjets_mc).Get(variable)
        stack.Draw()
        stack.GetHistogram().GetXaxis().SetTitle(xaxis)
        if xrange:
            stack.GetXaxis().SetRangeUser(xrange[0], xrange[1])
        if yrange:
            stack.SetMinimum(yrange[0])
            stack.SetMaximum(yrange[1])
            #stack.GetYaxis().SetRangeUser(yrange[0], yrange[1])
        stack.Draw()
        self.keep.append(stack)

        #print os.path.join(folder,variable)
        data_h = data.Get(variable)
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
        self.add_cms_blurb(self.sqrts)
        if logscale:
            self.pad.SetLogy()

        self.add_ratio_plot(data_h, stack, xrange, 0.2)


    def get_signal_views(self, iso_name, variable):
        folder    = 'h2Tau/os/'
        nbins     = self.get_view('Zjets*').Get('/'.join([folder,variable])).GetNbinsX()
        return self.make_folder_views(folder, nbins)
    
    def map_interesting_directories(self, selection_region, var, himt_region = 'MT70_120'):
        ret = {}
        for sign in ['os', 'ss']:
            ret[sign] = self.map_yields(
                os.path.join( selection_region, sign ),
                var
            )
        return ret

    def dump_selection_info(self, tauids, var):
        ret = {}
        for tau_id in tauids:
            ret[tau_id] = self.map_interesting_directories(
                tau_id, 
                var
            )
        return ret

if __name__ <> "__main__":
    sys.exit(0)

jobid = os.environ['jobid']
toPlot  = {
    'nvtx'           : { 'xaxis' : 'N_{vtx}'               , 'rebin' : 1 , 'leftside' : False},
    "m1Pt"           : { 'xaxis' : 'p_{#mu_{1} T} (GeV)'   , 'rebin' : 2 , 'leftside' : False},
    "m2Pt"           : { 'xaxis' : 'p_{#mu_{2} T} (GeV)'   , 'rebin' : 2 , 'leftside' : False},
    "m1AbsEta"       : { 'xaxis' : '|#eta|_{#mu_{1}} (GeV)', 'rebin' : 2 , 'leftside' : False, 'xrange' : (0,3)},
    "m2AbsEta"       : { 'xaxis' : '|#eta|_{#mu_{2}} (GeV)', 'rebin' : 2 , 'leftside' : False, 'xrange' : (0,3)},
    "m1_m2_Mass"     : { 'xaxis' : 'M_{#mu#mu} (GeV)'      , 'rebin' : 1 , 'leftside' : False, 'xrange' : (60,120)},
    "nvtx"           : { 'xaxis' : 'Number of vertices'    , 'rebin' : 1 , 'leftside' : False},
    "type1_pfMetEt"  : { 'xaxis' : 'Type 1 MET E_{T} (GeV)', 'rebin' : 2, 'leftside' : False, 'xrange' : (0,90)},
    "type1_pfMetPhi" : { 'xaxis' : 'Type 1 MET #phi'       , 'rebin' : 1, 'leftside' : False, 'xrange' : (-3.5,3.5)},
    "MET_Z_perp"     : { 'xaxis' : 'Type 1 MET perp. to Z (GeV)', 'rebin' : 2, 'leftside' : False, 'xrange' : (0,80)},
    "MET_Z_para"     : { 'xaxis' : 'Type 1 MET parallel. to Z (GeV)', 'rebin' : 2, 'leftside' : False, 'xrange' : (0,80)},
    }

plotter = TauEffPlotterMM()

print '\n\nPlotting MM\n\n'

folder = 'h2Tau/os'
plotter.set_subdir('signal')
for var, kwargs in toPlot.iteritems():
    ## plotter.plot_mc_vs_data(folder, var, **kwargs)
    ## plotter.save('mc_vs_data_%s_%s' % ('signal_region', var) )
    plotter.plot_with_estimate(folder, var, **kwargs)
    plotter.save('final_%s_%s' % ('signal_region', var)  )
    if var == "m1_m2_Mass":
        kwargs['logscale']=True
        kwargs['yrange']=[1,3*10**8]
        plotter.plot_with_estimate(folder, var, **kwargs)
        plotter.save('final_%s_%s_logscale' % ('signal_region', var) )
        del kwargs['logscale']
        del kwargs['yrange']

plotter.set_subdir('')
plotter.write_summary('','m1_m2_Mass')

yield_dump = plotter.dump_selection_info(['h2Tau'], 'm1_m2_Mass')
with open('results/%s/plots/mm/yield_dump.json' % jobid, 'w') as jfile:
    jfile.write(prettyjson.dumps(yield_dump) )
            
#Make QCD region plots
folder = folder.replace('os','ss')
plotter.set_subdir('qcd')
for var, kwargs in toPlot.iteritems():
    plotter.plot_mc_vs_data(folder, var, **kwargs)
    plotter.save('mc_vs_data_%s_%s' % ('qcd_region', var) )

#FIXME: _understand systamtic uncertainties:
#       _ask Evan for Zrecoil correction in MVA MET
#       _make uncertainties on Zrecoil correction --> propagate to WJets Ztautau QCD ecc...
#       _make #evts passing cuts #of MC events passing cuts (+ stat+sys)
            
            

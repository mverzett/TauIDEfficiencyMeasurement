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


class TauEffPlotterBase(Plotter):
    def __init__(self, channel):
        self.jobid = os.environ['jobid']
        jobid = self.jobid
        self.samples = [ os.path.split(i)[1].split('.')[0] for i in glob.glob('results/%s/TauEffZ%s/*.root' % (jobid, channel))]
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
        super(TauEffPlotterBase, self).__init__(files, lumifiles, self.outputdir, None)
        self.mc_samples = filter(lambda x: not x.startswith('data_'), self.samples)
        self.systematic = ''
        self.shape_systematics = []
        #expressed in %
        self.scale_systematics = {
            '*' : {
                'lumi' : 0.,
                },
            'Z + jets' : {
                'xsec' : 0.,
                }
            }

    def get_view(self, *args): #Is it against Liskov Substitution Principle? I don't care
        if self.systematic != '':
            return views.SubdirectoryView( super(TauEffPlotterBase, self).get_view(*args), self.systematic)
        else:
            return super(TauEffPlotterBase, self).get_view(*args)

    def set_subdir(self, folder):
        self.outputdir = '/'.join([self.base_out_dir, folder])

    def get_view_dir(self, pattern, rebin, folder):
        return views.SubdirectoryView(self.rebin_view(self.get_view(pattern), rebin), folder)

    def plot_mc_vs_data(self, folder, variable, rebin=1, xaxis='', leftside=True,
                        xrange=None, logscale=False, **kwargs): #Is it against Liskov Substitution Principle? I don't care
        super(TauEffPlotterBase, self).plot_mc_vs_data(folder, variable, rebin, xaxis, leftside, xrange)
        if logscale:
            self.canvas.SetLogy()
        self.add_cms_blurb(self.sqrts)

    def write_summary(self, iso_name, variable):
        ''' Write final cut-and-count table of events passing the selection '''
        sys_name  = self.systematic
        all_views = self.get_signal_views(iso_name, variable)
        sys_views = {}
        print self.shape_systematics
        for sys in self.shape_systematics:
            self.systematic = sys
            sys_views[sys] = self.get_signal_views(iso_name, variable)
        self.systematic = sys_name
        
        tab_form  = '%20s%20s%20s%20s\n'
        output    = tab_form % ('sample','# evts.','stat. err.','syst. err.')
        expected_evts = {
            'exp' : [],
            'stat': [],
            'sys' : [],
            }
        
        for name, view in all_views.iteritems():
            hist    = view.Get(variable)
            n_events= hist.GetBinContent(1)
            #compute nevents shape_systematics
            shape_sys_events = [sys_views[sys][name].Get(variable).GetBinContent(1) for sys in self.shape_systematics]
            #makes abs difference
            shape_sys_diffs  = [abs(i-n_events) for i in shape_sys_events]
            #sum in quad
            sys_err = 0.
            if len(shape_sys_diffs):
                sys_err = quad(*shape_sys_diffs)
            if name == 'data':
                output += tab_form % ('Bkg. Sum', '%.1f' % sum(expected_evts['exp' ]), '%.1f' % quad(*expected_evts['stat']), '%.1f' % quad(*expected_evts['sys' ]))
            else:
                expected_evts['exp' ].append(n_events)
                expected_evts['stat'].append(hist.GetBinError(1))
                expected_evts['sys' ].append(sys_err)
            output += tab_form % (name, '%.1f' % n_events, '%.1f' % hist.GetBinError(1), '%.1f' % sys_err)
            
        with open(os.path.join(self.outputdir,'summary_table_%s.raw_txt' % variable),'w') as out_file:
            out_file.write(output)


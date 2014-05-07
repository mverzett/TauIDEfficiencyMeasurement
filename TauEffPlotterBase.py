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
from FinalStateAnalysis.StatTools.quad import quad
from pdb import set_trace
#from FinalStateAnalysis.Utilities.shelve_wrapper  import make_shelf
import json

#from RecoLuminosity.LumiDB import argparse''' 
import sys
import os
import glob
import pprint
import ROOT
import fnmatch
import math
import systematics
import logging

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

def remove_name_entry(dictionary):
    return dict( [ i for i in dictionary.iteritems() if i[0] != 'name'] )

def get_histo_integral(histo):
    nbins = histo.GetNbinsX()
    hclone = histo.Clone()
    hclone.Rebin(nbins)
    return hclone.GetBinContent(1), hclone.GetBinError(1)

class TauEffPlotterBase(Plotter):
    def __init__(self, channel):
        self.jobid = os.environ['jobid']
        jobid = self.jobid
        self.samples = [ os.path.split(i)[1].split('.')[0] for i in glob.glob('results/%s/TauEffZ%s/*.root' % (jobid, channel))]
        self.sample_file = glob.glob('results/%s/TauEffZ%s/*.root' % (jobid, channel))[0] #keep one file name, you will need it
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
        super(TauEffPlotterBase, self).__init__(files, lumifiles, self.outputdir, blinder=None)
        self.mc_samples = filter(lambda x: not x.startswith('data_'), self.samples)
        self.zero_systematics_point = 'NOSYS' if channel=='MT' else ''
        self.systematic = ''
        self.zero_systematic = ''
        self.shape_systematics = []
        self.sample_mapping = {
            'WplusJets*' : 'wjets',
            'Zjets_M50'  : 'ztt'  ,
            'Zjets_ZToMuMu_M50' : 'zmm',
            'TTplusJets*' : 'ttbar',
            'WZ*' : 'wz',
            'WW*' : 'ww',
            'ZZ*' : 'zz',
        }

    def get_view(self, *args): #Is it against Liskov Substitution Principle? I don't care
        if self.systematic != '':
            toget = self.systematic
            if args[0] == 'data' or (not self.systematic.endswith('es_p')):
                toget = self.zero_systematic
            return views.SubdirectoryView( super(TauEffPlotterBase, self).get_view(*args), toget)
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

    def write_summary(self, iso_name, variable, shape_only=False):
        ''' Write final cut-and-count table of events passing the selection '''
        store     = {} 
        sys_name  = self.systematic
        all_views = self.get_signal_views(iso_name, variable)
        sys_views = {}
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
            stat_err= hist.GetBinError(1)
            #compute nevents shape_systematics
            shape_sys_events = [(sys, sys_views[sys][name].Get(variable).GetBinContent(1)) for sys in self.shape_systematics]
            #makes abs difference
            sys_errors       = [('sys_'+sys, abs(i-n_events)) for sys, i in shape_sys_events]
            sys_err          = quad(*[i for _, i in sys_errors])
                    
            if name == 'data':
                store['bkg_sum'] = {
                    'value' : sum(expected_evts['exp' ]),
                    'stat'  : quad(*expected_evts['stat']),
                    'sys'   : quad(*expected_evts['sys' ]),
                    }
                output    += tab_form % ('Bkg. Sum', '%.1f' % sum(expected_evts['exp' ]), '%.1f' % quad(*expected_evts['stat']), '%.1f' % quad(*expected_evts['sys' ]))
                sys_errors = []
                sys_err    = 0.
            else:
                expected_evts['exp' ].append(n_events)
                expected_evts['stat'].append(stat_err)
                expected_evts['sys' ].append(sys_err)

            entry = {
                'value' : n_events,
                'stat'  : stat_err,
                }
            entry.update(dict(sys_errors))
            store[name] = entry
            output += tab_form % (name, '%.1f' % n_events, '%.1f' % hist.GetBinError(1), '%.1f' % sys_err)

        with open(os.path.join(self.outputdir,'summary_table_%s.raw_txt' % variable),'w') as out_file:
            out_file.write(output)
        with open(os.path.join(self.outputdir,'summary_table_%s.json' % variable),'w') as out_file:
            out_file.write(json.dumps(store, indent=4, separators=(',', ': ')))

    def write_json_dump(self, variable):
        ''' Write final cut-and-count table of events passing the selection '''

        def GetContent(dir):
            #print dir
            tempList = dir.GetListOfKeys()
            retList = []
            for it in range(0,tempList.GetSize()):
               retList.append(tempList.At(it).ReadObj())
            return retList

        def MapDirStructure( directory, dirName, objectList ):
            dirContent = GetContent(directory)
            for entry in dirContent:
                if entry.InheritsFrom('TDirectory') or entry.InheritsFrom('TDirectoryFile'):
                    subdirName = os.path.join(dirName,entry.GetName())
                    MapDirStructure(entry, subdirName,objectList)
                elif entry.InheritsFrom('TH1'):
                    pathname = os.path.join(dirName,entry.GetName())
                    objectList.append(pathname)

        tfile  = ROOT.TFile.Open(self.sample_file)
        print self.sample_file
        histos = []
        MapDirStructure(tfile,'',histos) #find all the histograms
        tfile.Close()
        #take only the ones that we are interested in 
        histos = [i for i in histos if variable in i]
        #all the histograms have the same number of bins, store it
        nbins  = self.views['data']['view'].Get(histos[0]).GetNbinsX()
        #remove the trailing directory (systematics)
        if self.zero_systematics_point:
            print 'chopping away the first dir'
            histos = ['/'.join(i.split('/')[1:]) for i in histos]
        #remove duplicates
        histos = list( set( histos ) )

        #where to store all the info
        store     = {}

        #move to nosys region
        self.systematic = self.zero_systematics_point
        #get all the available views
        all_views = dict(
            [
                (name,
                 self.rebin_view(self.get_view(name),nbins)
                 ) for name in self.mc_samples+['data']
                 ]
            )
        sys_views = {}
        for sys in self.shape_systematics:
            #move to that systematic
            self.systematic = sys
            #store the views
            sys_views[sys] = dict(
                [
                    (name,
                     self.rebin_view(self.get_view(name),nbins)
                     ) for name in self.mc_samples+['data']
                     ]
                )

        #now get the bin contents
        for path in histos:
            rerion_name        = path[:path.rfind('/')]
            store[rerion_name] = {}
            #loop on MC/data samples
            for sample, view in all_views.iteritems():
                histogram     = view.Get(path)
                central_value = histogram.GetBinContent(1)
                stat_error    = histogram.GetBinError(1)
                #get all the central values sys shifted
                sys_values    = [sys_views[sys][sample].Get(path).GetBinContent(1) for sys in self.shape_systematics]
                #compute the difference
                sys_values    = [abs(sys-central_value) for name, sys in sys_values]
                sys_error     = quad(*sys_values) if sample != 'data' else 0.
                store[rerion_name][sample] = {
                    'val'  : central_value,
                    'stat' : stat_error,
                    'sys'  : sys_error,
                    }

        #write to file
        with open(os.path.join(self.outputdir,'full_content_dump.json'),'w') as out_file:
            out_file.write(json.dumps(store, indent=4, separators=(',', ': ')))
 
    def get_nicer_name(self, sample):
        for regex, nicename in self.sample_mapping.iteritems():
            if fnmatch.fnmatch(sample, regex):
                return nicename
        return sample

    def map_yields(self, path, var, systematics_hook=None):
        samples = self.mc_samples + ['data']
        fullpath = os.path.join(path, var)
        ret = {}
        for sample in samples:
            #get central values
            sample_name = self.get_nicer_name(sample)
            self.systematic = self.zero_systematics_point
            central_value, stat_err = get_histo_integral(
                self.get_view(sample).Get(fullpath)
            )
            ret[sample_name] = {
                'val'  : central_value,
                'stat' : stat_err,
                }
            
            #get sys shifts
            for sys_shift in self.shape_systematics:
                sys_name = systematics_hook(sys_shift) if systematics_hook else sys_shift
                if sys_name: #allows to filter systematics
                    self.systematic = sys_shift
                    sys_central, unused = get_histo_integral(
                        self.get_view(sample).Get(fullpath)
                    )
                    sys_err = sys_central - central_value 
                    ret[sample_name][sys_name] = sys_err
        return ret
        

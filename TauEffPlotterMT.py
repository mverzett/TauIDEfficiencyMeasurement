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
from TauEffPlotterBase import TauEffPlotterBase
import itertools
import sys
import os
import glob
import pprint
import ROOT
import fnmatch
import math
import logging

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)

class SystematicsView(object):
    def __init__(self, central_view, *systematics_views):
        self.central_view = central_view
        self.systematics_views = systematics_views

    @staticmethod
    def quad(*xs):
        return math.sqrt(sum(x*x for x in xs))
    
    def Get(self,path):
        central_hist = self.central_view.Get(path)
        systematics_hists = [i.Get(path) for i in self.systematics_views]
        ret_hist = central_hist.Clone()
        for binno in range(1, high_hist.GetNbinsX() + 1):
            error = quad(
                central_hist.GetBinError(binno),
                *[ abs(central_hist.GetBinContent(binno) - i.GetBinContent(binno)) for i in systematics_hists]
                )
            ret_hist.SetBinError(binno,error)

class TauEffPlotterMT(TauEffPlotterBase):
    def __init__(self):
        super(TauEffPlotterMT, self).__init__('MT')
        self.systematic = 'NOSYS'
        self.shape_systematics = ['mes_p','tes_p','jes_p','ues_p']
        #set the proper style
        self.views['Zjets_ZToMuMu_M50']['view'] = views.TitleView(
            views.StyleView(self.views['Zjets_ZToMuMu_M50']['view'], **data_styles['WW*']),
            'Z#rightarrow#mu#mu'
            )

    def compute_contributions(self, iso_name, himtregion):
        '''uses ABCD method to extrapolate the WJets and QCD yield in SS region. Regions:
        A - AntiIso objects LowMt
        B - AntiIso objects HihgMt
        C - Iso pass LowMt
        D - Iso pass HighMt'''
        data = self.get_view('data')
        wjet = self.get_view('WplusJets*')
        
        qcd_ss_noIso_loMt = data.Get('ss/QCD/LoMT/mPt').Integral()
        qcd_ss_noIso_hiMt = data.Get('ss/QCD/%s/mPt' % himtregion).Integral()
        qcd_ratio_lo_hi_mt= qcd_ss_noIso_loMt / qcd_ss_noIso_hiMt

        data_ss_Iso_loMt = data.Get('%s/ss/LoMT/mPt' % iso_name).Integral()
        data_ss_Iso_hiMt = data.Get('%s/ss/%s/mPt' % (iso_name, himtregion)).Integral()
        
        wjet_ss_Iso_loMt = wjet.Get('%s/ss/LoMT/mPt' % iso_name).Integral()
        wjet_ss_Iso_hiMt = wjet.Get('%s/ss/%s/mPt' % (iso_name, himtregion)).Integral()
        w_ratio_lo_hi_mt = wjet_ss_Iso_loMt / wjet_ss_Iso_hiMt

        print himtregion, data_ss_Iso_hiMt, w_ratio_lo_hi_mt, data_ss_Iso_loMt, (w_ratio_lo_hi_mt - qcd_ratio_lo_hi_mt)
        qcd_yield_hiMt = (data_ss_Iso_hiMt*w_ratio_lo_hi_mt - data_ss_Iso_loMt)/(w_ratio_lo_hi_mt - qcd_ratio_lo_hi_mt)
        qcd_yield_loMt = qcd_yield_hiMt*qcd_ratio_lo_hi_mt

        return {
            'LoMT' : {
                'yield' : qcd_yield_loMt,
                'view'  : views.TitleView(
                    views.StyleView(
                        views.ScaleView(
                            views.SubdirectoryView(
                                data,
                                'ss/QCD/LoMT'
                                ),
                                qcd_yield_loMt/qcd_ss_noIso_loMt
                                )
                        , **data_styles['QCD*']),
                        'QCD'),                        
                },
            himtregion : {
                'yield' : qcd_yield_hiMt,
                'view'  : views.TitleView(
                    views.StyleView(
                        views.ScaleView(
                            views.SubdirectoryView(
                                data,
                                'ss/QCD/HiMT'
                                ),
                                qcd_yield_loMt/qcd_ss_noIso_hiMt
                                )
                        , **data_styles['QCD*']),
                        'QCD'),
                },
            }
        pass
        
        
        
    def get_w_estimation(self, iso_name, mtregion, sign_region='os'):
        mt_h_name= 'mMtToPfMet_Ty1'
        #mtregion = 'HiMT'
        wjet_reg =  '/'.join((iso_name, sign_region, mtregion,''))

        #qcd contribution
        qcd_view = self.compute_contributions(iso_name, mtregion)[mtregion]['view']
        qcd_view = views.ScaleView(qcd_view, -1)
        
        #Scale factor
        wjets_mc = self.get_view('WplusJets*')
        wjets_hi = wjets_mc.Get(wjet_reg+mt_h_name) 
        wjets_lo = wjets_mc.Get(wjet_reg.replace(mtregion,'LoMT')+mt_h_name)
        w_scale  = wjets_lo.Integral() / wjets_hi.Integral()

        #MC ZJets view
        neg_zjet = views.ScaleView(self.get_view('Zjets*'),-1)
        neg_zjet = views.SubdirectoryView(neg_zjet,wjet_reg)

        #Data view
        data_view= views.SubdirectoryView(self.get_view('data'),wjet_reg)

        wjet_est = views.SumView(data_view, neg_zjet)
        wjet_est = views.SumView(wjet_est, qcd_view)
        wjet_est = views.ScaleView(wjet_est, w_scale)
        wjet_est = views.StyleView(wjet_est, **data_styles['WplusJets*'])
        wjet_est = views.TitleView(wjet_est, 'WplusJets')
        return wjet_est

    def get_qcd_estimation(self, iso_name, mtregion):
        qcd_reg   = '%s/ss/LoMT/' % iso_name

        #MC ZJets view
        neg_zjet = views.ScaleView(self.get_view('Zjets*'),-1)
        neg_zjet = views.SubdirectoryView(neg_zjet,qcd_reg)

        #DATA WJets estimation
        neg_wjet = views.ScaleView(self.get_w_estimation(iso_name, mtregion, 'ss'), -1)
        
        #Data view
        data_view= views.SubdirectoryView(self.get_view('data'),qcd_reg)

        qcd_est  = views.SumView(data_view, neg_zjet)
        qcd_est  = views.SumView(data_view, neg_wjet)
        qcd_est  = views.StyleView(qcd_est, **data_styles['QCD*'])
        qcd_est  = views.TitleView(qcd_est, 'QCD')
        return qcd_est

    def make_folder_views(self, folder, rebin, mtregion='HiMT'):
        iso_name    = folder.split('/')[0]
        zjets_mc    = self.get_view_dir('Zjets_M50'  , rebin, folder)
        data        = self.get_view_dir('data'       , rebin, folder)
        ttbar       = self.get_view_dir('TTplusJets*', rebin, folder)
        zjets_mm_mc = self.get_view_dir('Zjets_ZToMuMu_M50'  , rebin, folder)

        diboson   = views.TitleView(
            views.StyleView(
                views.SumView(
                    self.get_view_dir('WZ*' , rebin, folder),
                    self.get_view_dir('WW*' , rebin, folder),
                    self.get_view_dir('ZZ*' , rebin, folder)
                    ),
                **data_styles['WZ*']
            ),
            'diboson'
            )

        zjets_mc    = views.TitleView(zjets_mc, 'Z#rightarrow#tau#tau')
        zjets_mm_mc = views.TitleView(
            views.StyleView(zjets_mm_mc, **data_styles['WW*']),
            'Z#rightarrow#mu#mu'
            )
                
        #makes WJets Estimation view
        wjet_est  = self.rebin_view(self.get_w_estimation(iso_name, mtregion), rebin)
        
        #makes QCD Estimation view
        qcd_est   = self.rebin_view(self.get_qcd_estimation(iso_name, mtregion), rebin)
        return {
            'Z_tautau' : zjets_mc,
            'Z_mumu'   : zjets_mm_mc,
            'ttbar'    : ttbar,
            'diboson'  : diboson,
            'WJets'    : wjet_est,
            'QCD'      : qcd_est,
            'data'     : data,
            }

    def plot_wjets_region(self, folder, variable, rebin=1, xaxis='', leftside=True,
                        xrange=None, logscale=False, **kwargs):
        iso_name    = folder.split('/')[0]
        mtregion    = folder.split('/')[-1]
        zjets_mc    = self.get_view_dir('Zjets_M50'  , rebin, folder)
        data        = self.get_view_dir('data'       , rebin, folder)
        ttbar       = self.get_view_dir('TTplusJets*', rebin, folder)
        zjets_mm_mc = self.get_view_dir('Zjets_ZToMuMu_M50'  , rebin, folder)
        wjets_mc    = self.get_view_dir('WplusJets*',  rebin, folder)
        qcd_est     = self.rebin_view(self.compute_contributions(iso_name, mtregion)[mtregion]['view'], rebin)
        diboson     = views.TitleView(
            views.StyleView(
                views.SumView(
                    self.get_view_dir('WZ*' , rebin, folder),
                    self.get_view_dir('WW*' , rebin, folder),
                    self.get_view_dir('ZZ*' , rebin, folder)
                    ),
                **data_styles['WZ*']
            ),
            'diboson'
            )
        mc_stack = views.StackView(ttbar, diboson, qcd_est, zjets_mm_mc, zjets_mc, wjets_mc).Get(variable)
        mc_stack.Draw()
        mc_stack.GetHistogram().GetXaxis().SetTitle(xaxis)
        if xrange:
            mc_stack.GetXaxis().SetRangeUser(xrange[0], xrange[1])
            mc_stack.Draw()
        self.keep.append(mc_stack)
        data = self.rebin_view(self.get_view('data'), rebin).Get(os.path.join(folder,variable))
        data.Draw('same')
        self.keep.append(data)
        # Make sure we can see everything
        if data.GetMaximum() > mc_stack.GetMaximum():
            mc_stack.SetMaximum(1.2*data.GetMaximum())
        # Add legend
        self.add_legend([data, mc_stack], leftside, entries=7)


        if logscale:
            self.canvas.SetLogy()
        self.add_cms_blurb(self.sqrts)


    def plot_with_estimate(self, folder, variable, rebin=1, xaxis='', leftside=True, xrange=None, show_error=True):
        folder_views = self.make_folder_views(folder, rebin)
        zjets_mc     = folder_views['Z_tautau']
        zjets_mm_mc  = folder_views['Z_mumu'  ]
        ttbar        = folder_views['ttbar'   ]
        diboson      = folder_views['diboson' ]
        wjet_est     = folder_views['WJets'   ]
        qcd_est      = folder_views['QCD'     ]
        data         = folder_views['data'    ]
        
        stack     = views.StackView(ttbar, diboson, qcd_est, wjet_est,  zjets_mm_mc, zjets_mc).Get(variable)
        stack.Draw()
        stack.GetHistogram().GetXaxis().SetTitle(xaxis)
        if xrange:
            stack.GetXaxis().SetRange(xrange[0], xrange[1])
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

    def get_signal_views(self, iso_name, variable):
        folder    = '%s/os/LoMT' % iso_name
        nbins     = self.get_view('Zjets*').Get('/'.join([folder,variable])).GetNbinsX()
        return self.make_folder_views(folder, nbins)
        
        
if __name__ <> "__main__":
    sys.exit(0)

jobid = os.environ['jobid']
toPlot  = {
    "mPt"      : { 'xaxis' : 'p_{#mu T} (GeV)'        , 'rebin' : 10, 'leftside' : False},
    "tPt"      : { 'xaxis' : 'p_{#tau_{h} T} (GeV)'   , 'rebin' : 10, 'leftside' : False},
    "mAbsEta"  : { 'xaxis' : '|#eta|_{#mu} (GeV)'     , 'rebin' : 10, 'leftside' : False},
    "tAbsEta"  : { 'xaxis' : '|#eta|_{#tau_{h}} (GeV)', 'rebin' : 10, 'leftside' : False},
    "m_t_Mass" : { 'xaxis' : 'M_{#mu#tau_{h}} (GeV)'  , 'rebin' : 5 , 'leftside' : False},
    "nvtx"       : { 'xaxis' : 'Number of vertices'    , 'rebin' : 1 , 'leftside' : False},
    }

plotter = TauEffPlotterMT()

print '\n\nPlotting MT\n\n'

#pprint.pprint(plotter.views)

ids =  [
            'LooseIso'    ,
            'MediumIso'   ,
            'TightIso'    , 
            'LooseMVAIso' ,
            'MediumMVAIso',
            'TightMVAIso' ,
            'LooseIso3Hits',
            'LooseMVA2Iso',
            'MediumIso3Hits',
            'MediumMVA2Iso',
            'TightIso3Hits',
            'TightMVA2Iso',
            'VLooseIso',
    ]

#Make Ztt WJets and TTbar region plots
for iso in ids:
    region = 'LoMT'
    for var, kwargs in toPlot.iteritems():
        folder = '/'.join([iso,'os',region])
        plotter.set_subdir(os.path.join(iso, 'signal'))
        plotter.plot_mc_vs_data(folder, var, **kwargs)
        plotter.save('mc_vs_data_os_%s_%s_%s' % (iso, region, var) )
        print 'plot_with_estimate: ',folder, var
        plotter.plot_with_estimate(folder, var, **kwargs)
        plotter.save('final_data_os_%s_%s_%s' % (iso, region, var) )
    plotter.set_subdir(iso)
    plotter.write_summary(iso,'m_t_Mass')
    #plotter.write_summary(iso,'mAbsEta')

    for region, dirname in [
            ('HiMT', 'MT_Greater_40'),
            ('VHiMT', 'MT_Greater_70'),
            ('MT70_120', 'MT_70_120'),
            ]:
        for var, kwargs in toPlot.iteritems():
            plotter.set_subdir(os.path.join(iso, dirname))
            folder = '/'.join([iso,'os',region])
            plotter.plot_mc_vs_data(folder, var, **kwargs)
            plotter.save('mc_vs_data_os_%s_%s_%s' % (iso, region, var) )
            print 'plot_wjets_region: ',folder, var
            plotter.plot_wjets_region(folder, var, **kwargs)
            plotter.save('final_os_%s_%s_%s' % (iso, region, var) )


plotter.set_subdir('')
#plotter.write_json_dump('m_t_Mass')


#Make QCD region plots
for iso in ids:
    region = 'LoMT'
    for var, kwargs in toPlot.iteritems():
        folder = '/'.join([iso,'ss',region])
        plotter.set_subdir('/'.join([iso,'qcd']))
        plotter.plot_mc_vs_data(folder, var, **kwargs)
        plotter.save('mc_vs_data_os_%s_%s_%s' % (iso, region, var) )

integral_ss = sum([plotter.views['data']['view'].Get('NOSYS/ss/QCD/%s/mPt' % m).Integral() for m in ['LoMT']])
integral_os = sum([plotter.views['data']['view'].Get('NOSYS/os/QCD/%s/mPt' % m).Integral() for m in ['LoMT']])
ratio_ss_os = integral_ss/integral_os

print '(ss = %.1f +- %.1f) / (os = %.1f +- %.1f) = %.4f +- %.4f' % (integral_ss, math.sqrt(integral_ss), integral_os, math.sqrt(integral_os), ratio_ss_os, ratio_ss_os*quad(1./math.sqrt(integral_ss)+1./math.sqrt(integral_os)))


#FIXME: _understand systamtic uncertainties:
#       _ask Evan for Zrecoil correction in MVA MET
#       _make uncertainties on Zrecoil correction --> propagate to WJets Ztautau QCD ecc...
#       _make #evts passing cuts #of MC events passing cuts (+ stat+sys)
            
            
